#!/usr/bin/env python3
"""Run pairwise and germinal center zone analyses on the finalized tonsil profiles."""
from __future__ import annotations

import argparse
from collections import defaultdict
from copy import copy
import gzip
import hashlib
import json
import logging
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory
import time
from types import SimpleNamespace
from urllib.parse import unquote

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt, colors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy import sparse as sp, stats

from tonsil_processing import (
    HIGH_COVERAGE, LABEL, PROTOCOLS, ROOT, UNKNOWN, ad_module, digest, ensure_references,
    h5_read_elem, identity, label_values, load_study, log, output_lock, read_meta,
    require_files, save_h5ad, thread_budget, tools_required, write_json, write_table,
)


# Inputs and differential analysis utilities

MARKS = ('H3K27ac', 'H3K27me3')

SCOPES = ('merged', 'Tonsil1', 'Tonsil1_gt600')


def slug(text):
    # Hash prevents DC_CDC2 and DC-CDC2, for example, sharing a directory.
    return (re.sub(r'[^A-Za-z0-9]+','_',str(text)).strip('_') or 'label') + '_' + hashlib.sha256(str(text).encode()).hexdigest()[:8]


def aggregate(x,codes,n_units):
    codes=np.asarray(codes,dtype=int);valid=codes>=0
    membership=sp.csr_matrix((np.ones(int(valid.sum()),dtype=np.int64),(codes[valid],np.flatnonzero(valid))),shape=(n_units,x.shape[0]))
    return np.asarray((membership@x).toarray(),dtype=np.int64).T


def h5_rows(node,rows):
    """Read selected rows of a CSR H5AD without SciPy/AnnData backed-index quirks.

    Native fragment CSR can have duplicate indices with different lengths.
    Copy contiguous storage ranges verbatim; never canonicalize this data.
    """
    import h5py
    rows=np.asarray(rows,dtype=np.int64)
    if isinstance(node,h5py.Dataset):
        order=np.argsort(rows,kind='stable'); inverse=np.argsort(order,kind='stable')
        return np.asarray(node[rows[order],:])[inverse]
    encoding=node.attrs.get('encoding-type','')
    if isinstance(encoding,bytes):encoding=encoding.decode()
    shape=tuple(map(int,node.attrs['shape']))
    ptr=np.asarray(node['indptr'])
    if encoding=='csr_matrix':
        order=np.argsort(rows,kind='stable'); sorted_rows=rows[order]
        sizes=ptr[sorted_rows+1]-ptr[sorted_rows]
        target_ptr=np.concatenate([[0],np.cumsum(sizes,dtype=np.int64)])
        data=np.empty(int(target_ptr[-1]),dtype=node['data'].dtype)
        indices=np.empty(int(target_ptr[-1]),dtype=np.int64)
        boundaries=np.r_[0,np.flatnonzero(np.diff(sorted_rows)!=1)+1,len(rows)]
        for first,last in zip(boundaries[:-1],boundaries[1:]):
            if first==last:continue
            start,end=int(ptr[sorted_rows[first]]),int(ptr[sorted_rows[last-1]+1])
            dst=slice(int(target_ptr[first]),int(target_ptr[last]))
            data[dst]=node['data'][start:end];indices[dst]=node['indices'][start:end]
        result=sp.csr_matrix((data,indices,target_ptr),shape=(len(rows),shape[1]))
        return result[np.argsort(order,kind='stable'),:] if not np.array_equal(order,np.arange(len(rows))) else result
    if encoding=='csc_matrix':
        # RNA count matrices only: stream columns, retaining selected rows.
        mapping=np.full(shape[0],-1,dtype=np.int64);mapping[rows]=np.arange(len(rows))
        blocks=[]
        for first in range(0,shape[1],512):
            last=min(first+512,shape[1]);start,end=int(ptr[first]),int(ptr[last])
            source_rows=np.asarray(node['indices'][start:end]);target_rows=mapping[source_rows]
            cols=np.repeat(np.arange(last-first),np.diff(ptr[first:last+1]))
            values=np.asarray(node['data'][start:end]);keep=target_rows>=0
            blocks.append(sp.csr_matrix((values[keep],(target_rows[keep],cols[keep])),shape=(len(rows),last-first)))
        return sp.hstack(blocks,format='csr')
    raise ValueError(f'Unsupported sparse encoding: {encoding}')


def read_fragment_rows(path, key, rows):
    import h5py
    ad_module()
    with h5py.File(path,'r') as f:
        fragments = h5_rows(f['obsm'][key],rows)
        reference = h5_read_elem(f['uns']['reference_sequences'])
    return fragments,reference


def physical_ids(obs, sample=None):
    if 'merged_cell_id' in obs: values = obs.merged_cell_id.astype(str)
    elif obs.index.astype(str).str.contains('::', regex=False).all(): values = obs.index.astype(str)
    elif {'sample','source_cell_id'}.issubset(obs):
        values = obs['source_cell_id'].astype(str)
        values = values.where(values.str.contains('::', regex=False), obs['sample'].astype(str)+'::'+values)
    elif sample is not None:
        values = str(sample)+'::'+obs.index.astype(str)
    else: raise ValueError('Exact cell IDs required: merged_cell_id or sample::source_cell_id')
    idx = pd.Index(values, name='cell_id')
    if not idx.is_unique or not idx.str.contains('::', regex=False).all(): raise ValueError('Invalid or duplicate physical cell IDs')
    return idx


def raw_integer(x):
    from scipy import sparse as sp
    x = sp.csr_matrix(x)
    if x.data.size and (not np.isfinite(x.data).all() or (x.data < 0).any()
            or not np.equal(x.data, np.rint(x.data)).all()):
        raise ValueError('Raw nonnegative integer counts are required')
    return x.astype(np.int64)


def normalize_dataset(value):
    key=str(value).lower().replace('_','').replace('-','')
    aliases={'merged':'merged','tonsil1':'Tonsil1','tonsilgt600':'Tonsil1_gt600',
             'tonsil1gt600':'Tonsil1_gt600','gt600':'Tonsil1_gt600'}
    if key not in aliases:raise ValueError('Choose merged, tonsil1 or tonsilgt600')
    return aliases[key]


def comparison_key(first,second):
    def readable(label):
        return (re.sub(r'[^A-Za-z0-9]+','_',label).strip('_') or 'label')[:50]
    suffix=hashlib.sha256(json.dumps([first,second],ensure_ascii=True).encode()).hexdigest()[:12]
    return f'{readable(first)}_vs_{readable(second)}_{suffix}'


def title_for(scope):
    return {'merged':'Merged tonsils', 'Tonsil1':'Tonsil1', 'Tonsil1_gt600':'Tonsil1 >600'}[scope]


def object_path(final, publication, scope, modality):
    if scope == 'Tonsil1_gt600':
        suffix = 'counts' if modality == 'RNA' else 'tiles'
        return publication / scope / 'labelled_only/objects' / f'Tonsil1_{modality}_gt600_{suffix}.h5ad'
    return final / 'labelled_only/objects' / f'{modality}_merged.h5ad'


def chosen_rows(s, ids):
    local = s['obs'].index.get_indexer(ids)
    if (local < 0).any():
        raise ValueError(f'{s["path"]}: could not map selected physical cell IDs')
    return s['rows'][local]


def select(s, mask):
    return dict(s, obs=s['obs'].loc[mask].copy(), rows=s['rows'][np.asarray(mask, dtype=bool)])


def count_matrix(s, rows=None, layer='X'):
    import h5py
    rows = s['rows'] if rows is None else np.asarray(rows, dtype=np.int64)
    if not len(rows):
        raise ValueError(f'{s["path"]}: no selected count rows')
    with h5py.File(s['path'], 'r') as handle:
        node = handle['X'] if layer == 'X' else handle['layers'][layer]
        matrix = raw_integer(h5_rows(node, rows))
        var = h5_read_elem(handle['var'])
    if identity(s['path']) != s['stamp']:
        raise ValueError(f'{s["path"]}: changed while loading counts')
    return matrix, var


def rna_symbols(x, var, scope, *, sorted_genes=False):
    column = next((k for k in ('gene_symbol','gene_symbols','symbol','gene_name') if k in var), None)
    ids = pd.Series(var.index.astype(str), index=var.index, dtype='string')
    names = var[column].astype('string').fillna('').str.strip() if column else ids
    names = names.mask(names.str.lower().isin(UNKNOWN | {'.'}), ids)
    upper = names.str.upper()
    keep = ~upper.str.contains('MALAT1',regex=False) & ~upper.str.startswith(('MT-','MT_','MT.')) & ~upper.str.contains('ENSG',regex=False)
    cols = np.flatnonzero(keep.to_numpy())
    genes = pd.Index(pd.unique(names.iloc[cols]), name='gene')
    if sorted_genes:
        genes = genes.sort_values()
    if not len(genes):
        raise ValueError('No resolved, nonmitochondrial RNA features left for DGE')
    projection = sp.csr_matrix((np.ones(len(cols), dtype=np.int64), (cols,genes.get_indexer(names.iloc[cols]))),
                               shape=(len(var),len(genes)))
    totals = np.asarray(x.sum(axis=1)).ravel()
    if scope == 'Tonsil1_gt600' and not (totals > HIGH_COVERAGE).all():
        raise ValueError('Saved >600 RNA object contains a selected cell with ≤600 raw UMIs')
    return (x @ projection).tocsr(), genes, totals


def save_pdf(fig, path):
    from matplotlib import pyplot as plt
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.stem + '.tmp.pdf')
    try:
        fig.savefig(temporary, bbox_inches='tight', format='pdf')
        temporary.replace(path)
    finally:
        plt.close(fig)
        temporary.unlink(missing_ok=True)


def plotting_style():
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams.update({'pdf.fonttype':42, 'ps.fonttype':42, 'font.family':'DejaVu Sans',
                                'axes.spines.top':False, 'axes.spines.right':False})


def labels(obs):
    if LABEL not in obs:
        raise ValueError(f'Expected {LABEL!r}; use the stage 02 objects')
    return label_values(obs[LABEL])


def final_fragment_path(final, sample, mark):
    return final / 'donors' / sample / 'fragments' / f'{sample}_{mark}.final_pass.fragments.tsv.gz'


def gene_annotation():
    import snapatac2 as snap
    chosen = Path(snap.genome.hg38.annotation)
    require_files(chosen)
    return chosen


def check_r():
    tools_required('Rscript')
    subprocess.run(['Rscript','-e','stopifnot(requireNamespace("DESeq2", quietly=TRUE))'], check=True)


def fit_model(counts, meta, model, paired, fdr, settings=()):
    if not counts.index.is_unique or not counts.columns.equals(meta.index):
        raise ValueError('Invalid pseudobulk gene or replicate order')
    with TemporaryDirectory(prefix='tonsil-dge-') as name:
        folder = Path(name)
        counts.to_csv(folder/'counts.tsv.gz',sep='\t')
        meta.to_csv(folder/'units.tsv',sep='\t')
        result = subprocess.run(['Rscript',str(ROOT/'differential.R'),model,name,
                                 'paired' if paired else 'random',str(fdr),*map(str,settings)],
                                capture_output=True,text=True)
        result_file = folder/'DESeq2_results.tsv.gz'
        if result.returncode or not result_file.is_file():
            raise RuntimeError(f'DESeq2 {model} fit failed:\n{result.stderr.strip()}\n{result.stdout.strip()}')
        numeric = ('logFC','baseMean','logCPM','PValue','FDR','lfcSE','stat','dispersion')
        table = pd.read_csv(result_file,sep='\t',dtype={'gene':str},keep_default_na=False,
                            na_values={col:['NA'] for col in numeric})
        if not pd.Index(table.gene).equals(pd.Index(counts.index.astype(str))):
            raise ValueError('Unexpected gene identities or order in DESeq2 output')
        return table


def source(final, publication, scope, modality):
    path = object_path(final, publication, scope, modality)
    if not path.is_file():
        raise FileNotFoundError(f'Missing saved final {modality} object for {scope}: {path}')
    stamp = identity(path)
    original, keys, native_ref = read_meta(path)
    obs = original.copy()
    obs['_source_index'] = obs.index.astype(str)
    obs.index = physical_ids(obs, 'Tonsil1' if scope != 'merged' else None)
    prefixes = pd.Series(obs.index.str.split('::').str[0], index=obs.index)
    donors = obs['sample'].astype(str) if 'sample' in obs else prefixes
    if not donors.isin(PROTOCOLS).all() or not (donors == prefixes).all():
        raise ValueError(f'{path}: sample column and exact cell ID prefixes differ')
    if 'tagmentation_protocol' in obs and not obs.tagmentation_protocol.astype(str).eq(donors.map(PROTOCOLS)).all():
        raise ValueError(f'{path}: inconsistent tagmentation protocol provenance')
    obs['sample'] = obs['biological_donor'] = donors.to_numpy()
    obs[LABEL] = labels(obs).to_numpy()
    mask = np.ones(len(obs), dtype=bool) if scope == 'merged' else donors.eq('Tonsil1').to_numpy()
    if scope == 'Tonsil1_gt600':
        if not mask.all():
            raise ValueError(f'{path}: expected Tonsil1 only')
        if modality != 'RNA':
            metric = next((k for k in ('publication_gate_value','unique_fragments','n_unique_fragments',
                          'n_fragment','n_fragments','counts') if k in obs), None)
            if metric is None or not pd.to_numeric(obs[metric], errors='coerce').gt(HIGH_COVERAGE).all():
                raise ValueError(f'{path}: cannot verify strict >600 DNA gate')
    if identity(path) != stamp:
        raise ValueError(f'{path}: changed while loading metadata')
    return dict(path=path, stamp=stamp, keys=keys, native_ref=native_ref, obs=obs.loc[mask].copy(),
                rows=np.flatnonzero(mask))


# Stage 03: pairwise contrasts and shared gene signals

def fit_umap(rep, modality, args):
    import umap
    if len(rep) < 3:
        raise ValueError('At least three cells required to recalculate UMAP')
    metric = 'euclidean' if modality == 'RNA' else args.dna_metric
    return umap.UMAP(n_neighbors=min(args.neighbors, len(rep)-1), n_components=2,
        min_dist=args.min_dist, metric=metric, random_state=args.seed, n_jobs=1,
        init='random').fit_transform(rep)


def random_splits(x,obs,genes,target,args):
    condition=np.where(obs[LABEL].eq(target),'cluster','rest')
    units,codes=[],np.full(len(obs),-1,dtype=int)
    rng=np.random.default_rng(int(digest([args.seed,target])[:8],16))
    for name in ('rest','cluster'):
        idx=np.flatnonzero(condition==name)
        if len(idx)<args.min_cells_condition: raise ValueError('Too few cells for three disjoint random subsets')
        for rep,subset in enumerate(np.array_split(rng.permutation(idx),args.random_replicates),1):
            if not len(subset): continue
            unit=f'{name}__rep{rep}';codes[subset]=len(units)
            units.append(dict(unit=unit,condition=name,donor='pooled_donors',n_cells=len(subset),unit_type='random_cell_split'))
    meta=pd.DataFrame(units).set_index('unit')
    counts=pd.DataFrame(aggregate(x,codes,len(meta)),index=genes,columns=meta.index)
    return counts,meta


def donor_pairs(x,obs,genes,target,args):
    units,codes=[],np.full(len(obs),-1,dtype=int)
    for donor in sorted(obs.biological_donor.unique()):
        d=obs.biological_donor.eq(donor).to_numpy();a=obs[LABEL].eq(target).to_numpy()
        n_target=int((d&a).sum());n_rest=int((d&~a).sum())
        included=min(n_target,n_rest)>=args.min_cells_condition
        if not included: continue
        for name,mask in [('rest',d&~a),('cluster',d&a)]:
            codes[mask]=len(units)
            units.append(dict(unit=f'{donor}__{name}',condition=name,donor=donor,n_cells=int(mask.sum()),
                              protocol=PROTOCOLS[donor],unit_type='biological_donor_by_condition'))
    meta=pd.DataFrame(units,columns=['unit','condition','donor','n_cells','protocol','unit_type']).set_index('unit')
    counts=pd.DataFrame(aggregate(x,codes,len(meta)),index=genes,columns=meta.index)
    return counts,meta


def annotate_result(table,raw,obs,genes,method,args):
    target=obs[LABEL].eq(args.target_label).to_numpy()
    detection=pd.DataFrame(dict(gene=genes,
        pct_detected_target=np.asarray((raw[target]>0).mean(axis=0)).ravel(),
        pct_detected_rest=np.asarray((raw[~target]>0).mean(axis=0)).ravel()))
    table=table.merge(detection,on='gene',how='left',validate='one_to_one')
    table['pct_detected_reference']=table['pct_detected_rest']
    table['target_label']=args.target_label;table['reference_label']=args.reference_label
    table['contrast']=f'{args.target_label} / {args.reference_label}'
    table['n_target_cells']=int(target.sum());table['n_reference_cells']=int((~target).sum())
    table['n_model_target_cells']=int(target.sum());table['n_model_reference_cells']=int((~target).sum())
    table['method']=method
    table['analysis_unit']='biological donor; paired ~ donor + condition' if method=='donor_replicates' else (
        'pooled random cell subset across donors; exploratory' if args.scope=='merged' else
        'disjoint random cell subsets from one donor; exploratory')
    return table


def select_pair_markers(table,args):
    """Both directions use the same RNA-only FDR/effect/detection criteria."""
    finite=np.isfinite(table.logFC)&np.isfinite(table.FDR)
    selected=[]
    for group,direction,detection in ((args.target_label,1,'pct_detected_target'),
                                      (args.reference_label,-1,'pct_detected_rest')):
        eligible=table.loc[finite&table.FDR.le(args.fdr)&
            (direction*table.logFC).ge(args.min_lfc)&table[detection].ge(args.min_detection)].copy()
        eligible['_effect']=direction*eligible.logFC
        eligible=eligible.sort_values(['FDR','_effect','gene'],ascending=[True,False,True]).head(args.top_genes)
        eligible=eligible.drop(columns='_effect').assign(marker_for=group,
            direction='higher_in_target' if direction==1 else 'higher_in_reference',
            rank_within_direction=np.arange(1,len(eligible)+1))
        selected.append(eligible)
    return pd.concat(selected,ignore_index=True)


def normalize_rows(x, scale, totals=None):
    x = sp.csr_matrix(x, dtype=float)
    totals = np.asarray(x.sum(axis=1)).ravel() if totals is None else np.asarray(totals, dtype=float)
    if len(totals) != x.shape[0] or (totals < 0).any() or not np.isfinite(totals).all():
        raise ValueError('Invalid normalization denominator')
    weights = np.divide(scale, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
    return (sp.diags(weights) @ x).tocsr(), totals


def color_cap(values, percentile):
    values = np.asarray(values, dtype=float)
    positive = values[np.isfinite(values) & (values > 0)]
    return float(np.percentile(positive, percentile)) if len(positive) else 0.


def get_gene(signal, gene, ids):
    rows = signal['obs'].index.get_indexer(ids)
    if (rows < 0).any():
        raise ValueError('Signal/UMAP cells differ; exact IDs are required')
    if gene not in signal['genes']:
        return np.full(len(ids), np.nan)
    column = signal['genes'].get_loc(gene)
    # Slice columns first. Signals are CSC in the controller; this avoids a
    # complete sparse row copy for every gene and plasma state.
    return np.asarray(signal['x'][:, column].toarray()).ravel()[rows]


def gene_aliases(annotation):
    aliases={}
    opener=gzip.open if str(annotation).endswith('.gz') else open
    with opener(annotation,'rt') as handle:
        for line in handle:
            if line.startswith('#'):continue
            parts=line.rstrip('\n').split('\t')
            if len(parts)<9 or parts[2] not in {'gene','transcript','mRNA'}:continue
            attrs=dict(re.findall(r'([A-Za-z0-9_]+)\s+"([^"]*)"',parts[8]))
            if not attrs:
                attrs={k:unquote(v) for entry in parts[8].split(';') if '=' in entry for k,v in [entry.split('=',1)]}
            gid=attrs.get('gene_id',attrs.get('gene',attrs.get('ID',''))).removeprefix('gene:')
            name=attrs.get('gene_name',attrs.get('Name','') if parts[2]=='gene' else '')
            if gid and name:
                aliases[gid]=name;aliases[gid.split('.')[0]]=name
    return aliases


def pair_sources(sources, args):
    labels = [args.target_label, args.reference_label]
    reference = sources['RNA']['obs']
    counts = reference[LABEL].value_counts()
    small = [f'{label}: {counts.get(label, 0)}' for label in labels if counts.get(label, 0) < args.min_cells_group]
    if small:
        raise ValueError(f'Need {args.min_cells_group} RNA cells per label; found {", ".join(small)}. Use --list-labels.')
    reference = reference.loc[reference[LABEL].isin(labels), LABEL]
    result = {}
    for modality, s in sources.items():
        s = select(s, s['obs'].index.isin(reference.index))
        s['obs'][LABEL] = reference.reindex(s['obs'].index).to_numpy()
        result[modality] = s
    return result


def embedding(s, modality, args):
    if modality == 'RNA':
        candidates = ['X_pca_uncorrected','X_pca','X_pca_cohort']
    else:
        candidates = ['X_spectral_uncorrected','X_spectral','X_pca_uncorrected','X_pca']
        if args.scope == 'merged' and modality != 'H3K9me3':
            candidates = ['X_spectral_harmony','X_pca_harmony','X_harmony']
    key = next((key for key in candidates if key in s['keys']), None)
    if key is None:
        raise ValueError(f'{s["path"]}: missing {modality} representation; tried {candidates}')
    import h5py
    with h5py.File(s['path'], 'r') as handle:
        node = handle['obsm'][key]
        if isinstance(node, h5py.Dataset):
            rep = node[s['rows'], :min(args.n_dimensions, node.shape[1])]
        else:
            rep = h5_rows(node, s['rows'])[:, :args.n_dimensions].toarray()
    if identity(s['path']) != s['stamp'] or rep.ndim != 2 or len(rep) != len(s['obs']) or not np.isfinite(rep).all():
        raise ValueError(f'{s["path"]}: invalid embedding or source changed')
    return np.asarray(rep, dtype=float), key


def context_for(sources, args):
    native = {}
    for modality, s in sources.items():
        obs = s['obs']
        rep, key = embedding(s, modality, args)
        log(f'{modality}: {len(obs):,} matched cells, using {key}')
        xy = fit_umap(rep, modality, args) if len(obs) >= 3 else np.full((len(obs), 2), np.nan)
        if xy.shape != (len(obs), 2) or (len(obs) >= 3 and not np.isfinite(xy).all()):
            raise ValueError(f'{modality}: UMAP returned invalid coordinates')
        native[modality] = dict(umap=xy, limits=limits(xy), labels=obs[LABEL].to_numpy())
    reference = sources['RNA']['obs'].index
    return dict(reference=native['RNA'], native=native,
                mapping={m:s['obs'].index.get_indexer(reference) for m,s in sources.items()},
                palette={args.target_label:'#E69F00', args.reference_label:'#0072B2'})


def gene_activity(s, annotation):
    import h5py
    import snapatac2 as snap
    obs = s['obs']
    if obs.empty:
        return sp.csr_matrix((0,0)), pd.Index([], dtype=str)
    native = [k for k in ('fragment_paired','fragment_single','insertion') if k in s['keys']]
    if native:
        dna = ad_module().AnnData(X=sp.csr_matrix((len(obs),0)), obs=pd.DataFrame(index=obs.index))
        with h5py.File(s['path'], 'r') as handle:
            if not s['native_ref']:
                raise ValueError(f'{s["path"]}: native fragments lack a reference sequence catalog')
            dna.uns['reference_sequences'] = h5_read_elem(handle['uns']['reference_sequences'])
            for key in native:
                dna.obsm[key] = h5_rows(handle['obsm'][key], s['rows'])
        log(f'{s["path"].name}: counting gene activities from native fragments')
    else:
        x, var = count_matrix(s)
        if not all(re.fullmatch(r'[^:]+:\d+-\d+', str(v)) for v in var.index):
            raise ValueError(f'{s["path"]}: no native fragments or genomic tile IDs')
        dna = ad_module().AnnData(X=x, obs=pd.DataFrame(index=obs.index), var=var)
        log(f'{s["path"].name}: using the original tile approximation for gene activities')
    ga = snap.pp.make_gene_matrix(dna, str(annotation), upstream=1000, downstream=1000,
                                  include_gene_body=True, id_type='gene', use_x=not native)
    if not pd.Index(ga.obs_names.astype(str)).equals(obs.index) or identity(s['path']) != s['stamp']:
        raise ValueError(f'{s["path"]}: changed source or gene activity cell order')
    column = next((k for k in ('gene_symbol','gene_name','gene_names') if k in ga.var), None)
    if column is None:
        aliases = gene_aliases(annotation)
        names = [aliases.get(str(v), aliases.get(str(v).split('.')[0], str(v))) for v in ga.var_names]
    else:
        names = ga.var[column].astype(str).tolist()
    genes = pd.Index(pd.unique(np.asarray(names, dtype=object)), name='gene')
    assignment = sp.csr_matrix((np.ones(len(names), dtype=np.int64),(np.arange(len(names)),genes.get_indexer(names))),
                               shape=(len(names),len(genes)))
    return (sp.csr_matrix(ga.X) @ assignment).tocsr(), genes


def signals_for(x, genes, features, obs, totals=None):
    if totals is None:
        totals = np.asarray(x.sum(axis=1)).ravel()
    use = np.flatnonzero(genes.isin(features))
    normalized, _ = normalize_rows(x[:, use], 1e4, totals)
    return dict(x=normalized.tocsc(), genes=genes[use], obs=obs)


OTHER = '#d9dfe3'


UNAVAILABLE = '#8e969e'


def save_pairwise_pdf(fig, path, args):
    save_pdf(fig, Path(path).with_suffix('.pdf'))


def comparison_title(args):
    return f'{title_for(args.scope)} | {args.target_label} vs {args.reference_label}'


SIGNAL_MAP = LinearSegmentedColormap.from_list('gene_signal',
    ['#dfe4e8', '#fde0dd', '#f5a3b5', '#df5279', '#a5124b', '#680026'])


def limits(coords):
    coords=np.asarray(coords)
    if not len(coords) or not np.isfinite(coords).all():return ((-1.,1.),(-1.,1.))
    low=np.min(coords,axis=0);high=np.max(coords,axis=0)
    pad=np.maximum(high-low,1.)*.07
    return ((float(low[0]-pad[0]),float(high[0]+pad[0])),
            (float(low[1]-pad[1]),float(high[1]+pad[1])))


def shared_values(values,context):
    result={}
    for modality,positions in context['mapping'].items():
        available=positions>=0
        output=np.full(len(positions),np.nan,dtype=float)
        native=np.asarray(values[modality],dtype=float)
        if len(native)!=len(context['native'][modality]['umap']):
            raise ValueError('Signal length does not match the exact modality cell order')
        output[available]=native[positions[available]]
        result[modality]=output
    return result


def point_size(n, args):
    return args.point_size if args.point_size else max(.7, min(5., 22000/max(n,1)))


def human(value):
    number,suffix=(value/1e6,'M') if value>=1e6 else ((value/1e3,'k') if value>=1e3 else (value,''))
    if number>=10: label=f'{number:,.0f}'
    elif number>=.01: label=f'{number:.2f}'.rstrip('0').rstrip('.')
    elif number==0: label='0'
    else: label=f'{number:.2g}'
    return label+suffix


def volcano_highlights(table, genes, args):
    """Every specifically highlighted point is named, including all top hits.

    --annotate controls only the extra significant context labels; it never
    limits the selected top markers or requested examples present in the test.
    """
    d=table.loc[np.isfinite(table.logFC)&np.isfinite(table.FDR)].copy()
    required=set(map(str,genes))
    absent=required-set(d.gene.astype(str))
    if absent:raise ValueError(f'Selected volcano markers have no finite test result: {sorted(absent)}')
    up=d.FDR.le(args.fdr)&d.logFC.ge(args.min_lfc)
    down=d.FDR.le(args.fdr)&d.logFC.le(-args.min_lfc)
    context=d.loc[up|down].sort_values(['FDR','logFC'],ascending=[True,False]).head(args.annotate)
    if up.any() and down.any():
        context=pd.concat([d.loc[up].nsmallest((args.annotate+1)//2,'FDR'),
                           d.loc[down].nsmallest(args.annotate//2,'FDR')])
    extra=set(map(str,args.extra_genes))
    chosen=required|set(context.gene.astype(str))|extra
    result=d.loc[d.gene.astype(str).isin(chosen)].copy()
    result['highlight_reason']=[
        '; '.join(reason for condition,reason in [(str(g) in required,'selected marker'),
                   (str(g) in extra,'requested example'),(str(g) in set(context.gene.astype(str)),'significant context')]
                  if condition) for g in result.gene]
    return result.sort_values(['logFC','FDR','gene'],ascending=[False,True,True])


def label_positions(desired,lower=.035,upper=.965):
    """Monotone callout positions with a guaranteed gap and no random jitter."""
    desired=np.asarray(desired,dtype=float)
    if not len(desired):return desired
    gap=min(.037,(upper-lower)/max(len(desired)-1,1))
    result=np.clip(desired,lower,upper)
    for i in range(1,len(result)):result[i]=max(result[i],result[i-1]+gap)
    if result[-1]>upper:
        result[-1]=upper
        for i in range(len(result)-2,-1,-1):result[i]=min(result[i],result[i+1]-gap)
    return result


def volcano(table, target, subtitle, path, args, highlight_genes=(),reference=None):
    d=table.loc[np.isfinite(table.logFC)&np.isfinite(table.FDR)].copy()
    selected=volcano_highlights(table,highlight_genes,args)
    n_side=max(int(selected.logFC.ge(0).sum()),int(selected.logFC.lt(0).sum()))
    # Fixed callout columns keep even 30+ selected genes individually legible.
    # The height grows with the label count; every label has a leader line.
    fig,ax=plt.subplots(figsize=(7.8,max(5.6,1.8+.17*n_side)))
    fig.subplots_adjust(left=.20,right=.79,bottom=.17,top=.85)
    comparison=f'{target} vs {reference}' if reference else f'{target} vs other B cells'
    ax.set(xlabel=f'log2 fold change: {comparison}',ylabel='−log10 adjusted p-value',
           title=f'{comparison}\n{subtitle}')
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(labelsize=8)
    if d.empty:
        ax.text(.5,.5,'No genes with finite adjusted p-values',transform=ax.transAxes,ha='center')
    else:
        positive=d.FDR[d.FDR>0]
        floor=max(np.finfo(float).tiny,min(float(positive.min())*.1,args.fdr*.01)) if len(positive) else 1e-300
        d['y']=-np.log10(d.FDR.clip(lower=floor))
        up=d.FDR.le(args.fdr)&d.logFC.ge(args.min_lfc)
        down=d.FDR.le(args.fdr)&d.logFC.le(-args.min_lfc)
        colors=np.where(up,'#be3655',np.where(down,'#2879b0','#c3c9cf'))
        ax.scatter(d.logFC,d.y,c=colors,s=7,alpha=.65,lw=0,rasterized=False)
        ax.axhline(-np.log10(args.fdr),c='#777777',ls='--',lw=.7)
        for v in (-args.min_lfc,args.min_lfc): ax.axvline(v,c='#777777',ls='--',lw=.7)
        selected=selected.merge(d[['gene','y']],on='gene',how='left',validate='one_to_one')
        selected['highlight_class']=np.where(selected.FDR.le(args.fdr)&selected.logFC.ge(args.min_lfc),'Higher',
            np.where(selected.FDR.le(args.fdr)&selected.logFC.le(-args.min_lfc),'Lower','Other tested'))
        highlight_colors={'Higher':'#9e1638','Lower':'#155781','Other tested':'#59636e'}
        ax.scatter(selected.logFC,selected.y,s=23,c=selected.highlight_class.map(highlight_colors),
                   edgecolors='white',linewidths=.45,zorder=4,rasterized=False)
        ax.set_ylim(-max(float(d.y.max())*.03,.2),max(float(d.y.max())*1.08,1.))
        ymin,ymax=ax.get_ylim()
        left_labels=[]
        for right in (False,True):
            side=selected.loc[selected.logFC.ge(0).eq(right)].sort_values(['y','gene'])
            locations=label_positions((side.y.to_numpy()-ymin)/(ymax-ymin))
            for r,y in zip(side.itertuples(),locations):
                label=ax.annotate(str(r.gene),xy=(r.logFC,r.y),xycoords='data',
                    xytext=(1.035 if right else -.16,y),textcoords='axes fraction',
                    ha='left' if right else 'right',va='center',fontsize=7,fontstyle='italic',
                    color=highlight_colors[r.highlight_class],annotation_clip=False,
                    arrowprops=dict(arrowstyle='-',lw=.45,color='#8e969e',alpha=.85),zorder=5)
                if not right:left_labels.append(label)
        if left_labels:
            # Keep the axis title outside the callout column, and leave the
            # tick-label gutter free. Long symbols expand the saved page.
            fig.canvas.draw()
            renderer=fig.canvas.get_renderer()
            leftmost=min(t.get_window_extent(renderer).x0 for t in left_labels)
            xlabel=ax.transAxes.inverted().transform((leftmost-18*fig.dpi/72,0))[0]
            ax.yaxis.set_label_coords(xlabel,.5)
        handles=[Line2D([],[],marker='o',ls='',color=c,label=f'{name} ({n:,})',ms=5)
                 for c,name,n in [('#be3655','Higher',int(up.sum())),('#2879b0','Lower',int(down.sum())),
                                  ('#919aa3','Other tested',int((~(up|down)).sum()))]]
        ax.legend(handles=handles,frameon=False,fontsize=8,loc='upper center',
                  bbox_to_anchor=(.5,-.14),ncol=3)
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    save_pairwise_pdf(fig,path,args)


def clean(ax,geometry):
    ax.set_xlim(*geometry['limits'][0]);ax.set_ylim(*geometry['limits'][1])
    ax.set_aspect('equal',adjustable='box')
    ax.set_xticks([]);ax.set_yticks([])
    ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
    for spine in ax.spines.values():spine.set_visible(False)


def cell_types(ax,geometry,palette,args):
    labels=geometry['labels'];coords=geometry['umap']
    if len(coords)<3 or not np.isfinite(coords).all():
        clean(ax,geometry)
        ax.text(.5,.5,'Too few matched cells\nfor a native UMAP',ha='center',va='center',transform=ax.transAxes,fontsize=8)
        return
    order=np.random.default_rng(args.seed).permutation(len(coords))
    ax.scatter(coords[order,0],coords[order,1],c=[palette[x] for x in labels[order]],
               s=point_size(len(coords),args),linewidths=0,alpha=.85,rasterized=False)
    clean(ax,geometry)


def type_handles(labels,palette):
    return [Line2D([],[],marker='o',ls='',color=palette[label],label=label,ms=5)
            for label in sorted(set(labels))]


def unavailable_handle():
    return Line2D([],[],marker='o',ls='',markerfacecolor='none',markeredgecolor=UNAVAILABLE,
                  label='Signal unavailable',ms=4)


def prepare_signal(values,cap,transform):
    measured=np.isfinite(values)
    order=np.flatnonzero(measured)
    order=order[np.argsort(values[order],kind='stable')]
    colors=np.log1p(values[order]) if transform=='log1p' else values[order]
    vmax=(np.log1p(cap) if transform=='log1p' else cap) if cap>0 else 1.
    return measured,order,colors,vmax


def signal(ax,geometry,values,cap,args,transform,prepared=None):
    values=np.asarray(values,dtype=float);coords=geometry['umap']
    if len(values)!=len(coords) or np.isinf(values).any() or (values[np.isfinite(values)]<0).any():
        raise ValueError('Invalid signal values/cell order')
    if len(coords)<3 or not np.isfinite(coords).all():
        clean(ax,geometry)
        ax.text(.5,.5,'Too few matched cells\nfor a native UMAP',ha='center',va='center',transform=ax.transAxes,fontsize=8)
        return None,0
    measured,order,colors,vmax=prepare_signal(values,cap,transform) if prepared is None else prepared
    if (~measured).any():
        ax.scatter(coords[~measured,0],coords[~measured,1],s=point_size(len(coords),args),
            facecolors='none',edgecolors=UNAVAILABLE,linewidths=.3,alpha=.5,rasterized=False)
    artist=ax.scatter(coords[order,0],coords[order,1],c=colors,vmin=0,vmax=vmax,cmap=SIGNAL_MAP,
        s=point_size(len(coords),args),linewidths=0,rasterized=False)
    if not measured.any():
        ax.text(.5,.98,'Signal unavailable',ha='center',va='top',transform=ax.transAxes,fontsize=8)
    clean(ax,geometry)
    return artist,int(measured.sum())


def colorbar(fig,cax,artist,cap,unit,transform,available=True):
    if not available:cax.set_axis_off();return
    cb=fig.colorbar(artist,cax=cax,orientation='horizontal')
    if cb.solids is not None:
        cb.solids.set_rasterized(False);cb.solids.set_edgecolor('face')
    raw=[0.] if cap<=0 else [0.,float(np.expm1(np.log1p(cap)/2)) if transform=='log1p' else cap/2,cap]
    cb.set_ticks(np.log1p(raw) if transform=='log1p' else raw)
    cb.set_ticklabels(['0'] if cap<=0 else [human(v) for v in raw[:-1]]+['≥'+human(cap)])
    cb.ax.tick_params(labelsize=6,length=2)
    cb.set_label(unit+(' (log1p scale)' if transform=='log1p' else ''),fontsize=7,labelpad=1)


def grid(ncols):
    fig = plt.figure(figsize=(3.35*ncols, 4.6))
    gs = fig.add_gridspec(2,ncols,height_ratios=[1,.045],hspace=.18,wspace=.24)
    fig.subplots_adjust(left=.025,right=.985,top=.79,bottom=.23)
    return fig,gs


def pair_types(ax,geometry,labels,palette,args):
    coords=geometry['umap'];labels=np.asarray(labels,dtype=str)
    if len(coords)<3 or not np.isfinite(coords).all():
        clean(ax,dict(geometry,limits=((-1.,1.),(-1.,1.))))
        ax.text(.5,.5,'Too few matched cells\nfor a native UMAP',ha='center',va='center',transform=ax.transAxes,fontsize=8)
        return
    chosen=np.isin(labels,[args.target_label,args.reference_label])
    order=np.random.default_rng(args.seed).permutation(np.flatnonzero(chosen))
    ax.scatter(coords[order,0],coords[order,1],s=point_size(len(coords),args),
               c=[palette[x] for x in labels[order]],linewidths=0,alpha=.85,rasterized=False)
    bounds=np.asarray(limits(coords[chosen]) if chosen.any() else ((-1.,1.),(-1.,1.)))
    center=bounds.mean(axis=1);radius=float(np.max(bounds[:,1]-bounds[:,0]))/2
    square=tuple((float(c-radius),float(c+radius)) for c in center)
    focused=dict(geometry,limits=square)
    clean(ax,focused)
    if not chosen.any():ax.text(.5,.5,'No matched cells for this pair',ha='center',va='center',transform=ax.transAxes,fontsize=8)


def pair_handles(palette,args):
    return [Line2D([],[],marker='o',ls='',color=palette[x],label=x,ms=5)
            for x in (args.target_label,args.reference_label)]


def cell_overview(context,folder,args):
    mods=list(context['native'])
    fig,axes=plt.subplots(1,len(mods),figsize=(3.5*len(mods),4.1),squeeze=False)
    for ax,m in zip(axes[0],mods):
        pair_types(ax,context['native'][m],context['native'][m]['labels'],context['palette'],args)
        ax.set_title(f'{m} | UMAP')
    fig.legend(handles=pair_handles(context['palette'],args),loc='lower center',ncol=3,frameon=False,fontsize=8)
    fig.suptitle(comparison_title(args)+' | RNA-defined groups',fontsize=12)
    fig.tight_layout(rect=(0,.10,1,.93));save_pairwise_pdf(fig,folder/'cell_umaps',args)


def gene_views(gene,context,values,caps,folder,title,args):
    mapped=shared_values(values,context)
    reference=context['reference'];palette=context['palette'];modalities=list(context['native'])
    fig,gs=grid(1+len(modalities))
    ax=fig.add_subplot(gs[0,0]);cell_types(ax,reference,palette,args);ax.set_title('Cell types\nRNA reference')
    fig.add_subplot(gs[1,0]).set_axis_off()
    for j,m in enumerate(modalities,1):
        transform='log1p' if m=='RNA' else args.dna_color_transform
        prepared=prepare_signal(mapped[m],caps[m][gene],transform)
        ax=fig.add_subplot(gs[0,j]);artist,n=signal(ax,reference,mapped[m],caps[m][gene],args,transform,prepared)
        ax.set_title(f'{m}\n{n:,} / {len(reference["umap"]):,} cells available')
        colorbar(fig,fig.add_subplot(gs[1,j]),artist,caps[m][gene],
                 'RNA CP10K' if m=='RNA' else 'Gene-body CP10K',transform,n>0)
    handles=type_handles(reference['labels'],palette)+[unavailable_handle()]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.02),ncol=min(7,len(handles)),frameon=False,fontsize=8)
    fig.text(.5,.115,'RNA-reference UMAP; hollow points indicate unavailable signal. Measured zero is filled.',ha='center',fontsize=8,color='#56616b')
    fig.suptitle(f'{gene} | {title} | UMAP',fontsize=13,fontstyle='italic')
    save_pairwise_pdf(fig,folder/'shared_umap',args)


SUMMARY_COLUMNS = ('gene','modality','rna_label','reference_cells','measured_cells',
    'unavailable_cells','fraction_detected','mean_cp10k','median_cp10k',
    'mean_log1p_cp10k','labels_from','unit','values_capped')


def group_signal_summary(gene, mapped, reference_obs, args):
    """Legacy descriptive statistics: exclude unavailable values, retain true zeros."""
    rows = []
    for modality, values in mapped.items():
        values = np.asarray(values, dtype=float)
        if values.shape != (len(reference_obs),):
            raise ValueError(f'{modality}: signal/reference cell order differs')
        for label in (args.target_label,args.reference_label):
            group = reference_obs[LABEL].eq(label).to_numpy()
            measured = values[group]
            measured = measured[np.isfinite(measured)]
            rows.append(dict(gene=gene,modality=modality,rna_label=label,
                reference_cells=int(group.sum()),measured_cells=len(measured),
                unavailable_cells=int(group.sum())-len(measured),
                fraction_detected=float(np.mean(measured>0)) if len(measured) else np.nan,
                mean_cp10k=float(np.mean(measured)) if len(measured) else np.nan,
                median_cp10k=float(np.median(measured)) if len(measured) else np.nan,
                mean_log1p_cp10k=float(np.mean(np.log1p(measured))) if len(measured) else np.nan,
                labels_from='RNA lymphoid_final; exact physical cell matches',
                unit='CP10K',values_capped=False))
    return pd.DataFrame(rows,columns=SUMMARY_COLUMNS)


def pair_arguments(scope, first, second, cli):
    return SimpleNamespace(scope=scope,target_label=first,reference_label=second,
        min_cells_group=50,min_cells_condition=25,random_replicates=3,min_donors=3,seed=0,
        n_dimensions=30,neighbors=15,min_dist=.3,dna_metric='cosine',dna_color_transform='log1p',
        point_size=0,annotate=0,fdr=.001,min_lfc=.5,min_detection=.02,top_genes=30,
        extra_genes=cli.extra_genes)


def pair_methods(scope, cli):
    return list(dict.fromkeys(cli.merged_methods)) if scope=='merged' else ['random_subsets']


def pair_signals(sources, features, args, annotation, rna=None):
    if rna is None:
        raw, var = count_matrix(sources['RNA'])
        rna = rna_symbols(raw,var,args.scope)
    x, genes, totals = rna
    signals = {'RNA':signals_for(x,genes,features,sources['RNA']['obs'],totals)}
    for modality, s in sources.items():
        if modality=='RNA':
            continue
        matrix, symbols = gene_activity(s,annotation)
        # signals_for computes the denominator over ALL gene activities first.
        signals[modality] = signals_for(matrix,symbols,features,s['obs'])
    return signals


def write_pair_summary(features, signals, reference_obs, args, out):
    positions = {m:s['obs'].index.get_indexer(reference_obs.index) for m,s in signals.items()}
    tables = []
    for gene in sorted(features):
        mapped = {}
        for modality, s in signals.items():
            native = get_gene(s,gene,s['obs'].index)
            rows = positions[modality]
            values = np.full(len(reference_obs),np.nan)
            available = rows>=0
            values[available] = native[rows[available]]
            mapped[modality] = values
        tables.append(group_signal_summary(gene,mapped,reference_obs,args))
    table = pd.concat(tables,ignore_index=True) if tables else pd.DataFrame(columns=SUMMARY_COLUMNS)
    write_table(table,out/'gene_group_summaries.tsv')
    return table


def saved_pair_features(out, args, methods):
    features = set(args.extra_genes)
    numeric = ('logFC','FDR','pct_detected_target','pct_detected_rest')
    for method in methods:
        path = out/('donor_grouped' if method=='donor_replicates' else 'pseudo_grouped')/'DGE_results.tsv'
        require_files(path)
        table = pd.read_csv(path,sep='\t',dtype={'gene':str},keep_default_na=False,
                            na_values={key:['NA',''] for key in numeric})
        if not table.gene.is_unique:
            raise ValueError(f'{path}: duplicate genes')
        for key in numeric:
            table[key] = pd.to_numeric(table[key],errors='raise')
        for key, expected in (('target_label',args.target_label),('reference_label',args.reference_label)):
            if key in table and not table[key].eq(expected).all():
                raise ValueError(f'{path}: unexpected {key}')
        features.update(select_pair_markers(table,args).gene)
    expected = {slug(gene) for gene in features}
    plotted = {p.parent.name for p in (out/'genes').glob('*/shared_umap.pdf')}
    if expected!=plotted:
        raise ValueError(f'{out}: saved marker/extra-gene union differs from Shared UMAPs; '
                         f'unplotted={sorted(expected-plotted)}, unexplained_plots={sorted(plotted-expected)}. '
                         'Check --merged-methods and supply the original --extra-genes before writing summaries.')
    return features


def compare(config, scope, first, second, cli, annotation):
    results = config['output']
    final, publication = results/'final', results/'publication'
    args = pair_arguments(scope,first,second,cli)
    modalities = (*MARKS,'H3K9me3') if scope=='merged' else MARKS
    sources = pair_sources({m:source(final,publication,scope,m) for m in ('RNA',*modalities)},args)
    methods = pair_methods(scope,cli)
    obs = sources['RNA']['obs']
    if 'donor_replicates' in methods:
        coverage = obs.groupby(['sample',LABEL]).size().unstack(fill_value=0)
        keep = coverage.min(axis=1).ge(args.min_cells_condition)
        if keep.sum() < args.min_donors or (coverage.loc[keep].sum()<args.min_cells_group).any():
            raise ValueError(f'{scope} {first} / {second}: needs three donors with ≥25 cells per label:\n{coverage}')
    key = comparison_key(first,second)
    out = results/'pairwise'/scope/key
    with output_lock(results/'.work'/'pairwise'/scope/key):
        context = context_for(sources,args)
        cell_overview(context,out,args)
        raw, var = count_matrix(sources['RNA'])
        x, genes, totals = rna_symbols(raw,var,scope)
        del raw
        features = set(args.extra_genes)
        for method in methods:
            paired = method=='donor_replicates'
            if paired:
                counts, meta = donor_pairs(x,obs,genes,first,args)
                used = obs.biological_donor.isin(meta.donor)
                subtitle = f'{meta.donor.nunique()} paired biological donors'
            else:
                counts, meta = random_splits(x,obs,genes,first,args)
                used = obs[LABEL].notna()
                subtitle = 'Random cell subsets (exploratory)'
            table = fit_model(counts,meta,'pair',paired,args.fdr)
            table = annotate_result(table,x[used.to_numpy()],obs.loc[used],genes,method,args)
            markers = select_pair_markers(table,args)
            folder = out/('donor_grouped' if paired else 'pseudo_grouped')
            write_table(table,folder/'DGE_results.tsv')
            if not paired:
                write_table(markers,folder/'selected_markers.tsv')
            volcano(table,first,f'{title_for(scope)} | {subtitle}',folder/'volcano',args,
                    highlight_genes=list(markers.gene),reference=second)
            features.update(markers.gene)

        if features:
            signals = pair_signals(sources,features,args,annotation,(x,genes,totals))
            write_pair_summary(features,signals,obs,args,out)
            for gene in sorted(features):
                values = {m:get_gene(s,gene,s['obs'].index) for m,s in signals.items()}
                caps = {m:{gene:color_cap(v,99.)} for m,v in values.items()}
                gene_views(gene,context,values,caps,out/'genes'/slug(gene),comparison_title(args),args)
        else:
            write_pair_summary(features,{},obs,args,out)
    log(f'Pairwise: {out}')


def pair_requests(config, cli):
    both = cli.label_a is not None and cli.label_b is not None
    if ((cli.label_a is not None or cli.label_b is not None) and not both) or (both and cli.label_a==cli.label_b):
        raise ValueError('Supply two different labels or omit both to run the study pairs')
    if cli.dataset=='Tonsil1_gt600' and not both:
        raise ValueError('Choose explicit labels for Tonsil1_gt600')
    scopes = [cli.dataset] if cli.dataset else config['pairwise']['datasets']
    pairs = [(cli.label_a,cli.label_b)] if both else config['pairwise']['pairs']
    return [(scope,first,second) for scope in scopes for first,second in pairs]


def pairwise_stage(config, cli):
    requests = pair_requests(config,cli)

    thread_budget(config['threads'])
    check_r()
    plotting_style()
    annotation = gene_annotation()
    for scope,first,second in requests:
        compare(config,scope,first,second,cli,annotation)


def pairwise_summary_stage(config, cli):
    """Read finalized profiles and saved fits; never fit models or build embeddings."""
    plans = []
    for scope,first,second in pair_requests(config,cli):
        args = pair_arguments(scope,first,second,cli)
        out = config['output']/'pairwise'/scope/comparison_key(first,second)
        plans.append((args,out,saved_pair_features(out,args,pair_methods(scope,cli))))
    # Validate every requested gene set before starting expensive signal counting.
    thread_budget(config['threads'])
    annotation = gene_annotation()
    final, publication = config['output']/'final', config['output']/'publication'
    for args,out,features in plans:
        modalities = ('RNA',*MARKS,'H3K9me3') if args.scope=='merged' else ('RNA',*MARKS)
        with output_lock(config['output']/'.work'/'pairwise'/args.scope/out.name):
            # Check again under the comparison lock in case a producer changed fits.
            if features!=saved_pair_features(out,args,pair_methods(args.scope,cli)):
                raise ValueError(f'{out}: marker set changed during summary preparation')
            sources = pair_sources({m:source(final,publication,args.scope,m) for m in modalities},args)
            signals = pair_signals(sources,features,args,annotation) if features else {}
            table = write_pair_summary(features,signals,sources['RNA']['obs'],args,out)
        log(f'Pairwise summary: {out / "gene_group_summaries.tsv"} ({len(features)} genes, {len(table)} rows)')


# Stage 04: germinal center zone contrasts and effect sizes

DZ, LZ = 'B_GC_DZ', 'B_GC_LZ'


PSEUDO_SCOPE = 'merged_pseudo_subsets'


POOLED_GROUP = 'Merged'


def pseudobulk(x, obs, genes, scope, args):
    units = []; codes = np.full(len(obs),-1,dtype=int)
    if scope == 'merged':
        for donor in sorted(obs['sample'].unique()):
            for zone in (LZ,DZ):
                idx = np.flatnonzero(obs['sample'].eq(donor).to_numpy() & obs[LABEL].eq(zone).to_numpy())
                if len(idx)<args.min_cells: raise ValueError(f'{scope}: {donor}/{zone} has fewer than {args.min_cells} cells')
                codes[idx] = len(units)
                units.append(dict(unit=f'{donor}__{zone}',donor=donor,zone='DZ' if zone==DZ else 'LZ',
                    n_cells=len(idx),unit_type='paired biological donor'))
        if obs['sample'].nunique()<args.min_donors: raise ValueError('Too few complete biological donor pairs')
    else:
        for zone in (LZ,DZ):
            idx = np.flatnonzero(obs[LABEL].eq(zone).to_numpy())
            if len(idx)<max(args.min_cells,args.random_replicates): raise ValueError('Too few single-donor cells per zone')
            # Identical physical cells receive the same split across modalities
            # whenever possible, by stable ID order; counts never affect splits.
            rng = np.random.default_rng(int(digest([args.seed,scope,zone])[:8],16))
            idx = idx[np.argsort(obs.index[idx].astype(str))]
            for rep,subset in enumerate(np.array_split(rng.permutation(idx),args.random_replicates),1):
                codes[subset] = len(units)
                units.append(dict(unit=f'{zone}__split{rep}',donor='Tonsil1',zone='DZ' if zone==DZ else 'LZ',
                    n_cells=len(subset),unit_type='pseudo subsets'))
    if (codes<0).any(): raise ValueError('Unexpected cells outside DZ/LZ')
    meta = pd.DataFrame(units).set_index('unit')
    assignment = sp.csr_matrix((np.ones(len(obs),dtype=np.int64),(codes,np.arange(len(obs)))),shape=(len(meta),len(obs)))
    counts = pd.DataFrame((assignment@x).toarray().T,index=genes,columns=meta.index)
    return counts,meta


def eligible_donors(observations, args):
    rows = []; eligible = []
    donors = sorted(set.union(*(set(x['sample'].astype(str)) for x in observations.values())))
    for donor in donors:
        complete = True
        for modality,obs in observations.items():
            for zone in (DZ,LZ):
                n = int((obs['sample'].eq(donor) & obs[LABEL].eq(zone)).sum())
                rows.append(dict(donor=donor,modality=modality,zone=zone,n_cells=n,meets_minimum=n>=args.min_cells))
                complete &= n>=args.min_cells
        if complete: eligible.append(donor)
    frame = pd.DataFrame(rows); frame['included_in_all_assays'] = frame.donor.isin(eligible)
    return eligible,frame


def pooled_pseudobulk(x,obs,genes,args):
    if not obs.index.is_unique:raise ValueError('Merged pseudo subsets require unique physical cell IDs')
    codes=np.full(len(obs),-1,dtype=int);units=[]
    for zone in (LZ,DZ):
        rows=np.flatnonzero(obs[LABEL].eq(zone).to_numpy())
        if len(rows)<max(args.min_cells,args.random_replicates):
            raise ValueError(f'{PSEUDO_SCOPE}: too few pooled {zone} cells for nonempty pseudo subsets')
        rows=rows[np.argsort(obs.index[rows].astype(str))]
        rng=np.random.default_rng(int(digest([args.seed,PSEUDO_SCOPE,zone])[:8],16))
        for rep,subset in enumerate(np.array_split(rng.permutation(rows),args.random_replicates),1):
            codes[subset]=len(units)
            units.append(dict(unit=f'{zone}__split{rep}',donor=POOLED_GROUP,
                zone='DZ' if zone==DZ else 'LZ',n_cells=len(subset),unit_type='pseudo subsets'))
    if (codes<0).any():raise ValueError('Merged pseudo subsets contain cells outside DZ/LZ')
    meta=pd.DataFrame(units).set_index('unit')
    assignment=sp.csr_matrix((np.ones(len(obs),dtype=np.int64),(codes,np.arange(len(obs)))),
        shape=(len(meta),len(obs)))
    counts=pd.DataFrame((assignment@sp.csr_matrix(x)).toarray().T,index=genes,columns=meta.index)
    return counts,meta


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def union_intervals(table):
    """Merge true overlaps per gene. Adjacent distinct intervals remain separate."""
    rows = []
    for (gene, chrom), g in table.groupby(['gene', 'chrom'], sort=True):
        current = None
        for start, end in sorted(set(zip(g.start, g.end))):
            if current is not None and start < current[1]: current[1] = max(current[1], int(end))
            else:
                if current is not None: rows.append((chrom, *current, gene))
                current = [int(start), int(end)]
        if current is not None: rows.append((chrom, *current, gene))
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'gene'])


def attributes(text):
    # GFF3 key=value fields must bypass the GTF regex: long unquoted
    # attributes otherwise trigger quadratic backtracking on every exon.
    if re.match(r'^\s*[^=\s;]+=', text):
        return {key:unquote(value) for item in text.split(';') if '=' in item
            for key,value in [item.strip().split('=',1)]}
    return dict(re.findall(r'([^\s;]+)\s+"([^"]*)"',text))


def read_gene_bodies(path,flank=1000,chrom_sizes=None):
    """Prefer gene features per locus; transcripts/exons supply missing spans.

    Gene features include introns. Fallback transcript/exon bounds supply one
    continuous gene body, not an exon-only sum. GTF/GFF coordinates are 1-based
    closed; returned BED coordinates are 0-based half-open.
    """
    if flank<0:raise ValueError('Gene-body flank must be nonnegative')
    genes={};fallback={};names={};tx_gene={};exons=[]
    skipped_contigs=set();excluded=0
    with open_text(path) as handle:
        for number,line in enumerate(handle,1):
            if not line.strip() or line.startswith('#'):continue
            parts=line.rstrip('\n').split('\t')
            if len(parts)!=9:raise ValueError(f'{path}:{number}: expected GTF/GFF with nine columns')
            chrom,_,kind,start,end,_,strand,_,raw=parts
            if kind not in ('gene','transcript','mRNA','exon'):continue
            try:start,end=int(start)-1,int(end)
            except ValueError:raise ValueError(f'{path}:{number}: invalid gene coordinates')
            if start<0 or end<=start:raise ValueError(f'{path}:{number}: invalid gene interval')
            if chrom_sizes is not None:
                if chrom not in chrom_sizes:skipped_contigs.add(chrom);continue
                if end>chrom_sizes[chrom]:raise ValueError(f'{path}:{number}: gene exceeds chromosome length; check genome build')
            attrs=attributes(raw)
            gid=attrs.get('gene_id') or (attrs.get('ID') if kind=='gene' else attrs.get('Parent'))
            name=attrs.get('gene_name') or (attrs.get('Name') if kind=='gene' else None)
            if gid:gid=gid.removeprefix('gene:')
            if not gid and name:gid=name
            if not gid:continue
            if name:names[gid]=name
            if kind in ('transcript','mRNA') and attrs.get('ID'):
                tx_gene[attrs['ID']]=gid
                tx_gene[attrs['ID'].removeprefix('transcript:')]=gid
            if kind=='exon' and not attrs.get('gene_id'):
                exons.append((gid,chrom,strand,start,end));continue
            target=genes if kind=='gene' else fallback
            key=(gid,chrom,strand)
            if key in target:
                previous=target[key];target[key]=(min(previous[0],start),max(previous[1],end))
            else:target[key]=(start,end)
    for parent,chrom,strand,start,end in exons:
        gid=tx_gene.get(parent,tx_gene.get(parent.removeprefix('transcript:'),parent))
        key=(gid,chrom,strand)
        previous=fallback.get(key,(start,end));fallback[key]=(min(previous[0],start),max(previous[1],end))
    loci=dict(fallback);loci.update(genes)
    rows=[]
    for (gid,chrom,strand),(start,end) in sorted(loci.items()):
        name=names.get(gid,gid);upper=name.upper()
        if (not name.strip() or upper in ('NA','NAN','NONE','.') or 'ENSG' in upper
                or 'MALAT1' in upper or upper.startswith(('MT-','MT_','MT.'))):
            excluded+=1;continue
        left=max(0,start-flank);right=end+flank
        if chrom_sizes is not None:right=min(right,int(chrom_sizes[chrom]))
        rows.append((chrom,left,right,name,'gene_body',gid,strand,start,end,None))
    if not rows:raise ValueError(f'{path}: no named gene-body intervals on the selected reference chromosomes')
    table=pd.DataFrame(rows,columns=['chrom','start','end','gene','class','input_name',
        'strand','original_start','original_end','tss'])
    table.attrs['annotation_summary']=dict(source=str(Path(path).resolve()),flank_bp=flank,
        gene_features=len(genes),fallback_loci=len(set(fallback)-set(genes)),
        excluded_symbol_loci=excluded,ignored_contigs=sorted(skipped_contigs),
        coordinate_rule='Full gene start-1 through gene end, then +/- flank; introns retained')
    return table


def ensembl_id(value):
    value = str(value).removeprefix('gene:')
    match = re.fullmatch(r'(ENSG\d+)(?:\.\d+)?(_PAR_Y)?', value)
    return match.group(1)+(match.group(2) or '') if match else None


def annotation_symbols(path, targets):
    """Read only names for requested stable IDs; accept GTF and GFF3."""
    result = {target:set() for target in targets}
    with open_text(path) as handle:
        for number,line in enumerate(handle,1):
            if not line.strip() or line.startswith('#'): continue
            parts = line.rstrip('\n').split('\t')
            if len(parts)!=9: raise ValueError(f'{path}:{number}: expected GTF/GFF with nine columns')
            kind = parts[2]
            if kind not in ('gene','transcript','mRNA','exon'): continue
            attrs = attributes(parts[8])
            gene = attrs.get('gene_id') or (attrs.get('ID') if kind=='gene' else attrs.get('Parent'))
            stable = ensembl_id(gene) if gene else None
            if stable not in result: continue
            name = attrs.get('gene_name') or (attrs.get('Name') if kind=='gene' else None)
            if name and name.strip().lower() not in ('','na','nan','none','.') and ensembl_id(name) is None:
                result[stable].add(name.strip())
    return result


def resolve_targets(table, annotation):
    """Map stable Ensembl IDs to unique symbols; retain unresolved IDs."""
    ids = table.gene.map(ensembl_id)
    targets = set(ids.dropna())
    symbols = annotation_symbols(annotation,targets) if targets else {}
    names = []
    for raw,gene in zip(table.gene,ids):
        choices = symbols.get(gene,set())
        names.append(next(iter(choices)) if len(choices)==1 else (gene or raw))
    return table.assign(gene=names)



def enhancer_bed(path, chrom_sizes):
    rows = []
    with open_text(path) as handle:
        for n, line in enumerate(handle, 1):
            line = line.lstrip('\ufeff').strip()
            if not line or line.startswith(('#','track ','browser ')):
                continue
            fields = line.split()
            if not rows and len(fields) >= 4 and fields[0].lower() in ('chrom','chr','chromosome') and [x.lower() for x in fields[1:3]] == ['start','end']:
                continue
            if len(fields) < 4:
                raise ValueError(f'{path}:{n}: expected BED4 with the linked gene in the fourth column')
            chrom, start, end, gene = fields[:4]
            start, end = int(start), int(end)
            if chrom not in chrom_sizes or not (0 <= start < end <= chrom_sizes[chrom]):
                raise ValueError(f'{path}:{n}: invalid coordinates for the supplied reference')
            if gene.lower() in ('','nan','na','none','.'):
                raise ValueError(f'{path}:{n}: no linked gene')
            rows.append((chrom,start,end,gene,gene))
    if not rows:
        raise ValueError(f'{path}: no linked enhancers')
    return pd.DataFrame(rows, columns=['chrom','start','end','gene','input_name'])


class RegionPlan:
    def __init__(self, annotation, enhancers, chrom_sizes, flank=1000):
        body = read_gene_bodies(annotation,flank,chrom_sizes)
        linked = resolve_targets(enhancers,annotation)
        tables = {
            'gene_body':union_intervals(body),
            'combined':union_intervals(pd.concat([body,linked], ignore_index=True)),
        }
        coordinates = pd.concat(tables.values(), ignore_index=True).loc[:, ['chrom','start','end']]
        coordinates = coordinates.drop_duplicates().sort_values(['chrom','start','end'])
        self.ids = region_ids(coordinates)
        self.arms = {}
        for key, table in tables.items():
            genes = pd.Index(sorted(table.gene.unique()), name='gene')
            projection = sp.csr_matrix((np.ones(len(table), dtype=np.int64),
                (self.ids.get_indexer(region_ids(table)), genes.get_indexer(table.gene))),
                shape=(len(self.ids),len(genes)))
            self.arms[key] = (genes, projection)


def region_ids(table):
    return pd.Index([f'{r.chrom}:{r.start}-{r.end}' for r in table.itertuples(index=False)], name='region')


NATIVE = ('fragment_paired','fragment_single')


def barcode_candidates(obs, sample):
    """Literal final-object identifiers only; no barcode correction or fuzzy matching."""
    result = {}
    for cell, row in obs.iterrows():
        cell = str(cell)
        if not cell.startswith(sample+'::'): raise ValueError(f'{sample}: invalid physical cell ID {cell}')
        values = {cell, cell.split('::',1)[1]}
        for field in ('source_cell_id','barcode','_source_index','merged_cell_id'):
            value = row.get(field)
            if isinstance(value,str) and value.strip():
                values.add(value)
                if value.startswith(sample+'::'): values.add(value.split('::',1)[1])
        result[cell] = sorted(values)
    return result


def subset_valid(item):
    folder = Path(item['cache_dir'])
    try:
        meta = json.loads((folder/'SUBSET_COMPLETED.json').read_text())
        return (meta['signature']==item['subset_signature'] and
            meta['subset_file']==identity(folder/'selected.fragments.tsv.gz') and
            meta['mapping_file']==identity(folder/'selected_cells.tsv'))
    except (OSError,ValueError,KeyError): return False


def prepare_subset(item, args):
    folder = Path(item['cache_dir'])
    if identity(item['path'])!=item['source']: raise ValueError(f'{item["path"]}: fragment source changed after preflight; rerun with stable inputs')
    if subset_valid(item):
        log(f'{item["sample"]} {item["mark"]}: reusing verified selected-cell fragment TSV')
        return item
    with output_lock(folder):
        if subset_valid(item): return item
        candidates = item['barcode_candidates']; lookup = {}
        for cell, names in candidates.items():
            for name in names: lookup.setdefault(name,set()).add(cell)
        mapping = {}; records = dict.fromkeys(candidates,0)
        target = folder/'selected.fragments.tsv.gz'
        temporary = folder/'selected.fragments.tmp.tsv.gz'
        marker = folder/'SUBSET_COMPLETED.json'; marker.unlink(missing_ok=True)
        source = Path(item['path']); before = identity(source)
        log(f'{item["sample"]} {item["mark"]}: selecting {len(candidates):,} final cells from {source}')
        started = time.monotonic(); last = started; total = 0; retained = 0
        try:
            with gzip.open(source,'rt') as src, gzip.open(temporary,'wt',compresslevel=1) as dest:
                for total,line in enumerate(src,1):
                    if not line.strip() or line.startswith('#'): continue
                    fields = line.rstrip('\r\n').split('\t',4)
                    if len(fields)<4: raise ValueError(f'{source}:{total}: expected tab-delimited fragment BED')
                    barcode = fields[3]; cells = lookup.get(barcode)
                    if cells:
                        if len(cells)!=1: raise ValueError(f'{source}: ambiguous literal barcode {barcode!r} maps to multiple final cells')
                        cell = next(iter(cells)); previous = mapping.setdefault(cell,barcode)
                        if previous!=barcode:
                            raise ValueError(f'{source}: final cell {cell} matches two different source barcodes ({previous}, {barcode})')
                        dest.write(line if line.endswith('\n') else line+'\n')
                        records[cell] += 1; retained += 1
                    if total%1000000==0 and time.monotonic()-last>=30:
                        log(f'{item["sample"]} {item["mark"]}: scanned {total:,} rows; retained {retained:,}')
                        last = time.monotonic()
            missing = [cell for cell,n in records.items() if not n]
            if missing: raise ValueError(f'{source}: no fragment records for {len(missing)} final cells; examples: {missing[:5]}. Literal final-object IDs were checked; no cell was silently dropped.')
            if identity(source)!=before: raise ValueError(f'{source}: source changed while subsetting')
            temporary.replace(target)
            frame = pd.DataFrame([dict(cell_id=cell,source_barcode=mapping[cell],fragment_records=records[cell]) for cell in sorted(candidates)])
            frame.to_csv(folder/'selected_cells.tsv',sep='\t',index=False)
            write_json(marker,dict(signature=item['subset_signature'],source=before,
                subset_file=identity(target),mapping_file=identity(folder/'selected_cells.tsv'),
                selected_cells=len(frame),selected_records=retained,source_rows_scanned=total,
                seconds=time.monotonic()-started,rule='Original rows and multiplicity fields retained without deduplication'))
        finally: temporary.unlink(missing_ok=True)
        log(f'{item["sample"]} {item["mark"]}: cached {retained:,} records for all {len(candidates):,} cells')
    return item


def load_selected_native(item, obs, args):
    """Reuse a per-donor native cache, then attach current final labels in RAM."""
    import snapatac2 as snap
    ad = ad_module(); folder = Path(item['cache_dir'])
    if not subset_valid(item): prepare_subset(item,args)
    dest = folder/'selected.native.h5ad'; marker = folder/'IMPORT_COMPLETED.json'
    sig = digest(dict(subset=item['subset_signature'],reference=args.chrom_sizes_dict,
        snap=snap.__version__,import_rule='native_preserve_rows_v2'))
    def valid():
        try:
            info=json.loads(marker.read_text())
            return info['signature']==sig and info['file']==identity(dest)
        except (OSError,ValueError,KeyError): return False
    if not valid():
        with output_lock(folder):
            if not valid():
                mapping = pd.read_csv(folder/'selected_cells.tsv',sep='\t',dtype={'cell_id':str,'source_barcode':str})
                temporary = folder/'importing.tmp.h5ad'; temporary.unlink(missing_ok=True)
                try:
                    source = snap.pp.import_fragments(str(folder/'selected.fragments.tsv.gz'),
                        chrom_sizes=args.chrom_sizes_dict,file=str(temporary),min_num_fragments=0,
                        sorted_by_barcode=False,whitelist=mapping.source_barcode.tolist(),
                        n_jobs=args.threads,tempdir=str(folder))
                    source.close()
                    imported,keys,has_ref=read_meta(temporary)
                    keys=[k for k in ('fragment_paired','fragment_single') if k in keys]
                    if len(keys)!=1 or not has_ref: raise ValueError(f'{temporary}: unsupported imported fragment structure')
                    names=pd.Index(imported.index.astype(str)); rows=names.get_indexer(mapping.source_barcode)
                    if not names.is_unique or (rows<0).any(): raise ValueError('Fragment import omitted selected cells')
                    fragments,reference=read_fragment_rows(temporary,keys[0],rows)
                    totals=np.diff(fragments.indptr).astype(np.int64)
                    if (totals<=0).any(): raise ValueError('Selected cell has no imported fragments in the specified reference')
                    # No export metadata is copied into this shared fragment cache.
                    metadata=pd.DataFrame({'source_barcode':mapping.source_barcode.to_numpy(),
                        'total_native_fragments':totals},index=pd.Index(mapping.cell_id,name='cell_id'))
                    cached=ad.AnnData(X=sp.csr_matrix((len(mapping),0)),obs=metadata)
                    cached.obsm[keys[0]]=fragments;cached.uns['reference_sequences']=reference
                    cached.uns['signature']=sig;save_h5ad(cached,dest)
                    write_json(marker,dict(signature=sig,file=identity(dest),cells=len(mapping)))
                    log(f'{item["sample"]} {item["mark"]}: saved reusable native-fragment cache')
                finally: temporary.unlink(missing_ok=True)
    cached_obs,keys,_=read_meta(dest)
    rows=cached_obs.index.get_indexer(obs.index)
    if not cached_obs.index.is_unique or (rows<0).any(): raise ValueError('Shared native cache misses final physical cell IDs')
    key=next(k for k in keys if k in ('fragment_paired','fragment_single'))
    fragments,reference=read_fragment_rows(dest,key,rows)
    result=ad.AnnData(X=sp.csr_matrix((len(obs),0)),obs=pd.DataFrame({
        'sample':obs['sample'].astype(str).to_numpy(),
        'lymphoid_final':obs['lymphoid_final'].astype(str).to_numpy(),
        'merged_cell_id':obs.index.to_numpy(dtype=str),
        'total_native_fragments':np.diff(fragments.indptr).astype(np.int64)},index=obs.index.copy()))
    result.obsm[key]=fragments;result.uns['reference_sequences']=reference
    return result


def fragment_item(root, cache, mark, donor, obs):
    path = final_fragment_path(root, donor, mark)
    if not path.is_file():
        raise FileNotFoundError(f'{donor} {mark}: missing required stage 02 fragments: {path}')
    names = barcode_candidates(obs, donor)
    stamp = identity(path)
    signature = digest(dict(source=stamp, candidates=names, format='exact_literal_subset_v2'))
    return dict(sample=donor, mark=mark, path=str(path), source=stamp, barcode_candidates=names,
                subset_signature=signature, cache_dir=str(cache / donor / mark / signature))


def region_counts(native, plan, args):
    import snapatac2 as snap
    result = snap.pp.make_peak_matrix(native, use_rep=plan.ids.tolist(), use_x=False,
        counting_strategy='fragment', chunk_size=args.chunk_size, inplace=False)
    ids = pd.Index(result.var_names.astype(str))
    cells = pd.Index(result.obs_names.astype(str))
    rows, cols = cells.get_indexer(native.obs_names), ids.get_indexer(plan.ids)
    if (rows < 0).any() or (cols < 0).any() or not cells.is_unique or not ids.is_unique:
        raise ValueError('SnapATAC2 changed or omitted cell or interval identities')
    return raw_integer(result.X)[rows, :][:, cols]


def cohorts(final, publication, scopes, args):
    sources, dge, effect = {}, {}, {}
    for scope in scopes:
        all_sources = {m: source(final, publication, scope, m) for m in ('RNA', *MARKS)}
        zones = {m: select(s, s['obs'][LABEL].isin([DZ,LZ])) for m,s in all_sources.items()}
        donors, coverage = eligible_donors({m:s['obs'] for m,s in zones.items()}, args)
        required = args.min_donors if scope == 'merged' and 'donor_replicates' in args.merged_methods else 1
        if len(donors) < required:
            display = coverage.pivot(index=['donor','modality'], columns='zone', values='n_cells')
            raise ValueError(f'{scope}: need {required} donors with ≥{args.min_cells} cells in both zones in each modality; '
                             f'found {len(donors)}\n{display.to_string()}')
        dge[scope] = {m:select(s, s['obs']['sample'].isin(donors)) for m,s in zones.items()}
        sources[scope] = all_sources
        if scope == 'merged' and args.k9_dge:
            k9 = source(final, publication, scope, 'H3K9me3')
            k9 = select(k9, k9['obs'][LABEL].isin([DZ,LZ]))
            k9_donors, _ = eligible_donors({'H3K9me3':k9['obs']}, args)
            if len(k9_donors) >= args.min_k9_donors:
                all_sources['H3K9me3'] = k9
                dge[scope]['H3K9me3'] = select(k9, k9['obs']['sample'].isin(k9_donors))
            else:
                log(f'{scope}: skipping the separate H3K9me3 DGE cohort; only {len(k9_donors)} qualifying donors')
        matched = all_sources['RNA']['obs']
        matched = matched.loc[matched[LABEL].isin([DZ,'B_Cycling',LZ]) & matched['sample'].isin(donors)].copy()
        for mark in MARKS:
            matched = matched.loc[matched.index.isin(all_sources[mark]['obs'].index)]
        matched['_state'] = matched[LABEL]
        matched['_donor'] = matched['sample']
        effect[scope] = matched
        for m,s in dge[scope].items():
            log(f'{scope} {m}: {len(s["obs"]):,} retained DZ/LZ cells across {s["obs"]["sample"].nunique()} qualifying donor(s)')
        log(f'{scope}: {len(matched):,} exact matched GC cells for effect sizes (including cycling)')
    return sources, dge, effect


def fragment_inputs(sources, dge, effect, root, cache):
    result = {}
    for mark in sorted({m for states in dge.values() for m in states if m != 'RNA'}):
        views = {scope: d[mark]['obs'] for scope,d in dge.items() if mark in d}
        dge_ids = pd.Index(np.concatenate([v.index for v in views.values()])).unique()
        chosen = pd.concat([sources[scope][mark]['obs'].loc[obs.index.union(effect[scope].index, sort=False)]
                            if mark in MARKS else obs for scope,obs in views.items()])
        chosen = chosen.loc[~chosen.index.duplicated()].sort_index()
        records = []
        for donor, donor_obs in chosen.groupby('sample', sort=True):
            item = fragment_item(root, cache, mark, donor, donor_obs)
            records.append((item, donor_obs.loc[donor_obs.index.isin(dge_ids)]))
        result[mark] = records
    return result


def histone_count_matrix(items, plan, args):
    chunks, cells = [], []
    for item, obs in items:
        if obs.empty:
            prepare_subset(item,args)
            continue
        log(f'{item["sample"]} {item["mark"]}: counting exact fragment overlaps in {len(plan.ids):,} intervals')
        native = load_selected_native(item, obs, args)
        chunks.append(region_counts(native, plan, args))
        cells.extend(obs.index)
    if not chunks:
        raise ValueError('No histone counts were produced')
    return sp.vstack(chunks, format='csr'), pd.Index(cells)


def zone_counts(x, obs, genes, scope, method, args):
    if scope == 'merged' and method == 'random_subsets':
        return pooled_pseudobulk(x, obs, genes, args)
    return pseudobulk(x, obs, genes, scope, args)


def method_dir(out, scope, method):
    return out / scope / ('donor_grouped' if method == 'donor_replicates' else 'pseudo_grouped')


DZ_COLOR, LZ_COLOR, GRAY = '#A74470','#168C83','#C5CDD3'


def highlights(table,args,rna_genes=(),limit=None):
    d = table.loc[np.isfinite(table.logFC)&np.isfinite(table.FDR)].copy()
    n = args.labels_per_direction if limit is None else limit
    chosen = set(args.highlight_genes)|set(rna_genes)
    # Rank finite plotted results even when no gene passes significance.
    # Raw P breaks flat-FDR ties; effect magnitude and symbol make ties stable.
    ranked=d.assign(_label_p=d.PValue.fillna(np.inf) if 'PValue' in d else np.inf,
        _label_effect=d.logFC.abs())
    for sign in (1,-1):
        hits=ranked.loc[ranked.logFC.ge(0) if sign==1 else ranked.logFC.lt(0)]
        chosen.update(hits.sort_values(['FDR','_label_p','_label_effect','gene'],
            ascending=[True,True,False,True],kind='stable').head(n).gene.astype(str))
    return d.loc[d.gene.isin(chosen)].copy()


def draw_volcano(ax,table,title,args,selected,compact=False):
    d = table.loc[np.isfinite(table.logFC)&np.isfinite(table.FDR)].copy()
    up = (d.FDR<=args.fdr)&(d.logFC>=args.min_lfc)
    down = (d.FDR<=args.fdr)&(d.logFC<=-args.min_lfc)
    ax.set(xlabel='log2 fold change (DZ / LZ)',ylabel='-log10 FDR')
    ax.set_title(title,pad=11)
    if not len(d):
        status=table.attrs.get('fit_status',{})
        if status.get('status')=='not_testable':
            message='Not testable'
        else: message='No genes with finite adjusted P'
        ax.text(.5,.5,message,transform=ax.transAxes,ha='center',va='center'); return
    y = -np.log10(d.FDR.clip(lower=1e-300))
    ax.scatter(d.logFC,y,c=np.where(up,DZ_COLOR,np.where(down,LZ_COLOR,GRAY)),s=7 if compact else 9,
        alpha=.75,linewidths=0,rasterized=False)
    ax.axhline(-np.log10(args.fdr),color='#7D8892',lw=.6,ls='--')
    for v in (-args.min_lfc,args.min_lfc): ax.axvline(v,color='#9AA4AC',lw=.6,ls=':')
    ax.set_ylim(-max(.1,float(y.max())*.025),max(float(y.max())*1.12,-np.log10(args.fdr)*1.15,1))
    axismax = max(float(d.logFC.abs().max())*1.10,args.min_lfc*2,1)
    ax.set_xlim(-axismax,axismax)
    picked = selected.merge(pd.DataFrame({'gene':d.gene,'y':y}),on='gene',how='inner',validate='one_to_one')
    ax.scatter(picked.logFC,picked.y,s=20,facecolors='none',edgecolors='#29343D',lw=.55,zorder=4,rasterized=False)
    # All emphasized points receive a name. The callout columns grow in height
    # with label count rather than dropping names or silently overlapping them.
    for right in (False,True):
        side = picked.loc[picked.logFC.ge(0).eq(right)].sort_values(['y','gene'])
        if not len(side): continue
        bottom,top = ax.get_ylim()
        locations = (side.y.to_numpy()-bottom)/(top-bottom)
        gap = min(.07 if compact else .055,.9/max(1,len(locations)-1))
        locations = np.clip(locations,.035,.965)
        for i in range(1,len(locations)): locations[i] = max(locations[i],locations[i-1]+gap)
        if locations[-1]>.965:
            locations[-1]=.965
            for i in range(len(locations)-2,-1,-1): locations[i] = min(locations[i],locations[i+1]-gap)
        for row,position in zip(side.itertuples(),locations):
            ax.annotate(row.gene,xy=(row.logFC,row.y),xycoords='data',
                xytext=(1.035 if right else -.12,position),textcoords='axes fraction',
                ha='left' if right else 'right',va='center',fontsize=6.5 if compact else 8,
                fontstyle='italic',annotation_clip=False,
                arrowprops=dict(arrowstyle='-',color='#9BA4AC',lw=.45))
    if not compact and len(selected):
        # Keep the y-axis title outside the entire left callout column.
        ax.yaxis.set_label_coords(-.43,.5)


def volcano_job(job):
    table,title,path,kwargs,rna_genes = job
    from types import SimpleNamespace
    args = SimpleNamespace(**kwargs)
    selected = highlights(table,args,rna_genes)
    n = max(int(selected.logFC.ge(0).sum()),int(selected.logFC.lt(0).sum()))
    fig,ax = plt.subplots(figsize=((7.8 if len(selected) else 5.2),max(4.3,2.5+.24*n)))
    fig.subplots_adjust(left=.27 if len(selected) else .16,right=.77 if len(selected) else .96,bottom=.16,top=.88)
    draw_volcano(ax,table,title,args,selected)
    save_pdf(fig,path)
    return str(path)


STATES = ['B_GC_DZ', 'B_Cycling', 'B_GC_LZ']


CONTRASTS = [('B_GC_LZ','B_GC_DZ'), ('B_Cycling','B_GC_DZ'), ('B_GC_LZ','B_Cycling')]


LOG = logging.getLogger(__name__)


def logcpm(count,depth,prior):
    return np.log2((np.asarray(count)+prior)/(np.asarray(depth)+2*prior)*1e6)


def subset_estimates(obs,matrices,depths,catalogs,cfg,labels=None,contrasts=None):
    """Repeated disjoint partitions within donor/state, preserving modality pairs.

    Percentiles describe conditional sampling variability. They are NOT donor
    confidence intervals, and never certify biological equivalence.
    """
    labels=labels or STATES;contrasts=contrasts or CONTRASTS
    started=time.perf_counter();rng=np.random.default_rng(cfg['seed'])
    # Build exactly the original partitions, in the original RNG order. One
    # sparse membership matrix pools all partitions without repeated cell slices.
    pool_ids=[];draw_meta=[]
    donor_ids={donor:{s:np.flatnonzero((obs['_donor'].to_numpy()==donor)&
                    (obs['_state'].to_numpy()==s)) for s in labels}
               for donor in sorted(obs['_donor'].unique())}
    for repeat in range(cfg['subset_repeats']):
        for donor,state_ids in donor_ids.items():
            for high,low in contrasts:
                parts=min(cfg['subset_parts'],len(state_ids[high])//cfg['min_cells'],len(state_ids[low])//cfg['min_cells'])
                if parts<2: continue
                aa=np.array_split(rng.permutation(state_ids[high]),parts)
                bb=np.array_split(rng.permutation(state_ids[low]),parts)
                for a,b in zip(aa,bb):
                    pool_ids.extend((a,b));draw_meta.append((f'{high}_vs_{low}',donor))
    if not draw_meta:
        LOG.info('Random subsets: no donor/contrast has enough cells for two disjoint partitions')
        return {key:pd.DataFrame() for key in matrices}
    membership=sp.csr_matrix((np.ones(sum(map(len,pool_ids)),dtype=np.int64),
        (np.repeat(np.arange(len(pool_ids)),list(map(len,pool_ids))),np.concatenate(pool_ids))),
        shape=(len(pool_ids),len(obs)))
    del pool_ids
    groups={}
    for contrast in [f'{a}_vs_{b}' for a,b in contrasts]:
        ids=np.array([i for i,x in enumerate(draw_meta) if x[0]==contrast],dtype=int)
        if not len(ids): continue
        donors=sorted({draw_meta[i][1] for i in ids})
        by_donor=[np.array([i for i in ids if draw_meta[i][1]==d],dtype=int) for d in donors]
        groups[contrast]=(ids,by_donor)
    LOG.info('Random subsets: %s paired partitions, %s repeats; pooling in blocks of %s features',
             len(draw_meta),cfg['subset_repeats'],cfg.get('stats_chunk_features',2048))
    summaries={}
    for key,m in matrices.items():
        modality_start=time.perf_counter();n_features=m.shape[1]
        LOG.info('Random subsets %s: pooling %s features',key,n_features)
        # CSC makes each feature block a bounded slice. Retain only summaries,
        # rather than all partition-by-feature draws for all modalities.
        columns=m.tocsc();lib=np.asarray(membership@depths[key]).ravel()
        results={contrast:{name:np.empty(n_features,dtype=bool if name=='coverage' else float)
                 for name in ('log2FC','P10','P90','coverage','fraction_within_margin')}
                 for contrast in groups}
        chunk=int(cfg.get('stats_chunk_features',2048));last_log=time.perf_counter()
        for first in range(0,n_features,chunk):
            last=min(first+chunk,n_features)
            counts=(membership@columns[:,first:last]).toarray()
            ca,cb=counts[::2],counts[1::2]
            effects=logcpm(ca,lib[::2,None],cfg['prior_count'])-logcpm(cb,lib[1::2,None],cfg['prior_count'])
            total=ca+cb
            for contrast,(ids,by_donor) in groups.items():
                all_effect=effects[ids];result=results[contrast]
                # Equal donor weight and the original conditional percentiles.
                result['log2FC'][first:last]=np.median(np.stack([np.median(effects[ix],axis=0) for ix in by_donor]),axis=0)
                result['P10'][first:last],result['P90'][first:last]=np.quantile(all_effect,[.10,.90],axis=0)
                result['coverage'][first:last]=np.stack([total[ix].mean(axis=0)>=cfg['min_total_count'] for ix in by_donor]).all(axis=0)
                result['fraction_within_margin'][first:last]=np.mean(np.abs(all_effect)<cfg['histone_margin_log2'],axis=0)
            now=time.perf_counter()
            if now-last_log>=10 or last==n_features:
                LOG.info('Random subsets %s: %s/%s features (%.1fs)',key,last,n_features,now-modality_start)
                last_log=now
        tables=[]
        for contrast,(ids,by_donor) in groups.items():
            result=results[contrast];table=catalogs[key].copy();table['contrast']=contrast
            for name in ('log2FC','P10','P90'): table[name]=result[name]
            table['n_partitions']=len(ids);table['n_donors']=len(by_donor)
            for name in ('coverage','fraction_within_margin'): table[name]=result[name]
            table['inference']='conditional_random_subsets';tables.append(table)
        summaries[key]=pd.concat(tables,ignore_index=True) if tables else pd.DataFrame()
        del columns,counts,effects,total,results
    LOG.info('Random subsets finished in %.1fs',time.perf_counter()-started)
    return summaries


def display_gene_mask(table,cfg):
    """Filter displayed symbols, falling back to an ID only when unnamed.

    Canonical ENSG IDs of normally named genes must not exclude every human
    gene from the shortlist. Raw matrices/normalization denominators stay intact.
    """
    prefixes=tuple(cfg.get('exclude_gene_prefixes',['LINC','ENSG']))
    if not prefixes: return pd.Series(True,index=table.index)
    names=table.gene_name.astype('string').fillna('').str.strip()
    fallback=table.gene_id.astype('string').fillna('').str.strip()
    labels=names.mask(names.eq(''),fallback)
    return ~labels.str.startswith(prefixes,na=False)


def directional_gene_order(rna_stats,cfg):
    """RNA-only shortlist: signed effect ranking or lowest FDR per direction.

    A sign is an endpoint difference, not a requirement of monotonic expression
    along mDC. Significance is retained for review but does not filter the list.
    """
    columns=['gene_id','gene_name','RNA_direction','direction_rank','row_order',
             'RNA_LZ_vs_DZ_log2FC','RNA_LZ_vs_DZ_q_difference','RNA_LZ_vs_DZ_n_donors']
    if rna_stats.empty: return pd.DataFrame(columns=columns)
    state=rna_stats.loc[rna_stats.contrast.eq('B_GC_LZ_vs_B_GC_DZ')].copy()
    state=state.loc[np.isfinite(state.log2FC)&display_gene_mask(state,cfg)]
    if state.gene_id.duplicated().any(): raise ValueError('Expected one RNA LZ-versus-DZ estimate per gene')
    by_fdr=cfg.get('sort_genes_by_fdr',False)
    if by_fdr:
        if 'q_difference' not in state:state['q_difference']=np.nan
        state=state.loc[np.isfinite(state.q_difference)&state.q_difference.between(0,1)].copy()
        columns+=['RNA_gene_ranking','RNA_LZ_vs_DZ_q_source']
    groups=[];number=int(cfg.get('genes_per_direction',15))
    for direction,increasing in [('increase_DZ_to_LZ',True),('decrease_DZ_to_LZ',False)]:
        part=state.loc[state.log2FC.gt(0) if increasing else state.log2FC.lt(0)]
        if by_fdr:
            part=part.sort_values(['q_difference','log2FC','gene_id'],ascending=[True,not increasing,True],kind='mergesort').head(number)
        else:
            part=part.sort_values(['log2FC','gene_id'],ascending=[not increasing,True],kind='mergesort').head(number)
        order=part[['gene_id','gene_name']].copy();order['RNA_direction']=direction
        order['direction_rank']=np.arange(1,len(order)+1)
        order['RNA_LZ_vs_DZ_log2FC']=part.log2FC.to_numpy()
        order['RNA_LZ_vs_DZ_q_difference']=part.q_difference.to_numpy() if 'q_difference' in part else np.nan
        order['RNA_LZ_vs_DZ_n_donors']=part.n_donors.to_numpy() if 'n_donors' in part else np.nan
        if by_fdr:
            order['RNA_gene_ranking']='FDR'
            order['RNA_LZ_vs_DZ_q_source']=part.RNA_FDR_source.to_numpy() if 'RNA_FDR_source' in part else 'RNA_LZ_DZ_input_FDR'
        groups.append(order)
    result=pd.concat(groups,ignore_index=True);result['row_order']=np.arange(1,len(result)+1)
    return result.reindex(columns=columns)


def body_enhancer_effect_inputs(matrices,catalogs,order,cfg):
    """Select RNA genes and sum body/enhancer counts BEFORE normalization.

    Both source contexts count each fragment at most once. With unknown
    overlap, Var(B+E) <= 2*E[B+E] under the counting-noise model. Factor two
    is therefore a conservative floor; zero-default enhancer genes use one.
    """
    genes=order.gene_id.tolist();output={};metadata={}
    for key,matrix in matrices.items():
        catalog=catalogs[key].reset_index(drop=True).copy()
        catalog['_column']=np.arange(len(catalog))
        context='RNA' if key=='RNA' else 'gene_body'
        body=catalog.loc[catalog.context.eq(context)].set_index('gene_id')
        if not body.index.is_unique: raise ValueError('Duplicate gene/region count columns')
        body=body.reindex(genes)
        if body._column.isna().any(): raise ValueError(f'{key}: a selected gene lacks {context} counts')
        values=matrix[:,body._column.to_numpy(dtype=int)].tocsr()
        info=body.reset_index().drop(columns='_column')
        info['sampling_variance_factor']=1.
        info['signal_definition']='RNA counts' if key=='RNA' else 'padded_gene_body_count + linked_enhancer_count'
        if key!='RNA':
            enhancer=catalog.loc[catalog.context.eq('enhancer')].set_index('gene_id')
            if not enhancer.index.is_unique: raise ValueError('Duplicate enhancer count columns')
            enhancer=enhancer.reindex(genes)
            valid=enhancer._column.notna().to_numpy()
            unresolved=np.zeros(len(genes),bool)
            for source in (body,enhancer):
                if 'enhancer_mapping_complete' in source:
                    unresolved|=source.enhancer_mapping_complete.eq(False).to_numpy()
            valid&=~unresolved
            positions=np.flatnonzero(valid)
            if len(positions):
                selected=matrix[:,enhancer.loc[valid,'_column'].to_numpy(dtype=int)].tocsr()
                mapping=sp.csr_matrix((np.ones(len(positions),dtype=np.int64),(np.arange(len(positions)),positions)),
                                          shape=(len(positions),len(genes)))
                values=values+selected@mapping
            info['context']='gene_body_plus_enhancer'
            info['feature_id']=info.gene_id+'|gene_body_plus_enhancer'
            info['enhancer_zero_default']=~valid
            info['enhancer_status']=np.where(unresolved,'unresolved_mapping_zero',np.where(valid,'resolved','missing_mapping_zero'))
            info['sampling_variance_factor']=np.where(valid,2.,1.)
            info['gene_body_flank_bp']=cfg.get('gene_body_flank_bp',2000)
        output[key]=values;metadata[key]=info
    return output,metadata


CONTRAST = f'{LZ}_vs_{DZ}'


def gene_ids(var):
    column = next((x for x in ('gene_id','gene_ids') if x in var), None)
    values = var[column].astype('string').fillna('') if column else pd.Series(var.index, index=var.index, dtype='string')
    values = values.where(values.ne(''), pd.Series(var.index, index=var.index, dtype='string'))
    return pd.Index([ensembl_id(v) or str(v) for v in values], name='gene_id')


def effect_annotation(gtf, chrom_sizes, allowed, flank):
    body = read_gene_bodies(gtf, flank, chrom_sizes).rename(columns={'input_name':'gene_id','gene':'gene_name'})
    body['gene_id'] = body.gene_id.map(lambda v: ensembl_id(v) or v)
    body = body.loc[body.gene_id.isin(allowed)]
    conflicting = body.groupby('gene_id').agg(chroms=('chrom','nunique'), strands=('strand','nunique'))
    if conflicting.gt(1).any(axis=None):
        raise ValueError('Some selected genome IDs span multiple contigs or strands; use a matching primary-assembly annotation')
    return body.groupby('gene_id',as_index=False).agg(gene_name=('gene_name','first'), chrom=('chrom','first'),
        start=('start','min'), end=('end','max'))


def choose_order(dge, annotation, cfg):
    # DZ–LZ DESeq2 writes DZ / LZ. The effect figure's axis is LZ / DZ.
    symbols = annotation.loc[~annotation.gene_name.duplicated(keep=False), ['gene_id','gene_name']]
    rows = dge[['gene','logFC','FDR']].rename(columns={'gene':'gene_name'})
    ranked = rows.merge(symbols, on='gene_name', validate='one_to_one')
    ranked = ranked.assign(contrast=CONTRAST, log2FC=-ranked.logFC, q_difference=ranked.FDR,
                           RNA_FDR_source='matching_DZ_LZ_RNA_FDR')
    return directional_gene_order(ranked, cfg)


def region_features(annotation, enhancers, genes):
    body = annotation.set_index('gene_id').loc[genes].reset_index()
    symbols = annotation.groupby('gene_name').gene_id.agg(list).to_dict()
    exact_ids = set(annotation.gene_id)
    parts = {}
    ambiguous = set()
    for row in enhancers.itertuples(index=False):
        direct = ensembl_id(row.gene) or row.gene
        candidates = [direct] if direct in exact_ids else symbols.get(row.gene, [])
        if len(candidates) == 1:
            parts.setdefault(candidates[0], set()).add((row.chrom,row.start,row.end))
        elif len(candidates) > 1:
            ambiguous.update(candidates)
    catalog, intervals = [], []
    for row in body.itertuples(index=False):
        for context, components in (('gene_body',[(row.chrom,row.start,row.end)]), ('enhancer',sorted(parts.get(row.gene_id,())))):
            if not components:
                continue
            index = len(catalog)
            catalog.append(dict(feature_id=f'{row.gene_id}|{context}', gene_id=row.gene_id,
                gene_name=row.gene_name, context=context, enhancer_mapping_complete=row.gene_id not in ambiguous))
            intervals.extend((chrom,start,end,index) for chrom,start,end in components)
    return pd.DataFrame(catalog), pd.DataFrame(intervals, columns=['chrom','start','end','feature_index'])


def fragment_effect_counts(item, ids, chunk_size, catalog, intervals, chrom_sizes):
    from ncls import NCLS
    folder = Path(item['cache_dir'])
    mapping = pd.read_csv(folder/'selected_cells.tsv', sep='\t', dtype=str).set_index('cell_id')
    if not ids.isin(mapping.index).all() or not mapping.source_barcode.is_unique:
        raise ValueError(f'{item["path"]}: selected fragment rows do not have unique exact cell IDs')
    lookup = {value:j for j,value in enumerate(mapping.loc[ids,'source_barcode'])}
    trees = {chrom:NCLS(frame.start.to_numpy(np.int64), frame.end.to_numpy(np.int64),
                        frame.feature_index.to_numpy(np.int64)) for chrom,frame in intervals.groupby('chrom',sort=False)}
    depth = np.zeros(len(ids), dtype=np.int64)
    x = sp.csr_matrix((len(ids),len(catalog)), dtype=np.int64)
    chunks = {}
    def collect():
        nonlocal x
        out_r, out_c = [], []
        for chrom, items in chunks.items():
            block = np.asarray(items, dtype=np.int64)
            query, column = trees[chrom].all_overlaps_both(block[:,0],block[:,1],np.arange(len(block),dtype=np.int64))
            if len(query):
                hits = np.unique(np.column_stack([query,column]),axis=0)
                out_r.append(block[hits[:,0],2])
                out_c.append(hits[:,1])
        if out_r:
            r, c = np.concatenate(out_r), np.concatenate(out_c)
            x += sp.csr_matrix((np.ones(len(r),dtype=np.int64),(r,c)),shape=x.shape)
        chunks.clear()
    n = on_reference = 0
    with gzip.open(folder/'selected.fragments.tsv.gz','rt') as handle:
        for line in handle:
            chrom, start, end, barcode, *_ = line.rstrip('\r\n').split('\t')
            ix = lookup.get(barcode)
            if ix is None:
                continue
            start, end = int(start), int(end)
            if not 0 <= start < end:
                raise ValueError(f'{item["path"]}: invalid fragment interval')
            depth[ix] += 1
            n += 1
            on_reference += chrom in chrom_sizes
            if chrom in trees:
                chunks.setdefault(chrom,[]).append((start,end,ix))
            if n % chunk_size == 0:
                collect()
    collect()
    if (depth == 0).any() or on_reference < n / 2:
        raise ValueError(f'{item["path"]}: missing selected fragments or incompatible contigs')
    return x, depth


def paired_estimates(obs, matrices, depths, metadata, cfg):
    pairs = []
    for donor in sorted(obs._donor.unique()):
        d = obs._donor.eq(donor)
        lz = np.flatnonzero(d & obs._state.eq(LZ))
        dz = np.flatnonzero(d & obs._state.eq(DZ))
        if min(len(lz),len(dz)) >= cfg['min_cells']:
            pairs.append((donor, lz, dz))
    result = {}
    if not pairs:
        return result
    for key, x in matrices.items():
        lz = np.stack([np.asarray(x[a].sum(axis=0)).ravel() for _,a,_ in pairs])
        dz = np.stack([np.asarray(x[b].sum(axis=0)).ravel() for _,_,b in pairs])
        size_a = np.array([depths[key][a].sum() for _,a,_ in pairs])[:,None]
        size_b = np.array([depths[key][b].sum() for _,_,b in pairs])[:,None]
        prior = cfg['prior_count']
        delta = logcpm(lz,size_a,prior) - logcpm(dz,size_b,prior)
        keep = lz + dz >= cfg['min_total_count']
        n = keep.sum(axis=0)
        mean = np.divide(np.where(keep,delta,0).sum(axis=0),n,out=np.full(x.shape[1],np.nan),where=n>0)
        variance = np.divide((np.where(keep,delta-mean,0)**2).sum(axis=0),n-1,
                             out=np.full(x.shape[1],np.nan),where=n>1)
        factor = metadata[key].sampling_variance_factor.to_numpy()
        noise = factor*(1/(lz+prior)+1/(dz+prior))/(np.log(2)**2)
        se = np.maximum(np.sqrt(variance/np.maximum(n,1)),
                        np.sqrt(np.where(keep,noise,0).sum(axis=0))/np.maximum(n,1))
        width = stats.t.ppf(.975,np.maximum(n-1,1))*se
        supported = n >= cfg['min_donors']
        result[key] = metadata[key].assign(contrast=CONTRAST, log2FC=mean, n_donors=n, count_support=supported,
            CI95_low=np.where(supported,mean-width,np.nan), CI95_high=np.where(supported,mean+width,np.nan))
    return result


def estimates(obs, matrices, depths, catalogs, order, cfg, method, scope):
    if order.empty or obs.empty:
        return pd.DataFrame()
    matrices, metadata = body_enhancer_effect_inputs(matrices,catalogs,order,cfg)
    if method == 'donor_replicates':
        result = paired_estimates(obs,matrices,depths,metadata,cfg)
    else:
        if scope == 'merged':
            obs = obs.assign(_donor='Merged_pool')
        result = subset_estimates(obs,matrices,depths,metadata,cfg)
        result = {m:t.loc[t.contrast.eq(CONTRAST)].assign(count_support=lambda z:z.coverage.eq(True))
                  for m,t in result.items() if not t.empty}
    rows = [t.assign(modality=m) for m,t in result.items() if not t.empty]
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True).reindex(columns=['gene_id','contrast','modality','log2FC',
            'CI95_low','CI95_high','P10','P90','count_support'])


def render_integrated_effects(sources, cohort, dge, fragments, enhancers, annotation, args, out, cfg):
    raw = {}
    allowed = set()
    for scope, obs in cohort.items():
        if obs.empty:
            continue
        s = sources[scope]['RNA']
        counts, var = count_matrix(s, chosen_rows(s,obs.index))
        ids = gene_ids(var)
        totals = np.asarray(counts.sum(axis=1)).ravel()
        if scope == 'Tonsil1_gt600' and not (totals > 600).all():
            raise ValueError('Saved >600 RNA object contains a matched GC cell with ≤600 raw UMIs')
        raw[scope] = counts, ids, totals
        allowed.update(ids)
    genes = effect_annotation(annotation,args.chrom_sizes_dict,allowed,cfg['gene_body_flank_bp'])
    orders = {(scope,method):choose_order(table,genes.loc[genes.gene_id.isin(raw[scope][1])],cfg) if scope in raw
              else choose_order(table,genes.iloc[:0],cfg) for (scope,method),table in dge.items()}
    selected = sorted({g for order in orders.values() for g in order.gene_id})
    usable = {scope for (scope,_),order in orders.items() if not order.empty}
    all_cells = pd.Index([], dtype=str)
    for scope in cohort:
        if scope in usable:
            all_cells = all_cells.union(cohort[scope].index, sort=False)
    histone = {}
    if selected and len(all_cells):
        catalog, intervals = region_features(genes,enhancers,selected)
        for mark in MARKS:
            chunks, libs, cells = [], [], []
            for item,_ in fragments[mark]:
                ids = pd.Index([i for i in all_cells if i.startswith(item['sample'] + '::')])
                if not len(ids):
                    continue
                log(f'{item["sample"]} {mark}: counting selected features for the effect-size plot')
                x, depth = fragment_effect_counts(item,ids,50000,catalog,intervals,args.chrom_sizes_dict)
                chunks.append(x); libs.append(depth); cells.extend(ids)
            histone[mark] = sp.vstack(chunks,format='csr'), np.concatenate(libs), pd.Index(cells), catalog

    for scope, obs in cohort.items():
        pending = {m:o for (s,m),o in orders.items() if s == scope}
        matrices, depths, catalogs = {}, {}, {}
        if scope in usable:
            counts, ids, totals = raw[scope]
            needed = pd.Index(sorted({g for order in pending.values() for g in order.gene_id}))
            cols = needed.get_indexer(ids)
            valid = cols >= 0
            projection = sp.csr_matrix((np.ones(valid.sum(),dtype=np.int64),(np.flatnonzero(valid),cols[valid])),
                                       shape=(len(ids),len(needed)))
            matrices['RNA'] = (counts @ projection).tocsr()
            depths['RNA'] = totals
            catalogs['RNA'] = genes.set_index('gene_id').loc[needed,['gene_name']].rename_axis('gene_id').reset_index().assign(
                feature_id=lambda t:t.gene_id+'|RNA', context='RNA')
            for mark in MARKS:
                x, depth, index, catalog = histone[mark]
                rows = index.get_indexer(obs.index)
                if (rows < 0).any() or not index.is_unique:
                    raise ValueError(f'{scope}/{mark}: invalid cross-modal cell ordering')
                matrices[mark] = x[rows]
                depths[mark] = depth[rows]
                catalogs[mark] = catalog
        for method, order in pending.items():
            values = estimates(obs,matrices,depths,catalogs,order,cfg,method,scope)
            render_effects(values, order, scope, method, obs['sample'].nunique(), cfg,
                method_dir(out,scope,method)/'all_genes_effect_sizes.pdf')


MARK_COLORS={'RNA':'#344255','H3K27ac':'#8A4D9E','H3K27me3':'#D27F28','H3K9me3':'#278774'}


DIRECTION_LABELS={'increase_DZ_to_LZ':'RNA up: DZ to LZ','decrease_DZ_to_LZ':'RNA down: DZ to LZ'}


DIRECTION_COLORS={'increase_DZ_to_LZ':'#357C9A','decrease_DZ_to_LZ':'#B66A3C'}


def no_data(ax,text='Insufficient data'):
    ax.text(.5,.5,text,ha='center',va='center',transform=ax.transAxes,color='#667085',wrap=True)
    ax.set(xticks=[],yticks=[])


def no_genes_message(cfg):
    if cfg.get('sort_genes_by_fdr',False):return 'No genes with available RNA FDR pass the selection filters'
    return 'No eligible RNA genes'


def region_label(context,cfg):
    flank=cfg.get('gene_body_flank_bp',2000)
    body='gene body'+(f' +/- {flank/1000:g} kb' if flank else '')
    if context=='gene_body': return body
    if context=='gene_body_plus_enhancer': return body+' + enhancers'
    return context.replace('_',' ')


def effect_limit(effects,cfg):
    if effects.empty: return max(1.,cfg['histone_margin_log2']*1.25)
    values=effects[['log2FC','CI95_low','CI95_high','P10','P90']].to_numpy(dtype=float)
    finite=np.abs(values[np.isfinite(values)])
    return max(1.,float(finite.max())*1.10 if len(finite) else 1.,cfg['histone_margin_log2']*1.25)


def effect_panel(ax,effects,order,key,cfg,method,limit):
    title=f'{key} | LZ vs DZ'
    if key!='RNA': title+='\n'+region_label('gene_body_plus_enhancer',cfg)
    ax.set_title(title)
    if order.empty: return no_data(ax,no_genes_message(cfg))
    if effects.empty: return no_data(ax,'No LZ-versus-DZ estimates')
    table=effects.loc[effects.modality.eq(key)&effects.contrast.eq('B_GC_LZ_vs_B_GC_DZ')]
    if table.gene_id.duplicated().any(): raise ValueError('Effect-size panel requires one LZ/DZ estimate per gene')
    table=table.set_index('gene_id').reindex(order.gene_id)
    handles=[]
    if key!='RNA':
        margin=cfg['histone_margin_log2'];ax.axvspan(-margin,margin,color='#E4E8ED',zorder=0)
        handles.append(Patch(facecolor='#E4E8ED',label=f'Histone margin (+/-{margin:g})'))
    color=MARK_COLORS[key];hollow=False;missing=False
    for i,(_,row) in enumerate(table.iterrows()):
        def number(name):
            value=row.get(name,np.nan)
            return float(value) if pd.notna(value) else np.nan
        value=number('log2FC')
        if not np.isfinite(value):
            ax.text(.98,i,'NA',ha='right',va='center',transform=ax.get_yaxis_transform(),fontsize=7,color='#667085')
            missing=True;continue
        lo,hi=(number('CI95_low'),number('CI95_high')) if method=='donor_replicates' else (number('P10'),number('P90'))
        interval=np.isfinite(lo) and np.isfinite(hi)
        supported=interval if method=='donor_replicates' else row.get('count_support',False)==True
        if interval: ax.plot([lo,hi],[i,i],color=color,lw=1.2,ls='-' if supported else '--',rasterized=False,zorder=2)
        ax.plot(value,i,'o',color=color,markerfacecolor=color if supported else 'none',ms=4,rasterized=False,zorder=3)
        hollow|=not supported
    if 'RNA_direction' in order:
        groups=order.RNA_direction.to_numpy()
        for boundary in np.flatnonzero(groups[1:]!=groups[:-1])+1:
            ax.axhline(boundary-.5,color='#344054',lw=.8)
        if key=='RNA':
            handles.extend(Patch(facecolor=DIRECTION_COLORS[d],label=DIRECTION_LABELS[d]) for d in order.RNA_direction.unique())
    if hollow:
        handles.append(Line2D([],[],marker='o',markerfacecolor='none',color=color,ls='none',
            label='CI unavailable' if method=='donor_replicates' else 'Low count support'))
    if missing: handles.append(Line2D([],[],color='none',label='NA: not estimable'))
    interval_label='95% CI when supported' if method=='donor_replicates' else 'P10-P90 subset range'
    ax.axvline(0,color='#77808D',lw=.8,zorder=1)
    ax.set(yticks=np.arange(len(order)),yticklabels=order.gene_name.fillna(order.gene_id).tolist(),
           ylim=(len(order)-.5,-.5),xlim=(-limit,limit),xlabel='log2 FC (LZ / DZ)\n'+interval_label)
    ax.tick_params(axis='y',length=0,labelsize=8);ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.grid(axis='x',color='#D0D5DD',lw=.5,alpha=.5,zorder=0)
    if handles:
        ax.set_title(title,pad=25)
        ax.legend(handles=handles,frameon=False,ncol=min(3,len(handles)),fontsize=7,columnspacing=.7,
                  handlelength=1.2,loc='lower left',bbox_to_anchor=(0,1.005),borderaxespad=0)


def render_effects(values, order, scope, method, donors, cfg, path):
    fig, axes = plt.subplots(1,3,figsize=(20.4,max(5.,.24*len(order)+2.5)),layout='constrained')
    limit = effect_limit(values,cfg)
    for ax, key in zip(axes, ('RNA','H3K27ac','H3K27me3')):
        effect_panel(ax,values,order,key,cfg,method,limit)
    mode = 'Paired biological donors' if method == 'donor_replicates' else (
        'Pooled random subsets (exploratory)' if scope == 'merged' else 'Random subsets (exploratory)')
    fig.suptitle(f'{title_for(scope)} | {mode} | {donors} matched donor(s) | LZ vs DZ',
                 fontsize=14, fontweight='bold')
    save_pdf(fig,path)


def write_dge(counts, meta, scope, method, name, args, out):
    paired = method == 'donor_replicates'
    keep = counts.ge(args.min_count).sum(axis=1).ge(args.min_count_units)
    units = counts.columns[counts.loc[keep].sum(axis=0).eq(0)]
    reason = ''
    if keep.sum() < 2:
        reason = 'Fewer than two genes pass the pseudobulk count filter'
    elif len(units):
        reason = f'Empty pseudobulk after count filtering: {", ".join(units)}'
    if reason:
        table = pd.DataFrame(dict(gene=counts.index, logFC=np.nan, baseMean=np.nan, PValue=np.nan, FDR=np.nan,
                                  lfcSE=np.nan, stat=np.nan, dispersion=np.nan, count_filter_pass=keep.to_numpy()))
        table.attrs['fit_status'] = dict(status='not_testable', reason=reason)
    else:
        table = fit_model(counts, meta, 'zone', paired, args.fdr, (args.min_count,args.min_count_units))
    unit = 'Paired biological donors' if paired else (
        'Pooled random cell subsets (exploratory)' if scope == 'merged' else 'Random cell subsets (exploratory)')
    table = table.assign(contrast='DZ / LZ', scope=scope, measurement=name, analysis_unit=unit,
                         test_status='not_testable' if reason else 'fitted', test_reason=reason)
    root = method_dir(out,scope,method)
    write_table(table, root / f'{name}.tsv.gz')
    title = name.replace('_combined',' | Gene body + enhancers').replace('_gene_body',f' | Gene body ±{args.gene_body_flank_bp:,} bp')
    volcano_job((table, f'{title_for(scope)} | {title}\n{unit}', root / f'{name}_volcano.pdf', vars(args), ()))
    log(f'{scope}/{method}/{name}: {reason or "fit saved"}')
    return table



def dz_lz_stage(config, cli):
    scopes = list(dict.fromkeys(cli.scopes))
    merged_methods = tuple(dict.fromkeys(cli.merged_methods))
    results = config['output']
    final, publication = results/'final', results/'publication'
    out, cache = results/'dz_lz', results/'.work'/'dz_lz'
    spec,refs = config['dz_lz'],config['refs']
    ensure_references(config,'linked_enhancers','effect_annotation')
    thread_budget(config['threads'])

    import snapatac2 as snap
    args = SimpleNamespace(min_cells=25,min_donors=3,min_k9_donors=2,random_replicates=3,
        min_count=10,min_count_units=2,fdr=.001,min_lfc=.5,labels_per_direction=6,
        highlight_genes=cli.highlight_genes,seed=0,threads=config['threads'],chunk_size=512,
        merged_methods=merged_methods,k9_dge=spec['include_h3k9me3'],
        chrom_sizes_dict=dict(snap.genome.hg38.chrom_sizes),gene_body_flank_bp=spec['gene_body_flank_bp'])
    methods = defaultdict(lambda:('random_subsets',),merged=merged_methods)
    sources,dge,effect = cohorts(final,publication,scopes,args)
    fragments = fragment_inputs(sources,dge,effect,final,cache)
    annotation = gene_annotation()
    effect_annotation = refs['effect_annotation']
    enhancers = enhancer_bed(refs['linked_enhancers'],args.chrom_sizes_dict)
    plan = RegionPlan(annotation,enhancers,args.chrom_sizes_dict,args.gene_body_flank_bp)
    check_r()
    plotting_style()

    effect_cfg = dict(gene_body_flank_bp=spec['effect_gene_body_flank_bp'],min_cells=10,min_donors=3,
        min_total_count=20,prior_count=.5,histone_margin_log2=.5,subset_parts=3,subset_repeats=20,seed=73,
        genes_per_direction=spec['effect_genes_per_direction'],sort_genes_by_fdr=spec['effect_rank_by_fdr'])

    with output_lock(cache / 'run'):
        rna_dge = {}
        for scope, inputs in dge.items():
            s = inputs['RNA']
            raw, var = count_matrix(s)
            x, genes, _ = rna_symbols(raw,var,scope,sorted_genes=True)
            del raw
            for method in methods[scope]:
                counts, meta = zone_counts(x,s['obs'],genes,scope,method,args)
                rna_dge[scope,method] = write_dge(counts,meta,scope,method,'RNA',args,out)
        for mark, items in fragments.items():
            regions, cells = histone_count_matrix(items,plan,args)
            model = copy(args)
            if mark == 'H3K9me3':
                model.min_donors = args.min_k9_donors
            for arm,(genes,projection) in plan.arms.items():
                x = (regions @ projection).tocsr()
                for scope,inputs in dge.items():
                    if mark not in inputs:
                        continue
                    obs = inputs[mark]['obs']
                    rows = cells.get_indexer(obs.index)
                    if (rows < 0).any() or not cells.is_unique:
                        raise ValueError(f'{scope}/{mark}: interval counts cannot be mapped to the input cells')
                    for method in methods[scope]:
                        counts, meta = zone_counts(x[rows],obs,genes,scope,method,model)
                        write_dge(counts,meta,scope,method,f'{mark}_{arm}',args,out)
        render_integrated_effects(sources,effect,rna_dge,fragments,enhancers,effect_annotation,args,out,effect_cfg)
    log(f'DZ–LZ: {out}')


# Command line

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command')
    commands.add_parser('all', help='run the study pairwise and DZ–LZ analyses')
    modes = ('donor_replicates','random_subsets')
    for command,help_text in (
        ('pairwise','RNA differential expression, shared gene-signal UMAPs and gene-group summaries'),
        ('pairwise-summary','Only gene-group summaries from saved DGE tables and finalized profiles')):
        pair = commands.add_parser(command,help=help_text)
        pair.add_argument('--dataset',type=normalize_dataset)
        pair.add_argument('--label-a')
        pair.add_argument('--label-b')
        pair.add_argument('--extra-genes',nargs='*',default=[])
        pair.add_argument('--merged-methods',nargs='+',choices=modes,default=list(modes))
    zone = commands.add_parser('dz-lz',help='zone differential analyses and effect-size plots')
    zone.add_argument('--scopes',nargs='+',type=normalize_dataset,default=list(SCOPES))
    zone.add_argument('--merged-methods',nargs='+',choices=modes,default=list(modes))
    zone.add_argument('--highlight-genes',nargs='*',default=[])
    args = parser.parse_args(argv)
    config = load_study()
    if args.command == 'pairwise':
        pairwise_stage(config, args)
    elif args.command == 'pairwise-summary':
        pairwise_summary_stage(config, args)
    elif args.command == 'dz-lz':
        dz_lz_stage(config, args)
    else:
        require_files(config['refs']['linked_enhancers'])
        pairwise_stage(config, parser.parse_args(['pairwise']))
        dz_lz_stage(config, parser.parse_args(['dz-lz']))


if __name__ == '__main__':
    main()

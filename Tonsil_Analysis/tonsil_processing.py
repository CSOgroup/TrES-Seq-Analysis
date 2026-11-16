#!/usr/bin/env python3
"""Process TrESFlow outputs, transfer lymphoid labels and export final profiles and figures."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from copy import deepcopy
from fcntl import flock, LOCK_EX, LOCK_NB, LOCK_UN
import gc
import gzip
import hashlib
import importlib.metadata
import inspect
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
from tempfile import NamedTemporaryFile, TemporaryDirectory
from types import SimpleNamespace
from urllib.error import URLError
from urllib.request import Request, urlopen

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt, colors
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator, NullLocator
import numpy as np
import pandas as pd
from scipy import sparse as sp
from scipy.stats import gaussian_kde


# Study inputs and shared utilities

ROOT = Path(__file__).resolve().parent
STUDY = ROOT / 'study.json'
CURATION = ROOT / 'cell_curation.tsv.gz'
MODALITIES = ('RNA', 'H3K27ac', 'H3K27me3', 'H3K9me3')
LABEL = 'lymphoid_final'
UNLABELLED = 'Unlabelled'

PROTOCOLS = {'Tonsil1':'SingleTag', 'Tonsil2':'SingleTag', 'Tonsil3':'DualTag', 'Tonsil4':'DualTag'}

UNKNOWN = {'', 'nan', 'none', 'na', 'n/a', '<na>', 'unknown', 'unassigned', 'unlabelled', 'unlabeled'}
DD_SCORE = 'doubletdetection_doublet_score'
DD_STATUS = 'final_doubletdetection_score_status'
DD_SOURCE = 'final_doubletdetection_score_source'
CURATION_FLAGS = ('signal_exclude','fine_exclude','black_exclude','merged_exclude','cdc2_target','fine_review_available')

RNA_METRICS = {
    "min_umi": ("total_counts", "ge"), "max_umi": ("total_counts", "le"),
    "min_genes": ("n_genes_by_counts", "ge"), "max_genes": ("n_genes_by_counts", "le"),
    "max_mt": ("pct_counts_mt", "lt"), "max_hb": ("pct_counts_hb", "lt"),
    "max_rb": ("pct_counts_ribo", "lt"), "max_top50": ("pct_counts_in_top_50_genes", "le"),
}
DNA_METRICS = {
    "min_fragments": ("n_fragment", "ge"), "max_fragments": ("n_fragment", "le"),
    "max_mito": ("frac_mito", "lt"), "max_dup": ("frac_dup", "lt"),
    "max_blacklist": ("frac_blacklist", "lt"), "min_tsse": ("tsse", "ge"),
    "max_nucleosome": ("nucleosome_signal", "le"),
}
COUNT_KEYS = {"min_umi", "max_umi", "min_genes", "max_genes", "min_fragments", "max_fragments"}
PERCENT_KEYS = {"max_mt", "max_hb", "max_rb", "max_top50"}
FRACTION_KEYS = {"max_mito", "max_dup", "max_blacklist"}


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, default=str) + '\n')
    tmp.replace(path)


def ids_hash(ids):
    return hashlib.sha256(('\n'.join(sorted(map(str, ids))) + '\n').encode()).hexdigest()


def exact_obs(obs, sample=None):
    obs = obs.copy()
    if sample:
        obs['sample'] = sample
    required = {'sample', 'library', 'barcode'}
    if not required.issubset(obs):
        raise ValueError(f'Missing identity columns: {required - set(obs)}')
    source = obs.library.astype(str) + '_' + obs.barcode.astype(str)
    exact = obs['sample'].astype(str) + '::' + source
    if not (np.array_equal(obs.index.astype(str), source) or np.array_equal(obs.index.astype(str), exact)):
        raise ValueError('Cell IDs disagree with donor/library/barcode metadata')
    obs['source_cell_id'] = source.to_numpy()
    obs.index = pd.Index(exact, name='cell_id')
    if not obs.index.is_unique:
        raise ValueError('Duplicate full cell ID')
    obs['donor'] = obs['sample'].astype(str)
    obs['tagmentation_protocol'] = obs['sample'].astype(str).map(PROTOCOLS)
    if obs.tagmentation_protocol.isna().any():
        raise ValueError('Unknown donor')
    return obs


def integer_counts(x):
    values = x.data if sp.issparse(x) else np.asarray(x)
    if not np.isfinite(values).all() or (values < 0).any() or not np.allclose(values, np.rint(values), atol=1e-6, rtol=0):
        raise ValueError('Expected finite nonnegative raw integer counts')


def validate_quality(layer, label="quality_thresholds"):
    if not isinstance(layer, dict):
        raise ValueError(f"{label} must be a modality-to-threshold dictionary")
    for modality, limits in layer.items():
        if modality not in MODALITIES or not isinstance(limits, dict):
            raise ValueError(f"{label}: expected RNA or a supported histone mark with a threshold dictionary")
        metrics = RNA_METRICS if modality == "RNA" else DNA_METRICS
        for key, value in limits.items():
            if key not in metrics:
                raise ValueError(f"{label}/{modality}: unknown threshold {key!r}; choices: {sorted(metrics)}")
            if value is None:
                continue
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                raise ValueError(f"{label}/{modality}/{key}: use a finite nonnegative number or null")
            if key in COUNT_KEYS and int(value) != value:
                raise ValueError(f"{label}/{modality}/{key}: counts must be whole numbers")
            upper = 100 if key in PERCENT_KEYS else 1 if key in FRACTION_KEYS else None
            if upper is not None and value > upper:
                raise ValueError(f"{label}/{modality}/{key}: maximum allowed value is {upper}")
        for low, high in (("min_umi", "max_umi"), ("min_genes", "max_genes"),
                          ("min_fragments", "max_fragments")):
            if limits.get(low) is not None and limits.get(high) is not None and limits[low] > limits[high]:
                raise ValueError(f"{label}/{modality}: {low} exceeds {high}")
    return deepcopy(layer)


def quality_filter(frame, modality, limits):
    """Return pass flags and every failed threshold; enabled missing metrics fail loudly."""
    limits = validate_quality({modality: limits})[modality]
    metrics = RNA_METRICS if modality == "RNA" else DNA_METRICS
    failures = pd.DataFrame(index=frame.index)
    for key, cutoff in limits.items():
        if cutoff is None:
            continue
        column, comparison = metrics[key]
        if column not in frame:
            raise ValueError(f"{modality}/{key} needs metric {column!r}, which is unavailable")
        values = pd.to_numeric(frame[column], errors="raise")
        finite = values.notna() & values.ne(float("inf")) & values.ne(float("-inf"))
        failures[key] = ~(finite & getattr(values, comparison)(cutoff).fillna(False))
    if failures.empty:
        return pd.Series(True, index=frame.index), pd.Series("", index=frame.index, dtype=object)
    passed = ~failures.any(axis=1)
    reasons = failures.apply(lambda row: ";".join(row.index[row.to_numpy(bool)]), axis=1)
    return passed, reasons


def doubletdetection_pass(frame, cutoff, eligible=None):
    """Strict > cutoff removal; equality and negative infinity pass.

    Confirmed ambiguous scores stay NaN and are retained by this gate; other
    missing scores block only cells that would otherwise be retained. Nothing
    is converted to zero or a singlet call. Positive infinity fails the gate.
    """
    if cutoff is None:
        return pd.Series(True, index=frame.index)
    if DD_SCORE not in frame:
        raise ValueError(f"RNA raw DoubletDetection gate needs {DD_SCORE!r}; rerun initial doublet scoring if absent")
    scores = pd.to_numeric(frame[DD_SCORE], errors="raise")
    eligible = pd.Series(True, index=frame.index) if eligible is None else pd.Series(eligible, index=frame.index).astype(bool)
    confirmed = frame.get(DD_STATUS, pd.Series("", index=frame.index)).eq("unscored_ambiguous")
    unknown = scores.isna() & ~confirmed & eligible
    if unknown.any():
        raise ValueError(f"RNA raw DoubletDetection scores are missing without a confirmed ambiguous result for "
                         f"{int(unknown.sum())} otherwise-retained cells; examples: {frame.index[unknown].tolist()[:5]}. "
                         "Inspect this run's processing/<sample>/RNA_doublets.tsv.gz; no scores were invented.")
    return scores.le(cutoff).fillna(False) | (scores.isna() & (confirmed | ~eligible))


def identity(path):
    path=Path(path).resolve(); stat=path.stat()
    return dict(path=str(path),size=stat.st_size,mtime_ns=stat.st_mtime_ns)


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,default=str).encode()).hexdigest()


def ad_module():
    # Newer source files may use the optional HDF5 compression filters.
    try: import hdf5plugin  # noqa: F401
    except ImportError: pass
    import anndata as ad
    return ad


def h5_read_elem(node):
    from anndata.io import read_elem
    return read_elem(node)


def read_meta(path):
    """Reading backed AnnData also eagerly reads obsm; avoid that for inventory."""
    import h5py
    ad_module()
    with h5py.File(path,'r') as f:
        obs = h5_read_elem(f['obs'])
        keys = list(f['obsm']) if 'obsm' in f else []
        has_reference = 'uns' in f and 'reference_sequences' in f['uns']
    return obs,keys,has_reference


def save_h5ad(obj, path):
    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    temp = p.with_name(p.stem + '.tmp.h5ad')
    obj.write_h5ad(temp, compression='gzip'); temp.replace(p)


@contextmanager
def output_lock(folder):
    folder = Path(folder); folder.mkdir(parents=True, exist_ok=True)
    with (folder/'.analysis.lock').open('a') as f:
        try: flock(f, LOCK_EX | LOCK_NB)
        except BlockingIOError: raise RuntimeError(f'Another process is using {folder}; resume after it finishes')
        try: yield
        finally: flock(f, LOCK_UN)


def log(message):
    print(message, flush=True)


def thread_budget(n):
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','RAYON_NUM_THREADS'):
        os.environ[key] = str(n)
    from threadpoolctl import threadpool_limits
    threadpool_limits(limits=n)


def write_table(table, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    table.to_csv(temporary, sep='\t', index=False, compression='gzip' if path.suffix == '.gz' else None)
    temporary.replace(path)


def check_layout(out, *inputs):
    dest = out.resolve()
    for value in inputs:
        ref = value.resolve()
        if dest.is_relative_to(ref) or ref.is_relative_to(dest):
            raise ValueError(f'Output/cache directory overlaps an input directory: {dest} and {ref}')


HIGH_COVERAGE = 600


def load_study():
    config = read_json(STUDY)
    config['tresflow'] = (ROOT / config['tresflow_root']).resolve()
    config['output'] = (ROOT / config['output_root']).resolve()
    ref_dir = (ROOT / config['reference_dir']).resolve()
    config['refs'] = {k:(ref_dir / spec['file']).resolve()
                      for k,spec in config['references'].items()}
    sample_dirs = [config['tresflow'] / name for name in PROTOCOLS]
    check_layout(config['output'],ref_dir,*sample_dirs,*config['refs'].values())
    if len(config['samples']) != len(PROTOCOLS) or {
        s['name']:s['tagmentation_protocol'] for s in config['samples']
    } != PROTOCOLS:
        raise ValueError('Expected the four study samples with their recorded protocols')
    for sample in config['samples']:
        validate_quality(sample['quality'],sample['name'])
    return config


def require_files(*paths):
    missing = [str(p) for p in paths if p is None or not Path(p).is_file()]
    if missing:
        raise FileNotFoundError('Missing required input files:\n  ' + '\n  '.join(missing))


def tools_required(*names):
    for name in names:
        if not shutil.which(name):
            raise RuntimeError(f'{name} is missing from PATH. Activate the environment from environment.yml.')


def count_stamp(paths):
    return {name:identity(path) for name,path in paths.items()}


def label_values(values):
    values = values.astype('string').fillna(UNLABELLED).str.strip()
    return values.mask(values.str.lower().isin(UNKNOWN), UNLABELLED).astype(str)


def metadata(path):
    return read_meta(path)[0]


def load_curation():
    table = pd.read_csv(CURATION, sep='\t', low_memory=False)
    if table.duplicated(['cell_id','modality']).any():
        raise ValueError('Duplicate cell/modality in cell_curation.tsv.gz')
    if not table.modality.isin(MODALITIES).all() or any(table[key].dtype != bool for key in CURATION_FLAGS):
        raise ValueError('Invalid modality or flag in cell_curation.tsv.gz')
    for column in (DD_SCORE,'doubletdetection_call'):
        if column not in table:
            table[column] = np.nan
    return {m:table.loc[table.modality.eq(m)].set_index('cell_id') for m in MODALITIES}


def curation_for(curation, modality, ids):
    rows = curation[modality].reindex(ids)
    flags = rows.loc[:,list(CURATION_FLAGS)].eq(True)
    for column in (DD_SCORE, 'doubletdetection_call'):
        flags[column] = rows[column]
    return flags


def rna_paths(config, sample, library):
    folder = config['tresflow'] / sample / 'rna_align' / f'{library}.Solo.outGeneFull' / 'filtered'
    return tuple(folder / name for name in ('matrix.mtx.gz','barcodes.tsv.gz','features.tsv.gz'))


def bam_path(config, sample, library, mark):
    return config['tresflow'] / sample / 'dna_align' / f'{library}_{mark}_MarkedDup.bam'


def processing_paths(output, sample):
    name = sample['name']
    folder = output / 'processing' / name
    marks = {m for library in sample['libraries'].values() for m in library}
    return {'RNA':folder/'RNA_counts.h5ad',
            **{mark:folder/f'{name}_{mark}_candidates.h5ad' for mark in MODALITIES if mark in marks}}


def read_chrom_sizes(path):
    contents = pd.read_csv(path, sep='\t', header=None, usecols=[0,1])
    if contents[0].duplicated().any() or contents[1].le(0).any():
        raise ValueError(f'{path}: invalid chromosome sizes')
    return dict(zip(contents[0].astype(str), contents[1].astype(int)))


def safe_obs(obs):
    obs = obs.copy()
    obs.index = pd.Index(obs.index.to_numpy(dtype=str), dtype=object, name=obs.index.name)
    for key in obs:
        if isinstance(obs[key].dtype, pd.CategoricalDtype) and isinstance(obs[key].cat.categories.dtype, pd.StringDtype):
            categories = pd.Index(obs[key].cat.categories.to_numpy(dtype=str), dtype=object)
            obs[key] = pd.Categorical(obs[key].astype(object), categories=categories, ordered=obs[key].cat.ordered)
        if isinstance(obs[key].dtype, pd.StringDtype):
            obs[key] = obs[key].astype(object).where(obs[key].notna(), None)
    return obs


def save_object(data, path):
    integer_counts(data.X)
    data.obs = safe_obs(data.obs)
    save_h5ad(data,path)


def snap_obs(data):
    frame = data.obs.copy() if isinstance(data.obs,pd.DataFrame) else data.obs[:].to_pandas()
    frame.index = pd.Index(list(map(str,data.obs_names)), name='cell_id')
    return frame


def set_snap_obs(data, frame):
    data.obs = safe_obs(frame)
    data.obs_names = list(frame.index.astype(str))


def ensure_references(config, *keys):
    """Use existing files; fetch missing public references on demand."""
    for key in keys:
        path = config['refs'][key]
        if path.is_file() and path.stat().st_size:
            continue
        spec = config['references'][key]
        # Accept a locally supplied copy of the compressed blacklist as well.
        compressed = path.with_suffix(path.suffix + '.gz')
        if spec.get('gunzip') and compressed.is_file():
            path.parent.mkdir(parents=True,exist_ok=True)
            install_blacklist(compressed,path)
            continue
        if not spec.get('url'):
            require_files(path)
            raise ValueError(f'Empty reference file: {path}')
        download_reference(key,path,spec)


def install_blacklist(compressed, destination):
    with NamedTemporaryFile(dir=destination.parent,prefix=f'.{destination.name}.',delete=False) as temp:
        prepared = Path(temp.name)
    try:
        with gzip.open(compressed,'rb') as source, prepared.open('wb') as out:
            shutil.copyfileobj(source,out,1024*1024)
        if prepared.stat().st_size == 0:
            raise ValueError(f'Empty blacklist after decompression: {compressed}')
        prepared.replace(destination)
    finally:
        prepared.unlink(missing_ok=True)


def download_reference(name, path, spec):
    url = spec['url']
    if not url.startswith('https://'):
        raise ValueError(f'{name}: reference download URLs must use HTTPS')
    path.parent.mkdir(parents=True,exist_ok=True)
    log(f'Downloading {name} from {url}')
    with NamedTemporaryFile(dir=path.parent,prefix=f'.{path.name}.',delete=False) as temp:
        downloaded = Path(temp.name)
        try:
            size = 0
            checksum = hashlib.sha256() if spec.get('sha256') else None
            request = Request(url,headers={'User-Agent':'TrES-seq tonsil-analysis/1.0'})
            with urlopen(request,timeout=60) as response:
                length = response.headers.get('Content-Length')
                for block in iter(lambda: response.read(1024*1024),b''):
                    temp.write(block)
                    size += len(block)
                    if checksum is not None:
                        checksum.update(block)
            if not size or (length is not None and size != int(length)):
                raise ValueError(f'incomplete download ({size:,} bytes)')
            if checksum is not None and checksum.hexdigest() != spec['sha256']:
                raise ValueError('SHA-256 mismatch')
        except (OSError,URLError,ValueError) as exc:
            downloaded.unlink(missing_ok=True)
            raise RuntimeError(f'Could not download {name} from {url}: {exc}') from exc
        except BaseException:
            downloaded.unlink(missing_ok=True)
            raise
    try:
        if spec.get('gunzip'):
            install_blacklist(downloaded,path)
        else:
            with downloaded.open('rb') as handle:
                magic = handle.read(8)
            if path.suffix == '.h5ad' and magic != b'\x89HDF\r\n\x1a\n':
                raise ValueError('expected an HDF5 AnnData file')
            if path.suffix == '.gz' and not magic.startswith(b'\x1f\x8b'):
                raise ValueError('expected a gzip-compressed annotation')
            if name == 'chrom_sizes':
                read_chrom_sizes(downloaded)
            downloaded.replace(path)
    except (OSError,ValueError) as exc:
        raise RuntimeError(f'Invalid downloaded reference {name} ({url}): {exc}') from exc
    finally:
        downloaded.unlink(missing_ok=True)
    log(f'{name}: {path}')


# Stage 01: processing and quality control

def sample_metadata(obs, sample):
    obs = obs.copy()
    obs['sample'] = obs['donor'] = sample['name']
    obs['tagmentation_protocol'] = sample['tagmentation_protocol']
    obs['profiled_marks'] = obs['library'].map({k:','.join(v) for k,v in sample['libraries'].items()})
    obs['assay_group'] = obs.profiled_marks.str.replace('H3','',regex=False).str.replace(',','_',regex=False)
    return obs


def read_rna(config, sample, sc, ad):
    parts = []
    for library in sample['libraries']:
        directory = rna_paths(config,sample['name'],library)[0].parent
        log(f'{sample["name"]}/{library}: {directory}')
        part = sc.read_10x_mtx(directory,var_names='gene_symbols',make_unique=True)
        barcodes = part.obs_names.astype(str).str.removeprefix(library+'_')
        if barcodes.has_duplicates or not barcodes.str.fullmatch(r'\d{2}_[ACGT]{24}').all():
            raise ValueError(f'{directory}: expected unique SB_24nt cell barcodes')
        part.obs['barcode'] = barcodes.to_numpy()
        part.obs['library'] = library
        part.obs_names = pd.Index(library+'_'+barcodes,name='cell_id')
        parts.append(part)
    rna = ad.concat(parts,join='inner',merge='same',index_unique=None)
    integer_counts(rna.X)
    genes = rna.var_names.str.upper()
    rna.var['mt'] = genes.str.startswith('MT-')
    rna.var['ribo'] = genes.str.startswith(('RPS','RPL'))
    rna.var['hb'] = genes.str.match(r'^HB(?!P)')
    sc.pp.calculate_qc_metrics(rna,qc_vars=('mt','ribo','hb'),percent_top=(20,50,100),log1p=False,inplace=True)
    rna.obs = sample_metadata(rna.obs,sample)
    return rna


def score_rna(rna, sample, curation, threads):
    exact = exact_obs(rna.obs,sample['name'])
    technical = exact.index[exact.qc_pass]
    if len(technical) < 3:
        raise ValueError(f'{sample["name"]}: fewer than three RNA cells pass QC')
    supplied = curation['RNA']
    eligible = supplied.index.str.startswith(sample['name']+'::') & supplied.doubletdetection_call.notna()
    supplied = supplied.loc[eligible,[DD_SCORE,'doubletdetection_call']]
    missing = technical.difference(supplied.index)
    extras = supplied.index.difference(technical)
    scored = pd.to_numeric(supplied[DD_SCORE],errors='coerce')
    has_scores = (scored.notna() | supplied.doubletdetection_call.eq('ambiguous')).all()
    if not len(missing) and not len(extras) and has_scores:
        table = supplied.loc[technical].copy()
        source = 'cell_curation.tsv.gz'
        log(f'{sample["name"]}: reusing recorded DoubletDetection scores for {len(table):,} exact QC cells')
    else:
        log(f'{sample["name"]}: recomputing DoubletDetection scores '
            f'({len(missing):,} missing calls; {len(extras):,} extra calls; valid scores: {has_scores})')
        import doubletdetection
        # Sort exact IDs before fitting to keep the model row order deterministic.
        cells = technical.sort_values()
        rows = exact.index.get_indexer(cells)
        x = sp.csr_matrix(rna.X[rows])
        fit = doubletdetection.BoostClassifier(boost_rate=.25,n_components=min(30,min(x.shape)-1),
            n_top_var_genes=min(10000,x.shape[1]),clustering_algorithm='phenograph',n_iters=25,
            pseudocount=.1,random_state=0,verbose=True,standard_scaling=False,n_jobs=threads)
        calls = fit.fit(x).predict(p_thresh=1e-7,voter_thresh=.8)
        table = pd.DataFrame({DD_SCORE:fit.doublet_score(),
                             'doubletdetection_call':np.where(calls==0,'singlet',np.where(calls==1,'doublet','ambiguous'))},
                            index=cells)
        source = 'DoubletDetection'
    table.index.name = 'cell_id'
    matched = table.reindex(exact.index)
    scores = matched[DD_SCORE]
    calls = matched.doubletdetection_call
    status = np.where(scores.notna(),'scored',np.where(calls.eq('ambiguous'),'unscored_ambiguous','missing_unexplained'))
    rna.obs[DD_SCORE] = scores.to_numpy()
    rna.obs['doubletdetection_call'] = calls.fillna('not_tested').to_numpy()
    rna.obs[DD_STATUS] = status
    rna.obs[DD_SOURCE] = np.where(calls.notna(),source,'')
    table['source'] = source
    return table


def process_rna(config, sample, curation, path, sc, ad):
    rna = read_rna(config,sample,sc,ad)
    rna.obs['qc_pass'],rna.obs['quality_failure_reasons'] = quality_filter(rna.obs,'RNA',sample['quality']['RNA'])
    table = score_rna(rna,sample,curation,config['threads'])
    decisions = curation_for(curation,'RNA',exact_obs(rna.obs,sample['name']).index)
    keep = rna.obs.qc_pass.to_numpy() & ~decisions.signal_exclude.to_numpy()
    save_object(rna[keep].copy(),path)
    write_table(table.reset_index(),path.with_name('RNA_doublets.tsv.gz'))
    log(f'{sample["name"]}/RNA: {keep.sum():,} retained of {rna.n_obs:,} candidates')


def import_library(config, sample, library, mark, work, sizes, snap):
    fragments = work/f'{library}_{mark}.fragments.tsv.gz'
    bam = bam_path(config,sample['name'],library,mark)
    log(f'{sample["name"]}/{mark}: reading {bam.name}')
    snap.pp.make_fragment_file(str(bam),str(fragments),barcode_tag='CB',
                              shift_left=0,shift_right=0,min_mapq=30,tempdir=str(work))
    require_files(fragments)
    data = snap.pp.import_fragments(str(fragments),sizes,file=str(work/f'{library}_{mark}.h5ad'),
                                    min_num_fragments=100,sorted_by_barcode=True,tempdir=str(work))
    try:
        obs = snap_obs(data)
        raw = obs['barcode'] if 'barcode' in obs else obs.index
        barcode = pd.Index(raw.astype(str)).str.removeprefix(library+'_')
        if barcode.has_duplicates or not barcode.str.fullmatch(r'\d{2}_[ACGT]{24}').all():
            raise ValueError(f'{bam}: expected unique SB_24nt fragment barcodes')
        obs['barcode'] = barcode.to_numpy()
        obs['library'] = library
        obs['mark'] = mark
        obs.index = pd.Index(library+'_'+barcode,name='cell_id')
        set_snap_obs(data,sample_metadata(obs,sample))
        return data
    except Exception:
        data.close()
        raise


def process_dna(config, sample, mark, curation, path, sizes, snap):
    refs = config['refs']
    with TemporaryDirectory(prefix=mark+'-',dir=path.parent) as tmp:
        work = Path(tmp)
        parts = []
        try:
            for library,marks in sample['libraries'].items():
                if mark in marks:
                    parts.append(import_library(config,sample,library,mark,work,sizes,snap))
            destination = str(work/'all_cells.h5ad')
            data = parts[0].copy(destination) if len(parts)==1 else snap.concat(parts,file=destination,backend='hdf5')
        finally:
            for part in parts:
                part.close()
        try:
            snap.metrics.frip(data,{'frac_blacklist':str(refs['blacklist'])})
            snap.metrics.frag_size_distr(data)
            limits = sample['quality'][mark]
            if mark == 'H3K27ac' or limits.get('min_tsse') is not None:
                snap.metrics.tsse(data,str(refs['tss_annotation']))
            obs = snap_obs(data)
            passed, reasons = quality_filter(obs,mark,limits)
            obs['quality_failure_reasons'] = reasons
            decisions = curation_for(curation,mark,exact_obs(obs,sample['name']).index)
            keep = passed.to_numpy() & ~decisions.signal_exclude.to_numpy()
            if keep.sum() < 3:
                raise ValueError(f'{sample["name"]}/{mark}: fewer than three analysis candidates')
            temporary = work/'selected.h5ad'
            selected = data.subset(np.flatnonzero(keep),out=str(temporary),inplace=False)
            try:
                retained = obs.loc[keep].copy()
                retained['qc_pass'] = retained['analysis_pass'] = True
                set_snap_obs(selected,retained)
            finally:
                selected.close()
            temporary.replace(path)
            log(f'{sample["name"]}/{mark}: {keep.sum():,} retained of {len(obs):,} candidates')
        finally:
            data.close()


def processing_samples(config, sample=None):
    chosen = [s for s in config['samples'] if sample in (None,s['name'])]
    inputs = [CURATION]
    for donor in chosen:
        for library,marks in donor['libraries'].items():
            inputs.extend(rna_paths(config,donor['name'],library))
            inputs.extend(bam_path(config,donor['name'],library,m) for m in marks)
    require_files(*inputs)
    return chosen


def process_stage(config, sample=None):
    chosen = processing_samples(config,sample)
    ensure_references(config,'chrom_sizes','blacklist','tss_annotation')
    refs = config['refs']
    thread_budget(config['threads'])
    ad = ad_module()
    import scanpy as sc
    import snapatac2 as snap

    curation = load_curation()
    sizes = read_chrom_sizes(refs['chrom_sizes'])
    for sample in chosen:
        paths = processing_paths(config['output'],sample)
        paths['RNA'].parent.mkdir(parents=True,exist_ok=True)
        with output_lock(config['output']/'.work'/'processing'/sample['name']):
            process_rna(config,sample,curation,paths['RNA'],sc,ad)
            for mark in MODALITIES[1:]:
                if mark in paths:
                    process_dna(config,sample,mark,curation,paths[mark],sizes,snap)


# Study figures and agreement summaries

T1_MODALITIES = MODALITIES[:3]
SAMPLE_COLORS = dict(zip(PROTOCOLS, ('#0072B2', '#D55E00', '#009E73', '#CC79A7')))


def natural(value):
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r'(\d+)', str(value))]

def save_publication_pdf(fig, path, args):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    for ax in fig.axes:
        ax.tick_params(width=.6, length=2.5)
    # Resolution names contain decimal points (e.g. r1.5); do not replace .5 as a suffix.
    fig.savefig(path.parent / (path.name + '.pdf'), dpi=args.dpi, bbox_inches='tight',
        metadata={'Creator': 'Tonsil analysis', 'Title': path.name})
    plt.close(fig)

def primary_coord(ds):
    for name in ('local_' + ds.recommended, 'cohort_' + ds.recommended, 'cohort_uncorrected'):
        if name in ds.coords:
            return name, ds.coords[name]
    raise ValueError(f'{ds.modality}: no usable saved/recomputed UMAP')

def composition(ds, folder, args, palette):
    table = pd.crosstab(ds.obs[LABEL], ds.obs['sample'].astype(str)).reindex(sorted(set(ds.obs[LABEL]), key=natural))
    fraction = table.divide(table.sum(axis=0), axis=1) * 100
    dest = folder / 'tables' / ds.modality
    dest.mkdir(parents=True, exist_ok=True)
    table.to_csv(dest / 'cell_type_counts.tsv', sep='\t')
    fraction.to_csv(dest / 'cell_type_percentages.tsv', sep='\t')
    counts = table.sum(axis=1).sort_values()
    fig, ax = plt.subplots(figsize=(4, max(2.8, .17 * len(counts) + .7)))
    ax.barh(range(len(counts)), counts / counts.sum() * 100, color=[palette.get(v, '#BBBBBB') for v in counts.index], height=.75)
    for i, n in enumerate(counts):
        ax.text(n / counts.sum() * 100 + .2, i, f'n={n:,}', va='center', fontsize=6)
    ax.set_yticks(range(len(counts)), counts.index)
    ax.set(xlabel='Cells in this modality (%)', title=ds.modality + ' · final cell-type composition')
    ax.margins(x=.2)
    save_publication_pdf(fig, folder / 'composition' / f'{ds.modality}_cell_types', args)



def comparison_heatmap(counts, context, modality, key, ari, nmi, path, args,
                       *, coverage=None, score_n=None):
    """Plot saved labels versus clusters from this modality's own obs table."""
    fractions = counts.divide(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.)
    values = fractions.to_numpy(dtype=float)
    totals = counts.sum(axis=0)
    row_labels = [f'{name} ({int(n):,})' for name, n in counts.sum(axis=1).items()]
    column_labels = [f'{name} ({int(n):,})' for name, n in totals.items()]
    left = max(1.8, .053 * max(map(len, row_labels)) + .35)
    bottom, top, right = 1.6, 1.15, .72
    matrix_width, matrix_height = max(3.6, .36 * counts.shape[1]), max(2.4, .30 * counts.shape[0])
    width, height = left + matrix_width + right, bottom + matrix_height + top
    fig = plt.figure(figsize=(width, height))
    ax = fig.add_axes([left / width, bottom / height, matrix_width / width, matrix_height / height])
    artist = ax.imshow(values, cmap='coolwarm', vmin=0, vmax=1, aspect='auto', interpolation='nearest')
    ax.set_xticks(range(len(column_labels)), column_labels, rotation=45, ha='right', rotation_mode='anchor')
    ax.set_yticks(range(len(row_labels)), row_labels)
    ax.set(xlabel=f'{modality} {key}', ylabel='RNA lymphoid_final')
    ax.set_title(f'P({modality} cluster | lymphoid_final); each row sums to 1', fontsize=9, pad=13)
    ax.tick_params(axis='both', labelsize=8, width=.6, length=3)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(.65)
    for i, j in np.ndindex(values.shape):
        value = values[i, j]
        label = f'{value:.2f}'
        if value > 0 and label == '0.00':
            label = '<0.01'
        ax.text(j, i, label, ha='center', va='center', fontsize=6.5,
                color='white' if value >= .6 else 'black', zorder=2)
    cax = fig.add_axes([(left + matrix_width + .14) / width, bottom / height, .16 / width, matrix_height / height])
    fig.colorbar(artist, cax=cax, ticks=np.linspace(0, 1, 6))
    fmt = lambda value: f'{value:.3f}' if np.isfinite(value) else 'NA'
    n = int(counts.to_numpy().sum())
    score_n = n if score_n is None else int(score_n)
    scored = f'; labelled n={score_n:,}' if score_n != n else ''
    fig.suptitle(f'{context}\nn={n:,}{scored}; ARI={fmt(ari)}; NMI={fmt(nmi)}', x=.5, y=1 - .12 / height, fontsize=10)
    save_publication_pdf(fig, path, args)
    return fractions

def positioned_linked_coordinates(datasets):
    centers = {'RNA': (0., 0.), 'H3K27ac': (-3.0, 1.6), 'H3K27me3': (3.0, 1.6), 'H3K9me3': (0., -3.3)}
    positioned, recipes = {}, {}
    for modality, ds in datasets.items():
        if ds.scope in ('Tonsil1', 'Tonsil1_gt600'):
            basis = 'local_uncorrected'
            if basis not in ds.coords:
                raise ValueError(f'{modality}: Tonsil1 linked PDF requires the saved local_uncorrected UMAP')
            xy = ds.coords[basis]
        else:
            basis, xy = primary_coord(ds)
        xy = np.asarray(xy, dtype=float)
        if xy.shape != (len(ds.obs), 2) or not np.isfinite(xy).all():
            raise ValueError(f'{modality}: invalid coordinates for the linked PDF')
        midpoint = (xy.min(axis=0) + xy.max(axis=0)) / 2
        scale = 2.2 / max(float(np.ptp(xy, axis=0).max()), 1e-12)
        positioned[modality] = (xy - midpoint) * scale + np.asarray(centers[modality])
        recipes[modality] = dict(basis=basis, center=list(centers[modality]), scale=scale,
            original_midpoint=midpoint.tolist(), n_cells=len(ds.obs), recommended_method=ds.recommended,
            selected_method=basis.split('_', 1)[1])
    return positioned, recipes

def linked_umap(datasets, folder, args, palette):
    from matplotlib.collections import LineCollection
    names = [m for m in MODALITIES if m in datasets]
    if 'RNA' not in datasets:
        raise ValueError('The linked UMAP requires RNA to connect physical cells to histone marks')
    coordinates, recipes = positioned_linked_coordinates(datasets)
    pairs = []
    for left, right in (('RNA', mark) for mark in names if mark != 'RNA'):
        shared = datasets[left].obs.index.intersection(datasets[right].obs.index, sort=False)
        if len(shared):
            pairs.append((left, right, shared))
    available = sum(len(ids) for _, _, ids in pairs)
    selected = np.sort(np.random.default_rng(args.seed).choice(available, min(available, args.max_links), replace=False)) if available else np.array([], dtype=int)
    segments, link_colors, links = [], [], []
    offset = 0
    for left, right, shared in pairs:
        indices = selected[(selected >= offset) & (selected < offset + len(shared))] - offset
        chosen = shared[indices]
        offset += len(shared)
        a = coordinates[left][datasets[left].obs.index.get_indexer(chosen)]
        b = coordinates[right][datasets[right].obs.index.get_indexer(chosen)]
        for cell, xy_a, xy_b in zip(chosen, a, b):
            label = datasets['RNA'].obs.loc[cell, LABEL] if cell in datasets['RNA'].obs.index else datasets[left].obs.loc[cell, LABEL]
            if label == UNLABELLED:
                label = datasets[right].obs.loc[cell, LABEL]
            segments.append([xy_a, xy_b])
            link_colors.append(colors.to_rgba(palette.get(label, '#BBBBBB'), .045))
            links.append(dict(merged_cell_id=cell, left=left, right=right, lymphoid_final=label))
    has_k9 = 'H3K9me3' in datasets
    fig, ax = plt.subplots(figsize=(11.0, 10.5 if has_k9 else 6.6))
    fig.subplots_adjust(left=.025, right=.79, bottom=.035, top=.90)
    if segments:
        ax.add_collection(LineCollection(segments, colors=link_colors, linewidths=.3, zorder=1,
            rasterized=args.rasterize_points))
    for modality in names:
        ds, xy = datasets[modality], coordinates[modality]
        order = np.random.default_rng(args.seed).permutation(len(xy))
        ax.scatter(xy[order, 0], xy[order, 1], c=[palette.get(v, '#BDBDBD') for v in ds.obs[LABEL].iloc[order]],
            s=args.linked_point_size, marker='o', edgecolors='white', linewidths=.08, alpha=.96,
            antialiased=True, zorder=2, rasterized=args.rasterize_points)
        center = recipes[modality]['center']
        ax.text(center[0], float(xy[:, 1].max()) + .22, f'{modality}  (n={len(ds.obs):,})',
            ha='center', va='bottom', fontsize=10, fontweight='bold', zorder=3)
    scope = folder.parent.name.replace('_gt', ' >')
    scope = 'All tonsils' if scope == 'merged' else scope
    basis_title = 'local uncorrected' if folder.parent.name in ('Tonsil1', 'Tonsil1_gt600') else 'recommended DNA; uncorrected RNA'
    ax.set_title(f'{scope} | {folder.name.replace("_", " ")} | {basis_title} | lymphoid_final\n'
        f'{len(links):,}/{available:,} exact-cell RNA-histone links drawn', fontsize=10, pad=22)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(-4.25, 4.25)
    ax.set_ylim(-4.65 if has_k9 else -1.3, 3.1)
    ax.axis('off')
    labels = sorted(set().union(*(set(ds.obs[LABEL]) for ds in datasets.values())), key=natural)
    ncols = max(1, (len(labels) + 24) // 25)
    handles = [Line2D([], [], marker='o', linestyle='', markersize=4, color=palette.get(v, '#BBBBBB'),
        markeredgecolor='white', markeredgewidth=.15, label=v) for v in labels]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.025, .5), frameon=False,
        ncol=ncols, fontsize=7, handletextpad=.4, title='lymphoid_final')
    save_publication_pdf(fig, folder / 'linked_umap/linked_umap_recommended', args)
    pd.DataFrame(links, columns=['merged_cell_id', 'left', 'right', LABEL]).to_csv(folder / 'tables/linked_umap_exact_cells.tsv', sep='\t', index=False)
    write_json(folder / 'tables/linked_umap_layout.json', dict(modalities=recipes, links_available=available,
        links_drawn=len(links), link_cap=args.max_links, pairs='RNA-histone only; exact physical cell IDs',
        point_style=dict(area_pt2=args.linked_point_size, marker='o', edgecolor='white', linewidth_pt=.08, alpha=.96),
        transform='translation and uniform scaling for display only; no rotation, alignment fitting or UMAP refit'))
    fig, ax = plt.subplots(figsize=(2 * ncols, max(1, .14 * min(len(labels), 25) + .4)))
    ax.axis('off')
    fig.legend(handles=handles, loc='center', frameon=False, ncol=ncols, title='lymphoid_final')
    save_publication_pdf(fig, folder / 'linked_umap/lymphoid_final_legend', args)

COUNT_METRICS = {'umi', 'genes', 'fragments'}
LABELS = {'umi': 'UMIs per cell', 'genes': 'Detected genes per cell', 'fragments': 'Unique fragments per cell', 'frip': 'Fragments in peaks (%)'}
GATE_COLOR = '#546D59'
CUTOFF_COLOR = '#A16A32'


def guide_lines(samples, cutoffs, modality, metric, high_coverage):
    gate_visible = high_coverage and metric in {'umi','fragments'}
    lines = []
    for i,sample in enumerate(samples):
        for value in sorted(set(cutoffs[sample][modality].get(metric,[]))):
            if value > 0 and (not gate_visible or value >= HIGH_COVERAGE):
                lines.append((value,i-.40,i+.40,'QC cutoff'))
    if gate_visible:
        lines.append((float(HIGH_COVERAGE),-.43,len(samples)-.57,f'>{HIGH_COVERAGE} gate'))
    return lines


def number_label(value):
    if abs(value) < 1e-10:
        return "0"
    return f"{value:,.6f}".rstrip("0").rstrip(".")

def count_ticks(low, high, max_ticks):
    """Sparse original-unit 1/2/5 ticks, approximately even in log space."""
    high = max(high, 1.0)
    positive_low = max(low, min(1.0, high))
    first_power = int(math.floor(math.log10(positive_low))) - 1
    last_power = int(math.ceil(math.log10(high))) + 1
    pool = np.unique([m * 10.0**e for e in range(first_power, last_power + 1)
                      for m in (1, 2, 5)])
    lower = 0.0 if low <= 0 else float(pool[pool <= low][-1])
    upper = float(pool[pool >= high][0])
    candidates = pool[(pool >= lower) & (pool <= upper)]
    if lower == 0:
        candidates = np.r_[0.0, candidates[candidates >= 1]]
    if len(candidates) <= max_ticks:
        return candidates
    positions = np.log10(1 + candidates)
    selected = {0, len(candidates) - 1}
    for target in np.linspace(positions[0], positions[-1], max_ticks)[1:-1]:
        selected.add(int(np.argmin(abs(positions - target))))
    return candidates[sorted(selected)]

def axis_style(ax, values, guides, metric):
    is_count = metric in COUNT_METRICS
    use_log = is_count
    transform = (lambda x: np.log10(1.0 + np.asarray(x))) if use_log else np.asarray
    raw_range = np.r_[values, [g[0] for g in guides]]
    low, high = float(raw_range.min()), float(raw_range.max())
    if is_count and low < 0:
        raise ValueError(f"{metric}: negative counts found; verify saved metric units")
    if low == high:
        if high == 0:
            high = 1.0
        else:
            margin = max(abs(high) * 0.15, 0.1)
            low, high = low - margin, high + margin
    if use_log:
        ticks = count_ticks(low, high, 5)
    else:
        if low >= 0:
            low = 0.0
        nbins = 5 - 1
        while True:
            ticks = MaxNLocator(nbins=nbins, steps=[1, 2, 2.5, 5, 10]).tick_values(low, high)
            if metric == "frip" and 0 <= low <= high <= 100:
                ticks = ticks[(ticks >= 0) & (ticks <= 100)]
                if ticks[-1] < high:
                    ticks = np.r_[ticks, 100.0]
            if len(ticks) <= 5 or nbins <= 1:
                break
            nbins -= 1
    unit = 1000.0 if is_count and high >= 1000 else 1.0
    positions = transform(ticks)
    tick_labels = {float(pos): number_label(raw / unit) for pos, raw in zip(positions, ticks)}
    ax.yaxis.set_major_locator(FixedLocator(positions))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: tick_labels.get(float(y), "")))
    ax.yaxis.set_minor_locator(NullLocator())
    lower, upper = float(positions[0]), float(positions[-1])
    pad = max(upper - lower, 0.1) * 0.035
    ax.set_ylim(lower - pad * 0.2, upper + pad)
    ylabel = LABELS.get(metric, metric.replace("_", " ").capitalize())
    if unit == 1000:
        ylabel += r" ($\times$1,000)"
    ax.set_ylabel(ylabel, labelpad=9)
    return transform


def draw_violin(ax, values, x, color, transform):
    transformed = np.asarray(transform(values), dtype=float)
    if len(transformed) > 1 and np.ptp(transformed) > 1e-10:
        grid = np.linspace(transformed.min(), transformed.max(), 240)
        density = gaussian_kde(transformed, bw_method=.25)(grid)
        width = density / density.max() * 0.40
        ax.fill_betweenx(grid, x - width, x + width, color=color,
                         alpha=0.92, linewidth=0, zorder=2)
    else:
        ax.plot([x-0.18, x+0.18], [transformed[0]] * 2,
                color=color, linewidth=3, solid_capstyle="round", zorder=2)
    q25, median, q75 = np.quantile(values, [0.25, 0.5, 0.75])
    ax.plot([x, x], transform([q25, q75]), color="#333333", linewidth=1.05, zorder=4)
    ax.scatter([x], transform([median]), s=21, facecolor="white",
               edgecolor="#333333", linewidth=0.7, zorder=5)

def report_title(report):
    cohort, arm = report.parts[-2:]
    cohort = {"merged": "Merged tonsils", "Tonsil1_gt600": "Tonsil1 >600"}.get(cohort, cohort)
    arm = {"labelled_only": "labelled only", "full": "full"}[arm]
    return f"{cohort} | {arm}"

def qc_violin(frame, metric, modality, report, out, cutoffs):
    groups = sorted(frame['sample'].unique(), key=natural)
    guides = guide_lines(groups,cutoffs,modality,metric,report.parent.name=='Tonsil1_gt600')
    fig, ax = plt.subplots(figsize=(max(4.4,1.65+len(groups)*.85),3.75))
    transform = axis_style(ax,frame['value'].to_numpy(),guides,metric)
    for i,group in enumerate(groups):
        values = frame.loc[frame['sample']==group,'value'].to_numpy()
        draw_violin(ax,values,i,SAMPLE_COLORS.get(group,'#7C93AB'),transform)
    legend = {}
    for value, left, right, label in guides:
        color = GATE_COLOR if label == ">600 gate" else CUTOFF_COLOR
        ax.hlines(float(transform([value])[0]), left, right, color=color,
                  linewidth=0.85, linestyle=(0, (3, 3)), alpha=0.85, zorder=3)
        legend[label] = Line2D([0], [0], color=color, linewidth=0.9,
                              linestyle=(0, (3, 3)), label=label)
    ax.set_xlim(-0.65, len(groups)-0.35)
    ax.set_xticks(range(len(groups)),groups,ha="center",rotation_mode="anchor")
    ax.set_xlabel("Sample",labelpad=10)
    ax.set_title(f"{modality} | {report_title(report)}", loc="left", fontsize=11, pad=27 if legend else 12)
    if legend:
        ax.legend(handles=list(legend.values()), loc="lower left", bbox_to_anchor=(0, 1.005),
                  borderaxespad=0, frameon=False, ncol=len(legend), fontsize=8,
                  handlelength=2.1, columnspacing=1.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.8)
        ax.spines[side].set_color("#404040")
    ax.tick_params(axis="both", which="major", labelsize=9, width=0.7, length=3.5,
                   direction="out", color="#404040")
    ax.tick_params(axis="x", length=0, pad=8)
    ax.grid(False)
    fig.tight_layout(pad=1.1)
    name = f"{modality}_{metric}_by_sample.pdf"
    path = out / name
    fig.savefig(path, bbox_inches="tight", facecolor="white", metadata={"Creator": "Tonsil QC"})
    plt.close(fig)



def agreement(datasets, folder, args):
    from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
    dest=folder/'tables/cluster_agreement';dest.mkdir(parents=True,exist_ok=True)
    tri=sorted(set.intersection(*(set(datasets[m].obs.index) for m in T1_MODALITIES)))
    pd.Series(tri,name='cell_id').to_csv(folder/'tables/agreement_exact_cells.tsv.gz',sep='\t',index=False)
    reference=label_values(datasets['RNA'].obs[LABEL].reindex(tri))
    rows=[]; mappings=[]; coverage_rows=[]
    for modality in T1_MODALITIES:
        ds=datasets[modality]
        method=ds.recommended if ds.scope=='merged' else 'uncorrected'
        for resolution in (1,1.5,2,3):
            key=f'leiden_{method}_r{resolution:g}'
            if key not in ds.obs:raise KeyError(f'{modality}: required graph assignment {key} missing')
            raw=ds.obs[key].astype('string')
            axis=sorted(raw.dropna().unique(),key=natural)
            numbering={v:str(i) for i,v in enumerate(axis)}
            for v in axis:mappings.append(dict(modality=modality,cluster_key=key,original_cluster_id=v,publication_cluster_id=numbering[v],scope_cells=int(raw.eq(v).sum())))
            clustered=raw.reindex(tri).map(numbering)
            valid=clustered.notna()
            labels=reference.loc[valid]
            counts=pd.crosstab(labels,clustered.loc[valid]).reindex(columns=list(numbering.values()),fill_value=0)
            counts=counts.reindex(sorted(counts.index,key=lambda v:(v==UNLABELLED,natural(v))))
            labelled=valid & reference.ne(UNLABELLED)
            ari=adjusted_rand_score(reference.loc[labelled],clustered.loc[labelled]) if labelled.sum()>=2 else float('nan')
            nmi=normalized_mutual_info_score(reference.loc[labelled],clustered.loc[labelled]) if labelled.sum()>=2 else float('nan')
            assert counts.to_numpy().sum()==valid.sum()
            for original,num in numbering.items():
                coverage_rows.append(dict(modality=modality,cluster_key=key,original_cluster_id=original,cluster=num,total_scope_cells=int(raw.eq(original).sum()),tri_profiled_cells=int(clustered.eq(num).sum()),labelled_cells=int((clustered.eq(num)&labelled).sum())))
            stem=f'{modality}_{key}_lymphoid_final_all_clusters'
            counts.to_csv(dest/(stem+'_counts.tsv'),sep='\t')
            fractions=counts.div(counts.sum(axis=1).replace(0,np.nan),axis=0).fillna(0)
            fractions.to_csv(dest/(stem+'_row_proportions.tsv'),sep='\t')
            assert np.allclose(fractions.sum(axis=1),1)
            comparison_heatmap(counts,f'{ds.scope} | {folder.name} | exact tri-profiled cells',modality,key,ari,nmi,folder/'cluster_agreement'/(stem+'_row_proportions'),args,score_n=int(labelled.sum()))
            rows.append(dict(modality=modality,cluster_key=key,method=method,resolution=resolution,n_tri_profiled=len(tri),n_matched=len(tri),n_plotted=int(valid.sum()),n_scored=int(labelled.sum()),n_unlabelled=int((valid&~labelled).sum()),missing_cluster=int((~valid).sum()),ARI=ari,NMI=nmi,cell_ids_sha256=ids_hash(tri),label_source='exact retained RNA lymphoid_final'))
    pd.DataFrame(rows).to_csv(folder/'tables/cluster_agreement_summary.tsv',sep='\t',index=False)
    pd.DataFrame(mappings).to_csv(folder/'tables/cluster_id_mapping.tsv',sep='\t',index=False)
    pd.DataFrame(coverage_rows).to_csv(folder/'tables/cluster_agreement_coverage.tsv',sep='\t',index=False)
    return rows


def render_scope(objects, scope, population, config, publication):
    folder=publication/scope/population
    folder.mkdir(parents=True, exist_ok=True)
    (folder/'tables').mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':6.5,'axes.linewidth':.6,'lines.linewidth':.8,'pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'white','figure.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,'pdf.compression':6})
    palette=config['palette']
    args=SimpleNamespace(seed=0,dpi=600,max_links=10000,linked_point_size=3.2,rasterize_points=False)
    datasets={}
    for m,data in objects.items():
        obs=data.obs.copy()
        obs[LABEL]=label_values(obs[LABEL])
        ds=SimpleNamespace(modality=m,scope=scope,obs=obs,recommended=data.uns['integration']['recommended_method'],coords={})
        if scope=='merged':
            ds.coords['cohort_'+ds.recommended]=data.obsm['X_umap'].copy()
        else:
            ds.coords['local_uncorrected']=data.obsm['X_umap'].copy()
        obs.to_csv(folder/f'tables/{m}_cell_metadata.tsv.gz',sep='\t')
        datasets[m]=ds
    agreement(datasets,folder,args)
    composition(datasets['RNA'],folder,args,palette)
    from matplotlib_venn import venn3
    fig,ax=plt.subplots(figsize=(3.7,3.4))
    venn3([set(datasets[m].obs.index) for m in T1_MODALITIES],set_labels=T1_MODALITIES,ax=ax,set_colors=('#0072B2','#D55E00','#009E73'),alpha=.5)
    ax.set_title('Retained cells · exact physical IDs')
    save_publication_pdf(fig,folder/'overlap/venn_RNA_H3K27ac_H3K27me3',args)
    out=folder/'qc_violins_readable';out.mkdir(exist_ok=True)
    cutoffs=qc_cutoffs(config)
    for m,metrics in {'RNA':['genes','umi'],'H3K27ac':['fragments','frip'],'H3K27me3':['fragments','frip']}.items():
        for metric in metrics:
            column={'genes':'n_genes_by_counts','umi':'total_counts','fragments':'n_fragment','frip':'frip'}[metric]
            frame=pd.DataFrame(dict(sample=datasets[m].obs['sample'].astype(str),value=pd.to_numeric(datasets[m].obs[column])))
            if metric=='frip':frame['value']*=100
            with plt.rc_context({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'ps.fonttype':42,'axes.labelsize':11,'axes.unicode_minus':False,'figure.facecolor':'white','savefig.facecolor':'white'}):
                qc_violin(frame,metric,m,folder,out,cutoffs)
    linked_umap(datasets,folder,args,palette)


def qc_cutoffs(config):
    """Plot the configured thresholds with the same units as the violins."""
    names = {'umi':('min_umi','max_umi'),'genes':('max_genes',),'mt':('max_mt',),
             'hb':('max_hb',),'ribo':('max_rb',),'top50':('max_top50',),
             'fragments':('min_fragments','max_fragments')}
    return {s['name']:{m:{metric:[limits[k] for k in keys if limits.get(k) is not None]
                       for metric,keys in names.items()}
                       for m,limits in s['quality'].items()} for s in config['samples']}


# Stage 02: lymphoid labels and final profiles

SCANVI_LABEL_KEY = "labels"
BATCH_KEY = "dataset"
UNLABELED = "Unknown"
COUNTS_LAYER = "counts"
RANDOM_SEED = 0
N_HVG = 2_000
N_LATENT = 20
TRAIN_BATCH_SIZE = 128
TRAIN_DATASPLITTER_KWARGS = {"drop_last": True}
GENE_FALLBACK_TRIGGER = 1_000
MIN_SHARED_GENES = 500
MERGE_RULES = {
    "B_mem": ["B_mem", "B_IFN"],
    "B_naive": [
        "B_activated,4",
        "B_activated,3",
        "B_activated,0",
        "B_activated,2",
        "B_activated,1",
        "B_naive",
    ],
    "T_CD4+": ["T_CD4+", "T_TfR", "T_Treg"],
    "B_GC_LZ": ["B_GC_prePB"],
    "T_CD8+": ["T_CD8+_CD161+", "T_CD8+_cytotoxic", "T_CD8+_naive"],
    "FDC": ["VSMC", "Endo"],
    "Macrophages": ["Macrophages_M2", "Macrophages_M1"],
    "DC": ["DC_cDC1", "DC_cDC2", "DC_pDC", "DC_CCR7+"],
}


def matrix_is_counts(matrix) -> bool:
    values = matrix.data if sp.issparse(matrix) else np.asarray(matrix).ravel()
    if values.size == 0:
        return False
    return all(np.isfinite(part).all() and np.all(part >= 0)
               and np.allclose(part, np.rint(part), rtol=0, atol=1e-6)
               for part in (values[i:i+1_000_000] for i in range(0,values.size,1_000_000)))


def raw_count_anndata(source: ad.AnnData, name: str) -> ad.AnnData:
    """Return a compact AnnData whose X is verified raw counts."""
    candidates = []
    if COUNTS_LAYER in source.layers:
        candidates.append(("layers['counts']", source.layers[COUNTS_LAYER], source.var))
    candidates.append(("X", source.X, source.var))
    if source.raw is not None:
        candidates.append(("raw.X", source.raw.X, source.raw.var))

    for location, matrix, var in candidates:
        if not matrix_is_counts(matrix):
            continue
        counts = matrix.tocsr() if sp.issparse(matrix) else sp.csr_matrix(matrix)
        result = ad.AnnData(
            X=counts.astype(np.float32, copy=False),
            obs=source.obs.copy(),
            var=var.copy(),
        )
        result.obs_names = pd.Index(source.obs_names.astype(str), name="cell_id")
        result.var_names = pd.Index(var.index.astype(str), name="gene")
        if result.obs_names.has_duplicates:
            raise ValueError(f"{name} has duplicate cell identifiers")
        totals = np.asarray(result.X.sum(axis=1)).ravel()
        if np.any(totals <= 0):
            raise ValueError(
                f"{name} has {(totals <= 0).sum():,} zero-count cells; refusing "
                "to change the query/reference cell set silently"
            )
        print(f"{name}: raw counts from {location}; shape={result.shape}")
        return result

    raise ValueError(
        f"{name}: verified raw UMI counts are required; none were found in "
        "layers['counts'], X, or raw.X"
    )


def clean_gene_name(value: str) -> str:
    """Canonical feature name for the secondary gene match."""
    value = str(value).strip().upper()
    return re.sub(r"\.\d+$", "", value)


def collapse_duplicate_genes(adata: ad.AnnData, name: str) -> ad.AnnData:
    """Sum same-named feature columns so AnnData can be subset safely."""
    if not adata.var_names.has_duplicates:
        return adata
    names = pd.Index(adata.var_names.astype(str))
    codes, unique = pd.factorize(names, sort=False)
    projection = sp.csr_matrix(
        (
            np.ones(len(codes), dtype=np.float32),
            (np.arange(len(codes)), codes),
        ),
        shape=(len(codes), len(unique)),
    )
    collapsed = (sp.csr_matrix(adata.X) @ projection).tocsr()
    result = ad.AnnData(
        X=collapsed,
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=pd.Index(unique.astype(str), name="gene")),
    )
    print(
        f"{name}: summed {adata.n_vars - result.n_vars:,} duplicate gene columns"
    )
    return result


def collapse_clean_genes(adata: ad.AnnData, name: str) -> ad.AnnData:
    """Apply cleaned names and sum duplicate columns instead of dropping them."""
    cleaned = pd.Index([clean_gene_name(gene) for gene in adata.var_names])
    valid = cleaned.notna() & (cleaned != "") & (cleaned != "NAN")
    cleaned = cleaned[valid]
    matrix = adata.X[:, np.asarray(valid)]

    codes, unique = pd.factorize(cleaned, sort=False)
    projection = sp.csr_matrix(
        (
            np.ones(len(codes), dtype=np.float32),
            (np.arange(len(codes)), codes),
        ),
        shape=(len(codes), len(unique)),
    )
    collapsed = (sp.csr_matrix(matrix) @ projection).tocsr()
    result = ad.AnnData(
        X=collapsed,
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=pd.Index(unique.astype(str), name="gene")),
    )
    print(
        f"{name}: cleaned {adata.n_vars:,} genes to {result.n_vars:,} unique genes "
        f"({adata.n_vars - result.n_vars:,} duplicate/invalid columns collapsed)"
    )
    return result


def harmonize_genes(
    reference: ad.AnnData,
    query: ad.AnnData,
    atlas_name: str,
) -> tuple[ad.AnnData, ad.AnnData, str]:
    """Match shared features by exact name, then normalized names if needed."""
    reference = collapse_duplicate_genes(reference, f"{atlas_name} reference")
    query = collapse_duplicate_genes(query, "query")
    exact = reference.var_names.intersection(query.var_names, sort=False)
    method = "exact"
    if len(exact) < GENE_FALLBACK_TRIGGER:
        print(
            f"{atlas_name}: {len(exact):,} exact gene matches; matching normalized names"
        )
        reference = collapse_clean_genes(reference, f"{atlas_name} reference")
        query = collapse_clean_genes(query, "query")
        shared = reference.var_names.intersection(query.var_names, sort=False)
        method = "uppercase_and_strip_version"
    else:
        shared = exact

    if len(shared) < MIN_SHARED_GENES:
        raise ValueError(
            f"{atlas_name}: only {len(shared):,} genes overlap after {method} "
            f"matching; expected at least {MIN_SHARED_GENES:,}"
        )
    reference = reference[:, shared].copy()
    query = query[:, shared].copy()
    if not reference.var_names.equals(query.var_names):
        raise RuntimeError(f"{atlas_name}: reference/query gene order differs")
    print(f"{atlas_name}: {len(shared):,} shared genes ({method} matching)")
    return reference, query, method


def probability_frame(values, index: pd.Index) -> pd.DataFrame:
    if isinstance(values, pd.DataFrame):
        result = values.copy()
        result.index = index
        return result
    array = np.asarray(values)
    if array.ndim != 2:
        raise TypeError("scANVI predict(soft=True) did not return a 2-D table")
    raise TypeError(
        "scANVI predict(soft=True) returned an unlabeled array, so label names "
        "cannot be assigned safely. Upgrade/verify scvi-tools 1.5.1."
    )


def graph_umap(data, rep, method, settings, modality):
    spec = settings["rna"] if modality == "RNA" else settings["dna"]
    key = "neighbors_" + method
    sc.pp.neighbors(data, use_rep=rep, n_neighbors=min(spec["n_neighbors"], data.n_obs-1),
                    metric="euclidean", random_state=settings["seed"], key_added=key)
    sc.tl.umap(data, neighbors_key=key, random_state=settings["seed"], min_dist=spec["min_dist"])
    data.obsm["X_umap_" + method] = np.asarray(data.obsm["X_umap"], dtype=np.float32).copy()
    for resolution in settings["leiden_resolutions"]:
        sc.tl.leiden(data, neighbors_key=key, resolution=resolution, random_state=settings["seed"],
                     flavor="leidenalg", key_added=f"leiden_{method}_r{resolution:g}")


def embed_rna(data, settings):
    spec = settings["rna"]
    # Full-gene raw counts stay in data.X. Only the temporary HVG object is scaled.
    work = data.copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    sc.pp.highly_variable_genes(work, n_top_genes=min(spec["n_hvg"], work.n_vars),
                              flavor="seurat", batch_key="sample")
    hvg = work.var.highly_variable.to_numpy()
    if hvg.sum() < 3:
        raise ValueError("Too few variable RNA genes")
    data.var["highly_variable"] = hvg
    data.var["highly_variable_nbatches"] = work.var["highly_variable_nbatches"].to_numpy()
    work = work[:, hvg].copy()
    sc.pp.scale(work, max_value=10)
    n_components = min(spec["n_components"], work.n_obs-1, work.n_vars-1)
    sc.tl.pca(work, n_comps=n_components, random_state=settings["seed"], svd_solver="arpack")
    data.obsm["X_pca"] = work.obsm["X_pca"].astype(np.float32)
    data.varm["PCs"] = np.zeros((data.n_vars, n_components), dtype=np.float32)
    data.varm["PCs"][hvg] = work.varm["PCs"]
    data.uns["pca"] = work.uns["pca"]
    del work
    graph_umap(data, "X_pca", "uncorrected", settings, "RNA")
    data.uns["integration"] = dict(recommended_method="uncorrected", batch_corrected=False,
                                    label_smoothing="Selected per-donor labels are retained; no new cohort smoothing",
                                    normalization="Temporary log1p(counts/total*10000), sample-aware HVG selection, scaling, PCA; X remains raw counts")
    choose_default(data, "uncorrected")
    return data


def embed_spectral(data, settings, threads, modality):
    limit = settings["dna"]["features"][modality]
    kwargs = dict(n_features=min(limit, data.n_vars), max_iter=1)
    if "n_jobs" in inspect.signature(snap.pp.select_features).parameters:
        kwargs["n_jobs"] = threads
    snap.pp.select_features(data, **kwargs)
    chromosomes = data.var_names.to_series().str.split(":").str[0]
    excluded = chromosomes.isin(settings["dna"].get("exclude_chromosomes", ["chrM", "MT", "M"])).to_numpy()
    data.var.loc[excluded, "selected"] = False
    selected = np.asarray(data.var["selected"], dtype=bool)
    n_components = min(settings["dna"]["n_components"], data.n_obs - 1, int(selected.sum()) - 1)
    if n_components < 2:
        raise ValueError(f"{modality}: too few informative windows")
    snap.tl.spectral(data, n_comps=n_components, features="selected", random_state=settings["seed"],
                     distance_metric="cosine", num_threads=threads)
    data.obsm["X_spectral_uncorrected"] = np.asarray(data.obsm["X_spectral"], dtype=np.float32).copy()
    graph_umap(data, "X_spectral_uncorrected", "uncorrected", settings, modality)
    return data


def harmonize_dna(data, settings, threads):
    key = 'X_spectral_harmony'
    options = dict(settings['dna']['harmony'])
    import harmonypy
    options["random_state"] = settings["seed"]
    # Harmony's K-means default is invalid on tiny subsets; this also limits
    # the number of centroids in a reproducible, explicit way on real runs.
    options.setdefault("nclust", min(100, max(2, round(data.n_obs / 30))))
    # harmonypy 2.1 owns its thread pool; OMP_NUM_THREADS does not limit it.
    if "ncores" in inspect.signature(harmonypy.run_harmony).parameters:
        options.setdefault("ncores", threads)
    source = np.asarray(data.obsm["X_spectral_uncorrected"], dtype=np.float64)
    batches = data.obs[[settings["batch_key"]]].copy()
    if batches.isna().any().any():
        raise ValueError("Harmony batch labels contain missing values")
    # Call the backend directly: released SnapATAC2 wrappers can transpose
    # modern cells-by-components Z_corr output before AnnData validates it.
    output = harmonypy.run_harmony(source.copy(), batches,
                                  [settings["batch_key"]], **options)
    value = np.asarray(output.Z_corr)
    if value.shape != source.shape and value.ndim == 2 and value.T.shape == source.shape:
        value = value.T
    if value.shape != source.shape or not np.isfinite(value).all():
        raise ValueError(f"Harmony returned {value.shape}; expected finite cells-by-components {source.shape}")
    data.obsm[key] = value.astype(np.float32)
    data.uns["harmony_backend_json"] = json.dumps(dict(
        entrypoint="harmonypy.run_harmony", version=importlib.metadata.version("harmonypy"),
        output_orientation="cells_by_components", parameters=options), sort_keys=True)

def choose_default(data, method):
    data.obsm["X_umap"] = data.obsm["X_umap_" + method].copy()
    if "X_spectral_" + method in data.obsm:
        data.obsm["X_integrated"] = data.obsm["X_spectral_" + method].copy()
    key = "neighbors_" + method
    data.obsp["connectivities"] = data.obsp[key + "_connectivities"].copy()
    data.obsp["distances"] = data.obsp[key + "_distances"].copy()
    data.uns["neighbors"] = dict(data.uns[key], connectivities_key="connectivities", distances_key="distances")


def peak_paths(state, args, mark):
    name = f"{args.sample}_{mark}"
    directory = state["output"] / "peaks" / mark
    broad = state["config"]["dna_processing"]["peak_calling"][mark]["broad"]
    return dict(directory=directory, fragments=state["output"] / "fragments" / f"{name}.final_pass.fragments.tsv.gz",
        raw=directory / f"{name}_peaks.{'broadPeak' if broad else 'narrowPeak'}",
        filtered=directory / f"{name}_peaks.blacklist_filtered.bed",
        report=directory / f"{name}_peaks.xls", auxiliary=directory / (f"{name}_peaks.gappedPeak" if broad else f"{name}_summits.bed"))


def peak_command(state, args, mark, paths):
    config = state["config"]["dna_processing"]
    params = config["peak_calling"][mark]
    command = ["macs3", "callpeak", "-t", str(paths["fragments"]), "-f", "FRAG", "-g", config["genome_size"],
               "-n", f"{args.sample}_{mark}", "--outdir", str(paths["directory"]), "-q", str(params["qvalue"]),
               "--max-count", "1", "--nolambda", "--nomodel"]
    if params["broad"]:
        command.extend(["--broad", "--broad-cutoff", str(params["broad_cutoff"])])
    return command


def finalize_dna_fragments(state, args, mark, snap):
    destination = state["output"] / "objects" / f"{args.sample}_{mark}_analysis_pass.h5ad"
    paths = peak_paths(state, args, mark)
    phase = mark + "_fragments_and_peaks"
    outputs = [destination, paths["fragments"], paths["raw"], paths["filtered"], paths["report"], paths["auxiliary"]]
    if phase_complete(state, phase, outputs):
        print(f"Reusing completed {mark} final fragments/new peaks", flush=True)
        return destination, paths["filtered"]
    destination.parent.mkdir(parents=True, exist_ok=True)
    paths["directory"].mkdir(parents=True, exist_ok=True)
    paths["fragments"].parent.mkdir(parents=True, exist_ok=True)
    for path in outputs:
        path.unlink(missing_ok=True)
    source = snap.read(str(state["dna_sources"][mark]), backed="r")
    selected = None
    try:
        source_ids = snap_obs(source).index
        expected = state["retained"][mark].index
        positions = np.flatnonzero(source_ids.isin(expected))
        if len(positions) != len(expected):
            raise ValueError(f"{mark}: source cell membership changed")
        selected = source.subset(positions, out=str(destination), inplace=False)
        obs = state["retained"][mark].loc[snap_obs(selected).index].copy()
        obs["analysis_pass"] = True
        obs["final_qc_pass"] = True
        rna = state['retained']['RNA']
        for key in state['label_columns']:
            obs[key] = rna[key].reindex(obs.index)
        obs['direct_rna_match'] = obs.index.isin(rna.index)
        set_snap_obs(selected, obs)
        selected.uns["final_qc_provenance"] = dict(branch_signature=state["signature"], source=str(state["dna_sources"][mark]),
                                                  labels="exact matches to retained RNA")
        exported = snap.ex.export_fragments(selected, groupby=["cells"] * selected.n_obs, selections=["cells"],
                                           ids=list(obs.index), out_dir=str(paths["fragments"].parent),
                                           prefix=f"{args.sample}_{mark}.", suffix=".final_pass.fragments.tsv.gz")
        generated = Path(exported["cells"])
        if generated.resolve() != paths["fragments"].resolve():
            if not generated.resolve().is_relative_to(state["output"]):
                raise ValueError("Fragment export escaped the output folder")
            generated.replace(paths["fragments"])
        command = peak_command(state, args, mark, paths)
        print("+ " + " ".join(command), flush=True)
        subprocess.run(command, check=True)
        with paths["filtered"].open("w") as handle:
            subprocess.run(["bedtools", "intersect", "-v", "-a", str(paths["raw"]), "-b", str(state["blacklist"])],
                           stdout=handle, check=True)
        require_files(*outputs)
        if paths["filtered"].stat().st_size == 0:
            raise RuntimeError(f"No {mark} peaks remain after blacklist filtering")
        snap.metrics.frip(selected, {"frip": str(paths["filtered"])})
        final_obs = snap_obs(selected)
        if not np.isfinite(pd.to_numeric(final_obs["frip"], errors="raise")).all():
            raise ValueError(f"Missing/nonfinite new FRiP for {mark}")
    finally:
        if selected is not None:
            selected.close()
        source.close()
    finish_phase(state, phase, outputs, n_cells=len(state["retained"][mark]), peak_command=command)
    return destination, paths["filtered"]


def count_dna(state, args, mark, representation, native, peaks):
    size = state['config']['dna_processing']['tile_size'][mark]
    folder = state['output'] / 'objects'
    path = folder / f'{args.sample}_{mark}_{representation}.h5ad'
    phase = f'{mark}_{representation}'
    if phase_complete(state,phase,[path]):
        return path
    path.unlink(missing_ok=True)
    source = snap.read(str(native),backed='r')
    data = None
    try:
        obs = snap_obs(source)
        if representation == 'tiles':
            data = snap.pp.add_tile_matrix(source,bin_size=size,inplace=False,file=str(path))
        else:
            data = snap.pp.make_peak_matrix(source,peak_file=str(peaks),inplace=False,file=str(path))
        target = pd.Index(list(map(str,data.obs_names)))
        if not target.is_unique or set(target) != set(obs.index):
            raise ValueError(f'{mark}: feature counting changed cell membership')
        set_snap_obs(data,obs.loc[target])
        if representation == 'peaks':
            snap.pp.select_features(data,n_features=min(100000,data.n_vars),max_iter=1)
            n = min(30,data.n_obs-1,int(np.asarray(data.var['selected']).sum())-1)
            if n < 2:
                raise ValueError(f'{mark}: too few selected peaks or cells')
            snap.tl.spectral(data,n_comps=n,features='selected',random_state=0,
                             distance_metric='cosine',num_threads=args.threads)
            snap.tl.umap(data,use_rep='X_spectral',random_state=0)
    finally:
        if data is not None:
            data.close()
        source.close()
    finish_phase(state,phase,[path])
    return path


def phase_file(state,phase):
    return state['output']/'completed'/f'{phase}.json'


def phase_complete(state,phase,paths):
    marker = phase_file(state,phase)
    if not marker.exists():
        return False
    old = read_json(marker)
    return old.get('signature')==state['signature'] and old.get('files')==[identity(p) for p in paths if Path(p).is_file()] and all(Path(p).is_file() for p in paths)


def finish_phase(state,phase,paths,**details):
    write_json(phase_file(state,phase),dict(signature=state['signature'],files=[identity(p) for p in paths],**details))


def stable_features(data,modality):
    if modality=='RNA':
        key=next((x for x in ('gene_ids','gene_id','feature_id') if x in data.var),None)
        if key is None:raise ValueError('Stable RNA gene IDs are required')
        data.var['gene_symbol']=data.var_names.astype(str)
        data.var['gene_id']=data.var[key].astype(str)
        data.var_names=pd.Index(data.var.gene_id,name='gene_id_index')
    else:
        positions=data.var_names.to_series().str.extract(r'^(.*):(\d+)-(\d+)$')
        size=1000 if modality=='H3K27ac' else 10000
        if positions.isna().any().any():raise ValueError('Invalid common tile IDs')
        lo,hi=positions[1].astype(int),positions[2].astype(int)
        if ((lo%size!=0)|(hi<=lo)|(hi-lo>size)).any():raise ValueError('Incompatible common tile grid')
    if not data.var_names.is_unique:raise ValueError('Duplicate stable features')
    return data


def merge_counts(parts,modality):
    parts=[stable_features(x.copy(),modality) for x in parts]
    first=parts[0].var_names
    if any(set(x.var_names)!=set(first) for x in parts[1:]):
        raise ValueError('Feature universes differ; missing measurements cannot be filled with zero')
    data=ad.concat([x[:,first] for x in parts],join='inner',merge='first',index_unique=None)
    data.X=sp.csr_matrix(data.X)
    data.obs=exact_obs(data.obs)
    data.obs[LABEL]=label_values(data.obs[LABEL])
    data.obs['merged_total_counts']=np.asarray(data.X.sum(axis=1)).ravel()
    data.obs['merged_detected_features']=np.diff(data.X.indptr)
    # Raw single-donor graph blocks do not define a merged neighbor graph.
    data.obsp.clear();data.obsm.clear();data.uns.clear()
    integer_counts(data.X)
    return data


def fit_merged(data,modality,settings):
    if modality=='RNA':
        return embed_rna(data,settings)
    embed_spectral(data,settings,THREADS,modality)
    method='uncorrected'
    levels=data.obs.tagmentation_protocol.nunique()
    if modality in ('H3K27ac','H3K27me3') and levels==2:
        harmonize_dna(data,settings,THREADS)
        graph_umap(data,'X_spectral_harmony','harmony',settings,modality)
        method='harmony'
    data.uns['integration']=dict(recommended_method=method,batch_key='tagmentation_protocol',batch_corrected=method=='harmony')
    choose_default(data,method)
    return data


def fit_local(data,modality,settings):
    # Local UMAPs use the uncorrected RNA PCA or DNA spectral space.
    basis='X_pca' if modality=='RNA' else 'X_spectral_uncorrected'
    local=ad.AnnData(X=sp.csr_matrix((data.n_obs,0)),obs=data.obs.copy())
    local.obsm['X_saved']=data.obsm[basis].copy()
    sc.pp.neighbors(local,use_rep='X_saved',n_neighbors=min(30,local.n_obs-1),metric='euclidean' if modality=='RNA' else 'cosine',random_state=0)
    sc.tl.umap(local,random_state=0,min_dist=.3)
    data.obsm['X_umap_uncorrected']=local.obsm['X_umap'].copy()
    data.obsm['X_umap']=local.obsm['X_umap'].copy()
    data.obsp.clear()
    for key in list(data.uns):
        if key=='neighbors' or key.startswith(('neighbors_','harmony_')):del data.uns[key]
    for key in ('connectivities','distances'):
        data.obsp[key]=local.obsp[key].copy()
        data.obsp['neighbors_uncorrected_'+key]=local.obsp[key].copy()
    data.uns['neighbors']=local.uns['neighbors'].copy()
    data.uns['neighbors_uncorrected']=dict(local.uns['neighbors'],connectivities_key='neighbors_uncorrected_connectivities',distances_key='neighbors_uncorrected_distances')
    for resolution in (1,1.5,2,3):
        key=f'leiden_uncorrected_r{resolution:g}'
        sc.tl.leiden(local,resolution=resolution,random_state=0,flavor='leidenalg',key_added=key)
        data.obs[key]=local.obs[key].copy()
    for key in list(data.obsm):
        if 'harmony' in key or key=='X_integrated':del data.obsm[key]
    for key in list(data.obs):
        if key.startswith('leiden_harmony'):del data.obs[key]
    data.uns['integration']=dict(recommended_method='uncorrected',batch_corrected=False,source_basis=basis)
    return data


def filter_reference(reference, label_column, excluded):
    if label_column not in reference.obs:
        raise ValueError(f'Lymphoid reference is missing obs[{label_column!r}]')
    labels = reference.obs[label_column].astype('string').str.strip()
    omitted = {'','nan','na','unknown','unassigned','tumor'} | {label.strip().casefold() for label in excluded}
    keep = labels.notna() & ~labels.str.casefold().isin(omitted)
    if not keep.any():
        raise ValueError('No lymphoid reference cells after label filtering')
    reference = raw_count_anndata(reference[keep.to_numpy()],'lymphoid reference')
    reference.obs = pd.DataFrame({SCANVI_LABEL_KEY:labels.loc[reference.obs_names].astype(str)},index=reference.obs_names)
    return reference


def load_reference(config):
    ensure_references(config,'lymphoid')
    path = config['refs']['lymphoid']
    global ad
    ad = ad_module()
    unfiltered = ad.read_h5ad(path)
    try:
        reference = filter_reference(unfiltered,config['label_column'],config['exclude_reference_labels'])
    except ValueError as exc:
        raise ValueError(
            f'{path}: {exc}. scVI needs raw UMI counts after excluding Tumor. '
            'If this GEO file contains only processed expression, provide a count-bearing H5AD '
            'and set references.lymphoid.file in study.json (remove its url for a local override).'
        ) from exc
    log(f'{path.name}: using {reference.n_obs:,} of {unfiltered.n_obs:,} cells as the lymphoid reference')
    return reference


def transfer(reference, query, folder):
    ref, q, _ = harmonize_genes(reference,raw_count_anndata(query,'query'),'lymphoid reference')
    ref.obs[BATCH_KEY] = 'ref'
    q.obs[BATCH_KEY] = 'query'
    q.obs[SCANVI_LABEL_KEY] = UNLABELED
    q.obs['_cell_id'] = q.obs_names.astype(str)
    combined = ad.concat([ref,q],label='source',keys=['ref','query'],join='inner',
                         index_unique='__',merge='same')
    combined.layers[COUNTS_LAYER] = combined.X.copy()
    sc.pp.normalize_total(combined,target_sum=1e4)
    sc.pp.log1p(combined)
    sc.pp.highly_variable_genes(combined,flavor='seurat_v3',n_top_genes=min(N_HVG,combined.n_vars),
                               layer=COUNTS_LAYER,batch_key=BATCH_KEY,subset=True)
    scvi.settings.seed = RANDOM_SEED
    scvi.model.SCVI.setup_anndata(combined,layer=COUNTS_LAYER,batch_key=BATCH_KEY)
    model = scvi.model.SCVI(combined,n_latent=N_LATENT)
    training = dict(batch_size=TRAIN_BATCH_SIZE,datasplitter_kwargs=TRAIN_DATASPLITTER_KWARGS.copy(),
                    accelerator='cpu',devices=1,default_root_dir=str(folder))
    model.train(**training)
    scanvi = scvi.model.SCANVI.from_scvi_model(model,unlabeled_category=UNLABELED,
                                              labels_key=SCANVI_LABEL_KEY)
    scanvi.train(max_epochs=20,n_samples_per_label=100,**training)
    probabilities = probability_frame(scanvi.predict(combined,soft=True),combined.obs_names)
    rows = combined.obs['source'].eq('query')
    probabilities = probabilities.loc[rows].set_axis(combined.obs.loc[rows,'_cell_id'].astype(str)).reindex(query.obs_names)
    if probabilities.isna().all(axis=1).any() or any(str(x).strip().casefold()=='tumor' for x in probabilities.columns):
        raise ValueError('Missing query predictions or unexpected Tumor class')
    predictions = pd.DataFrame({'predicted_labels':probabilities.idxmax(axis=1),'label_confidence':probabilities.max(axis=1)})
    predictions.to_csv(folder/'lymphoid_predictions.tsv.gz',sep='\t',index_label='cell_id')
    return predictions


def smooth_labels(obs, graph, clusters, options):
    labels = obs.lymphoid_atlas_merged_labels
    order = list(pd.Categorical(labels).categories.astype(str))
    codes = pd.Categorical(labels,categories=order).codes.astype(int)
    valid = codes >= 0
    if not valid.any() or clusters.isna().any():
        raise ValueError('Label smoothing requires labelled cells and complete Leiden clusters')
    binary = sp.csr_matrix(graph).copy()
    binary.setdiag(0)
    binary.eliminate_zeros()
    binary.sum_duplicates()
    binary.data[:] = 1
    onehot = sp.csr_matrix((np.ones(valid.sum()),(np.flatnonzero(valid),codes[valid])),
                           shape=(len(obs),len(order)))
    votes = (binary @ onehot).toarray()
    labeled_degree = votes.sum(axis=1)
    total_degree = np.asarray(binary.sum(axis=1)).ravel()
    own_votes = votes[np.arange(len(obs)),codes]
    _, cluster = np.unique(clusters.astype(str),return_inverse=True)
    size = np.bincount(cluster)
    by_cluster = np.zeros((len(size),len(order)),dtype=int)
    np.add.at(by_cluster,(cluster[valid],codes[valid]),1)
    majority = by_cluster.argmax(axis=1)[cluster]
    purity = by_cluster.max(axis=1)[cluster] / size[cluster]
    local = np.divide(votes[np.arange(len(obs)),majority],labeled_degree,out=np.zeros(len(obs)),where=labeled_degree>0)
    coverage = np.divide(labeled_degree,total_degree,out=np.zeros(len(obs)),where=total_degree>0)
    minimum = options['min_neighbor_majority']
    eligible = (valid & (codes != majority) & (size[cluster] >= options['min_cluster_cells'])
                & (purity >= options['min_cluster_majority']) & (labeled_degree >= options['min_labeled_neighbors'])
                & (coverage >= minimum) & (local >= minimum) & (own_votes <= options['max_same_label_neighbors']))
    rank = np.lexsort((obs.index.to_numpy(dtype=str),-purity,-local))
    budget = int(np.floor(len(obs)*options['max_changed_fraction']))
    changed = np.zeros(len(obs),dtype=bool)
    changed[rank[eligible[rank]][:budget]] = True
    codes[changed] = majority[changed]
    return [order[i] if i>=0 else None for i in codes], changed


def donor_rna(query, reference, sample, folder, curation, config):
    prediction = transfer(reference,query,folder)
    query.obs['lymphoid_atlas_predicted_labels'] = prediction.predicted_labels.reindex(query.obs_names)
    query.obs['lymphoid_atlas_label_confidence'] = prediction.label_confidence.reindex(query.obs_names)
    merged = {label:label_type for label_type,labels in MERGE_RULES.items() for label in labels}
    query.obs['lymphoid_atlas_merged_labels'] = query.obs.lymphoid_atlas_predicted_labels.map(lambda x:merged.get(x,x))

    exact = exact_obs(query.obs,sample['name'])
    passed, _ = quality_filter(exact,'RNA',sample['quality']['RNA'])
    if not passed.all():
        raise ValueError('Stage 01 RNA violates configured technical thresholds')
    decisions = curation_for(curation,'RNA',exact.index)
    eligible = ~decisions.fine_exclude
    keep = eligible & doubletdetection_pass(exact,config['rna_doublets']['max_score'],eligible)
    counts = query[keep.to_numpy()].copy()
    work = counts.copy()
    spec = config['rna_processing']
    sc.pp.filter_genes(work,min_cells=spec['min_cells_per_gene'])
    sc.pp.highly_variable_genes(work,flavor=spec['hvg_flavor'],n_top_genes=min(spec['n_hvg'],work.n_vars),
                               span=spec['hvg_span'],batch_key=spec['hvg_batch_key'])
    sc.pp.normalize_total(work,target_sum=spec['target_sum'])
    sc.pp.log1p(work)
    work = work[:,work.var.highly_variable.to_numpy()].copy()
    sc.pp.scale(work,max_value=spec['scale_max_value'])
    pcs = min(spec['n_pcs'],work.n_obs-1,work.n_vars-1)
    solver = spec['pca_solver']
    if sp.issparse(work.X) and solver == 'randomized':
        solver = 'arpack'
    sc.pp.pca(work,n_comps=pcs,svd_solver=solver,random_state=spec['seed'])
    sc.pp.neighbors(work,n_neighbors=min(spec['n_neighbors'],work.n_obs-1),n_pcs=pcs,
                    metric=spec['metric'],random_state=spec['seed'])
    sc.tl.umap(work,min_dist=spec['min_dist'],spread=spec['spread'],random_state=spec['seed'])
    options = sample['smoothing']
    key = f'leiden_final_fine_r{options["fine_resolution"]:g}'
    sc.tl.leiden(work,resolution=options['fine_resolution'],random_state=0,flavor='leidenalg',key_added=key)
    labels, changed = smooth_labels(work.obs,work.obsp['connectivities'],work.obs[key],options)

    counts.obs = work.obs.copy()
    counts.obs['lymphoid_raw'] = work.obs.lymphoid_atlas_predicted_labels.to_numpy()
    counts.obs['lymphoid_original'] = work.obs.lymphoid_atlas_merged_labels.to_numpy()
    counts.obs[LABEL] = labels
    counts.obs['lymphoid_smoothing_changed'] = changed
    counts.obs['lymphoid_smoothing_parameters'] = json.dumps(options,sort_keys=True)
    counts.obsm['X_pca_donor'] = work.obsm['X_pca'].copy()
    counts.obsm['X_umap_donor'] = work.obsm['X_umap'].copy()
    counts.obsp['donor_connectivities'] = work.obsp['connectivities'].copy()
    counts.obs = exact_obs(counts.obs,sample['name'])
    decisions = curation_for(curation,'RNA',counts.obs_names)
    counts.obs.loc[decisions.cdc2_target,LABEL] = 'DC_CDC2'
    counts.obs['fine_review_available'] = decisions.fine_review_available
    counts = counts[~(decisions.black_exclude | decisions.merged_exclude).to_numpy()].copy()
    return counts


def donor_dna(sample, mark, source, rna, curation, config, folder, signature):
    source_obs = metadata(source)
    technical, _ = quality_filter(source_obs,mark,sample['quality'][mark])
    if not technical.all():
        raise ValueError('Stage 01 histone profiles violate configured technical thresholds')
    decisions = curation_for(curation,mark,exact_obs(source_obs,sample['name']).index)
    keep = ~(decisions.fine_exclude | decisions.black_exclude | decisions.merged_exclude).to_numpy()
    retained = source_obs.loc[keep].copy()
    ref = rna.obs.copy()
    ref.index = ref.source_cell_id.astype(str)
    keys = [LABEL,'lymphoid_raw','lymphoid_original','lymphoid_atlas_predicted_labels',
            'lymphoid_atlas_merged_labels','lymphoid_atlas_label_confidence']
    state = dict(output=folder,signature=signature,dna_sources={mark:source},retained={mark:retained,'RNA':ref},
                 label_columns=keys,blacklist=config['refs']['blacklist'],config={'dna_processing':config['dna_processing']})
    args = SimpleNamespace(sample=sample['name'],threads=config['threads'])
    native, peaks = finalize_dna_fragments(state,args,mark,snap)
    output = {'native_fragments':native,'fragment_tsv':peak_paths(state,args,mark)['fragments'],'peaks_bed':peaks}
    for representation in ('tiles','peaks'):
        path = count_dna(state,args,mark,representation,native,peaks)
        data = ad.read_h5ad(path)
        data.obs = exact_obs(data.obs,sample['name'])
        data.obs[LABEL] = label_values(data.obs[LABEL])
        data.obs['fine_review_available'] = decisions.fine_review_available.reindex(data.obs_names).eq(True)
        output[representation] = data
    return output


def export_population(counts, donor_data, population, config, final, publication):
    root = final/population
    palette = config['palette']
    write_json(root/'palette.json',palette)
    settings = config['embedding']
    merged = {}
    for mark, full in counts.items():
        keep = label_values(full.obs[LABEL]).ne(UNLABELLED)
        data = full.copy() if population=='full' else full[keep.to_numpy()].copy()
        data = fit_merged(data,mark,settings)
        data.obs[LABEL] = pd.Categorical(label_values(data.obs[LABEL]))
        data.uns[LABEL+'_colors'] = np.array([palette.get(x,'#000000') for x in data.obs[LABEL].cat.categories])
        save_object(data,root/f'objects/{mark}_merged.h5ad')
        merged[mark] = data
        for sample in config['samples']:
            name = sample['name']
            mask = data.obs['sample'].astype(str).eq(name)
            if not mask.any():
                continue
            part = data[mask.to_numpy()].copy()
            part.obs['merged_cell_id'] = part.obs_names.astype(str)
            part.obs_names = part.obs.source_cell_id.astype(str)
            folder = root/'per_sample'/name
            suffix = 'counts' if mark=='RNA' else 'tiles'
            save_object(part,folder/f'{name}_{mark}_final_{suffix}.h5ad')
            part.obs.to_csv(folder/f'{name}_{mark}_cell_metadata.tsv.gz',sep='\t',index_label='cell_id')
            if mark!='RNA':
                aux = donor_data[name][mark]
                positions = aux['peaks'].obs_names.get_indexer(data.obs_names[mask])
                if (positions<0).any():
                    raise ValueError('Final peaks are missing retained profiles')
                peak = aux['peaks'][positions].copy()
                peak.obs = part.obs.copy()
                save_object(peak,folder/f'{name}_{mark}_final_peaks.h5ad')
                write_json(folder/f'{name}_{mark}_fragment_sources.json',
                    {k:aux[k].relative_to(final) for k in ('native_fragments','fragment_tsv','peaks_bed')})
    render_scope(merged,'merged',population,config,publication)

    for scope in ('Tonsil1','Tonsil1_gt600'):
        local = {}
        for mark, data in merged.items():
            keep = data.obs['sample'].astype(str).eq('Tonsil1').to_numpy()
            if not keep.any():
                continue
            depth = np.asarray(data.X.sum(axis=1)).ravel() if mark=='RNA' else data.obs.n_fragment.to_numpy()
            if scope.endswith('gt600'):
                keep &= depth>HIGH_COVERAGE
            part = data[keep].copy()
            part.obs['publication_gate_value'] = depth[keep]
            part = fit_local(part,mark,settings)
            suffix = 'counts' if mark=='RNA' else 'tiles'
            tag = 'gt600_' if scope.endswith('gt600') else ''
            save_object(part,publication/scope/population/'objects'/f'Tonsil1_{mark}_{tag}{suffix}.h5ad')
            local[mark] = part
        render_scope(local,scope,population,config,publication)


def build_final_results(config, sources, final, publication, signature, reference):
    curation = load_curation()
    counts = {m:[] for m in MODALITIES}
    donors = {}
    for sample in config['samples']:
        name = sample['name']
        folder = final/'donors'/name
        folder.mkdir(parents=True,exist_ok=True)
        destination = folder/'RNA_counts.h5ad'
        state = dict(output=folder,signature=signature)
        if phase_complete(state,'rna',[destination]):
            rna = ad.read_h5ad(destination)
        else:
            query = ad.read_h5ad(sources[name]['RNA'])
            rna = donor_rna(query,reference,sample,folder,curation,config)
            save_object(rna,destination)
            finish_phase(state,'rna',[destination])
            del query
            gc.collect()
        counts['RNA'].append(rna)
        donors[name] = {}
        for mark, source in sources[name].items():
            if mark=='RNA':
                continue
            aux = donor_dna(sample,mark,source,rna,curation,config,folder,signature)
            donors[name][mark] = aux
            counts[mark].append(aux['tiles'])
    del reference
    merged = {m:merge_counts(parts,m) for m,parts in counts.items() if parts}
    for population in ('full','labelled_only'):
        export_population(merged,donors,population,config,final,publication)


def label_stage(config, reference=None):
    global ad, sc, snap, scvi, THREADS
    THREADS = config['threads']
    thread_budget(THREADS)
    sources = {s['name']:processing_paths(config['output'],s) for s in config['samples']}
    require_files(*(p for group in sources.values() for p in group.values()),CURATION)
    ensure_references(config,'blacklist','lymphoid')
    if reference is None:
        reference = load_reference(config)
    tools_required('macs3','bedtools')
    ad = ad_module()
    import scanpy as sc
    import snapatac2 as snap
    import scvi
    import torch
    torch.set_num_threads(THREADS)

    final = config['output']/'final'
    publication = config['output']/'publication'
    inputs = {name:count_stamp(paths) for name,paths in sources.items()}
    script = Path(__file__).resolve()
    signature = digest(dict(config=identity(STUDY),curation=identity(CURATION),
        reference=identity(config['refs']['lymphoid']),blacklist=identity(config['refs']['blacklist']),inputs=inputs,
        code=digest(script.read_text())))
    with output_lock(config['output']/'.work'/'final'):
        plan = final/'run.json'
        if plan.exists() and read_json(plan)['signature'] != signature:
            raise ValueError('Stage 02 inputs or settings changed. Choose a new output_root in study.json.')
        write_json(plan,dict(signature=signature,reference_accession='GSE269695',
                             reference=identity(config['refs']['lymphoid']),inputs=inputs))
        build_final_results(config,sources,final,publication,signature,reference)


# Command line

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command')
    commands.add_parser('all', help='run processing, label transfer and study figures')
    refs = commands.add_parser('references', help='download public references and validate lymphoid counts')
    processing_refs = ('lymphoid', 'chrom_sizes', 'blacklist', 'tss_annotation')
    refs.add_argument('--only', nargs='+', choices=(*processing_refs, 'effect_annotation'))
    process = commands.add_parser('process', help='build RNA and histone profiles and apply QC')
    process.add_argument('--sample', choices=tuple(PROTOCOLS))
    commands.add_parser('label', help='transfer lymphoid labels and make final profiles and study figures')
    args = parser.parse_args(argv)
    config = load_study()
    if args.command == 'references':
        keys = args.only or processing_refs
        for key in dict.fromkeys(keys):
            if key != 'lymphoid':
                ensure_references(config, key)
        if 'lymphoid' in keys:
            load_reference(config)
    elif args.command == 'process':
        process_stage(config, args.sample)
    elif args.command == 'label':
        label_stage(config)
    else:
        processing_samples(config)
        tools_required('macs3', 'bedtools')
        reference = load_reference(config)
        process_stage(config)
        label_stage(config, reference)


if __name__ == '__main__':
    main()

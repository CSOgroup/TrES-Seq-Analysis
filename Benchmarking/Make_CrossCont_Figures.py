#!/usr/bin/env python3
"""
Requires MACS3 3.0.5, deepTools 3.5.6, pysam, numpy, pyBigWig, UCSC bedGraph tools,
liftOver, GNU sort, awk, SciPy
"""
from pathlib import Path

# Editable paths and source metadata.
# Data roots
DATA_ROOT = Path('/mnt/dataFast/ahrmad/ScH2R_202409/FullBench')
NAS_ROOT = Path('/mnt/dataNas/ahrmad')

# Reference and output directories
REFERENCE_DIR = Path(__file__).resolve().parent / "references"
OUTPUT_ROOT = Path(__file__).resolve().parent / 'crosscont_figures'
SCRATCH_ROOT = Path('/tmp')

# RNA sublibraries for published scMTR-seq Day 0 assignment.
MTR_RNA_ACCESSIONS = {f'S{i}': f'GSM{8529551 + i}' for i in range(1, 9)}
MTR_RNA_FILES = {'matrix': 'matrix.mtx.gz', 'barcodes': 'barcodes.tsv.gz', 'features': 'features.tsv.gz'}

# Input paths
INPUTS = {
    'MTRseq_H3K27ac': NAS_ROOT / 'MTRseq/dna_align/MTRseq1_A_H3K27ac_NoDup.bam',
    'MTRseq_H3K27ac_index': NAS_ROOT / 'MTRseq/dna_align/MTRseq1_A_H3K27ac_NoDup.bam.bai',
    'MTRseq_H3K27me3': NAS_ROOT / 'MTRseq/dna_align/MTRseq1_A_H3K27me3_NoDup.bam',
    'MTRseq_H3K27me3_index': NAS_ROOT / 'MTRseq/dna_align/MTRseq1_A_H3K27me3_NoDup.bam.bai',
    'TrES_H3K27ac': NAS_ROOT / 'TrES_K422/dna_align/Ori_K422_H3K27ac_NoDup.bam',
    'TrES_H3K27ac_index': NAS_ROOT / 'TrES_K422/dna_align/Ori_K422_H3K27ac_NoDup.bam.bai',
    'TrES_H3K27me3': NAS_ROOT / 'TrES_K422/dna_align/Ori_K422_H3K27me3_NoDup.bam',
    'TrES_H3K27me3_index': NAS_ROOT / 'TrES_K422/dna_align/Ori_K422_H3K27me3_NoDup.bam.bai',
    'uCoTargetX_OverNight_H3K27ac': NAS_ROOT / 'uCOTx_OverNight/dna_align/uCOTx_OvN_A_H3K27ac_NoDup.bam',
    'uCoTargetX_OverNight_H3K27ac_index': NAS_ROOT / 'uCOTx_OverNight/dna_align/uCOTx_OvN_A_H3K27ac_NoDup.bam.bai',
    'uCoTargetX_OverNight_H3K27me3': NAS_ROOT / 'uCOTx_OverNight/dna_align/uCOTx_OvN_A_H3K27me3_NoDup.bam',
    'uCoTargetX_OverNight_H3K27me3_index': NAS_ROOT / 'uCOTx_OverNight/dna_align/uCOTx_OvN_A_H3K27me3_NoDup.bam.bai',
    'InHouse_NanoCNT_H3K27ac': DATA_ROOT / 'InHouse_NanoCNT/ScKDMA_H3K27ac.tsv.gz',
    'InHouse_NanoCNT_H3K27me3': DATA_ROOT / 'InHouse_NanoCNT/ScKDMA_H3K27me3.tsv.gz',
    'public_NanoCNT_ac': DATA_ROOT / 'NanoCNT/GSE198467_H3K27ac_fragments.tsv.gz',
    'public_NanoCNT_me3': DATA_ROOT / 'NanoCNT/GSE198467_H3K27me3_fragments.tsv.gz',
    'public_uCoTargetX_H3K27ac': DATA_ROOT / 'uCOTx_published_BWs/GSM8148409_K562.H3K27ac.bw',
    'public_uCoTargetX_H3K27me3': DATA_ROOT / 'uCOTx_published_BWs/GSM8148410_K562.H3K27me3.bw',
    'MTR_S1_H3K27ac': DATA_ROOT / 'MTR/Technical/dna/GSM8529544_HIST_S1_H3K27ac.fragments.tsv.gz',
    'MTR_S2_H3K27ac': DATA_ROOT / 'MTR/Technical/dna/GSM8529545_HIST_S2_H3K27ac.fragments.tsv.gz',
    'MTR_S3_H3K27ac': DATA_ROOT / 'MTR/Technical/dna/GSM8529546_HIST_S3_H3K27ac.fragments.tsv.gz',
    'MTR_S4_H3K27ac': DATA_ROOT / 'MTR/Technical/dna/GSM8529547_HIST_S4_H3K27ac.fragments.tsv.gz',
    'MTR_S5_H3K27ac': DATA_ROOT / 'MTR/Technical/dna/GSM8529548_HIST_S5_H3K27ac.fragments.tsv.gz',
    'MTR_S6_H3K27ac': DATA_ROOT / 'MTR/Technical/dna/GSM8529549_HIST_S6_H3K27ac.fragments.tsv.gz',
    'MTR_S7_H3K27ac': DATA_ROOT / 'MTR/Technical/dna/GSM8529550_HIST_S7_H3K27ac.fragments.tsv.gz',
    'MTR_S8_H3K27ac': DATA_ROOT / 'MTR/Technical/dna/GSM8529551_HIST_S8_H3K27ac.fragments.tsv.gz',
    'MTR_S1_H3K27me3': DATA_ROOT / 'MTR/Technical/dna/GSM8529544_HIST_S1_H3K27me3.fragments.tsv.gz',
    'MTR_S2_H3K27me3': DATA_ROOT / 'MTR/Technical/dna/GSM8529545_HIST_S2_H3K27me3.fragments.tsv.gz',
    'MTR_S3_H3K27me3': DATA_ROOT / 'MTR/Technical/dna/GSM8529546_HIST_S3_H3K27me3.fragments.tsv.gz',
    'MTR_S4_H3K27me3': DATA_ROOT / 'MTR/Technical/dna/GSM8529547_HIST_S4_H3K27me3.fragments.tsv.gz',
    'MTR_S5_H3K27me3': DATA_ROOT / 'MTR/Technical/dna/GSM8529548_HIST_S5_H3K27me3.fragments.tsv.gz',
    'MTR_S6_H3K27me3': DATA_ROOT / 'MTR/Technical/dna/GSM8529549_HIST_S6_H3K27me3.fragments.tsv.gz',
    'MTR_S7_H3K27me3': DATA_ROOT / 'MTR/Technical/dna/GSM8529550_HIST_S7_H3K27me3.fragments.tsv.gz',
    'MTR_S8_H3K27me3': DATA_ROOT / 'MTR/Technical/dna/GSM8529551_HIST_S8_H3K27me3.fragments.tsv.gz',
    **{f'RT_{lib}_{kind}': DATA_ROOT / f'MTR/Technical/rna/{sample}_Endo.{lib}/{name}'
       for lib, sample in MTR_RNA_ACCESSIONS.items() for kind, name in MTR_RNA_FILES.items()},
    'hg19_sizes': REFERENCE_DIR / 'hg19.chrom.sizes',
    'hg19_blacklist': REFERENCE_DIR / 'hg19-blacklist.v2.bed.gz',
    'hg38_sizes': REFERENCE_DIR / 'hg38.chrom.sizes',
    'hg38_blacklist': REFERENCE_DIR / 'hg38-blacklist.v2.bed.gz',
    'mm10_sizes': REFERENCE_DIR / 'mm10.chrom.sizes',
    'mm10_blacklist': REFERENCE_DIR / 'mm10-blacklist.v2.bed.gz',
    'hg38_to_hg19': REFERENCE_DIR / 'hg38ToHg19.over.chain.gz',
}

# File URLs or accessions
SOURCES = {
    'MTRseq_H3K27ac': 'GSE324511',
    'MTRseq_H3K27ac_index': 'GSE324511',
    'MTRseq_H3K27me3': 'GSE324511',
    'MTRseq_H3K27me3_index': 'GSE324511',
    'TrES_H3K27ac': 'GSE324511',
    'TrES_H3K27ac_index': 'GSE324511',
    'TrES_H3K27me3': 'GSE324511',
    'TrES_H3K27me3_index': 'GSE324511',
    'uCoTargetX_OverNight_H3K27ac': 'GSE324511',
    'uCoTargetX_OverNight_H3K27ac_index': 'GSE324511',
    'uCoTargetX_OverNight_H3K27me3': 'GSE324511',
    'uCoTargetX_OverNight_H3K27me3_index': 'GSE324511',
    'InHouse_NanoCNT_H3K27ac': 'GSE324511',
    'InHouse_NanoCNT_H3K27me3': 'GSE324511',
    'public_NanoCNT_ac': 'GSE198467',
    'public_NanoCNT_me3': 'GSE198467',
    'public_uCoTargetX_H3K27ac': 'GSM8148409',
    'public_uCoTargetX_H3K27me3': 'GSM8148410',
    'MTR_S1_H3K27ac': 'GSM8529544',
    'MTR_S2_H3K27ac': 'GSM8529545',
    'MTR_S3_H3K27ac': 'GSM8529546',
    'MTR_S4_H3K27ac': 'GSM8529547',
    'MTR_S5_H3K27ac': 'GSM8529548',
    'MTR_S6_H3K27ac': 'GSM8529549',
    'MTR_S7_H3K27ac': 'GSM8529550',
    'MTR_S8_H3K27ac': 'GSM8529551',
    'MTR_S1_H3K27me3': 'GSM8529544',
    'MTR_S2_H3K27me3': 'GSM8529545',
    'MTR_S3_H3K27me3': 'GSM8529546',
    'MTR_S4_H3K27me3': 'GSM8529547',
    'MTR_S5_H3K27me3': 'GSM8529548',
    'MTR_S6_H3K27me3': 'GSM8529549',
    'MTR_S7_H3K27me3': 'GSM8529550',
    'MTR_S8_H3K27me3': 'GSM8529551',
    **{f'RT_{lib}_{kind}': sample
       for lib, sample in MTR_RNA_ACCESSIONS.items() for kind in MTR_RNA_FILES},
    'hg19_sizes': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes',
    'hg19_blacklist': 'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz',
    'hg38_sizes': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
    'hg38_blacklist': 'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz',
    'mm10_sizes': 'https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes',
    'mm10_blacklist': 'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz',
    'hg38_to_hg19': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz',
}
# End editable block.


import argparse
from collections import Counter
from functools import lru_cache
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile

for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from GEO_Preparation.merge_libraries import fragment_path

MARKS = ('H3K27ac', 'H3K27me3')
# Aggregate published tracks retain the approved signal-threshold exception.
PUBLISHED_UCOTARGETX_SIGNAL_CUTOFF = 10
DATASETS = {
    'MTRseq': ('hg38', 'bam', {m: f'MTRseq_{m}' for m in MARKS}),
    'MTR_published': ('hg38', 'fragments', {}),
    'uCoTargetX_overnight': ('hg38', 'bam', {m: f'uCoTargetX_OverNight_{m}' for m in MARKS}),
    'uCoTargetX_published': ('hg19', 'bigwig', {m: f'public_uCoTargetX_{m}' for m in MARKS}),
    'NanoCNT_inhouse': ('hg38', 'fragments', {m: f'InHouse_NanoCNT_{m}' for m in MARKS}),
    'NanoCNT_published': ('mm10', 'fragments', dict(zip(MARKS, ('public_NanoCNT_ac', 'public_NanoCNT_me3')))),
    'TrES': ('hg38', 'bam', {m: f'TrES_{m}' for m in MARKS}),
}
GEO_DATASETS = {
    'MTRseq': 'MTRseq', 'uCoTargetX_overnight': 'uCoTargetX_OverNight',
    'NanoCNT_inhouse': 'InHouse_NanoCNT', 'TrES': 'TrES',
}
FLANK, BIN_SIZE, QVALUE = 5000, 100, 1e-4
LABELS = {
    'MTRseq': 'MTR-seq (1 h)', 'MTR_published': 'scMTR-seq (H9 D0)',
    'uCoTargetX_overnight': 'uCoTargetX (overnight)', 'uCoTargetX_published': 'uCoTargetX (K562)',
    'NanoCNT_inhouse': 'Nano-CUT&Tag (K422)', 'NanoCNT_published': 'Nano-CUT&Tag (mouse)', 'TrES': 'TrES',
}
STYLE = '''font.family: DejaVu Sans
font.size: 9
pdf.fonttype: 42
ps.fonttype: 42
svg.fonttype: none
axes.spines.top: False
axes.spines.right: False
axes.linewidth: 0.7
axes.labelsize: 9
axes.titlesize: 10
xtick.labelsize: 8
ytick.labelsize: 8
lines.linewidth: 1.8
figure.facecolor: white
savefig.facecolor: white
'''


def assemble_mtr(mark, work, day0):
    """Pool published fragments for the Day 0 cells assigned from GEO RNA counts."""
    work = Path(work)
    target = work / f'MTR_{mark}.fragments.tsv'
    if target.is_file(): return target
    if not day0: raise ValueError('No Day 0 cells in the GEO RNA matrices')
    cells = work / 'MTR_day0_ids.txt'
    if not cells.is_file(): cells.write_text('\n'.join(sorted(day0)) + '\n')
    pending = target.with_suffix('.pending')
    # Keep all five fields, then restore cell/chromosome/coordinate order within each sublibrary.
    program = 'NR==FNR{keep[$1]=1;next} ($4 in keep){print $0}'
    with pending.open('wb') as out:
        for i in range(1, 9):
            raw = work / f'MTR_{mark}_S{i}.selected.tsv'
            with gzip.open(input_path(f'MTR_S{i}_{mark}'), 'rb') as source, raw.open('wb') as filtered:
                process = subprocess.Popen(['awk', '-F', '\t', program, str(cells), '-'],
                                           stdin=subprocess.PIPE, stdout=filtered)
                try:
                    shutil.copyfileobj(source, process.stdin, 1024*1024)
                    process.stdin.close()
                    if process.wait(): raise RuntimeError('MTR assembly cell filtering failed')
                finally:
                    if process.poll() is None: process.kill(); process.wait()
            subprocess.run(['sort', '-s', '-t', '\t', '-k4,4', '-k1,1', '-k2,2n', '-k3,3n',
                            '-S', '1G', '--parallel=2', '-T', str(work), str(raw)],
                           stdout=out, check=True, env=dict(os.environ, LC_ALL='C'))
            raw.unlink()
    pending.replace(target)
    return target


def input_path(key, replacements=None):
    path = INPUTS[key] if replacements is None else replacements.get(key, INPUTS[key])
    return path.expanduser().resolve()


def source_path(dataset, mark, work, replacements=None):
    if dataset == 'MTR_published': return assemble_mtr(mark, work, mtr_day0_ids())
    return input_path(DATASETS[dataset][2][mark], replacements)


def textopen(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def run(command, directory, env=None, stdout=None):
    with open(directory / 'commands.log', 'a') as log:
        log.write(repr(list(map(str, command))) + '\n'); log.flush()
        subprocess.run(list(map(str, command)), check=True, cwd=directory,
                       stdout=stdout if stdout is not None else log, stderr=log, env=env)


@lru_cache(None)
def references(genome):
    allowed = {f'chr{i}' for i in range(1, 20 if genome == 'mm10' else 23)} | {'chrX', 'chrY'}
    sizes = dict(sorted((c, int(n)) for c, n in
                        (line.split() for line in input_path(f'{genome}_sizes').read_text().splitlines())
                        if c in allowed))
    blocks = {c: [] for c in sizes}
    with textopen(input_path(f'{genome}_blacklist')) as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#', 'track', 'browser')): continue
            c, a, b, *_ = line.split(); c = 'chr' + c.removeprefix('chr')
            if c in sizes:
                a, b = max(0, int(a)), min(int(b), sizes[c])
                if a < b: blocks[c].append((a, b))
    for c, rows in blocks.items():
        merged = []
        for a, b in sorted(rows):
            if merged and a <= merged[-1][1]: merged[-1] = (merged[-1][0], max(b, merged[-1][1]))
            else: merged.append((a, b))
        blocks[c] = np.asarray(merged, dtype=np.int64).reshape(-1, 2)
    return sizes, blocks


def overlaps(blocks, c, a, b):
    rows = blocks[c]; i = np.searchsorted(rows[:, 1], a, side='right')
    return i < len(rows) and rows[i, 0] < b


def eligible(coordinate, genome):
    c, a, b = coordinate; sizes, blocks = references(genome)
    # The run checked full peak +/-5 kb, then deepTools centered on its midpoint.
    return (c in sizes and 0 <= a < b <= sizes[c] and a >= FLANK and b + FLANK <= sizes[c]
            and not overlaps(blocks, c, a - FLANK, b + FLANK))


@lru_cache(maxsize=1)
def mtr_day0_ids():
    from scipy.io import mmread

    barcode_pattern = re.compile(r'(?:(GSM\d+|S[1-8])[_:])?(\d{2}:\d{2}:\d{2}):(\d{2})(?:-1)?')
    day0 = set()
    for lib, sample in MTR_RNA_ACCESSIONS.items():
        with textopen(input_path(f'RT_{lib}_barcodes')) as handle:
            barcodes = [line.rstrip('\r\n').split('\t')[0] for line in handle]
        with textopen(input_path(f'RT_{lib}_features')) as handle:
            features = [line.rstrip().split('\t') for line in handle]
        with gzip.open(input_path(f'RT_{lib}_matrix'), 'rb') as handle:
            matrix = mmread(handle).tocsr()
        if matrix.shape != (len(features), len(barcodes)) or len(barcodes) != len(set(barcodes)):
            raise ValueError(f'Invalid RNA matrix dimensions or barcodes for {sample}')
        if not np.isfinite(matrix.data).all() or np.any(matrix.data < 0) or np.any(matrix.data != np.rint(matrix.data)):
            raise ValueError(f'Expected nonnegative integer RNA counts for {sample}')
        if features and len(features[0]) >= 3:
            expressed = np.array([bool(re.search(r'Gene Expression|RNA', row[2], re.I)) for row in features])
            if not expressed.any():
                raise ValueError(f'No RNA features for {sample}')
            matrix = matrix[expressed]
        counts = np.asarray(matrix.sum(axis=0, dtype=np.int64)).ravel()

        # Each physical cell retains its sublibrary. RT01 must hold a strict majority.
        total, rt01 = Counter(), Counter()
        for barcode, count in zip(barcodes, counts):
            match = barcode_pattern.fullmatch(barcode.strip())
            if match is None or match[1] not in (None, lib, sample):
                raise ValueError(f'Invalid RNA barcode in {sample}: {barcode}')
            cell, rt = match[2], match[3]
            if count:
                total[cell] += int(count)
                if rt == '01':
                    rt01[cell] += int(count)
        day0.update(f'{cell}:{lib}' for cell, count in rt01.items() if 2 * count > total[cell])
        del matrix
    return day0


def bam_rows(path, sizes):
    import pysam

    # Source BAMs are already filtered/deduplicated. One read1 represents one fragment.
    with pysam.AlignmentFile(str(path), 'rb') as bam, pysam.AlignmentFile(str(path), 'rb') as mates:
        lengths = {'chr' + c.removeprefix('chr'): n for c, n in zip(bam.references, bam.lengths)}
        if lengths.get('chr1') != sizes['chr1'] or any(lengths[c] != n for c, n in sizes.items() if c in lengths):
            raise ValueError(f'BAM assembly mismatch: {path}')
        for read in bam.fetch(until_eof=True):
            if not read.is_read1: continue
            c = 'chr' + (read.reference_name or 'unmapped').removeprefix('chr')
            if c not in sizes: continue
            a = min(read.reference_start, read.next_reference_start); b = a + abs(read.template_length)
            valid = read.reference_id == read.next_reference_id and 0 <= a < b <= sizes[c]
            if not valid:
                if (not read.is_paired or read.is_unmapped or read.mate_is_unmapped or a < 0
                        or read.reference_id != read.next_reference_id or read.is_secondary or read.is_supplementary): continue
                b = -1
                if read.has_tag('MC'):
                    cigar = read.get_tag('MC'); parts = re.findall(r'([1-9][0-9]*)([MIDNSHP=X])', cigar)
                    length = sum(int(n) for n, op in parts if op in 'MDN=X')
                    if ''.join(n + op for n, op in parts) == cigar and length > 0 and read.reference_end:
                        b = max(read.reference_end, read.next_reference_start + length)
                if not 0 <= a < b <= sizes[c]:
                    try: mate = mates.mate(read)
                    except ValueError: continue
                    if (mate.is_unmapped or mate.is_secondary or mate.is_supplementary or not mate.is_read2
                            or mate.reference_id != read.reference_id or mate.reference_start != read.next_reference_start
                            or mate.reference_end is None or read.reference_end is None
                            or (read.has_tag('RG') and mate.has_tag('RG') and read.get_tag('RG') != mate.get_tag('RG'))): continue
                    a, b = min(read.reference_start, mate.reference_start), max(read.reference_end, mate.reference_end)
                    if not 0 <= a < b <= sizes[c]: continue
            cell = read.get_tag('CB')
            if not cell or any(x.isspace() for x in cell): raise ValueError('Invalid CB tag')
            yield c, a, b, cell


def fragment_rows(path, sizes):
    with textopen(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#', 'track', 'browser')): continue
            c, a, b, cell, *_ = line.split(); c = 'chr' + c.removeprefix('chr')
            if c not in sizes: continue
            a, b = int(a), int(b)
            if not 0 <= a < b <= sizes[c]: raise ValueError(f'Invalid source interval: {line.strip()}')
            yield c, a, b, cell


def prepare_fragments(dataset, mark, directory, replacements=None):
    genome = DATASETS[dataset][0]
    sizes, blocks = references(genome)
    source = source_path(dataset, mark, directory.parent, replacements)
    records = bam_rows(source, sizes) if source.suffix == '.bam' else fragment_rows(source, sizes)
    counts = Counter(); raw = directory / 'canonical.tsv'; selected = directory / 'selected.bedpe'
    with open(raw, 'w') as out:
        for c, a, b, cell in records:
            if dataset == 'MTR_published' and not re.fullmatch(r'\d{2}:\d{2}:\d{2}:S[1-8]', cell):
                raise ValueError('Published MTR requires complete BC1:BC2:BC3:S# fragment identities')
            counts[cell] += 1; out.write(f'{c}\t{a}\t{b}\t{cell}\n')
    minimum = (200 if mark == MARKS[0] else 300) if dataset in ('MTRseq', 'MTR_published') else 500
    keep = {c for c, n in counts.items() if minimum <= n <= 15000}
    n = 0
    with open(raw) as handle, open(selected, 'w') as out:
        for line in handle:
            c, a, b, cell = line.split(); a, b = int(a), int(b)
            if cell in keep and not overlaps(blocks, c, a, b):
                out.write(f'{c}\t{a}\t{b}\n'); n += 1
    raw.unlink()
    if not n: raise ValueError(f'Empty cohort: {dataset}/{mark}')
    print(f'{dataset}/{mark}: {len(keep):,} cells; {n:,} retained fragments', flush=True)
    run(['macs3', 'callpeak', '-t', selected, '-f', 'BEDPE', '--keep-dup', 'all', '--nolambda',
         '-g', 'mm' if genome == 'mm10' else 'hs', '-q', QVALUE, '-B', '--SPMR', '--outdir', directory, '-n', 'own'], directory)
    ordered = directory / 'pileup.sorted.bdg'; sizes_file = directory / 'canonical.chrom.sizes'
    sizes_file.write_text(''.join(f'{c}\t{n}\n' for c, n in sizes.items()))
    with open(ordered, 'w') as out:
        run(['sort', '-k1,1', '-k2,2n', '-S', '4G', '-T', directory, directory / 'own_treat_pileup.bdg'],
            directory, env=dict(os.environ, LC_ALL='C'), stdout=out)
    signal = directory / 'signal.spmr.bw'
    run(['bedGraphToBigWig', ordered, sizes_file, signal], directory)
    return directory / 'own_peaks.narrowPeak', signal


def prepare_bigwig(dataset, mark, directory):
    import pyBigWig

    sizes, _ = references('hg19'); signal = source_path(dataset, mark, directory.parent)
    with pyBigWig.open(str(signal)) as bw:
        lengths = bw.chroms()
        if lengths.get('chr1') != sizes['chr1'] or any(lengths[c] != n for c, n in sizes.items() if c in lengths):
            raise ValueError('Published bigWig assembly mismatch')
    raw, continuous = directory / 'original.bdg', directory / 'continuous.bdg'
    run(['bigWigToBedGraph', signal, raw], directory)
    ends = {c: 0 for c in sizes}
    with open(raw) as handle, open(continuous, 'w') as out:
        for line in handle:
            c, a, b, value = line.split()
            if c not in sizes: continue
            a, b, value = int(a), int(b), float(value)
            if not (ends[c] <= a < b <= sizes[c]) or not np.isfinite(value) or value < 0:
                raise ValueError(f'Invalid published signal: {line.strip()}')
            if a > ends[c]: out.write(f'{c}\t{ends[c]}\t{a}\t0\n')
            out.write(line); ends[c] = b
        for c, b in ends.items():
            if b < sizes[c]: out.write(f'{c}\t{b}\t{sizes[c]}\t0\n')
    ordered = directory / 'sorted.bdg'
    with open(ordered, 'w') as out:
        run(['sort', '-k1,1', '-k2,2n', '-S', '4G', '-T', directory, continuous],
            directory, env=dict(os.environ, LC_ALL='C'), stdout=out)
    calls = directory / 'own.score_peaks.narrowPeak'
    run(['macs3', 'bdgpeakcall', '-i', ordered, '-o', calls, '-c', PUBLISHED_UCOTARGETX_SIGNAL_CUTOFF,
         '-l', '200', '-g', '150'], directory)
    return calls, signal


def peak_intervals(dataset, callfile, directory):
    genome = DATASETS[dataset][0]; coordinates = []
    with open(callfile) as handle:
        for line in handle:
            if not line.strip() or line.startswith(('#', 'track', 'browser')): continue
            fields = line.split(); c, a, b = fields[:3]
            row = ('chr' + c.removeprefix('chr'), int(a), int(b))
            if DATASETS[dataset][1] != 'bigwig' and float(fields[8]) < -np.log10(QVALUE):
                raise ValueError(f'Peak fails q={QVALUE}: {row}')
            if eligible(row, genome): coordinates.append(row)
    # These two reference datasets had a shared hg19 eligibility filter even in own plots.
    if dataset in ('TrES', 'NanoCNT_inhouse'):
        source, mapped = directory / 'liftover_input.bed', directory / 'mapped.bed'
        source.write_text(''.join(f'{c}\t{a}\t{b}\tp{i}\n' for i, (c, a, b) in enumerate(coordinates)))
        run(['liftOver', '-minMatch=0.95', source, input_path('hg38_to_hg19'),
             mapped, directory / 'unmapped.bed'], directory)
        rows = [line.split()[:4] for line in mapped.read_text().splitlines() if line.strip()]
        copies = Counter(row[3] for row in rows)
        keep = {int(ident[1:]) for c, a, b, ident in rows if copies[ident] == 1 and eligible((c, int(a), int(b)), 'hg19')}
        coordinates = [row for i, row in enumerate(coordinates) if i in keep]
    if not coordinates: raise ValueError(f'No eligible peaks: {dataset}')
    return sorted(coordinates)


def write_bed(path, rows):
    path.write_text(''.join(f'{c}\t{a}\t{b}\n' for c, a, b in rows))


def read_bed(path):
    with open(path) as handle:
        return [(p[0], int(p[1]), int(p[2])) for line in handle
                if line.strip() and not line.startswith(('#', 'track', 'browser')) for p in [line.split()]]


def plot(dataset, mark, matrix, count, directory):
    pdf, table, style = directory / 'profile.pdf', directory / 'profile.tsv', directory / 'matplotlibrc'
    style.write_text(STYLE)
    units = 'Provided bigWig units' if DATASETS[dataset][1] == 'bigwig' else 'Coverage / million fragments'
    run(['plotProfile', '-m', matrix, '-out', pdf, '--plotFileFormat', 'pdf', '--perGroup',
         '--averageType', 'mean', '--plotType', 'lines', '--colors', '#168A91', '#B94D70', '--samplesLabel', *MARKS,
         '--regionsLabel', f'{mark} | Own peaks | n={count:,}', '--plotTitle', LABELS[dataset] + '\n' + units,
         '--yAxisLabel', 'Mean signal', '--legendLocation', 'best', '--plotWidth', '10', '--plotHeight', '8',
         '--refPointLabel', 'Distance from peak center', '--outFileNameData', table], directory,
        env=dict(os.environ, MPLCONFIGDIR=str(directory / 'mpl'), MATPLOTLIBRC=str(style)))
    return pdf, table


def compute(dataset, mark, calls, signals, directory, threads):
    rows = peak_intervals(dataset, calls, directory)
    # A called peak must overlap signal in at least one fragment channel (the run's zero filter).
    if DATASETS[dataset][1] != 'bigwig':
        tracks = [pyBigWig.open(str(p)) for p in signals]
        try:
            rows = [r for r in rows if any((bw.stats(*r, type='max', exact=True)[0] or 0) > 0 for bw in tracks)]
        finally:
            for bw in tracks: bw.close()
    bed, matrix = directory / 'peaks.bed', directory / 'matrix.gz'; write_bed(bed, rows)
    run(['computeMatrix', 'reference-point', '--referencePoint', 'center', '-S', *signals, '-R', bed,
         '-b', FLANK, '-a', FLANK, '--binSize', BIN_SIZE, '--averageTypeBins', 'mean', '--missingDataAsZero',
         '--sortRegions', 'keep', '--samplesLabel', *MARKS, '-p', threads, '-o', matrix,
         '--outFileSortedRegions', directory / 'matrix_regions.bed'], directory)
    if Counter(read_bed(directory / 'matrix_regions.bed')) != Counter(rows): raise ValueError('deepTools changed the peak set')
    return bed, matrix, len(rows)


def check_inputs(datasets=DATASETS, replacements=None):
    if set(INPUTS) != set(SOURCES): raise ValueError('INPUTS and SOURCES keys differ')
    required = set()
    for dataset in datasets:
        genome, _, marks = DATASETS[dataset]
        required.update((f'{genome}_sizes', f'{genome}_blacklist'))
        required.update(marks.values())
        required.update(f'{key}_index' for key in marks.values() if input_path(key, replacements).suffix == '.bam')
        if dataset in ('TrES', 'NanoCNT_inhouse'):
            required.update(('hg19_sizes', 'hg19_blacklist', 'hg38_to_hg19'))
        if dataset == 'MTR_published':
            required.update(key for key in INPUTS if key.startswith(('MTR_S', 'RT_S')))
    paths = [input_path(k, replacements) for k in INPUTS if k in required]
    missing = [str(p) for p in paths if not p.is_file()]
    if missing: raise FileNotFoundError('Missing required input files:\n' + '\n'.join(missing))
    if len(set(paths)) != len(paths): raise ValueError('INPUTS must list each file only once')
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--datasets', nargs='+', choices=DATASETS, default=list(DATASETS))
    parser.add_argument('--dataset-root', type=Path,
                        help='Root of the merged in-house dataset folders')
    parser.add_argument('--check-inputs', action='store_true', help='Check required paths and stop')
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()
    if args.threads < 1: parser.error('--threads must be positive')
    replacements = {}
    if args.dataset_root is not None:
        root = args.dataset_root.expanduser().resolve()
        for dataset, name in GEO_DATASETS.items():
            for mark, key in DATASETS[dataset][2].items():
                replacements[key] = fragment_path(root, name, mark)
    paths = check_inputs(args.datasets, replacements)
    if args.check_inputs:
        print(f'All {len(paths)} required files exist.'); return
    if OUTPUT_ROOT.exists(): raise FileExistsError(f'Refusing to overwrite output root: {OUTPUT_ROOT}')
    SCRATCH_ROOT.mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix='run_', dir=SCRATCH_ROOT)).resolve()
    os.environ['TMPDIR'] = str(scratch)
    print(f'Scratch (including commands, matrices and profile tables): {scratch}', flush=True)
    OUTPUT_ROOT.mkdir(parents=True)
    for dataset in args.datasets:
        jobs = {}
        for mark in MARKS:
            directory = scratch / f'{dataset}_{mark}'; directory.mkdir()
            if DATASETS[dataset][1] == 'bigwig':
                jobs[mark] = prepare_bigwig(dataset, mark, directory)
            else:
                jobs[mark] = prepare_fragments(dataset, mark, directory, replacements)
        for mark in MARKS:
            directory = scratch / f'{dataset}_{mark}'
            bed, matrix, count = compute(dataset, mark, jobs[mark][0], [jobs[m][1] for m in MARKS], directory, args.threads)
            pdf, table = plot(dataset, mark, matrix, count, directory)
            stem = f'{dataset}_own_{mark}_profile'
            shutil.copyfile(bed, OUTPUT_ROOT / f'{dataset}_own_{mark}_peaks.bed')
            shutil.copyfile(pdf, OUTPUT_ROOT / f'{stem}.pdf')
            print(f'{stem}.pdf: {count:,} peaks', flush=True)


if __name__ == '__main__':
    main()

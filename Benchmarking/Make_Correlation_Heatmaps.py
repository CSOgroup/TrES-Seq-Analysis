#!/usr/bin/env python3
"""Two panel-k heatmaps: 10 kb bigWig bins and Spearman correlation.
Requires deepTools 3.5.5, numpy, pyBigWig and pysam. Existing bigWigs retain provided units.
"""
from pathlib import Path

# Editable paths and source metadata.
# Data roots
DATA_ROOT = Path('/mnt/dataFast/ahrmad/ScH2R_202409/FullBench')
BACKUP_ROOT = Path('/mnt/svraw1/Ahrmad/')
CORR_ROOT = BACKUP_ROOT / 'CORR'

# Reference and output directories
REFERENCE_DIR = Path(__file__).resolve().parents[1] / 'references'
OUTPUT_ROOT = Path(__file__).resolve().parent / 'panel_k_figures'
SCRATCH_ROOT = Path('/tmp')

# Input paths
INPUTS = {
    'K9_inhouse_ChIP': CORR_ROOT / 'ALL_marks/In_ChIP_H3K9me3.bw',
    'K9_ENCODE_ChIP': CORR_ROOT / 'ALL_marks/ENC_ChIP_H3K9me3.bw',
    'K9_TrES_dual': CORR_ROOT / 'ALL_marks/H2R_DT_Sc_H3K9me3.bw',
    'K27ac_TrES_rep1': CORR_ROOT / 'ALL_marks/H2R_Sc_H3K27ac_1.bw',
    'K27ac_TrES_rep2': CORR_ROOT / 'ALL_marks/H2R_Sc_H3K27ac_2.bw',
    'K27ac_nanoCT': CORR_ROOT / 'ALL_marks/NanoCutTag_H3K27ac.bw',
    'K27me3_TrES_rep1': CORR_ROOT / 'H3K27me3/H2R_SingT_Sn09_H3K27me3.bw',
    'K27me3_TrES_rep2': CORR_ROOT / 'ALL_marks/H2R_Sc_H3K27me3.bw',
    'K27me3_TrES_dual': CORR_ROOT / 'ALL_marks/H2R_DT_Sc_H3K27me3.bw',
    'K27me3_nanoCT': CORR_ROOT / 'ALL_marks/NanoCutTag_H3K27me3.bw',
    'RNA_TrES_rep1_plus': CORR_ROOT / 'RNA/H2R_Sn_1_RNAplus.bw',
    'RNA_TrES_rep1_minus': CORR_ROOT / 'RNA/H2R_Sn_1_RNAminus.bw',
    'RNA_TrES_rep2_plus': CORR_ROOT / 'RNA/H2R_Sn_2_RNAplus.bw',
    'RNA_TrES_rep2_minus': CORR_ROOT / 'RNA/H2R_Sn_2_RNAminus.bw',
    'RNA_TrES_dual_plus': CORR_ROOT / 'RNA/H2R_SnBulk_K9_RNAplus.bw',
    'RNA_TrES_dual_minus': CORR_ROOT / 'RNA/H2R_SnBulk_K9_RNAminus.bw',
    'RNA_10x_plus': CORR_ROOT / 'RNA/10x_Sc_RNAplus.bw',
    'RNA_10x_minus': CORR_ROOT / 'RNA/10x_Sc_RNAminus.bw',
    'RNA_ENCODE_plus': CORR_ROOT / 'RNA/ENCODE_Bulk_RNAplus.bw',
    'RNA_ENCODE_minus': CORR_ROOT / 'RNA/ENCODE_Bulk_RNAminus.bw',
    'hg38_sizes': REFERENCE_DIR / 'hg38.chrom.sizes',
    'hg38_blacklist': REFERENCE_DIR / 'hg38-blacklist.v2.bed.gz',
}
# To replace an input with a verified hg38 BAM, change its INPUTS path, add its
# explicit index key to INPUTS/SOURCES, and declare (index_key, RNA strand filter).
# RNA strand filters must be verified for the library; use forward or reverse.
BAM_TRACKS = {}

SOURCES = {
    'K9_inhouse_ChIP': 'GSE324511', 'K9_ENCODE_ChIP': 'https://www.encodeproject.org/experiments/ENCSR978BQI/', 'K9_TrES_dual': 'GSE324511',
    'K27ac_TrES_rep1': 'GSE324511', 'K27ac_TrES_rep2': 'GSE324511', 'K27ac_nanoCT': 'GSE324511',
    'K27me3_TrES_rep1': 'GSE324511', 'K27me3_TrES_rep2': 'GSE324511', 'K27me3_TrES_dual': 'GSE324511', 'K27me3_nanoCT': 'GSE324511',
    'RNA_TrES_rep1_plus': 'GSE324511', 'RNA_TrES_rep1_minus': 'GSE324511',
    'RNA_TrES_rep2_plus': 'GSE324511', 'RNA_TrES_rep2_minus': 'GSE324511',
    'RNA_TrES_dual_plus': 'GSE324511', 'RNA_TrES_dual_minus': 'GSE324511',
    'RNA_10x_plus': 'GSE324511', 'RNA_10x_minus': 'GSE324511', 'RNA_ENCODE_plus': 'https://www.encodeproject.org/experiments/ENCSR314LXG/', 'RNA_ENCODE_minus': 'https://www.encodeproject.org/experiments/ENCSR314LXG/',
    'hg38_sizes': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
    'hg38_blacklist': 'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz',
}
# End editable block.

import argparse
import csv
import gzip
import importlib.metadata
import os
import subprocess
import sys
import tempfile
import numpy as np
import pyBigWig
import pysam

HISTONE_TRACKS = [
    ('K9_inhouse_ChIP', 'H3K9me3 | ChIP-seq (in-house)'),
    ('K9_ENCODE_ChIP', 'H3K9me3 | ChIP-seq (ENCODE)'),
    ('K9_TrES_dual', 'H3K9me3 | TrES dual tag'),
    ('K27ac_TrES_rep1', 'H3K27ac | TrES single tag rep1'),
    ('K27ac_TrES_rep2', 'H3K27ac | TrES single tag rep2'),
    ('K27ac_nanoCT', 'H3K27ac | nano-CT'),
    ('K27me3_TrES_rep1', 'H3K27me3 | TrES single tag rep1'),
    ('K27me3_TrES_rep2', 'H3K27me3 | TrES single tag rep2'),
    ('K27me3_TrES_dual', 'H3K27me3 | TrES dual tag'),
    ('K27me3_nanoCT', 'H3K27me3 | nano-CT'),
]
RNA_TRACKS = [
    ('RNA_TrES_rep1_plus', 'TrES single tag rep1 (+)'),
    ('RNA_TrES_rep1_minus', 'TrES single tag rep1 (-)'),
    ('RNA_TrES_rep2_plus', 'TrES single tag rep2 (+)'),
    ('RNA_TrES_rep2_minus', 'TrES single tag rep2 (-)'),
    ('RNA_TrES_dual_plus', 'TrES dual tag (+)'),
    ('RNA_TrES_dual_minus', 'TrES dual tag (-)'),
    ('RNA_10x_plus', '10x scRNA-seq (+)'), ('RNA_10x_minus', '10x scRNA-seq (-)'),
    ('RNA_ENCODE_plus', 'ENCODE bulk RNA (+)'), ('RNA_ENCODE_minus', 'ENCODE bulk RNA (-)'),
]
BIN_SIZE, COVERAGE_BIN_SIZE, EFFECTIVE_GENOME_SIZE = 10000, 20, 2913022398


def sizes():
    return {c: int(n) for c, n in (line.split() for line in INPUTS['hg38_sizes'].read_text().splitlines())}


def check_inputs():
    if set(INPUTS) != set(SOURCES): raise ValueError('INPUTS/SOURCES keys differ')
    missing = [str(p) for p in INPUTS.values() if not p.is_file()]
    if missing: raise FileNotFoundError('\n'.join(missing))
    expected = sizes()
    for key, _ in HISTONE_TRACKS + RNA_TRACKS:
        if key in BAM_TRACKS:
            index, strand = BAM_TRACKS[key]
            if key.startswith('RNA_') and strand not in ('forward', 'reverse'):
                raise ValueError(f'{key}: specify the verified deepTools RNA strand convention')
            with pysam.AlignmentFile(str(INPUTS[key]), 'rb', index_filename=str(INPUTS[index])) as bam:
                if not bam.has_index() or bam.get_reference_length('chr1') != expected['chr1']:
                    raise ValueError(f'{key}: indexed hg38 BAM required')
        else:
            with pyBigWig.open(str(INPUTS[key])) as bw:
                chroms = bw.chroms()
                if chroms.get('chr1') != expected['chr1'] or any(n != expected[c] for c, n in chroms.items() if c in expected):
                    raise ValueError(f'{key}: bigWig assembly mismatch')
    plus, minus = INPUTS['RNA_TrES_dual_plus'], INPUTS['RNA_TrES_dual_minus']
    if (plus.parent != minus.parent or plus.name.removesuffix('RNAplus.bw') != 'H2R_SnBulk_K9_'
            or minus.name.removesuffix('RNAminus.bw') != 'H2R_SnBulk_K9_'):
        raise ValueError('Dual-tag RNA inputs must be the paired H2R_SnBulk_K9 plus/minus tracks')
    with pyBigWig.open(str(plus)) as a, pyBigWig.open(str(minus)) as b:
        if a.chroms() != b.chroms(): raise ValueError('K9 strand tracks have different reference headers')
    return expected


def input_rows():
    return [(k, str(p), SOURCES[k],
             'BAM: regenerated RPGC, 20 bp' if k in BAM_TRACKS else
             'bigWig: generation/normalization/bin/strand provenance unverified' if p.suffix in ('.bw', '.bigWig') else '')
            for k, p in INPUTS.items()]


def run(command, scratch):
    with (scratch / 'commands.log').open('a') as log:
        log.write(repr(list(map(str, command))) + '\n'); log.flush()
        subprocess.run(list(map(str, command)), check=True, stdout=log, stderr=log,
                       env=dict(os.environ, TMPDIR=str(scratch), MPLCONFIGDIR=str(scratch / 'mpl')))


def coverage(key, scratch, threads, region=None):
    if key not in BAM_TRACKS: return INPUTS[key]
    index, strand = BAM_TRACKS[key]
    # Custom index names are passed through a scratch BAM symlink.
    bam = scratch / (key + '.bam'); bam.symlink_to(INPUTS[key].resolve())
    Path(str(bam) + '.bai').symlink_to(INPUTS[index].resolve())
    output = scratch / (key + '.bw')
    command = ['bamCoverage', '-b', bam, '-o', output, '--normalizeUsing', 'RPGC',
               '--effectiveGenomeSize', EFFECTIVE_GENOME_SIZE, '--binSize', COVERAGE_BIN_SIZE, '-p', threads]
    if strand: command += ['--filterRNAstrand', strand]
    if region: command += ['--region', region]
    run(command, scratch)
    return output


def make_heatmap(name, tracks, scratch, output, threads, region=None):
    signals = [coverage(key, scratch, threads, region) for key, _ in tracks]
    matrix = output / f'panel_k_{name}_bins.npz'
    command = ['multiBigwigSummary', 'bins', '-b', *signals, '--labels', *[label for _, label in tracks],
               '--binSize', BIN_SIZE, '--numberOfProcessors', threads, '-o', matrix]
    # The supplied DNA commands used a blacklist; the RNA commands did not.
    if name == 'histone':
        blacklist = scratch / 'blacklist.bed'
        with gzip.open(INPUTS['hg38_blacklist'], 'rt') as source: blacklist.write_text(source.read())
        command += ['--blackListFileName', blacklist]
    if region: command += ['--region', region]
    run(command, scratch)
    values = np.load(matrix)
    if values['matrix'].shape[1] != len(tracks) or list(values['labels']) != [label for _, label in tracks]:
        raise ValueError('deepTools changed the track order')
    run(['plotCorrelation', '-in', matrix, '-c', 'spearman', '-p', 'heatmap',
         '-o', output / f'panel_k_{name}_heatmap.pdf', '--plotNumbers', '--colorMap', 'jet',
         '--zMin', '-1', '--zMax', '1', '--outFileCorMatrix', output / f'panel_k_{name}_spearman.tsv'], scratch)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-inputs', action='store_true')
    parser.add_argument('--list-inputs', action='store_true')
    parser.add_argument('--smoke', action='store_true', help='Both heatmaps over chr8:127000000:129000000 in scratch')
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()
    check_inputs()
    if args.list_inputs:
        writer = csv.writer(sys.stdout, delimiter='\t'); writer.writerow(['key', 'path', 'source', 'generation']); writer.writerows(input_rows()); return
    if args.check_inputs: print(f'All {len(INPUTS)} input files exist; hg38 track headers match.'); return
    if args.threads < 1: parser.error('--threads must be positive')
    if importlib.metadata.version('deepTools') != '3.5.5': raise RuntimeError('Use deepTools 3.5.5')
    SCRATCH_ROOT.mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix='panel_k_', dir=SCRATCH_ROOT))
    output = scratch / 'figures' if args.smoke else OUTPUT_ROOT
    output.mkdir(parents=True, exist_ok=False)
    with (output / 'track_inputs.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t'); writer.writerow(['key', 'path', 'source', 'generation']); writer.writerows(input_rows())
    print('FLAG: provided bigWig generation settings and exact public-file identities are unverified.', file=sys.stderr)
    region = 'chr8:127000000:129000000' if args.smoke else None
    make_heatmap('histone', HISTONE_TRACKS, scratch, output, args.threads, region)
    make_heatmap('rna', RNA_TRACKS, scratch, output, args.threads, region)
    print(f'Wrote both heatmaps: {output}\nCommands/intermediates: {scratch}')


if __name__ == '__main__': main()

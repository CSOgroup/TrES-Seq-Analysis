##########################################
##########################################
# GENERATE BW FOR RNA 10x
##########################################
##########################################

#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import glob

import scanpy as sc


# ----------------- helpers ----------------- #

def run_cmd(cmd, stdout=None):
    """
    Helper to run a shell command and fail loudly if it breaks.
    """
    print(">> CMD:", " ".join(cmd), file=sys.stderr)
    try:
        subprocess.run(cmd, check=True, stdout=stdout)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: command failed with exit code {e.returncode}", file=sys.stderr)
        sys.exit(e.returncode)


def infer_sample_name(args):
    """
    If --sample-name is not given, infer it from the h5ad basename (no extension).
    """
    if args.sample_name is not None:
        return args.sample_name

    base = os.path.basename(args.h5ad)
    if base.endswith(".h5ad"):
        base = base[:-5]
    print(f">> Inferring sample name from h5ad: {base}", file=sys.stderr)
    return base


# ----------------- core steps ----------------- #

def parse_args():
    ap = argparse.ArgumentParser(
        description="Use all cells from RNA h5ad, subset BAM by CB tag, and generate stranded/unstranded bigWigs with STAR."
    )

    # I/O
    ap.add_argument(
        "--h5ad",
        required=True,
        help="Input RNA .h5ad file (e.g. ../../ScRNA_Natalya_Post.h5ad)",
    )
    ap.add_argument(
        "--bam-in",
        required=True,
        help="Input BAM file with CB tags (e.g. 10x_ScRNA.bam)",
    )
    ap.add_argument(
        "--outdir",
        required=True,
        help="Output directory",
    )
    ap.add_argument(
        "--sample-name",
        required=False,
        default=None,
        help="Sample name for output files (default: basename of h5ad).",
    )

    # External tools / references
    ap.add_argument(
        "--genome-dir",
        required=True,
        help="STAR genomeDir",
    )
    ap.add_argument(
        "--chrom-sizes",
        required=True,
        help="chrom.sizes file for bedGraphToBigWig",
    )
    ap.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Threads for samtools & STAR",
    )

    # Explicit cell list path (optional)
    ap.add_argument(
        "--cell-list",
        default=None,
        help=(
            "Optional path for output list of cell barcodes (one per line). "
            "These must match the CB tag values in the BAM. "
            "If not set, a default is used in outdir."
        ),
    )

    # Normalization for STAR output
    ap.add_argument(
        "--wig-norm",
        default="RPM",
        choices=["None", "RPM", "RPKM"],
        help="Normalization for STAR --outWigNorm (default: RPM, good for correlations).",
    )

    return ap.parse_args()


def write_all_cells_from_h5ad(args, sample_name):
    """
    Take *all* cell barcodes from adata.obs.index and write them to a text file.
    These must match the CB tags in the BAM.
    Returns (cell_list_path, n_cells).
    """
    os.makedirs(args.outdir, exist_ok=True)

    print(f">> Reading h5ad: {args.h5ad}", file=sys.stderr)
    adata = sc.read_h5ad(args.h5ad)
    barcodes = list(adata.obs.index)
    n_cells = len(barcodes)

    print(f">> Number of cells in h5ad (used as-is): {n_cells}", file=sys.stderr)

    if n_cells == 0:
        print("ERROR: h5ad has zero cells. Aborting.", file=sys.stderr)
        sys.exit(1)

    # Where to write the barcodes list
    if args.cell_list is not None:
        cell_list_path = args.cell_list
    else:
        cell_list_path = os.path.join(
            args.outdir,
            f"{sample_name}_cell_barcodes.lst",
        )

    print(f">> Writing cell barcode list to: {cell_list_path}", file=sys.stderr)
    with open(cell_list_path, "w") as f:
        for bc in barcodes:
            f.write(f"{bc}\n")

    return cell_list_path, n_cells


def subset_bam(args, cell_list_path, sample_name):
    """
    Use samtools to subset BAM based on CB tag values.
    Uses: samtools view --tag-file CB:<file>
    Only alignments with CB tags listed in <file> are kept.
    Returns path to filtered BAM.
    """
    out_bam = os.path.join(args.outdir, f"{sample_name}_Filt.bam")

    cmd_view = [
        "samtools",
        "view",
        "-@",
        str(args.threads),
        "--with-header",
        "--tag-file",
        f"CB:{cell_list_path}",
        "-o",
        out_bam,
        args.bam_in,
    ]
    run_cmd(cmd_view)

    cmd_index = ["samtools", "index", "-@", str(args.threads), out_bam]
    run_cmd(cmd_index)

    return out_bam


def run_star_and_bws(args, filt_bam, sample_name):
    """
    Run STAR in inputAlignmentsFromBAM mode and convert bedGraph → bigWig.
    Produces stranded and unstranded bigWigs with chosen normalization.
    """
    common = [
        "STAR",
        "--runMode",
        "inputAlignmentsFromBAM",
        "--runThreadN",
        str(args.threads),
        "--genomeDir",
        args.genome_dir,
        "--inputBAMfile",
        filt_bam,
        "--outWigType",
        "bedGraph",
        "--outWigNorm",
        args.wig_norm,
        "--outWigReferencesPrefix",
        "chr",
    ]

    stranded_prefix = os.path.join(args.outdir, f"{sample_name}.stranded_")
    unstranded_prefix = os.path.join(args.outdir, f"{sample_name}.unstranded_")

    # Stranded
    cmd_stranded = common + [
        "--outWigStrand",
        "Stranded",
        "--outFileNamePrefix",
        stranded_prefix,
    ]
    run_cmd(cmd_stranded)

    # Unstranded
    cmd_unstranded = common + [
        "--outWigStrand",
        "Unstranded",
        "--outFileNamePrefix",
        unstranded_prefix,
    ]
    run_cmd(cmd_unstranded)

    # bedGraph → bigWig
    def convert_prefix_to_bw(prefix):
        # STAR names: <prefix>Signal.Unique.str1.out.bg etc.
        pattern = prefix + "Signal.Unique.str" + "*.bg"
        for bg in glob.glob(pattern):
            print(f">> Converting {bg} to bigWig", file=sys.stderr)

            if bg.endswith(".bg"):
                sorted_bg = bg[:-3] + ".sorted.bg"
            else:
                sorted_bg = bg + ".sorted.bg"

            bw = sorted_bg.replace(".sorted.bg", ".bw")

            # sort bedGraph
            with open(sorted_bg, "w") as out_f:
                run_cmd(["sort", "-k1,1", "-k2,2n", bg], stdout=out_f)

            # bedGraphToBigWig
            run_cmd(["bedGraphToBigWig", sorted_bg, args.chrom_sizes, bw])

    convert_prefix_to_bw(stranded_prefix)
    convert_prefix_to_bw(unstranded_prefix)

    # cleanup bedGraphs (.bg and .sorted.bg)
    print(">> Cleaning up bedGraph files", file=sys.stderr)
    for bg in glob.glob(os.path.join(args.outdir, "*.bg")):
        try:
            os.remove(bg)
        except OSError as e:
            print(f"Warning: could not remove {bg}: {e}", file=sys.stderr)


# ----------------- main ----------------- #

def main():
    args = parse_args()
    sample_name = infer_sample_name(args)

    cell_list_path, n_cells = write_all_cells_from_h5ad(args, sample_name)
    print(
        f">> {n_cells} cell barcodes written to: {cell_list_path}",
        file=sys.stderr,
    )

    filt_bam = subset_bam(args, cell_list_path, sample_name)
    print(f">> Filtered BAM: {filt_bam}", file=sys.stderr)

    run_star_and_bws(args, filt_bam, sample_name)
    print(">> Done.", file=sys.stderr)


if __name__ == "__main__":
    main()





python pl.py \
  --h5ad ../../ScRNA_Natalya_Post.h5ad \
  --bam-in 10x_ScRNA.bam \
  --outdir RNA_tracks_Natalya \
  --sample-name 10x_ScRNA \
  --genome-dir /mnt/dataFast/ahrmad/GRCh38_TrES/star \
  --chrom-sizes /mnt/dataFast/ahrmad/hg38.chrom.sizes \
  --threads 64


##########################################
##########################################
# GENERATE BW FOR RNA TrES
##########################################
##########################################

#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import glob

import scanpy as sc


def run_cmd(cmd, stdout=None):
    """
    Helper to run a shell command and fail loudly if it breaks.
    """
    print(">> CMD:", " ".join(cmd), file=sys.stderr)
    try:
        subprocess.run(cmd, check=True, stdout=stdout)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: command failed with exit code {e.returncode}", file=sys.stderr)
        sys.exit(e.returncode)


def parse_args():
    ap = argparse.ArgumentParser(
        description="Filter RNA h5ad, subset BAM and generate stranded/unstranded bigWigs with STAR."
    )

    # I/O
    ap.add_argument("--h5ad", required=True, help="Input RNA .h5ad file")
    ap.add_argument("--bam-in", required=True, help="Input UMI-collapsed BAM")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument(
        "--sample-name",
        required=False,
        default=None,
        help=(
            "Sample name (used in file names). If not given, "
            "it will be inferred from the h5ad basename."
        ),
    )

    # External tools / references
    ap.add_argument("--genome-dir", required=True, help="STAR genomeDir")
    ap.add_argument(
        "--chrom-sizes",
        required=True,
        help="chrom.sizes file for bedGraphToBigWig",
    )
    ap.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Threads for samtools & STAR",
    )

    # QC thresholds (defaults from your example)
    ap.add_argument("--rna-min-umi", type=int, default=1000)
    ap.add_argument("--rna-max-umi", type=int, default=25000)
    ap.add_argument("--rna-max-mt", type=float, default=10.0)      # %
    ap.add_argument("--rna-max-hb", type=float, default=5.0)       # %
    ap.add_argument("--rna-max-top50", type=float, default=30.0)   # %

    # Column names in .obs (override if your AnnData is different)
    ap.add_argument("--obs-total-counts", default="total_counts")
    ap.add_argument("--obs-pct-mt", default="pct_counts_mt")
    ap.add_argument("--obs-pct-hb", default="pct_counts_hb")
    ap.add_argument("--obs-pct-top50", default="pct_counts_in_top_50_genes")

    # Explicit cell list path (optional)
    ap.add_argument(
        "--cell-list",
        default=None,
        help=(
            "Optional path for output list of filtered (rotated) cell names "
            "(one per line). If not set, a default is used in outdir based "
            "on the sample name."
        ),
    )

    # Normalization for STAR output (default RPM for cross-sample correlation)
    ap.add_argument(
        "--wig-norm",
        default="RPM",
        choices=["None", "RPM", "RPKM"],
        help=(
            "Normalization for STAR --outWigNorm. Default: RPM, recommended for "
            "between-sample correlations."
        ),
    )

    return ap.parse_args()


def infer_sample_name(args):
    if args.sample_name is not None:
        return args.sample_name

    # Use h5ad basename without extension as sample name
    base = os.path.basename(args.h5ad)
    if base.endswith(".h5ad"):
        base = base[:-5]
    print(f">> Inferring sample name from h5ad: {base}", file=sys.stderr)
    return base


def rotate_barcode(bc: str) -> str:
    """
    Take the first 4 characters of the barcode and put them at the end.
    """
    return bc[4:] + bc[:4]


def filter_h5ad_and_write_cells(args, sample_name):
    """
    Filter the AnnData by QC metrics, rotate barcodes, and write the passing cell IDs
    to a text file. Returns (cell_list_path, n_pass).
    """
    os.makedirs(args.outdir, exist_ok=True)

    print(f">> Reading h5ad: {args.h5ad}", file=sys.stderr)
    adata = sc.read_h5ad(args.h5ad)
    obs = adata.obs

    needed_cols = [
        args.obs_total_counts,
        args.obs_pct_mt,
        args.obs_pct_hb,
        args.obs_pct_top50,
    ]
    for col in needed_cols:
        if col not in obs.columns:
            print(
                f"ERROR: required column '{col}' not found in adata.obs",
                file=sys.stderr,
            )
            print(f"Available columns: {list(obs.columns)}", file=sys.stderr)
            sys.exit(1)

    total = obs[args.obs_total_counts]
    mt = obs[args.obs_pct_mt]
    hb = obs[args.obs_pct_hb]
    top50 = obs[args.obs_pct_top50]

    mask = (
        (total >= args.rna_min_umi)
        & (total <= args.rna_max_umi)
        & (mt <= args.rna_max_mt)
        & (hb <= args.rna_max_hb)
        & (top50 <= args.rna_max_top50)
    )

    passing_cells = obs.index[mask].tolist()
    n_pass = len(passing_cells)
    n_total = adata.n_obs

    print(f">> Cells passing filters: {n_pass} / {n_total}", file=sys.stderr)

    if n_pass == 0:
        print("ERROR: no cells passed the filters. Aborting.", file=sys.stderr)
        sys.exit(1)

    # Rotate barcodes: take last 4 characters and move them to the front
    rotated_cells = [rotate_barcode(bc) for bc in passing_cells]

    # Where to write the rotated barcodes list
    if args.cell_list is not None:
        cell_list_path = args.cell_list
    else:
        # No Sc_CMD prefix anymore; just use sample_name directly
        cell_list_path = os.path.join(
            args.outdir,
            f"{sample_name}_rna_CountFiltered_cell_names.lst",
        )

    print(f">> Writing rotated cell list to: {cell_list_path}", file=sys.stderr)
    with open(cell_list_path, "w") as f:
        for bc in rotated_cells:
            f.write(f"{bc}\n")

    return cell_list_path, n_pass


def subset_bam(args, cell_list_path, sample_name):
    """
    Use samtools to subset BAM based on read groups / (rotated) cell IDs.
    Returns path to filtered BAM.
    """
    # Use sample_name in the BAM name, but otherwise keep naming simple
    out_bam = os.path.join(args.outdir, f"{sample_name}_rna_Filt.bam")

    cmd_view = [
        "samtools",
        "view",
        "-@",
        str(args.threads),
        "--with-header",
        "--read-group-file",
        cell_list_path,
        "-o",
        out_bam,
        args.bam_in,
    ]
    run_cmd(cmd_view)

    cmd_index = ["samtools", "index", "-@", str(args.threads), out_bam]
    run_cmd(cmd_index)

    return out_bam


def run_star_and_bws(args, filt_bam, sample_name):
    """
    Run STAR in inputAlignmentsFromBAM mode and convert bedGraph → bigWig.
    Uses RPM by default, which is suitable for between-sample correlation.
    """
    # STAR common options
    common = [
        "STAR",
        "--runMode",
        "inputAlignmentsFromBAM",
        "--runThreadN",
        str(args.threads),
        "--genomeDir",
        args.genome_dir,
        "--inputBAMfile",
        filt_bam,
        "--outWigType",
        "bedGraph",
        "--outWigNorm",
        args.wig_norm,
        "--outWigReferencesPrefix",
        "chr",
    ]

    stranded_prefix = os.path.join(args.outdir, f"{sample_name}.stranded_")
    unstranded_prefix = os.path.join(args.outdir, f"{sample_name}.unstranded_")

    # Stranded
    cmd_stranded = common + [
        "--outWigStrand",
        "Stranded",
        "--outFileNamePrefix",
        stranded_prefix,
    ]
    run_cmd(cmd_stranded)

    # Unstranded
    cmd_unstranded = common + [
        "--outWigStrand",
        "Unstranded",
        "--outFileNamePrefix",
        unstranded_prefix,
    ]
    run_cmd(cmd_unstranded)

    # bedGraph → bigWig
    def convert_prefix_to_bw(prefix):
        # Matches patterns like:
        #   {outdir}/{sample}.stranded_Signal.Unique.str1.bg
        pattern = prefix + "Signal.Unique.str" + "*.bg"
        for bg in glob.glob(pattern):
            print(f">> Converting {bg} to bigWig", file=sys.stderr)

            if bg.endswith(".bg"):
                sorted_bg = bg[:-3] + ".sorted.bg"
                bw = bg[:-3] + ".bw"
            else:
                sorted_bg = bg + ".sorted.bg"
                bw = bg + ".bw"

            # sort bedGraph
            with open(sorted_bg, "w") as out_f:
                run_cmd(["sort", "-k1,1", "-k2,2n", bg], stdout=out_f)

            # bedGraphToBigWig
            run_cmd(["bedGraphToBigWig", sorted_bg, args.chrom_sizes, bw])

    convert_prefix_to_bw(stranded_prefix)
    convert_prefix_to_bw(unstranded_prefix)

    # cleanup bedGraphs (.bg and .sorted.bg)
    print(">> Cleaning up bedGraph files", file=sys.stderr)
    for bg in glob.glob(os.path.join(args.outdir, "*.bg")):
        try:
            os.remove(bg)
        except OSError as e:
            print(f"Warning: could not remove {bg}: {e}", file=sys.stderr)


def main():
    args = parse_args()
    sample_name = infer_sample_name(args)

    cell_list_path, n_pass = filter_h5ad_and_write_cells(args, sample_name)
    print(
        f">> {n_pass} cells passed filters; cell list: {cell_list_path}",
        file=sys.stderr,
    )

    filt_bam = subset_bam(args, cell_list_path, sample_name)
    print(f">> Filtered BAM: {filt_bam}", file=sys.stderr)

    run_star_and_bws(args, filt_bam, sample_name)
    print(">> Done.", file=sys.stderr)


if __name__ == "__main__":
    main()



python make_rna_bw.py \
  --h5ad /mnt/dataFast/ahrmad/ScH2R_202409/FullBench/TrES/H2RRNA_Human.h5ad \
  --bam-in /mnt/dataFast/ahrmad/ScH2R_202409/RNA2/H2RRNA_Human.filtered_cells.bam \
  --outdir rna_bw_TrES \
  --sample-name TrES \
  --genome-dir /mnt/dataFast/ahrmad/GRCh38_TrES/star \
  --chrom-sizes /mnt/dataFast/ahrmad/hg38.chrom.sizes \
  --threads 64 \
  --rna-min-umi 1000 \
  --rna-max-umi 25000 \
  --rna-max-mt 10 \
  --rna-max-hb 5 \
  --rna-max-top50 30


##########################################
##########################################
# GENERATE BW FOR ENCODE RNA
##########################################
##########################################

#!/usr/bin/env python3
import os
import sys
import subprocess
import glob

# ==================== CONFIG ====================

THREADS = 64
GENOME_DIR = "/mnt/dataFast/ahrmad/GRCh38_TrES/star"
CHROM_SIZES = "/mnt/dataFast/ahrmad/hg38.chrom.sizes"
OUTDIR_BASE = "ENCODE_Bulk_RNA"

# ENCODE bulk RNA replicates
REPLICATES = [
    {
        "sample_name": "ENCODE_Bulk_rep1",
        "bam_name": "ENCODE_Bulk_RNA_rep1.bam",
        "url": "https://www.encodeproject.org/files/ENCFF325HWX/@@download/ENCFF325HWX.bam",
    },
    {
        "sample_name": "ENCODE_Bulk_rep2",
        "bam_name": "ENCODE_Bulk_RNA_rep2.bam",
        "url": "https://www.encodeproject.org/files/ENCFF936COB/@@download/ENCFF936COB.bam",
    },
]

# wig normalization: RPM is good for correlations
WIG_NORM = "RPM"   # can be "None", "RPM", "RPKM"


# ==================== HELPERS ====================

def run_cmd(cmd, stdout=None):
    """Run command and fail loudly if it breaks."""
    print("Running:", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True, stdout=stdout)


def download_bam_if_needed(bam_path, url):
    """Download BAM if it does not already exist."""
    if os.path.exists(bam_path):
        print(f"[SKIP] BAM already exists: {bam_path}", file=sys.stderr)
        return

    os.makedirs(os.path.dirname(bam_path), exist_ok=True)
    print(f"[DOWNLOAD] {bam_path} from {url}", file=sys.stderr)
    run_cmd(["wget", "-q", "--show-progress", "-O", bam_path, url])
    print(f"[OK] Downloaded {bam_path}", file=sys.stderr)


def index_bam_if_needed(bam_path):
    """Index BAM with samtools if index is missing."""
    bai_path = bam_path + ".bai"
    if os.path.exists(bai_path):
        print(f"[SKIP] BAM index already exists: {bai_path}", file=sys.stderr)
        return
    print(f"[INDEX] {bam_path}", file=sys.stderr)
    run_cmd(["samtools", "index", "-@", str(THREADS), bam_path])
    print(f"[OK] Indexed {bam_path}", file=sys.stderr)


def bigwigs_exist(outdir, sample_name):
    """
    Check if main bigWigs already exist for this sample.
    If stranded str1 + unstranded exist, we consider it done.
    """
    stranded_bw1 = os.path.join(
        outdir, f"{sample_name}.stranded_Signal.Unique.str1.out.bw"
    )
    unstranded_bw1 = os.path.join(
        outdir, f"{sample_name}.unstranded_Signal.Unique.str1.out.bw"
    )

    if os.path.exists(stranded_bw1) and os.path.exists(unstranded_bw1):
        print(f"[SKIP] bigWigs already found for {sample_name}", file=sys.stderr)
        return True
    return False


def convert_prefix_to_bw(prefix, chrom_sizes):
    """
    Convert STAR bedGraphs to bigWig for a given prefix.

    Only process unsorted *.out.bg files:
      <prefix>Signal.Unique.strX.out.bg

    For each bg:
      sort -> .sorted.bg
      bedGraphToBigWig -> .bw
      then remove both bg + sorted.bg

    Skips empty bedGraphs and already-existing bigWigs.
    """
    pattern = prefix + "Signal.Unique.str*.out.bg"
    bg_files = sorted(glob.glob(pattern))

    if not bg_files:
        print(f"[INFO] No bedGraph files matching {pattern}", file=sys.stderr)
        return

    for bg in bg_files:
        # Derive output names
        if not bg.endswith(".bg"):
            # Should not happen, but guard anyway
            print(f"[WARN] Unexpected bedGraph name (no .bg extension): {bg}", file=sys.stderr)
            continue

        sorted_bg = bg[:-3] + ".sorted.bg"
        bw = bg[:-3] + ".bw"

        # If bw already exists and is non-empty, skip conversion for this bg
        if os.path.exists(bw) and os.path.getsize(bw) > 0:
            print(f"[SKIP] bigWig already exists for {bg}: {bw}", file=sys.stderr)
            continue

        # Skip empty bg (STAR sometimes can output 0-byte files for a strand)
        if os.path.getsize(bg) == 0:
            print(f"[SKIP] Empty bedGraph, not converting: {bg}", file=sys.stderr)
            # Clean up the empty file to avoid future confusion
            try:
                os.remove(bg)
            except OSError:
                pass
            continue

        print(f">> Converting {bg} -> {bw}", file=sys.stderr)

        # sort bedGraph
        with open(sorted_bg, "w") as out_f:
            run_cmd(["sort", "-k1,1", "-k2,2n", bg], stdout=out_f)

        # bedGraphToBigWig
        run_cmd(["bedGraphToBigWig", sorted_bg, chrom_sizes, bw])

        # cleanup intermediates
        for f in (bg, sorted_bg):
            try:
                os.remove(f)
            except OSError as e:
                print(f"[WARN] could not remove {f}: {e}", file=sys.stderr)


def run_star_and_bws(bam_path, sample_name):
    """
    Run STAR in inputAlignmentsFromBAM mode and convert bedGraph → bigWig.

    Generates stranded and unstranded RPM-normalized tracks, with idempotent behavior.
    """
    outdir = os.path.join(OUTDIR_BASE, sample_name)
    os.makedirs(outdir, exist_ok=True)

    # If bigWigs already exist, skip STAR entirely
    if bigwigs_exist(outdir, sample_name):
        return

    print(f"[STAR] Processing sample {sample_name}", file=sys.stderr)
    print(f"       BAM: {bam_path}", file=sys.stderr)
    print(f"       Outdir: {outdir}", file=sys.stderr)

    common = [
        "STAR",
        "--runMode", "inputAlignmentsFromBAM",
        "--runThreadN", str(THREADS),
        "--genomeDir", GENOME_DIR,
        "--inputBAMfile", bam_path,
        "--outWigType", "bedGraph",
        "--outWigNorm", WIG_NORM,
        "--outWigReferencesPrefix", "chr",
    ]

    stranded_prefix = os.path.join(outdir, f"{sample_name}.stranded_")
    unstranded_prefix = os.path.join(outdir, f"{sample_name}.unstranded_")

    # Stranded
    cmd_stranded = common + [
        "--outWigStrand", "Stranded",
        "--outFileNamePrefix", stranded_prefix,
    ]
    run_cmd(cmd_stranded)

    # Unstranded
    cmd_unstranded = common + [
        "--outWigStrand", "Unstranded",
        "--outFileNamePrefix", unstranded_prefix,
    ]
    run_cmd(cmd_unstranded)

    # Convert bedGraphs to bigWigs
    convert_prefix_to_bw(stranded_prefix, CHROM_SIZES)
    convert_prefix_to_bw(unstranded_prefix, CHROM_SIZES)

    print(f"[OK] Finished STAR + bigWig generation for {sample_name}", file=sys.stderr)


# ==================== MAIN ====================

def main():
    os.makedirs(OUTDIR_BASE, exist_ok=True)

    for rep in REPLICATES:
        sample_name = rep["sample_name"]
        bam_path = os.path.join(OUTDIR_BASE, rep["bam_name"])
        url = rep["url"]

        print("\n=== Processing replicate:", sample_name, "===\n", file=sys.stderr)

        # 1) Download BAM if needed
        download_bam_if_needed(bam_path, url)

        # 2) Index BAM if needed
        index_bam_if_needed(bam_path)

        # 3) Run STAR + bigWig generation (idempotent)
        run_star_and_bws(bam_path, sample_name)

    # Summary
    print("\n=== Processing Summary ===", file=sys.stderr)
    print("Downloaded/used BAM files:", file=sys.stderr)
    for rep in REPLICATES:
        print(f"  - {os.path.join(OUTDIR_BASE, rep['bam_name'])}", file=sys.stderr)
    print("\nOutput directories:", file=sys.stderr)
    for rep in REPLICATES:
        outdir = os.path.join(OUTDIR_BASE, rep["sample_name"])
        print(f"  - {outdir}", file=sys.stderr)
        if os.path.isdir(outdir):
            bws = sorted(glob.glob(os.path.join(outdir, "*.bw")))
            for bw in bws:
                print(f"    - {os.path.basename(bw)}", file=sys.stderr)


if __name__ == "__main__":
    main()

##########################################
##########################################
# GENERATE BW FOR RNA K9 TrES (2025/10)
##########################################
##########################################
#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import glob


def run_cmd(cmd, stdout=None):
    """Run a shell command and fail loudly if it breaks."""
    print("Running:", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True, stdout=stdout)


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate stranded/unstranded bigWigs from single-cell RNA BAM (RG-tagged)."
    )
    p.add_argument("--bam-in", required=True, help="Input BAM (e.g. Sc_K9r_UMI.bam)")
    p.add_argument(
        "--cell-list",
        required=True,
        help="Text file with one RG ID per line (e.g. Sc_K9_rna_CountFiltered_cell_names.lst)",
    )
    p.add_argument(
        "--outdir",
        required=True,
        help="Output directory for filtered BAM + bigWigs",
    )
    p.add_argument(
        "--sample-name",
        required=False,
        default=None,
        help="Sample name for outputs (default: basename of BAM without .bam)",
    )
    p.add_argument(
        "--genome-dir",
        required=True,
        help="STAR genomeDir",
    )
    p.add_argument(
        "--chrom-sizes",
        required=True,
        help="chrom.sizes file for bedGraphToBigWig",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Threads for samtools & STAR (default: 8)",
    )
    p.add_argument(
        "--wig-norm",
        default="RPM",
        choices=["None", "RPM", "RPKM"],
        help="STAR --outWigNorm (default: RPM, good for correlation)",
    )
    return p.parse_args()


def infer_sample_name(args):
    if args.sample_name:
        return args.sample_name
    base = os.path.basename(args.bam_in)
    if base.endswith(".bam"):
        base = base[:-4]
    return base


def filter_bam_by_rg(args, sample_name):
    """
    Filter BAM using RG IDs from cell-list with samtools:
      samtools view -R cell_list ...
    Returns path to filtered BAM.
    """
    os.makedirs(args.outdir, exist_ok=True)

    out_bam = os.path.join(args.outdir, f"{sample_name}_Filt.bam")
    out_bai = out_bam + ".bai"

    # If filtered BAM already exists, reuse it
    if os.path.exists(out_bam):
        print(f"[SKIP] Filtered BAM already exists: {out_bam}", file=sys.stderr)
        # Ensure it is indexed
        if not os.path.exists(out_bai):
            print(f"[INDEX] Indexing existing filtered BAM: {out_bam}", file=sys.stderr)
            run_cmd(["samtools", "index", "-@", str(args.threads), out_bam])
        return out_bam

    print(f"[FILTER] Creating filtered BAM: {out_bam}", file=sys.stderr)
    run_cmd(
        [
            "samtools",
            "view",
            "-@",
            str(args.threads),
            "--with-header",
            "-R",
            args.cell_list,
            "-o",
            out_bam,
            args.bam_in,
        ]
    )

    print(f"[INDEX] Indexing filtered BAM: {out_bam}", file=sys.stderr)
    run_cmd(["samtools", "index", "-@", str(args.threads), out_bam])

    return out_bam


def bigwigs_exist(outdir, sample_name):
    """Check if main bigWigs already exist for this sample."""
    stranded_bw1 = os.path.join(
        outdir, f"{sample_name}.stranded_Signal.Unique.str1.out.bw"
    )
    unstranded_bw1 = os.path.join(
        outdir, f"{sample_name}.unstranded_Signal.Unique.str1.out.bw"
    )
    return os.path.exists(stranded_bw1) and os.path.exists(unstranded_bw1)


def convert_prefix_to_bw(prefix, chrom_sizes):
    """
    Convert STAR bedGraphs to bigWig for a given prefix.

    Only process unsorted *.out.bg files:
      <prefix>Signal.Unique.strX.out.bg
    """
    pattern = prefix + "Signal.Unique.str*.out.bg"
    bg_files = sorted(glob.glob(pattern))

    if not bg_files:
        print(f"[INFO] No bedGraph files matching {pattern}", file=sys.stderr)
        return

    for bg in bg_files:
        if not bg.endswith(".bg"):
            print(f"[WARN] Unexpected bedGraph name (no .bg): {bg}", file=sys.stderr)
            continue

        sorted_bg = bg[:-3] + ".sorted.bg"
        bw = bg[:-3] + ".bw"

        # Skip if bw already exists and is non-empty
        if os.path.exists(bw) and os.path.getsize(bw) > 0:
            print(f"[SKIP] bigWig already exists for {bg}: {bw}", file=sys.stderr)
            continue

        # Skip empty bedGraphs
        if os.path.getsize(bg) == 0:
            print(f"[SKIP] Empty bedGraph, not converting: {bg}", file=sys.stderr)
            try:
                os.remove(bg)
            except OSError:
                pass
            continue

        print(f">> Converting {bg} -> {bw}", file=sys.stderr)

        # sort bedGraph
        with open(sorted_bg, "w") as out_f:
            run_cmd(["sort", "-k1,1", "-k2,2n", bg], stdout=out_f)

        # bedGraphToBigWig
        run_cmd(["bedGraphToBigWig", sorted_bg, chrom_sizes, bw])

        # cleanup intermediates
        for f in (bg, sorted_bg):
            try:
                os.remove(f)
            except OSError as e:
                print(f"[WARN] could not remove {f}: {e}", file=sys.stderr)


def run_star_and_bws(args, filt_bam, sample_name):
    """
    Run STAR in inputAlignmentsFromBAM mode and convert bedGraph → bigWig.
    """
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # If bigWigs already exist, skip STAR entirely
    if bigwigs_exist(outdir, sample_name):
        print(f"[SKIP] bigWigs already exist for {sample_name}", file=sys.stderr)
        return

    print(f"[STAR] Sample: {sample_name}", file=sys.stderr)
    print(f"       BAM: {filt_bam}", file=sys.stderr)
    print(f"       Outdir: {outdir}", file=sys.stderr)

    common = [
        "STAR",
        "--runMode", "inputAlignmentsFromBAM",
        "--runThreadN", str(args.threads),
        "--genomeDir", args.genome_dir,
        "--inputBAMfile", filt_bam,
        "--outWigType", "bedGraph",
        "--outWigNorm", args.wig_norm,
        "--outWigReferencesPrefix", "chr",
    ]

    stranded_prefix = os.path.join(outdir, f"{sample_name}.stranded_")
    unstranded_prefix = os.path.join(outdir, f"{sample_name}.unstranded_")

    # Stranded
    cmd_stranded = common + [
        "--outWigStrand", "Stranded",
        "--outFileNamePrefix", stranded_prefix,
    ]
    run_cmd(cmd_stranded)

    # Unstranded
    cmd_unstranded = common + [
        "--outWigStrand", "Unstranded",
        "--outFileNamePrefix", unstranded_prefix,
    ]
    run_cmd(cmd_unstranded)

    # Convert bedGraphs to bigWigs
    convert_prefix_to_bw(stranded_prefix, args.chrom_sizes)
    convert_prefix_to_bw(unstranded_prefix, args.chrom_sizes)

    print(f"[OK] Finished STAR + bigWig generation for {sample_name}", file=sys.stderr)


def main():
    args = parse_args()
    sample_name = infer_sample_name(args)

    # 1) Filter BAM by RG IDs
    filt_bam = filter_bam_by_rg(args, sample_name)

    # 2) Run STAR + convert to bigWigs
    run_star_and_bws(args, filt_bam, sample_name)

    print(">> Done.", file=sys.stderr)


if __name__ == "__main__":
    main()



python ll.py \
  --bam-in Sc_K9r_UMI.bam \
  --cell-list Sc_K9_rna_CountFiltered_cell_names.lst \
  --outdir Sc_K9r_BWs \
  --sample-name Sc_K9r \
  --genome-dir /mnt/dataFast/ahrmad/GRCh38_TrES/star \
  --chrom-sizes /mnt/dataFast/ahrmad/hg38.chrom.sizes \
  --threads 64
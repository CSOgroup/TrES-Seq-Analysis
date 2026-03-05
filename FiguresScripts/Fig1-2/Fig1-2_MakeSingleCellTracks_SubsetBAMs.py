#!/usr/bin/env python3

import argparse
import os
import sys
import pysam

# =========================
# Part 1: FX tagging logic
# =========================

# Priority of tags to use as cell identifier.
# If a read has CB, we use that; otherwise we fall back to RG.
CELL_TAG_PRIORITY = ["CB", "RG"]


def get_cell_identifier(read):
    """
    Return the cell identifier string for a read, or None if not found.
    """
    for tag in CELL_TAG_PRIORITY:
        try:
            val = read.get_tag(tag)
        except KeyError:
            continue
        if val is not None and val != "":
            return val
    return None


def collect_cells_for_bam(bam_path):
    """
    Scan a single BAM and collect all unique cell identifiers (CB/RG).
    Returns a sorted list of cell IDs.
    """
    sys.stderr.write(f"[INFO] Scanning cells in BAM: {bam_path}\n")
    cells = set()
    with pysam.AlignmentFile(bam_path, "rb") as infile:
        for read in infile:
            cid = get_cell_identifier(read)
            if cid is not None:
                cells.add(cid)
    sorted_cells = sorted(cells)
    sys.stderr.write(f"[INFO]   Found {len(sorted_cells)} unique cells in {bam_path}\n")
    return sorted_cells


def build_fx_mapping_for_bam(sorted_cells):
    """
    For one BAM, build mapping:
        cell_identifier -> FX integer ID (1..N)
    """
    return {cid: idx + 1 for idx, cid in enumerate(sorted_cells)}


def write_fx_mapping_tsv_for_bam(bam_path, sorted_cells, cell_to_fx, suffix=".FX.map.tsv"):
    """
    Write BAM-specific mapping:
        FX_ID <tab> CELL_IDENTIFIER
    """
    base = os.path.splitext(bam_path)[0]
    out_path = base + suffix
    sys.stderr.write(f"[INFO] Writing FX mapping for {bam_path} to: {out_path}\n")
    with open(out_path, "w") as out:
        out.write("#FX_ID\tCELL_IDENTIFIER\n")
        for cid in sorted_cells:
            out.write(f"{cell_to_fx[cid]}\t{cid}\n")


def add_fx_to_bam(in_bam, out_bam, cell_to_fx):
    """
    Read a BAM and write a new BAM where we add FX tag
    per read, using the BAM-specific cell_to_fx mapping.
    """
    sys.stderr.write(f"[INFO] Adding FX tags: {in_bam} -> {out_bam}\n")
    with pysam.AlignmentFile(in_bam, "rb") as infile:
        with pysam.AlignmentFile(out_bam, "wb", header=infile.header) as outfile:
            for read in infile:
                cid = get_cell_identifier(read)
                if cid is not None:
                    fx_id = cell_to_fx.get(cid)
                    if fx_id is not None:
                        # FX as a 32-bit integer
                        read.set_tag("FX", fx_id, value_type="i")
                outfile.write(read)

    sys.stderr.write(f"[INFO] Indexing {out_bam}\n")
    pysam.index(out_bam)


def run_add_fx(args):
    """
    This function preserves the original 'main' logic of the FX-tagging script,
    but takes an argparse 'args' object instead of parsing sys.argv itself.
    """
    for bam_path in args.bams:
        if not os.path.exists(bam_path):
            sys.stderr.write(f"[ERROR] BAM not found: {bam_path}\n")
            sys.exit(1)

    for bam_path in args.bams:
        # 1) Collect cells for this BAM
        sorted_cells = collect_cells_for_bam(bam_path)

        # 2) Build BAM-specific FX mapping
        cell_to_fx = build_fx_mapping_for_bam(sorted_cells)

        # 3) Write mapping TSV for this BAM
        write_fx_mapping_tsv_for_bam(
            bam_path,
            sorted_cells,
            cell_to_fx,
            suffix=args.map_suffix
        )

        # 4) Write FX-tagged BAM
        base, ext = os.path.splitext(bam_path)
        if ext.lower() == ".bam":
            out_bam = base + args.suffix
        else:
            out_bam = bam_path + args.suffix

        add_fx_to_bam(bam_path, out_bam, cell_to_fx)

    sys.stderr.write("[INFO] Done.\n")


# =========================
# Part 2: Region subsetting
# =========================

# Regions of interest: (name, chrom, start, end) 1-based inclusive
REGIONS = [
    ("MYC",   "chr8", 127723000, 127753000),
    ("BCL6",  "chr3", 187713000, 187755000),
    ("IRF4",  "chr6",      336448,     436683),
    ("PRDM1", "chr6", 106062900, 106133200),
]


def subset_region(bam_path, region_name, chrom, start_1based, end_1based, suffix_template):
    """
    Write a BAM containing only reads overlapping the given region.
    No filtering on FX or anything else.
    """
    base, ext = os.path.splitext(bam_path)
    if ext.lower() != ".bam":
        base = bam_path  # weird extension, just use full name as base

    suffix = suffix_template.format(gene=region_name)
    out_bam_path = base + suffix

    sys.stderr.write(
        f"[INFO] Subsetting {bam_path} to region {region_name} "
        f"{chrom}:{start_1based}-{end_1based} -> {out_bam_path}\n"
    )

    n_in = 0
    n_out = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam_in:
        with pysam.AlignmentFile(out_bam_path, "wb", header=bam_in.header) as bam_out:
            for read in bam_in.fetch(chrom, start_1based - 1, end_1based):
                n_in += 1
                bam_out.write(read)
                n_out += 1

    sys.stderr.write(
        f"[INFO]   Reads overlapping region written: {n_out}\n"
    )

    sys.stderr.write(f"[INFO] Indexing {out_bam_path}\n")
    pysam.index(out_bam_path)


def run_subset_regions(args):
    """
    This function preserves the original 'main' logic of the region-subsetting script,
    but takes an argparse 'args' object instead of parsing sys.argv itself.
    """
    for bam_path in args.bams:
        if not os.path.exists(bam_path):
            sys.stderr.write(f"[ERROR] BAM not found: {bam_path}\n")
            sys.exit(1)

    for bam_path in args.bams:
        sys.stderr.write(f"[INFO] Processing BAM: {bam_path}\n")
        for region_name, chrom, start, end in REGIONS:
            subset_region(
                bam_path=bam_path,
                region_name=region_name,
                chrom=chrom,
                start_1based=start,
                end_1based=end,
                suffix_template=args.suffix_template,
            )

    sys.stderr.write("[INFO] Done.\n")


# =========================
# Combined CLI entry point
# =========================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Combined tool to:\n"
            "  1) Add an FX tag to BAMs for per-BAM cell IDs, and/or\n"
            "  2) Subset BAMs by predefined regions (MYC, BCL6, IRF4, PRDM1).\n\n"
            "Use subcommands 'add-fx' or 'subset-regions'."
        )
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand: add-fx (original first script)
    parser_fx = subparsers.add_parser(
        "add-fx",
        help="Add FX tags to BAMs for per-BAM cell IDs."
    )
    parser_fx.add_argument(
        "bams",
        nargs="+",
        help="Input BAM files to process."
    )
    parser_fx.add_argument(
        "--suffix",
        default=".FX.bam",
        help="Suffix for output BAMs (default: .FX.bam)."
    )
    parser_fx.add_argument(
        "--map-suffix",
        default=".FX.map.tsv",
        help="Suffix for per-BAM FX mapping tables (default: .FX.map.tsv)."
    )
    parser_fx.set_defaults(func=run_add_fx)

    # Subcommand: subset-regions (original second script)
    parser_sub = subparsers.add_parser(
        "subset-regions",
        help=(
            "Subset BAMs by predefined regions (MYC, BCL6, IRF4, PRDM1), "
            "without any cell/FX filtering."
        )
    )
    parser_sub.add_argument(
        "bams",
        nargs="+",
        help="Input BAM files (ideally already with FX tags, but not required)."
    )
    parser_sub.add_argument(
        "--suffix-template",
        default=".{gene}.bam",
        help=(
            "Suffix template for output BAMs. "
            "Use {gene} as placeholder for region name. "
            "Default: .{gene}.bam (e.g. sample.FX.MYC.bam)"
        ),
    )
    parser_sub.set_defaults(func=run_subset_regions)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

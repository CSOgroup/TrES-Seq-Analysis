import os
import re
import sys
import glob
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, plot, from_indicators
import polars as pl

pl.Config.set_fmt_str_lengths(100)

# ---------------------------
# Helpers tied to your new runs
# ---------------------------

SUBSAMPLE_BCS = ["GCAT", "CGAT", "GACT", "AGCT", "CAGT", "ACGT", "GCTA", "CGTA"]

def sb_slice_for_subsample(sample_name: str):
    """
    Given a sample name (e.g., 'Sc_ID_DMSO', 'Sc_SKD_K422_96H_Treatment', 'Sc_TKr_K422_8D_Recovery'),
    return the 3-barcode slice that Split_Reads2.codon used for that subsample.
    """
    name = sample_name
    # ISA
    if "Human" in name:
        return SUBSAMPLE_BCS[:-1]
    if "Mouse" in name:
        return SUBSAMPLE_BCS[-1:]

    # Default: no filtering (should not happen with your current cohorts)
    return None

def tag_records_for_sample(sample_name: str, out_folder: str) -> str:
    """
    Return the merged Tag_Records file path that corresponds to the sample prefix:

    """
    if sample_name.startswith("H2RD"):
        return f"{out_folder}/Tag_Records_H2RD.tsv.gz"
    if sample_name.startswith("H2RR"):
        return f"{out_folder}/Tag_Records_H2RR.tsv.gz"
    raise ValueError(f"Unrecognized sample prefix for Tag_Records: {sample_name}")

def good_bc_path_for_sample(sample_name: str, modality: str, out_folder: str, mark_for_dna="H3K27ac") -> str:
    """
    Build the reads-per-barcode file path used in the QC step.
    For DNA, we use mark_for_dna='H3K27ac' as the anchor and then also look up the me3 file inside make_upsetplot.
    For RNA, use ALL as in your pipeline.
    """
    if modality == "DNA":
        return f"{out_folder}/{sample_name}_{mark_for_dna}_ProperPairedMapped_reads_per_barcode.tsv"
    elif modality == "RNA":
        return f"{out_folder}/{sample_name}_ProperPairedMapped_reads_per_barcode.tsv"
    else:
        raise ValueError("modality must be 'DNA' or 'RNA'")

def bam_path_for_sample(sample_name: str, modality: str, out_folder: str, mark_for_dna="H3K27ac") -> str:
    """
    Build the main BAM path for plotting.
    - DNA: use MarkDuplicates output (MarkedDup.bam) anchored on H3K27ac; we’ll join me3 inside make_upsetplot.
    - RNA: use _ALL.bam (and also join the UMI bam inside make_upsetplot).
    """
    if modality == "DNA":
        return f"{out_folder}/{sample_name}_{mark_for_dna}_MarkedDup.bam"
    elif modality == "RNA":
        return f"{out_folder}/{sample_name}_ALL.bam"
    else:
        raise ValueError("modality must be 'DNA' or 'RNA'")

# ---------------------------
# Core functions (unchanged logic except filtering)
# ---------------------------

def create_alignment_df(bam_file_path, RNAinGenes=False):
    unique_query = set()
    primary_reads_data = []

    with pysam.AlignmentFile(bam_file_path, "rb", threads=12) as bamfile:
        for count, read in enumerate(bamfile.fetch()):
            if (not read.is_unmapped and read.is_paired and read.query_name not in unique_query):
                unique_query.add(read.query_name)
                if not RNAinGenes:
                    primary_reads_data.append([read.query_name, not read.is_unmapped, not read.is_duplicate])
                else:
                    xs_tag = read.get_tag("XS")
                    is_xs_valid = False if xs_tag == "Unassigned_NoFeatures" else True
                    primary_reads_data.append([read.query_name, not read.is_unmapped, not read.is_duplicate, is_xs_valid])

            if (count + 1) % 100000 == 0:
                sys.stdout.write(f"\rProcessed {(count + 1) / 1e6:.1f} million queries...")
                sys.stdout.flush()

    sys.stdout.write(f"\rDone processing {(count + 1)} queries.\n")
    sys.stdout.flush()

    sys.stdout.write('Preparing Polars Query DataFrame...\n')
    sys.stdout.flush()

    if RNAinGenes:
        schema = [('Read_Name_BAM', pl.String),
                  ('Mapped', pl.Boolean),
                  ('Unique', pl.Boolean),
                  ('inGenes', pl.Boolean)]
    else:
        schema = [('Read_Name_BAM', pl.String),
                  ('Mapped', pl.Boolean),
                  ('Unique', pl.Boolean)]

    reads_df = pl.DataFrame(primary_reads_data, schema=schema, orient="row")

    sys.stdout.write('Integrating Alignment Results...\n')
    sys.stdout.flush()

    return reads_df


def make_upsetplot(
        Tag_records_path,
        bam_file_path,
        good_bc_path,
        select_modality,
        sample,
        out_folder
):
    print(f'Load Tag Records: {Tag_records_path}')
    df = pl.read_csv(Tag_records_path, separator='\t', has_header=False)

    print('Loading Done. Cleaning...')
    first_row = list(df.row(0))[1:]
    column_prefixes = [s[:2] for s in list(first_row)]
    df.columns = ["Record_Name"] + column_prefixes

    df = df.with_columns((
        pl.col(col).str.slice(5)
        for col in column_prefixes
    ))
    df = df.drop(["L1", "L2", "L3"])  # remove ligation columns

    # Determine SB barcodes for this sample (use both 'A' and 'T' prefixes; your merged Tag_Records contain both)
    sb_slice = sb_slice_for_subsample(sample)
    if sb_slice is not None:
        bcs_with_prefix = [p + bc for bc in sb_slice for p in ("", "")]
        df = df.filter(pl.col("SB").is_in(set(bcs_with_prefix)))

    CheckinGenes = (select_modality == "RNA")

    if 'DNA' in select_modality:
        df = df.rename({"MO": "Modality"})
        filename = f'{sample}_UpsetPlot.pdf'
        CheckinGenes = False

    print('Cleaning Done.\n')
    print(df.head())

    cntBC = pl.read_csv(good_bc_path, separator='\t', has_header=False, new_columns=["MappedNoFilt", "bc"])
    cntBC = cntBC.with_columns(pl.lit(True).alias("MappedNoFilt"))
    cntBC = cntBC.unique(subset=["bc"])

    df = df.join(
        cntBC,
        left_on='CB',
        right_on='bc',
        how='left',
        validate='m:1',
        coalesce=True
    )

    print(f"\nCheck Alignment Status {bam_file_path}...")
    df = df.join(
        create_alignment_df(bam_file_path, RNAinGenes=CheckinGenes),
        left_on='Record_Name',
        right_on='Read_Name_BAM',
        how='left',
        validate='1:1',
        coalesce=True
    )

    if select_modality == 'RNA':
        filename = f'{sample}_UpsetPlot.pdf'
        print(f"\nCheck UMI Status {bam_file_path.replace('_ALL', '_UMI')}...")
        df = df.join(
            create_alignment_df(f"{bam_file_path.replace('_ALL', '_UMI')}"),
            left_on='Record_Name',
            right_on='Read_Name_BAM',
            how='left',
            validate='1:1',
            coalesce=True,
            suffix='_UMI'
        )
        df = df.drop(["Unique", "Unique_UMI"], strict=False)
        df = df.rename({"Mapped_UMI": "UMI"})
        df = df.drop(["Mapped_UMI", "UM"], strict=False)

    if select_modality == 'DNA':
        # Anchor is H3K27ac in bam_file_path; also pull H3K27me3
        print(f'\nCheck Alignment Status {bam_file_path.replace("H3K27ac", "H3K27me3")}...')
        df = df.join(
            create_alignment_df(bam_file_path.replace("H3K27ac", "H3K27me3")),
            left_on='Record_Name',
            right_on='Read_Name_BAM',
            how='left',
            validate='1:1',
            coalesce=True,
            suffix='_me3'
        )
        df = df.with_columns(pl.coalesce(["Mapped", "Mapped_me3"]).alias("Mapped"))
        df = df.with_columns(pl.coalesce(["Unique", "Unique_me3"]).alias("Unique"))
        df = df.drop("Mapped_me3")
        df = df.drop("Unique_me3")

        cntBC_me3 = pl.read_csv(good_bc_path.replace("H3K27ac", "H3K27me3"),
                                separator='\t', has_header=False, new_columns=["MappedNoFilt", "bc"])
        cntBC_me3 = cntBC_me3.with_columns(pl.lit(True).alias("MappedNoFilt"))

        df = df.join(
            cntBC_me3,
            left_on='CB',
            right_on='bc',
            how='left',
            validate='m:1',
            coalesce=True,
            suffix='_me3'
        )
        df = df.with_columns(pl.coalesce(["MappedNoFilt", "MappedNoFilt_me3"]).alias("MappedNoFilt"))
        df = df.drop("MappedNoFilt_me3")

    print("Done.")
    print("Preparing Boolean DF.")

    df = df.with_columns([
        pl.when(pl.col(col) == "NoMatch").then(False).otherwise(pl.col(col)).alias(col)
        for col in df.columns
    ])

    if ('DNA' in select_modality):
        df = df.with_columns(pl.col(['Mapped', 'Unique']).fill_null(False))
    else:
        df = df.with_columns(pl.col(['Mapped']).fill_null(False))

    df = df.with_columns([
        pl.when(pl.col(col) != False).then(True).otherwise(pl.col(col)).alias(col)
        for col in df.columns
    ])

    df = df.cast(pl.UInt8).cast(pl.Boolean)

    df = df.rename({"Mapped": "CB>50"})
    df = df.rename({"MappedNoFilt": "Mapped"})
    df = df.rename(lambda column_name: column_name + ' +')
    df = df.rename({"Record_Name +": "Reads", "SB +": "Sample +"})

    print('Creating Figure')
    pd_df = df.to_pandas()
    fig = plt.figure(figsize=(10, 12))
    plot(
        from_indicators(pd_df.columns, data=pd_df),
        sort_by='-degree',
        show_percentages=True,
        min_subset_size="0.5%",
        facecolor="darkblue",
        fig=fig, element_size=None
    )

    print('Saving Figure')
    plt.suptitle(sample.replace("_", " "), fontsize=16, fontweight='bold')
    fig.savefig(f'{out_folder}/{filename}')
    plt.close()
    return 0


def make_barplots(subdir_path, plot_dir):
    """
    Improved: label bars by what barcodes mean (mark / experimental condition),
    not the raw barcode strings. Also writes a labeled TSV per figure.

    Inputs:
      subdir_path: e.g. ${outdir}
      plot_dir:    where to write PNGs/TSVs
    """

    import os, re, glob
    import pandas as pd
    import matplotlib.pyplot as plt

    # ----- cohort mapping helpers (match your Split_Reads2.codon) -----

    SUBSAMPLE_BCS = ["GCAT", "CGAT", "GACT", "AGCT", "CAGT", "ACGT", "GCTA", "CGTA", "GTCA", "TGCA", "CTGA", "TCGA"]

    # ISA: 4 groups (3 SBs each)
    ISA_GROUPS = {
        "DMSO": SUBSAMPLE_BCS[0:3],
        "4D_Treatment": SUBSAMPLE_BCS[3:6],
        "4D_Recovery": SUBSAMPLE_BCS[6:9],
        "8D_Recovery": SUBSAMPLE_BCS[9:12],
    }

    # DIV1 (SKD / SKr): 3 groups
    DIV1_GROUPS = {
        "K422_DMSO": SUBSAMPLE_BCS[0:3],
        "K422_96H_Treatment": SUBSAMPLE_BCS[3:6],
        "K422_48H_Recovery": SUBSAMPLE_BCS[6:9],
    }

    # DIV2 (TKD / TKr): 3 groups
    DIV2_GROUPS = {
        "K422_96H_Recovery": SUBSAMPLE_BCS[0:3],
        "K422_8D_Recovery": SUBSAMPLE_BCS[3:6],
        "K422_96H_Treatment_Frozen": SUBSAMPLE_BCS[6:9],
    }

    # Modality mapping
    MO_MAP = {
        "AGGCTATA": "H3K27me3",
        "GCCTCTAT": "H3K27ac",
    }

    def detect_cohort_groups(sample_name: str):
        """Pick the right SB→group mapping dict based on sample prefix in the stats filename."""
        if sample_name.startswith(("Sc_ID", "Sc_Ir")):
            return ISA_GROUPS
        if sample_name.startswith(("Sc_SKD", "Sc_SKr")):
            return DIV1_GROUPS
        if sample_name.startswith(("Sc_TKD", "Sc_TKr")):
            return DIV2_GROUPS
        # Fallback: no relabeling
        return {}

    def label_sb(df_counts: pd.DataFrame, sample_name: str) -> pd.DataFrame:
        """
        Map SB barcodes to group labels (DMSO/… or K422_*).
        Handles raw 4bp, and optionally 'A'/'T'+4bp prefixed strings.
        Returns a two-column DataFrame: ['Label','Value']
        """
        groups = detect_cohort_groups(sample_name)
        if not groups:
            # unknown cohort: keep as-is
            return df_counts.rename(columns={"String": "Label", "Value": "Value"})

        # Build reverse lookup from barcode -> label (include optional A/T prefixed variants)
        bc_to_label = {}
        for label, bcs in groups.items():
            for bc in bcs:
                bc_to_label[bc] = label
                bc_to_label["A" + bc] = label
                bc_to_label["T" + bc] = label

        # Split out NoMatch and map the others
        nomatch = df_counts[df_counts["String"] == "NoMatch"]
        other = df_counts[df_counts["String"] != "NoMatch"].copy()
        other["Label"] = other["String"].map(bc_to_label).fillna("Unknown")

        # Aggregate by Label and add NoMatch back
        agg = other.groupby("Label", as_index=False)["Value"].sum()
        if not nomatch.empty:
            agg = pd.concat([agg, pd.DataFrame({"Label": ["NoMatch"], "Value": [int(nomatch["Value"].sum())]})],
                            ignore_index=True)
        return agg[["Label", "Value"]]

    def label_mo(df_counts: pd.DataFrame) -> pd.DataFrame:
        """Map MO 8-mers to histone marks; keep NoMatch."""
        nomatch = df_counts[df_counts["String"] == "NoMatch"]
        other = df_counts[df_counts["String"] != "NoMatch"].copy()
        other["Label"] = other["String"].map(MO_MAP).fillna("Unknown")
        agg = other.groupby("Label", as_index=False)["Value"].sum()
        if not nomatch.empty:
            agg = pd.concat([agg, pd.DataFrame({"Label": ["NoMatch"], "Value": [int(nomatch["Value"].sum())]})],
                            ignore_index=True)
        return agg[["Label", "Value"]]

    # ----- 1) Demux histogram (unchanged except default colors) -----
    demux_files = glob.glob(f'{subdir_path}/stats/Barcode_Statistics*.tsv')
    if not demux_files:
        raise FileNotFoundError("No TSV files with 'Barcode_Statistics' in their name found.")

    for file_path in demux_files:
        basen = os.path.splitext(os.path.basename(file_path))[0]
        parts = basen.split('_')
        samp = '_'.join(parts[-4:-1])  # e.g., Sc_ID11_S1
        L_num = parts[-1]              # e.g., CB/MO/SB etc
        print(f'{samp}_{L_num}_Demux.png')

        df = pd.read_csv(file_path, sep='\t', header=None, names=['Category', 'Count', 'Percentage'])
        df['Count'] = df['Count'].astype(str).str.replace(',', '', regex=False).astype(int)

        possible_categories = {
            'reads': 'Total Reads',
            'bc_reads_with_0_mismatches': '0',
            'bc_reads_with_1_mismatches': '1',
            'bc_reads_with_2_mismatches': '2',
            'bc_reads_with_3_mismatches': '3',
            'reads_without_bc': 'No bc'
        }

        labels, values = [], []
        for category, label in possible_categories.items():
            match = df.loc[df['Category'] == category, 'Count']
            if not match.empty:
                labels.append(label)
                values.append(int(match.iloc[0]))

        plt.figure(figsize=(10, 6))
        plt.bar(labels, values)
        plt.xlabel('Number of Mismatches in Barcode', fontsize=14)
        plt.ylabel('Number of Reads', fontsize=14)
        plt.title(f'{samp} {L_num}', fontsize=18)
        plt.tight_layout()
        plt.savefig(f'{plot_dir}/{samp}_{L_num}_Demux.png')
        plt.close()

    # ----- 2) Reads per barcode (NOW labeled by meaning) -----
    rpb_files = glob.glob(f'{subdir_path}/stats/Reads_Per_Barcode_*.tsv')
    rpb_files = [file for file in rpb_files if re.search(r'_(MO|SB)\.tsv$', file)]
    if not rpb_files:
        raise FileNotFoundError("No TSV files with 'Reads_Per_Barcode_*_(MO|SB).tsv' found.")

    for file_path in rpb_files:
        basen = os.path.splitext(os.path.basename(file_path))[0]
        parts = basen.split('_')
        # expected patterns like: Reads_Per_Barcode_Sc_ID11_S1_MO
        samp = '_'.join(parts[-4:-1])  # Sc_ID11_S1
        L_num = parts[-1]              # MO or SB
        print(f'{samp}_{L_num}_ReadsPerBC_labeled.png')

        df = pd.read_csv(file_path, sep='\t', header=None, names=['Value', 'String'])

        # Map to labels
        if L_num == 'MO':
            labeled = label_mo(df)
            title = 'Modality'
        else:  # 'SB'
            labeled = label_sb(df, samp)
            title = 'Sample (Condition)'

        # Plot
        plt.figure(figsize=(12, 8))
        plt.bar(labeled['Label'], labeled['Value'])
        plt.xlabel('')
        plt.ylabel('Number of Reads')
        plt.title(f'{samp} {title}', fontsize=18)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(f'{plot_dir}/{samp}_{L_num}_ReadsPerBC_labeled.png')
        plt.close()



if __name__ == "__main__":
    out_folder = '/mnt/dataFast/ahrmad/ScH2R_202409/processed'

    DNA = ["H2RD_Human", "H2RD_Mouse"]
    RNA = ["H2RR_Human", "H2RR_Mouse"]

    # --- DNA ---
    for sample in DNA:
        select_modality = "DNA"
        Tag_records_path = tag_records_for_sample(sample, out_folder)
        bam_file_path = bam_path_for_sample(sample, select_modality, out_folder, mark_for_dna="H3K27ac")
        good_bc_path = good_bc_path_for_sample(sample, select_modality, out_folder, mark_for_dna="H3K27ac")
        make_upsetplot(Tag_records_path, bam_file_path, good_bc_path, select_modality, sample, out_folder)

    # --- RNA ---
    for sample in RNA:
        select_modality = "RNA"
        Tag_records_path = tag_records_for_sample(sample, out_folder)
        bam_file_path = bam_path_for_sample(sample, select_modality, out_folder)
        good_bc_path = good_bc_path_for_sample(sample, select_modality, out_folder)
        make_upsetplot(Tag_records_path, bam_file_path, good_bc_path, select_modality, sample, out_folder)

    # Optional: barplots summary (demux + reads_per_barcode)
    make_barplots(out_folder, out_folder)

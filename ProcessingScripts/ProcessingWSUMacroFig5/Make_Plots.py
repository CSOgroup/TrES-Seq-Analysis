import os
import re
import sys
import glob
import pysam
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from upsetplot import UpSet, plot, from_indicators
import polars as pl
pl.Config.set_fmt_str_lengths(100)


def create_alignment_df(bam_file_path, RNAinGenes=False):
    # Set to store unique read names
    unique_query = set()
    primary_reads_data = []

    # Open the BAM file
    with pysam.AlignmentFile(bam_file_path, "rb", threads=12) as bamfile:
        
        for count, read in enumerate(bamfile.fetch()):
            
            # Check if the read is mapped (not unmapped), is paired, and ensure it's not a duplicate
            if (not read.is_unmapped and read.is_paired and read.query_name not in unique_query):
                unique_query.add(read.query_name)  # Add to the set
                
                if not RNAinGenes:
                    # If RNAinGenes is False, just add the read name, mapped, and duplicate status
                    primary_reads_data.append([read.query_name, not read.is_unmapped, not read.is_duplicate])
                
                # If RNAinGenes is True, check for the XS tag
                elif RNAinGenes:
                    xs_tag = read.get_tag("XS")
                    # Output False if the XS tag content is "Unassigned_NoFeatures", otherwise True
                    is_xs_valid = False if xs_tag == "Unassigned_NoFeatures" else True
                    
                    # Add the read name, mapped, duplicate status, and the XS boolean to primary_reads_data
                    primary_reads_data.append([read.query_name, not read.is_unmapped, not read.is_duplicate, is_xs_valid])

            # Update progress every 100000 reads, overwriting the previous progress
            if (count + 1) % 100000 == 0:
                sys.stdout.write(f"\rProcessed {(count + 1) / 1e6:.1f} million queries...")
                sys.stdout.flush()

    # Final progress message once all reads are processed
    sys.stdout.write(f"\rDone processing {(count + 1)} queries.\n")
    sys.stdout.flush()

    sys.stdout.write('Preparing Polars Query DataFrame...\n')
    sys.stdout.flush()

    # Define the schema dynamically based on RNAinGenes
    if RNAinGenes:
        schema = [('Read_Name_BAM', pl.String), 
                  ('Mapped', pl.Boolean), 
                  ('Unique', pl.Boolean),
                  ('inGenes', pl.Boolean)]
    else:
        schema = [('Read_Name_BAM', pl.String), 
                  ('Mapped', pl.Boolean), 
                  ('Unique', pl.Boolean)]
    
    # Create a polars DataFrame from the list
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
    # Load the CSV file into a Polars DataFrame
    df = pl.read_csv(Tag_records_path, separator='\t', has_header=False)

    print(f'Loading Done. Cleaning...')

    # Extract the first row to determine column names
    first_row = list(df.row(0))[1:]
    column_prefixes = [s[:2] for s in list(first_row)]
    df.columns = ["Record_Name"] + column_prefixes
    
    # Extract only the element part from the formatted data for each existing column
    df = df.with_columns((
            pl.col(col).str.slice(5)
            for col in column_prefixes
        ))
    df = df.drop(["L1","L2","L3"]) #remove ligaation columns
    #df = df.drop("UM", strict=False)
    CheckinGenes = True
    # All subsample barcodes
    subsample_bcs = ["GCAT", "CGAT", "GACT", "AGCT", "CAGT", "ACGT", "GCTA", "CGTA", "GTCA", "TGCA", "CTGA", "TCGA"]

    # For WSU samples, add "A" or "T" before all barcodes (both options)
    # For K422 and A20, add only "A" before all barcodes

    if "WSU_Untreated" in bam_file_path:
        bcs = subsample_bcs[0:3]
        bcs_with_prefix = [p + bc for bc in bcs for p in ("A", "T")]
        df = df.filter(pl.col("SB").is_in(set(bcs_with_prefix)))
    elif "WSU_CD47_Inhibitor" in bam_file_path:
        bcs = subsample_bcs[3:6]
        bcs_with_prefix = [p + bc for bc in bcs for p in ("A", "T")]
        df = df.filter(pl.col("SB").is_in(set(bcs_with_prefix)))
    elif "WSU_DMSO" in bam_file_path:
        bcs = subsample_bcs[6:9]
        bcs_with_prefix = [p + bc for bc in bcs for p in ("A", "T")]
        df = df.filter(pl.col("SB").is_in(set(bcs_with_prefix)))
    elif "WSU_PB1" in bam_file_path:
        bcs = subsample_bcs[9:12]
        bcs_with_prefix = [p + bc for bc in bcs for p in ("A", "T")]
        df = df.filter(pl.col("SB").is_in(set(bcs_with_prefix)))
    elif "K422_K9" in bam_file_path:
        bcs = subsample_bcs[0:4]
        bcs_with_prefix = ["A" + bc for bc in bcs]
        df = df.filter(pl.col("SB").is_in(set(bcs_with_prefix)))
    elif "K422_K36" in bam_file_path:
        bcs = subsample_bcs[4:8]
        bcs_with_prefix = ["A" + bc for bc in bcs]
        df = df.filter(pl.col("SB").is_in(set(bcs_with_prefix)))
    elif "A20_K9" in bam_file_path:
        bcs = subsample_bcs[8:12]
        bcs_with_prefix = ["A" + bc for bc in bcs]
        df = df.filter(pl.col("SB").is_in(set(bcs_with_prefix)))

    if 'DNA' in select_modality:
        df = df.rename({"MO":"Modality"})
        filename = f'{sample}_UpsetPlot.pdf'
        CheckinGenes= False

    print(f'Cleaning Done.\n')
    print(df.head())

    cntBC = pl.read_csv(good_bc_path, separator='\t', has_header=False, new_columns=["MappedNoFilt", "bc"])
    cntBC = cntBC.with_columns(pl.lit(True).alias("MappedNoFilt"))
    cntBC = cntBC.unique(subset=["bc"])


    df = df.join(cntBC,
    left_on='CB',
    right_on='bc',
    how='left',
    validate='m:1',
    coalesce=True
    )

    # Display the filtered first column
    print(f"\nCheck Alignment Status {bam_file_path}...")

    # Merge alignment_df with df based on the read names
    df = df.join(create_alignment_df(bam_file_path,RNAinGenes=CheckinGenes),
     left_on='Record_Name',
     right_on='Read_Name_BAM',
     how='left',
     validate='1:1',
     coalesce=True
     )
    
    if select_modality == 'RNA':
        filename = f'{sample}_UpsetPlot.pdf'
        
        print(f"\nCheck UMI Status {bam_file_path.replace('_ALL', '_UMI')}...")
        df = df.join(create_alignment_df(f"{bam_file_path.replace('_ALL', '_UMI')}"),
        left_on='Record_Name',
        right_on='Read_Name_BAM',
        how='left',
        validate='1:1',
        coalesce=True,
        suffix='_UMI'
        )

        df = df.drop(["Unique","Unique_UMI"], strict=False)
        df = df.rename({"Mapped_UMI": "UMI"})
        df = df.drop(["Mapped_UMI","UM"], strict=False)

    if select_modality == 'DNA':
        mark = 'H3K27ac'
        if "_K36" in sample:
            mark = 'H3K36me3'
        elif "_K9" in sample:
            mark = 'H3K9me3'

        print(f'\nCheck Alignment Status {bam_file_path.replace(mark, "H3K27me3")}...')
        df = df.join(create_alignment_df(bam_file_path.replace(mark, "H3K27me3")),
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

        cntBC = pl.read_csv(good_bc_path.replace(mark, "H3K27me3"), separator='\t', has_header=False, new_columns=["MappedNoFilt", "bc"])
        cntBC = cntBC.with_columns(pl.lit(True).alias("MappedNoFilt"))

        df = df.join(cntBC,
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
    # Replace "NoMatch" with False
    df = df.with_columns([
        pl.when(pl.col(col) == "NoMatch")
        .then(False)
        .otherwise(pl.col(col))
        .alias(col) 
        for col in df.columns
    ])
    
    # Fill the reads not found in df with 0s
    if ('DNA' in select_modality):
        df = df.with_columns(
            pl.col(['Mapped','Unique'])
            .fill_null(False)
            )
    else:
        df = df.with_columns(
            pl.col(['Mapped'])
            .fill_null(False)
            )

    # Set non-"NoMatch" elements to True
    df = df.with_columns([
        pl.when(pl.col(col) != False).
        then(True)
        .otherwise(pl.col(col))
        .alias(col)
        for col in df.columns
    ])

    # Cast df as boolean
    df = df.cast(pl.UInt8).cast(pl.Boolean)

    # Rename Columns for Upset Plot rows
    df = df.rename({"Mapped": "CB>50"})
    df = df.rename({"MappedNoFilt": "Mapped"})
    df = df.rename(lambda column_name: column_name + ' +')
    df = df.rename({"Record_Name +": "Reads", "SB +": "Sample +"})

    print('Creating Figure')
    # Create the upset plot
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
    plt.suptitle(sample.replace("_"," "), fontsize=16, fontweight='bold')
    fig.savefig(f'{out_folder}/{filename}')
    plt.close()
    
    return 0


def make_barplots(subdir_path,plot_dir):
    # Get all TSV files with "reads_per_barcode" in their names
    tsv_files = glob.glob(f'{subdir_path}/stats/Barcode_Statistics*.tsv')

    # Check if any files are found
    if not tsv_files:
        raise FileNotFoundError("No TSV files with 'demux_stats' in their name found.")

    # Process each file
    for file_path in tsv_files:
        basen = os.path.splitext(os.path.basename(file_path))[0]
        parts = basen.split('_')
        samp = '_'.join(parts[-4:-1])  # Extract the last four parts and join the first three
        L_num = parts[-1]             # Extract the last part
        print(f'{samp}_{L_num}_Demux.png')
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"The file {file_path} does not exist.")

        # Read the TSV file using pandas
        df = pd.read_csv(file_path, sep='\t', header=None, names=['Category', 'Count', 'Percentage'])

        # Convert 'Count' column to integer, removing commas if present
        df['Count'] = df['Count'].astype(str).str.replace(',', '').astype(int)

        # Define possible categories
        possible_categories = {
            'reads':'Total Reads',
            'bc_reads_with_0_mismatches': '0',
            'bc_reads_with_1_mismatches': '1',
            'bc_reads_with_2_mismatches': '2',
            'bc_reads_with_3_mismatches': '3',
            'reads_without_bc': 'No bc'
        }

        # Extracting values for the histogram
        labels = []
        values = []

        for category, label in possible_categories.items():
            if df['Category'].str.contains(category).any():
                count = df[df['Category'] == category]['Count'].values[0]
                labels.append(label)
                values.append(count)

        # Plotting the histogram
        plt.figure(figsize=(10, 6))
        plt.bar(labels, values, color='skyblue')
        plt.xlabel('Number of Mismatches in Barcode',fontsize=14)
        plt.ylabel('Number of Reads',fontsize=14)
        plt.title(f'{samp} {L_num}',fontsize=18)
        plt.tight_layout()
        plt.savefig(f'{plot_dir}/{samp}_{L_num}_Demux.png')
        plt.close()
    # Get all TSV files with "reads_per_barcode" in their names
    tsv_files = glob.glob(f'{subdir_path}/stats/Reads_Per_Barcode_*.tsv')
    tsv_files = [file for file in tsv_files if re.search(r'_(MO|SB)\.tsv$', file)]
    # Check if any files are found
    if not tsv_files:
        raise FileNotFoundError("No TSV files with 'reads_per_barcode' in their name found.")

    # Process each file
    for file_path in tsv_files:
        # Define barcode mappings
        sb_bcs = ["GCAT", "CGAT", "GACT", "AGCT", "CAGT", "ACGT", "GCTA", "CGTA", "GTCA", "TGCA", "CTGA", "TCGA"]

        # Define experiment mappings based on Split_Reads.codon context
        exp1_map = {
            sb_bcs[0]: 'WSU_Untreated',
            sb_bcs[1]: 'WSU_Untreated',
            sb_bcs[2]: 'WSU_Untreated',
            sb_bcs[3]: 'WSU_CD47_Inhibitor',
            sb_bcs[4]: 'WSU_CD47_Inhibitor',
            sb_bcs[5]: 'WSU_CD47_Inhibitor',
            sb_bcs[6]: 'WSU_DMSO',
            sb_bcs[7]: 'WSU_DMSO',
            sb_bcs[8]: 'WSU_DMSO',
            sb_bcs[9]: 'WSU_PB1',
            sb_bcs[10]: 'WSU_PB1',
            sb_bcs[11]: 'WSU_PB1'
        }

        exp2_map = {
            sb_bcs[0]: 'K422_K9',
            sb_bcs[1]: 'K422_K9',
            sb_bcs[2]: 'K422_K9',
            sb_bcs[3]: 'K422_K9',
            sb_bcs[4]: 'K422_K36',
            sb_bcs[5]: 'K422_K36',
            sb_bcs[6]: 'K422_K36',
            sb_bcs[7]: 'K422_K36',
            sb_bcs[8]: 'A20_K9',
            sb_bcs[9]: 'A20_K9',
            sb_bcs[10]: 'A20_K9',
            sb_bcs[11]: 'A20_K9'
        }


        basen = os.path.splitext(os.path.basename(file_path))[0]
        parts = basen.split('_')
        samp = '_'.join(parts[-4:-1])  # Extract the last four parts and join the first three
        L_num = parts[-1]             # Extract the last part
        print(f'{samp}_{L_num}_ReadsPerBC.png')
        df = pd.read_csv(file_path, sep='\t', header=None, names=['Value', 'String'])

        if L_num == 'SB':
            tit = 'Sample Barcode'
            # Map barcodes to sample names based on experiment
            df_nomatch = df[df['String'] == 'NoMatch']
            df_other = df[df['String'] != 'NoMatch']
            
            # Determine which experiment mapping to use based on sample name
            if 'MWD' in samp or 'MWr' in samp:
                mo_bcs = {'AGGCTATA': 'H3K27me3', 'GCCTCTAT': 'H3K27ac'}

                bc_map = exp1_map
                #df_other.loc[:, 'String'] = df_other['String'].str[1:]
            elif 'K9' in samp:
                mo_bcs = {'AGGCTATA': 'H3K27me3', 'GCCTCTAT': 'H3K9me3'}
                bc_map = exp2_map
            if 'K36' in samp:
                mo_bcs = {'AGGCTATA': 'H3K27me3', 'GCCTCTAT': 'H3K36me3'}
                bc_map = exp2_map
            
            # Map barcodes to their corresponding sample names using .loc to avoid SettingWithCopyWarning
            df_other.loc[:, 'String'] = df_other['String'].map(bc_map)
            # Group by mapped names and sum values
            df_other = df_other.groupby('String')['Value'].sum().reset_index()
            df = pd.concat([df_other, df_nomatch])

        elif L_num == 'MO':
            tit = 'Modality Barcode'
            # Map molecular barcodes to their names
            df_nomatch = df[df['String'] == 'NoMatch']
            df_other = df[df['String'] != 'NoMatch']
            df_other.loc[:, 'String'] = df_other['String'].map(mo_bcs)
            df_other = df_other.groupby('String')['Value'].sum().reset_index()
            df = pd.concat([df_other, df_nomatch])

        # Calculate the sum of all values
        total_reads = df['Value'].sum()

        # Create a histogram
        plt.figure(figsize=(12, 8))
        plt.bar(df['String'], df['Value'], label='Values')

        # Add the total reads bar at the end
        #plt.bar('totalreads', total_reads, label='Total Reads')

        # Add labels and title
        plt.xlabel('')
        plt.ylabel('Number of Reads')
        plt.title(f'{samp} {tit}', fontsize=18)
        plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
        #plt.legend()

        # Show the plot
        plt.tight_layout()
        plt.savefig(f'{plot_dir}/{samp}_{L_num}_ReadsPerBC.png')
        plt.close()
    plt.close()

    
    return 0


import polars as pl
import matplotlib.pyplot as plt

def process_sample(df, sample, subsample_bcs):
    # Determine which barcodes to use
    if "K422_K9" in sample:
        bcs = subsample_bcs[0:4]
    elif "K422_K36" in sample:
        bcs = subsample_bcs[4:8]
    elif "A20_K9" in sample:
        bcs = subsample_bcs[8:12]
    else:
        return  # Unrecognized sample

    bcs_with_prefix = ["A" + bc for bc in bcs]
    df1 = df.filter(pl.col("SB").is_in(set(bcs_with_prefix)))

    # Count total and matching reads
    total_reads = df1.height
    columns = ["CB", "L1", "L2", "L3", "MO", "SB"]
    match_counts_dict = {
        col: df1.filter(pl.col(col) != "NoMatch").height
        for col in columns
    }

    # Plot results
    ordered_labels = ["Total Reads", "SB", "MO", "L1", "L2", "L3", "CB"]
    ordered_values = [total_reads] + [match_counts_dict[col] for col in ordered_labels[1:]]

    plt.figure(figsize=(10, 6))
    plt.bar(ordered_labels, ordered_values)
    plt.ylabel("Number of Reads")
    plt.title(f"{sample} - Reads with Proper Barcodes")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{sample}_read_barcodes_summary.pdf")
    plt.close()

def process_all_samples(filepath, sample_list, subsample_bcs):
    df = pl.read_csv(filepath, separator='\t', has_header=False)

    # Extract and clean column names
    first_row = list(df.row(0))[1:]
    column_prefixes = [s[:2] for s in first_row]
    df.columns = ["Record_Name"] + column_prefixes

    # Trim barcode prefixes
    df = df.with_columns((
        pl.col(col).str.slice(5) for col in column_prefixes
    ))

    # Process individual samples
    for sample in sample_list:
        process_sample(df, sample, subsample_bcs)

    # Global summary
    total_reads = df.height
    columns = ["CB", "L1", "L2", "L3", "MO", "SB"]
    match_counts_dict = {
        col: df.filter(pl.col(col) != "NoMatch").height
        for col in columns
    }

    ordered_labels = ["Total Reads", "SB", "MO", "L1", "L2", "L3", "CB"]
    ordered_values = [total_reads] + [match_counts_dict[col] for col in ordered_labels[1:]]

    plt.figure(figsize=(10, 6))
    plt.bar(ordered_labels, ordered_values)
    plt.ylabel("Number of Reads")
    plt.title("Reads with Proper Barcodes - All Samples")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig("ALL_Sc_K9D45_S5_read_barcodes_summary.pdf")
    plt.close()


if __name__ == "__main__":
    
    # DNA samples
    sample_list_DNA1 = ["Sc_K9D45_S5_K422_K9","Sc_K9D45_S5_K422_K36","Sc_K9D45_S5_A20_K9"]

    sample_list_DNA2 = [
        "Sc_MWD_WSU_Untreated",
        "Sc_MWD_WSU_CD47_Inhibitor",
        "Sc_MWD_WSU_DMSO",
        "Sc_MWD_WSU_PB1",

    ]

    sample_list_DNA3 = [
        "Sc_K9D45",
    ]
    sample_list_DNA4 = [
        "Sc_MWD",
    ]
    
    # RNA samples
    sample_list_RNA1 = ["Sc_K9r46_S6_K422_K9","Sc_K9r46_S6_K422_K36","Sc_K9r46_S6_A20_K9"]

    sample_list_RNA2 = [
        "Sc_MWr_WSU_Untreated",
        "Sc_MWr_WSU_CD47_Inhibitor",
        "Sc_MWr_WSU_DMSO",
        "Sc_MWr_WSU_PB1",

    ]

    sample_list_RNA3 = [
        "Sc_K9r46",
    ]
    sample_list_RNA4 = [
        "Sc_MWr",
    ]
    
    out_folder='/mnt/dataFast/ahrmad/triseq_202506/processed'
    
    # Barplots
    make_barplots(out_folder,out_folder)

    for sample in sample_list_DNA2:
        Tag_records_path = f'{out_folder}/Tag_Records_Sc_MWD.tsv.gz'
        bam_file_path = f'{out_folder}/{sample}_H3K27ac_MarkedDup.bam'
        good_bc_path = f'{out_folder}/{sample}_H3K27ac_ProperPairedMapped_reads_per_barcode.tsv'
        select_modality = "DNA"
        make_upsetplot(Tag_records_path,bam_file_path,good_bc_path,select_modality,sample,out_folder)
    
    for sample in sample_list_RNA2:
        Tag_records_path = f'{out_folder}/Tag_Records_Sc_MWr.tsv.gz'
        bam_file_path = f'{out_folder}/{sample}_ALL.bam'
        good_bc_path = f'{out_folder}/{sample}_ProperPairedMapped_reads_per_barcode.tsv'
        select_modality = "RNA"
        make_upsetplot(Tag_records_path,bam_file_path,good_bc_path,select_modality,sample,out_folder)
    

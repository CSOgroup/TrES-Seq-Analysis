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

    if select_modality == 'H3K27me3':
        df = df.filter(pl.col("MO") == "AGGCTATA")
    elif select_modality == 'H3K27ac':
        df = df.filter(pl.col("MO") == "GCCTCTAT") 

    if 'DNA' in select_modality:
        df = df.rename({"MO":"Modality"})
        filename = f'{sample}_UpsetPlot.pdf'
        CheckinGenes= False
    
    df = df.filter(pl.col("SB") != "CGTA")

    if 'H3K' in select_modality:
        df = df.rename({"MO":f"{select_modality}"})
        filename = f'{sample}_{select_modality}_UpsetPlot.pdf'
        CheckinGenes= False

    print(f'Cleaning Done.\n')
    print(df.head())

    cntBC = pl.read_csv(good_bc_path, separator='\t', has_header=False, new_columns=["MappedNoFilt", "bc"])
    cntBC = cntBC.with_columns(pl.lit(True).alias("MappedNoFilt"))

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
        
        print(f"\nCheck UMI Status {bam_file_path.replace('_ALL.bam', '_UMI.bam')}...")
        df = df.join(create_alignment_df(f"{bam_file_path.replace('_ALL.bam', '_UMI.bam')}"),
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

if __name__ == "__main__":
    
    # DNA samples
    sample_list_DNA1 = ["H2RD1_S1"]
    
    # RNA samples
    sample_list_RNA1 = ["Sc_K9r"]
    
    mod_list = ["H3K27ac","H3K27me3"]
    out_folder='/mnt/dataFast/ahrmad/ScH2R_202409/processed'
    
    # Barplots
    #make_barplots(out_folder,out_folder)
    for sample in sample_list_DNA1:
        Tag_records_path = f'{out_folder}/Tag_Records_H2RD1_S1.tsv'
        bam_file_path = f'{out_folder}/{sample}_H3K27ac_MarkedDup.bam'
        good_bc_path = f'{out_folder}/{sample}_H3K27ac_ProperPairedMapped_reads_per_barcode.tsv'
        select_modality = "H3K27ac"
        make_upsetplot(Tag_records_path,bam_file_path,good_bc_path,select_modality,sample,out_folder)
        

    for sample in sample_list_DNA1:
        Tag_records_path = f'{out_folder}/Tag_Records_H2RD1_S1.tsv'
        bam_file_path = f'{out_folder}/{sample}_H3K27me3_MarkedDup.bam'
        good_bc_path = f'{out_folder}/{sample}_H3K27me3_ProperPairedMapped_reads_per_barcode.tsv'
        select_modality = "H3K27me3"
        make_upsetplot(Tag_records_path,bam_file_path,good_bc_path,select_modality,sample,out_folder)
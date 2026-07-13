#!/usr/bin/env python

import sys
import argparse
import zipfile
import pandas as pd
import json
import re


def parse_args():
    parser = argparse.ArgumentParser(
        prog='update_assay_table',
        description='Update MethylSeq assay table from ISA.zip with processed data file information.')
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--assay_suffix', required=True, action='store', 
                          help='Specifies the assay type for the dataset (e.g. _GLMethylSeq, _GLRNAMethylSeq).')
    required.add_argument('--runsheet', required=True, 
                          help='Runsheet')
    required.add_argument('--glds_accession', required=True, 
                          help='GLDS accession number (e.g. GLDS-123)')
    isa_file_input = required.add_mutually_exclusive_group(required=True)
    isa_file_input.add_argument('--isa_zip', action='store', default='',
                                help='Appropriate ISA file for the dataset (a zip archive, providing this instead of '
                                     'an assay table directly will attempt to extract the correct assay table given '
                                     'the provided assay and technology types.)')
    isa_file_input.add_argument('--assay_table', action='store', default='',
                                help='Assay table for the dataset provided directly instead of extracting from an ISA '
                                     'zip file')
    parser.add_argument('--is_deduplicated', action=argparse.BooleanOptionalAction,
                        help='Provide this flag if the data was deduplicated.')
    parser.add_argument('--read_counts_from_multiqc', action='store',
                        help='A raw multiqc data file to provide raw read counts.')
    return parser.parse_args()


tty_colors = {
    'green': '\033[0;32m%s\033[0m',
    'yellow': '\033[0;33m%s\033[0m',
    'red': '\033[0;31m%s\033[0m'
}


def color_text(text, color='green'):
    """
    Colors text for output in terminal

    Args:
        text (str): input text
        color (str): a valid tty color ('red', 'yellow', or 'green')

    Returns:
        str: colored text

    """
    if sys.stdout.isatty():
        return tty_colors[color] % text
    else:
        return text


def report_failure_and_exit(message, color="red"):
    """
    Reports a failure and exits with status '1'.

    Args:
        message (str): Error message to report.
        color (str): Color in which to render the error message, default = 'red'
    """
    print("")
    print(color_text(f"Error: {message}", color))
    print("\nAssay table update failed.\n")

    sys.exit(1)


def report_warning(message, color="yellow"):
    """
    Reports are warning message.

    Args:
        message (str): Error message to report
        color (str): Color in which to render the error message, default = 'yellow'
    """
    print("")
    print(color_text(f"Warning: {message}", color))


def load_runsheet(runsheet_file):
    """
    Load the runsheet as a pandas.DataFrame.

    Args:
        runsheet_file (PathLike[str]): a file containing the assay runsheet used to generate the processed data

    Returns:
        pandas.DataFrame: sample information from the runsheet
    """
    try:
        runsheet_df = pd.read_csv(runsheet_file)
        print(f"Runsheet has {len(runsheet_df)} rows and {len(runsheet_df.columns)} columns")
        return runsheet_df
    except Exception as e:
        report_warning(f"Cannot read runsheet, proceeding without it: {e}")
        return None


def get_runsheet_sample_name_map(runsheet_df, assay_sample_names):
    """
    Generates a mapping of sample names in the assay table to the samplenames in the runsheet

    Args:
        runsheet_df (pandas.DataFrame): runsheet sample information
        assay_sample_names (list): sample names from the assay table

    Returns:
        dict: sample name mapping
    """
    sample_name_map = {}
    if 'Sample Name' in runsheet_df.columns:
        # Check for 'Original Sample Name' column to map between assay table and runsheet
        if 'Original Sample Name' in runsheet_df.columns:
            for _, row in runsheet_df.iterrows():
                orig_name = row['Original Sample Name']
                rs_name = row['Sample Name']
                if orig_name in assay_sample_names:
                    sample_name_map[orig_name] = rs_name
    return sample_name_map


def is_paired_end_data(runsheet_df):
    """
    Determine if this is paired-end data based on information the runsheet.
    
    Args:
        runsheet_df (pandas.DataFrame): runsheet sample information

    Returns:
        bool: True if paired-end, False if single-end
    """
    if runsheet_df is None:
        print("Warning: No runsheet provided, assuming single-end data")
        return False

    # Check if there's a paired_end column
    if 'paired_end' in runsheet_df.columns:
        paired_end_values = runsheet_df['paired_end'].unique()
        if len(paired_end_values) == 1:
            # If all values are the same, use that
            value = paired_end_values[0]
            # Handle different types of values (string or boolean)
            if isinstance(value, bool):
                return value
            elif isinstance(value, str):
                return value.lower() == 'true'
            else:
                # Try to convert to string if it's not a boolean or string
                return str(value).lower() == 'true'

    # Check based on R1/R2 file presence
    if any(col for col in runsheet_df.columns if col.endswith('_R2.fastq.gz')):
        print("Detected paired-end data based on _R2 files in runsheet")
        return True
    else:
        print("Assuming single-end data (no _R2 files in runsheet)")
        return False


def get_assay_table_from_isa(isa_file, assay, technology):
    """
    tries to find an assay table in an ISA zip file that matches the type expected for the provided assay

    Args:
        isa_file (PathLike[str]): path to ISA zip file
        assay (str): valid OSDR methylseq assay type ('MethylSeq' or 'RNAMethylSeq')
        technology (str): valid OSDR methylseq technology type

    Returns:
        pandas.DataFrame: assay table from extracted from ISA zip
    """

    if assay == 'MethylSeq':
        valid_measurement = 'DNA methylation profiling'
    else:
        valid_measurement = 'RNA methylation profiling'

    zip_file = zipfile.ZipFile(isa_file)
    isa_files = zip_file.namelist()

    # Parse investigation file to build STUDY ASSAYS table
    study_assays_table = {}
    study_assays_section = False

    with zip_file.open('i_Investigation.txt', 'r') as f:
        for line in f:
            line = line.decode('utf-8').strip()

            # Track STUDY ASSAYS section
            if line == 'STUDY ASSAYS':
                study_assays_section = True
                continue
            elif study_assays_section and not line:
                study_assays_section = False
                continue

            # Extract data from section
            if study_assays_section and line:
                parts = line.split('\t')
                if parts and parts[0]:
                    key = parts[0]
                    values = [v.strip() for v in parts[1:] if v.strip()]
                    study_assays_table[key] = values
    # Check if we have all required keys
    required_keys = ['Study Assay Measurement Type', 'Study Assay Technology Type', 'Study Assay File Name']
    if not all(key in study_assays_table for key in required_keys):
        report_failure_and_exit("Missing required keys in STUDY ASSAYS section")

    # Get the values from the table
    measurement_types = study_assays_table['Study Assay Measurement Type']
    technology_types = study_assays_table['Study Assay Technology Type']
    file_names = study_assays_table['Study Assay File Name']

    # Ensure all lists have equal length
    if not (len(measurement_types) == len(technology_types) == len(file_names)):
        report_failure_and_exit("Measurement types, technology types, and file names have different lengths")

    # Find matching assay file
    matched_file = ""
    for i in range(len(measurement_types)):
        if (measurement_types[i].lower() == valid_measurement.lower() and
                technology_types[i].lower() == technology.lower()):
            matched_file = file_names[i]
            break

    if not matched_file:
        report_failure_and_exit(f"No assay file matched for {assay}. "
                                f"Measurement types: {measurement_types}, Technology types: {technology_types}")
    elif matched_file not in isa_files:
        # Load the matched assay file
        report_failure_and_exit(f"Matched assay file doesn't exist in ISA zip: {matched_file}")
    else:
        return pd.read_csv(zip_file.open(matched_file), sep='\t'), matched_file

    return pd.DataFrame(), matched_file


def fastqc_info_from_multiqc_json(multiqc_zip):
    zip_file = zipfile.ZipFile(multiqc_zip)

    # Find the multiqc_data.json file in the archive
    json_filename = None
    for name in zip_file.namelist():
        if name.endswith('multiqc_data.json'):
            json_filename = name
            break
    
    if json_filename is None:
        raise FileNotFoundError("multiqc_data.json not found in the archive")
    
    json_data = json.load(zip_file.open(json_filename))

    # figure out which index in the json report tables is associated with FastQC
    # assume the first index is FastQC if nothing is named FastQC
    fastqc_index = 0
    for i, key in enumerate(json_data['report_data_sources'].keys()):
        if key == 'FastQC':
            fastqc_index = i

    return pd.DataFrame.from_dict(json_data['report_general_stats_data'][fastqc_index], orient='index')


def get_read_count_from_df(sample_name: str, assay_suffix: str, fastqc_info_df: pd.DataFrame, is_paired_end: bool = True):
    """
    Retrieves read count from a pandas.DataFrame of the raw multiqc results

    Args:
        sample_name (str): name of the sample to get read count for
        assay_suffix (str): assay suffix to be added to sample names when looking up in multiqc report (e.g. "_GLMethylSeq")
        fastqc_info_df (pandas.DataFrame): FastQC general stats data from the multiqc.json file
        is_paired_end (bool): specifies if the data is paired-end (determines file naming)

    Returns:
        number: read count from fastqc info

    """
    full_sample_name = str(sample_name) + assay_suffix
    
    if is_paired_end:
        # Try both naming conventions
        key = full_sample_name + ' Read 1'
        if key not in fastqc_info_df.index:
            key = full_sample_name + '_R1'
        return round(fastqc_info_df.at[key, 'total_sequences'])
    else:
        return round(fastqc_info_df.at[full_sample_name, 'total_sequences'])
    

def add_read_counts(df, multiqc_zip, assay_suffix, sample_name_map, assay_sample_names, is_paired_end=True):
    """
    Add the read counts columns to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        multiqc_zip (PathLike[str]): multiqc data.zip file from raw fastq multiqc processing step
        assay_suffix (str): assay suffix to be added to sample names when looking up in multiqc report (e.g. "_GLMethylSeq")
        sample_name_map (dict): Dictionary mapping assay table sample names to runsheet sample names
        assay_sample_names (list): List of sample names in the assay table
        is_paired_end (bool): specifies if the data is paired-end

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Read Depth]"

    fastqc_info_df = fastqc_info_from_multiqc_json(multiqc_zip)

    values = []
    for assay_sample in assay_sample_names:
        # Use mapped name if available, otherwise use the assay table name
        sample = sample_name_map.get(assay_sample, assay_sample)
        values.append(get_read_count_from_df(sample, assay_suffix, fastqc_info_df, is_paired_end))

    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
    else:
        print(f"Updating column: {column_name}")
    df[column_name] = values

    return df


def add_parameter_column(df, column_name, value, glds_prefix=None):
    """
    Add a parameter column to the dataframe if it doesn't exist already. Use the same value for all rows in the table.
    
    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each value (only for values that are filenames)
        column_name (str): parameter column name to add (e.g., "Parameter Value[Entry]")
        value (str): value to set for all rows in the table

    Returns:
        pandas.DataFrame: update assay table
    """
    # Apply prefix to value if provided
    if glds_prefix and isinstance(value, str):
        # Check if value already has the prefix
        if not value.startswith(glds_prefix):
            prefixed_value = f"{glds_prefix}{value}"
        else:
            prefixed_value = value
    else:
        prefixed_value = value

    if column_name not in df.columns:
        print(f"Adding new column: {column_name}")
        df[column_name] = prefixed_value
    else:
        print(f"Column {column_name} already exists, updating values")
        df[column_name] = prefixed_value

    return df


def add_fastq_data_column(df, glds_prefix, assay_suffix, sample_name_map, assay_sample_names, fastq_type='trimmed',
                          is_paired_end=True):
    """
    Add the Trimmed Sequence Data column and the optional Merged Sequence Data column to the dataframe.
    
    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename
        assay_sample_names (list): sample names found in assay table
        sample_name_map (dict): maps from sample names in assay table to sample names used in files (from runsheet)
        is_paired_end (bool): specifies if the data is paired-end (determines file naming)
        fastq_type (str): 'trimmed' or 'raw'

    Returns:
        pandas.DataFrame: updated assay table
    """
    if fastq_type == 'trimmed':
        column_name = "Parameter Value[Trimmed Sequence Data]"
    elif fastq_type == 'raw':
        column_name = "Parameter Value[Merged Sequence Data]"
    else:
        print(f"Unrecognized fastq type: '{fastq_type}'. Must be one of ['raw', 'trimmed']")
        sys.exit(1)

    # Generate file paths using the appropriate sample names
    values = []
    for assay_sample in assay_sample_names:
        # Use mapped name if available, otherwise use the assay table name
        sample = sample_name_map.get(assay_sample, assay_sample)

        if is_paired_end:
            # For paired-end data, create entries with both R1 and R2 files, comma-separated without spaces
            values.append(
                f"{glds_prefix}{sample}{assay_suffix}_R1_{fastq_type}.fastq.gz,{glds_prefix}{sample}{assay_suffix}_R2_{fastq_type}.fastq.gz")
        else:
            # For single-end data
            values.append(f"{glds_prefix}{sample}{assay_suffix}_{fastq_type}.fastq.gz")

    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values

    return df


def add_trimming_reports_column(df, glds_prefix, assay_suffix, assay_sample_names, sample_name_map, is_paired_end=True):
    """
    Add the Trimming Reports column to the dataframe.
    
    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename
        assay_sample_names (list): sample names found in assay table
        sample_name_map (dict): maps from sample names in assay table to sample names used in files (from runsheet)
        is_paired_end (bool): specifies if the data is paired-end (determines file naming)

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Trimmed Sequence Data/Trimming Reports]"

    # Generate file paths using the appropriate sample names
    values = []
    for assay_sample in assay_sample_names:
        # Use mapped name if available, otherwise use the assay table name
        sample = sample_name_map.get(assay_sample, assay_sample)

        if is_paired_end:
            # For paired-end data, create report entries for both R1 and R2 files
            values.append(f"{glds_prefix}{sample}_R1{assay_suffix}_trimming_report.txt,"
                          f"{glds_prefix}{sample}_R2{assay_suffix}_trimming_report.txt")
        else:
            # For single-end data
            values.append(f"{glds_prefix}{sample}{assay_suffix}_trimming_report.txt")
    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values

    return df


def add_aligned_sequence_data_column(df, glds_prefix, assay_suffix, assay_sample_names, sample_name_map,
                                     is_deduplicated=False, is_rna=False):
    """
    Add the aligned sequence data column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename
        assay_sample_names (list): sample names found in assay table
        sample_name_map (dict): maps from sample names in assay table to sample names used in files (from runsheet)
        is_deduplicated (bool): specifies if the data was deduplicated (determines file naming)
        is_rna (bool): specifies if RNA methylseq data (determines file naming)

    """
    # Title Case column name

    column_name = "Parameter Value[Aligned Sequence Data]"

    if is_rna:
        bam_suffix = '_hisat2'
    else:
        bam_suffix = '_bt2'

    # Generate file paths using the appropriate sample names
    values = []
    for assay_sample in assay_sample_names:
        # Use mapped name if available, otherwise use the assay table name
        sample = sample_name_map.get(assay_sample, assay_sample)

        # Add BAM and BAI files
        if is_deduplicated:
            # include deduplicated BAM/BAI files
            values.append(f"{glds_prefix}{sample}{assay_suffix}_bismark{bam_suffix}_sorted.bam,"
                          f"{glds_prefix}{sample}{assay_suffix}_bismark{bam_suffix}_sorted.bam.bai,"
                          f"{glds_prefix}{sample}{assay_suffix}_bismark{bam_suffix}.deduplicated_sorted.bam,"
                          f"{glds_prefix}{sample}{assay_suffix}_bismark{bam_suffix}.deduplicated_sorted.bam.bai")
        else:
            values.append(f"{glds_prefix}{sample}{assay_suffix}_bismark{bam_suffix}_sorted.bam,"
                          f"{glds_prefix}{sample}{assay_suffix}_bismark{bam_suffix}_sorted.bam.bai")

    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values

    return df


def add_alignment_reports_column(df, glds_prefix, assay_suffix, assay_sample_names, sample_name_map,
                                 is_rna=False, is_paired_end=True):
    """
    Add the Alignment Reports column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename
        assay_sample_names (list): sample names found in assay table
        sample_name_map (dict): maps from sample names in assay table to sample names used in files (from runsheet)
        is_rna (bool): specifies if RNA methylseq data (determines file naming)
        is_paired_end (bool): specifies if the data is paired-end (determines file naming)

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Aligned Sequence Data/Alignment Reports]"

    if is_rna:
        bam_suffix = '_hisat2'
    else:
        bam_suffix = '_bt2'

    if is_paired_end:
        report_suffix = bam_suffix + '_PE'
    else:
        report_suffix = bam_suffix + '_SE'

    # Generate file paths using the appropriate sample names
    values = []
    for assay_sample in assay_sample_names:
        # Use mapped name if available, otherwise use the assay table name
        sample = sample_name_map.get(assay_sample, assay_sample)

        # Add the qualimap
        values.append(f"{glds_prefix}{sample}{assay_suffix}_bismark{bam_suffix}_qualimap.zip,"
                      f"{glds_prefix}{sample}{assay_suffix}_bismark{bam_suffix}.nucleotide_stats.txt,"
                      f"{glds_prefix}{sample}{assay_suffix}_bismark{report_suffix}_report.txt")

    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values

    return df


def add_deduplication_report_column(df, glds_prefix, assay_suffix, assay_sample_names=None, sample_name_map=None,
                                    is_rna=False, is_paired_end=True):
    """
    Add the Deduplication Reports column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename
        assay_sample_names (list): sample names found in assay table
        sample_name_map (dict): maps from sample names in assay table to sample names used in files (from runsheet)
        is_rna (bool): specifies if RNA methylseq data (determines file naming)
        is_paired_end (bool): specifies if the data is paired-end (determines file naming)

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Aligned Sequence Data/Deduplication Reports]"

    if is_rna:
        bam_suffix = '_hisat2'
    else:
        bam_suffix = '_bt2'

    if is_paired_end:
        bam_suffix += '_pe'

    if assay_sample_names is None:
        print("Warning: Could not find Sample Name column in assay table")
        # If no sample column, just use placeholder values
        values = [f"{glds_prefix}sample{i + 1}{assay_suffix}_bismark{bam_suffix}.deduplication_report.txt" for i in range(len(df))]
    else:
        # Generate file paths using the appropriate sample names
        values = []
        for assay_sample in assay_sample_names:
            # Use mapped name if available, otherwise use the assay table name
            sample = sample_name_map.get(assay_sample, assay_sample)

            # Create just the log file
            values.append(f"{glds_prefix}{sample}{assay_suffix}_bismark{bam_suffix}.deduplication_report.txt")

    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values

    return df


def add_methylation_call_data_column(df, glds_prefix, assay_suffix, assay_sample_names, sample_name_map, is_paired_end=True,
                                     is_deduplicated=False, is_rna=False):
    """
    Add the Methylation Call Data column to the dataframe.
    
    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename
        assay_sample_names (list): sample names found in assay table
        sample_name_map (dict): maps from sample names in assay table to sample names used in files (from runsheet)
        is_paired_end (bool): specifies if the data is paired-end (determines file naming)
        is_deduplicated (bool): specifies if the data was deduplicated (determines file naming)
        is_rna (bool): specifies if RNA methylseq data (determines file naming)

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Methylation Call Data]"

    if is_rna:
        file_suffix = '_hisat2'
    else:
        file_suffix = '_bt2'

    if is_paired_end:
        file_suffix += '_pe'

    if is_deduplicated:
        file_suffix += '.deduplicated'

    # Generate file paths using the appropriate sample names
    values = []
    for assay_sample in assay_sample_names:
        # Use mapped name if available, otherwise use the assay table name
        sample = sample_name_map.get(assay_sample, assay_sample)
        values.append(f"{glds_prefix}CHG_context_{sample}{assay_suffix}_bismark{file_suffix}.txt.gz,"
                      f"{glds_prefix}CHH_context_{sample}{assay_suffix}_bismark{file_suffix}.txt.gz,"
                      f"{glds_prefix}CpG_context_{sample}{assay_suffix}_bismark{file_suffix}.txt.gz,"
                      f"{glds_prefix}{sample}{assay_suffix}_bismark{file_suffix}.bedGraph.gz,"
                      f"{glds_prefix}{sample}{assay_suffix}_bismark{file_suffix}.bismark.cov.gz,"
                      f"{glds_prefix}{sample}{assay_suffix}_bismark{file_suffix}.CpG_report.txt.gz")

    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values

    return df


def add_methylation_call_reports_column(df, glds_prefix, assay_suffix, assay_sample_names, sample_name_map, is_paired_end=True,
                                        is_deduplicated=False, is_rna=False):
    """
    Add the Methylation Call Reports column to the dataframe.
    
    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename
        assay_sample_names (list): sample names found in assay table
        sample_name_map (dict): maps from sample names in assay table to sample names used in files (from runsheet)
        is_paired_end (bool): specifies if the data is paired-end (determines file naming)
        is_deduplicated (bool): specifies if the data was deduplicated (determines file naming)
        is_rna (bool): specifies if RNA methylseq data (determines file naming)

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Methylation Call Data/Methylation Call Reports]"

    if is_rna:
        report_suffix = '_hisat2'
    else:
        report_suffix = '_bt2'

    if is_paired_end:
        report_html_suffix = report_suffix + '_PE'
        report_suffix += '_pe'
    else:
        report_html_suffix = report_suffix + '_SE'

    if is_deduplicated:
        report_suffix += '.deduplicated'

    # Generate file paths using the appropriate sample names
    values = []
    for assay_sample in assay_sample_names:
        # Use mapped name if available, otherwise use the assay table name
        sample = sample_name_map.get(assay_sample, assay_sample)
        values.append(f"{glds_prefix}{sample}{assay_suffix}_bismark{report_suffix}.cytosine_context_summary.txt,"
                      f"{glds_prefix}{sample}{assay_suffix}_bismark{report_suffix}.M-bias.txt,"
                      f"{glds_prefix}{sample}{assay_suffix}_bismark{report_suffix}_splitting_report.txt,"
                      f"{glds_prefix}{sample}{assay_suffix}_bismark{report_html_suffix}_report.html")

    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values

    return df


def add_methylation_call_summary_reports_column(df, glds_prefix, assay_suffix):
    """
    Add the methylation call summary reports data column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Methylation Call Data/Methylation Call Reports/Methylation Call Summary Reports]"

    # Create the alignment multiqc report filename - same for all samples
    summary_reports = (f"{glds_prefix}bismark_summary_report{assay_suffix}.html,"
                       f"{glds_prefix}bismark_summary_report{assay_suffix}.txt")

    # Add the column to the dataframe with the same value for all rows
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = summary_reports
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = summary_reports

    return df


def add_multiqc_reports_column(df, glds_prefix, assay_suffix, multiqc_type):
    """
    Add the raw/trimmed/aligned sequence MultiQC reports data column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename
        multiqc_type (str): specifies if the processing step for which these multiqc reports were generated

    Returns:
        pandas.DataFrame: updated assay table
    """
    if multiqc_type == 'trimmed':
        column_name = "Parameter Value[Trimmed Sequence Data/MultiQC Reports]"
    elif multiqc_type == 'raw':
        column_name = "Parameter Value[Merged Sequence Data/MultiQC Reports]"
    elif multiqc_type == 'align_and_bismark':
        column_name = "Parameter Value[Methylation Call Data/MultiQC Reports]"
    else:
        print(f"Unrecognized fastq type: '{multiqc_type}'. Must be one of ['raw', 'trimmed', 'align_and_bismark']")
        sys.exit(1)

    # Create the alignment multiqc report filename - same for all samples
    multiqc_report = (f"{glds_prefix}{multiqc_type}_multiqc{assay_suffix}.html,"
                      f"{glds_prefix}{multiqc_type}_multiqc{assay_suffix}_data.zip")

    # Add the column to the dataframe with the same value for all rows
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = multiqc_report
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = multiqc_report

    return df


def add_differential_methylation_column(df, glds_prefix, assay_suffix):
    """
    Add the Differential Methylation Analysis Data column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename
        assay_suffix (str): a suffix to add to the end of each filename

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Differential Methylation Analysis Data]"

    # Create the differential expression filenames - same for all samples
    de_files = [
        f"{glds_prefix}SampleTable{assay_suffix}.csv",
        f"{glds_prefix}contrasts{assay_suffix}.csv",
        f"{glds_prefix}differential_methylation_bases{assay_suffix}.csv",
        f"{glds_prefix}differential_methylation_tiles{assay_suffix}.csv"
    ]

    # Join the files with commas
    combined_files = ','.join(de_files)

    # Add the column to the dataframe with the same value for all rows
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = combined_files
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = combined_files

    return df


def clean_comma_space(df):
    """
    Remove spaces after commas in all string columns of the dataframe.

    Args:
        df (pandas.DataFrame): current assay table

    Returns:
        pandas.DataFrame: updated assay table
    """
    # Loop through all columns in the dataframe
    for col in df.columns:
        # Only process string (object) columns
        if df[col].dtype == 'object':
            # Replace comma-space with just comma
            df[col] = df[col].str.replace(', ', ',', regex=False)

    print("Removed spaces after commas in all string columns")
    return df


def clean_column_names(df):
    """Clean column names by removing any .# suffixes pandas adds to duplicates.
    
    Args:
        df: The DataFrame to clean column names
        
    Returns:
        The DataFrame with cleaned column names
    """
    # Create a mapping of old_name -> new_name (without .# suffix)
    name_mapping = {}
    for col in df.columns:
        # Use regex to match column names with .digits suffix
        if re.search(r'\.\d+$', col):
            # Remove the .# suffix
            base_name = re.sub(r'\.\d+$', '', col)
            name_mapping[col] = base_name
    
    # Rename columns using the mapping if any found
    if name_mapping:
        print(f"Cleaning {len(name_mapping)} column names by removing .# suffixes:")
        for old_name, new_name in name_mapping.items():
            print(f"  - {old_name} -> {new_name}")
        df = df.rename(columns=name_mapping)
    
    return df


def main():
    args = parse_args()

    # Find and parse the runsheet
    runsheet_df = load_runsheet(args.runsheet)

    # Extract assay type from assay_suffix and technology from runsheet
    assay = args.assay_suffix.replace("_GL", "", 1)

    unique_tech = runsheet_df['Study Assay Technology Type'].unique()
    if len(unique_tech) != 1:
        raise ValueError(f"Runsheet contains multiple technologies: {unique_tech}")
    technology = unique_tech[0]

    is_rna = False
    if technology == "Whole Genome Bisulfite Sequencing":
        resource_category = 'wgbs'
    elif technology == "Whole Transcriptome Bisulfite Sequencing":
        resource_category = 'wtbs'
        is_rna = True
    else:
        if assay == 'MethylSeq':
            resource_category = 'RRBS'
        else:
            resource_category = 'tRRBS'
            is_rna = True

    # Create GLDS prefix for all filenames
    glds_id = args.glds_accession.upper()
    glds_prefix = f"{glds_id}_G{resource_category}_"

    # Find Methyl-Seq assay file and get its contents
    if args.isa_zip != '':
        print(f"Extracting assay table from {args.isa_zip} based on assay type: {assay} "
              f"and technology type: {technology}")
        assay_df, assay_filename = get_assay_table_from_isa(args.isa_zip, assay, technology)
        print(f"Original assay table has {len(assay_df)} rows and {len(assay_df.columns)} columns")
    else:
        print(f"Reading user-provided assay table from file: {args.assay_table}")
        assay_filename = args.assay_table
        assay_df = pd.read_csv(open(assay_filename), sep='\t')
        print(f"Original assay table has {len(assay_df)} rows and {len(assay_df.columns)} columns")

    # Determine if paired-end from runsheet
    is_paired_end = is_paired_end_data(runsheet_df)
    print(f"Data is {'paired-end' if is_paired_end else 'single-end'} based on runsheet")

    # Create a mapping from assay table sample names to runsheet sample names, exit if no sample column found
    sample_col = next((col for col in assay_df.columns if 'Sample Name' in col), None)
    if sample_col is None:
        report_failure_and_exit(f"Could not find 'Sample Name' column in assay table '{assay_filename}'")

    assay_sample_names = assay_df[sample_col].tolist()
    if runsheet_df is not None and 'Sample Name' in runsheet_df.columns:
        sample_name_map = get_runsheet_sample_name_map(runsheet_df, assay_sample_names)
    else:
        sample_name_map = {}

    is_deduplicated = args.is_deduplicated
    if 'Reduced-Representation' in resource_category and is_deduplicated:
        report_failure_and_exit(f"Data was reported to be deduplicated, but the technology type is '{technology}'."
                                "Reduced-representation data should not be deduplicated.")

    # Process and save assay file
    try:

        # Following the original workflow order:

        # 0. Merged Sequence Data column
        assay_df = add_fastq_data_column(assay_df, glds_prefix, args.assay_suffix, fastq_type='raw', is_paired_end=is_paired_end, 
                                         assay_sample_names=assay_sample_names, sample_name_map=sample_name_map)

        # Merged Sequence Data/MultiQC Reports column
        assay_df = add_multiqc_reports_column(assay_df, glds_prefix, args.assay_suffix, multiqc_type='raw')

        # Read counts
        assay_df = add_read_counts(assay_df, args.read_counts_from_multiqc, args.assay_suffix, sample_name_map=sample_name_map, 
                                   assay_sample_names=assay_sample_names, is_paired_end=is_paired_end)
        assay_df = add_parameter_column(assay_df, column_name='Unit', value='reads')

        # 1. Trimmed Sequence Data column
        assay_df = add_fastq_data_column(assay_df, glds_prefix, args.assay_suffix, fastq_type='trimmed', is_paired_end=is_paired_end,
                                         sample_name_map=sample_name_map, assay_sample_names=assay_sample_names)

        # 2. Trimmed Sequence Data/MultiQC Reports column
        assay_df = add_multiqc_reports_column(assay_df, glds_prefix, args.assay_suffix, multiqc_type='trimmed')

        # 3. Trimmed Sequence Data/Trimming Reports column
        assay_df = add_trimming_reports_column(assay_df, glds_prefix, args.assay_suffix, assay_sample_names=assay_sample_names,
                                               sample_name_map=sample_name_map, is_paired_end=is_paired_end)

        # 4. Aligned Sequence Data column
        assay_df = add_aligned_sequence_data_column(assay_df, glds_prefix, args.assay_suffix, assay_sample_names=assay_sample_names,
                                                    sample_name_map=sample_name_map, is_deduplicated=is_deduplicated,
                                                    is_rna=is_rna)

        # 5. Aligned Sequence Data/Alignment Reports
        assay_df = add_alignment_reports_column(assay_df, glds_prefix, args.assay_suffix, is_rna=is_rna, is_paired_end=is_paired_end,
                                                assay_sample_names=assay_sample_names, sample_name_map=sample_name_map)

        # 6. Aligned Sequence Data/Deduplication Reports
        if is_deduplicated:
            assay_df = add_deduplication_report_column(assay_df, glds_prefix, args.assay_suffix, assay_sample_names=assay_sample_names,
                                                       sample_name_map=sample_name_map, is_rna=is_rna,
                                                       is_paired_end=is_paired_end)

        # 7. Methylation Call Data column
        assay_df = add_methylation_call_data_column(assay_df, glds_prefix, args.assay_suffix, assay_sample_names=assay_sample_names,
                                                    sample_name_map=sample_name_map, is_paired_end=is_paired_end,
                                                    is_deduplicated=is_deduplicated, is_rna=is_rna)

        # 8. Methylation Call Data/Methylation Call Reports
        assay_df = add_methylation_call_reports_column(assay_df, glds_prefix, args.assay_suffix, assay_sample_names=assay_sample_names,
                                                       sample_name_map=sample_name_map, is_paired_end=is_paired_end,
                                                       is_deduplicated=is_deduplicated, is_rna=is_rna)

        # 9. Methylation Call Data/Methylation Call Reports/Methylation Call Summary Reports
        assay_df = add_methylation_call_summary_reports_column(assay_df, glds_prefix, args.assay_suffix)

        # 10. Methylation Call Data/MultiQC Reports
        assay_df = add_multiqc_reports_column(assay_df, glds_prefix, args.assay_suffix, multiqc_type='align_and_bismark')

        # 11. Differential Methylation Analysis Data column
        assay_df = add_differential_methylation_column(assay_df, glds_prefix, args.assay_suffix)

        # Clean comma-space in all string columns
        assay_df = clean_comma_space(assay_df)

        # Clean column names by removing any .# suffixes pandas adds
        assay_df = clean_column_names(assay_df)

        # Use the filename we found in extract_and_find_assay
        orig_filename = assay_filename

        # Create both original and modified output files
        # Original file (preserving the original name)
        assay_df.to_csv(orig_filename, sep='\t', index=False)
        print(f"Original assay table saved as: {orig_filename}")

        # Modified file with GLDS prefix
        if not orig_filename.startswith(glds_prefix):
            mod_filename = f"{glds_prefix}{orig_filename}"
        else:
            mod_filename = orig_filename

        assay_df.to_csv(mod_filename, sep='\t', index=False)
        print(f"Modified assay table saved as: {mod_filename}")

    except Exception as e:
        report_failure_and_exit(str(e))


if __name__ == "__main__":
    main()
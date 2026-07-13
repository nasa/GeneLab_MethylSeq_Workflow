#!/usr/bin/env python

"""
Validate GeneLab MethylSeq processed datasets.
"""

import os
import sys
import argparse
import textwrap
import zipfile
import csv
import json
from statistics import mean, median, stdev

parser = argparse.ArgumentParser(description="Validate GeneLab MethylSeq processed datasets.")

required = parser.add_argument_group('required arguments')

required.add_argument("-a", "--assay_suffix", choices=['_GLMethylSeq', '_GLRNAMethylSeq'],
                      help="Specifies which datatype (assay) is to be validated. Parsed from assay_suffix in ch_meta", action="store", required=True)
required.add_argument("-g", "--GLDS_ID", 
                      help='GLDS ID (e.g. "GLDS-47")', action="store", required=True)
required.add_argument("-s", "--runsheet", 
                      help="Runsheet with samplenames in the first column", action="store", required=True)
required.add_argument("--outdir", 
                      help="Directory of generated output files to be validated.", action="store", required=True)
parser.add_argument("--output", "--output_name", default="",
                    help='Name of output file ex. "GLDS-47-methylseq-validation.log"')
parser.add_argument("-p", "--output_prefix", 
                    help="Output additional file prefix if there is one", action="store", default="")
parser.add_argument("--single_ended", 
                    help="Add this flag if data are single-end sequencing. Parsed from paired_end in ch_meta", action="store_true")
parser.add_argument("--deduped", 
                    help="Add this flag if the data have been deduplicated. Parsed from deduplication in ch_meta", action="store_true")
parser.add_argument("--include_raw_fastq", 
                    help="Add this flag if we should check existence of the raw fastq files.", action="store_true")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_args()


######################### Aesthetic functions #################################

# setting some colors
tty_colors = {
    'green': '\033[0;32m%s\033[0m',
    'yellow': '\033[0;33m%s\033[0m',
    'red': '\033[0;31m%s\033[0m'
}


def color_text(text, color='green'):
    if sys.stdout.isatty():
        return tty_colors[color] % text
    else:
        return text


def wprint(text):
    """ print wrapper """

    print(textwrap.fill(text, width=120, initial_indent="  ",
                        subsequent_indent="  ", break_on_hyphens=False))


def report_failure(message, color="yellow", log=True):
    print("")
    wprint(color_text(message, color))
    print("\nValidation failed.\n")

    if log:
        with open(validation_log, "a") as log:
            log.write(f"{message}\nValidation failed.\n\n")

    sys.exit(1)


############################ Main functions ###################################

def setup_log():
    with open(validation_log, "w") as log:
        log.write(f"Performing baseline {args.assay_suffix.replace('_GL', '', 1)} V+V as per {V_V_guidelines_link}\n\n")
        command_run = " ".join(sys.argv)
        log.write(f"Validation program executed as:\n    {command_run}\n\n")


def append_message_to_log(message, one_return=False):
    with open(validation_log, "a") as log:

        if one_return:
            log.write(f"{message}\n")
        else:
            log.write(f"{message}\n\n")


def add_to_log_table(col1, col2):
    with open(validation_log, "a") as log:
        log.write(f"{col1}\t{col2}\n")


def report_success():
    print("")
    wprint(color_text("Validation has completed successfully :)", "green"))
    print(f"\n  Log written to: '{validation_log}'\n")

    with open(validation_log, "a") as log:
        log.write("   -----------------------------------------------------------------------------\n")
        log.write("                         Validation completed successfully." + "\n")
        log.write("   -----------------------------------------------------------------------------\n")


def report_present(message, color="green", log=True):
    print("")
    wprint(color_text(message, color))

    if log:
        with open(validation_log, "a") as log:
            log.write(f"{message}\n")
       

def check_expected_directories(expected):
    """ checks expected directories exist """

    for directory in expected:
        if not os.path.isdir(directory):

            report_failure(f"The directory '{directory}' was expected but not found.")
        else:
            report_present(f"    - directory '{directory}' found.")


def read_samples(file_path):
    """ reading unique sample names from runsheet into list """
    sample_names = []

    with open(file_path) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if reader.line_num > 1:
                sample_names.append(row[0])

    return sample_names


def fastqc_info_from_multiqc_json(multiqc_dir, multiqc_zip, count_key="total_sequences"):
    zip_file = zipfile.ZipFile(os.path.join(multiqc_dir, multiqc_zip))

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
    # assume the first index is FastQC if nothing is named FastQC (case of align_and_bismark_multiqc report)
    fastqc_index = 0
    fastqc_key = list(json_data['report_data_sources'].keys())[0]  # default to first key

    for i, key in enumerate(json_data['report_data_sources'].keys()):
        if key == 'FastQC':
            fastqc_index = i
            fastqc_key = key
            break
    
    # Get the first subsection regardless of its name (all_sections, genome_results, ...)
    first_subsection = list(json_data['report_data_sources'][fastqc_key].keys())[0]

    return (json_data['report_data_sources'][fastqc_key][first_subsection].keys(),
            [int(d.get(count_key)) for d in json_data['report_general_stats_data'][fastqc_index].values()])


def check_multiqc_samples(sample_names, file_prefixes_in_multiqc, per_read=True):
    """ ensures all expected samples are present in MultiQC outputs
    per_read = True  -> expect sample_R1 (FastQC) and sample_R2 (if paired-end)
    per_read = False -> expect sample only (alignment/Bismark)
    """

    if per_read:
        if not args.single_ended:
            R1_suffix = raw_R1_suffix.split(".")[0].replace("_raw", "")
            R2_suffix = raw_R2_suffix.split(".")[0].replace("_raw", "")

            for sample in sample_names:
                if sample + args.assay_suffix + R1_suffix not in file_prefixes_in_multiqc:
                    report_failure(f"The multiqc output is missing the expected '{sample + args.assay_suffix + R1_suffix}' entry.")
                if sample + args.assay_suffix + R2_suffix not in file_prefixes_in_multiqc:
                    report_failure(f"The multiqc output is missing the expected '{sample + args.assay_suffix + R2_suffix}' entry.")

        else:
            suffix = raw_suffix.split(".")[0].replace("_raw", "")

            for sample in sample_names:
                if sample + args.assay_suffix + suffix not in file_prefixes_in_multiqc:
                    report_failure(f"The multiqc output is missing the expected '{sample + args.assay_suffix + suffix}' entry.")
    else:
        for sample in sample_names:
            if sample + args.assay_suffix not in file_prefixes_in_multiqc:
                report_failure(f"The multiqc output is missing the expected '{sample + args.assay_suffix}' entry.")


def check_for_file_and_contents(file_path):
    """ used by various functions """

    if not os.path.exists(file_path):
        report_failure(f"The expected file '{file_path}' does not exist.")
    if not os.path.getsize(file_path) > 0:
       report_failure(f"The file '{file_path}' is empty.")


def check_fastq_files(sample_names):
    """ makes sure all expected read fastq files exist and hold something """

    for sample in sample_names:

        # if paired-end
        if not args.single_ended:

            # raw
            if args.include_raw_fastq:
                check_for_file_and_contents(os.path.join(raw_dir, sample + args.assay_suffix + raw_R1_suffix))
                check_for_file_and_contents(os.path.join(raw_dir, sample + args.assay_suffix + raw_R2_suffix))

            # quality filtered
            check_for_file_and_contents(os.path.join(trimmed_dir, sample + args.assay_suffix + trimmed_R1_suffix))
            check_for_file_and_contents(os.path.join(trimmed_dir, sample + args.assay_suffix + trimmed_R2_suffix))

        # if single-end
        else:

            # raw
            if args.include_raw_fastq:
                check_for_file_and_contents(os.path.join(raw_dir, sample + args.assay_suffix + raw_suffix))

            # trimmed
            check_for_file_and_contents(os.path.join(trimmed_dir, sample + args.assay_suffix + trimmed_suffix))
            
            
def get_files_in_dir(dir_path):
    return [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]


def check_dir_for_file_patterns(directory, expected_suffixes, sample_names, expected_prefixes=None, rna=False):
    # files directly under top directory
    all_files = get_files_in_dir(directory)

    # add files from one-level-down subdirectories
    for d in os.listdir(directory):
        sub = os.path.join(directory, d)
        
        if os.path.isdir(sub):
            all_files.extend(get_files_in_dir(sub))

    if rna:
        bismark_suffix = "_bismark_hisat2"
    else:
        bismark_suffix = "_bismark_bt2"

    for sample_ID in sample_names:        
        curr_files = [file for file in all_files if file.startswith(sample_ID)]

        # ---- suffix checks ----
        for suffix in expected_suffixes:
            if not any(f.endswith(suffix) for f in curr_files):
                report_failure(
                    f"An expected file ending with '{suffix}' for sample '{sample_ID}' "
                    f"was not found in the {directory} directory."
                )

        # ---- prefix checks ----
        if expected_prefixes is not None:
            for prefix in expected_prefixes:
                curr_prefix = f"{prefix}_context_{sample_ID + args.assay_suffix + bismark_suffix}"

                if not any(f.startswith(curr_prefix) for f in all_files):
                    report_failure(
                        f"An expected file beginning with '{curr_prefix}' for sample '{sample_ID}' "
                        f"was not found in the {directory} directory."
                    )


def check_dir_for_files(expected_files, directory):
    # getting all files in target dir
    all_files = get_files_in_dir(directory)

    for file in expected_files:

        if file not in all_files:
            report_failure(f"The expected file '{file}' was not found in the {directory} directory.")


def gen_stats(list_of_ints):
    """ returns min, max, mean, median of input integer list """

    min_val = min(list_of_ints)
    max_val = max(list_of_ints)

    mean_val = round(mean(list_of_ints), 2)
    median_val = int(median(list_of_ints))
    stdev_val = round(stdev(list_of_ints), 2)

    return min_val, max_val, mean_val, median_val, stdev_val


def get_read_count_stats(read_counts_list, read_type="Raw"):
    """ Summarizes read counts (min, max, mean, median) """

    (min_val, max_val, mean_val, median_val, stdev_val) = gen_stats(read_counts_list)

    print(f"\n  {read_type} read count summary:")
    print(f"    {'Min:':<10}{min_val:>0}")
    print(f"    {'Max:':<10}{max_val:>0}")
    print(f"    {'Mean:':<10}{mean_val:>0}")
    print(f"    {'Median:':<10}{median_val:>0}")
    print(f"    {'Std Dev:':<10}{stdev_val:>0}")
    print("")

    with open(validation_log, "a") as log:
        log.write(f"\n  {read_type} read count summary:")
        log.write(f"\n    {'Min:':<10}{min_val:>0}")
        log.write(f"\n    {'Max:':<10}{max_val:>0}")
        log.write(f"\n    {'Mean:':<10}{mean_val:>0}")
        log.write(f"\n    {'Median:':<10}{median_val:>0}")
        log.write(f"\n    {'Std Dev:':<10}{stdev_val:>0}")


# variable setup #

# universal settings

V_V_guidelines_link = "GeneLab's validation and verification protocol" # Dummy replacement of link


# top level folders
genelab_dir = f"{args.outdir}/GeneLab"
raw_dir = f"{args.outdir}/Merged_Sequence_Data"
trimmed_dir = f"{args.outdir}/Trimmed_Sequence_Data"
bismark_alignments_dir = f"{args.outdir}/Aligned_Sequence_Data"
bismark_meth_calls_dir = f"{args.outdir}/Methylation_Call_Data"
methylkit_dir = f"{args.outdir}/Differential_Methylation_Analysis_Data"

# subfolders
raw_multiqc_dir = os.path.join(raw_dir, "MultiQC_Reports")

trimming_reports_dir = os.path.join(trimmed_dir, "Trimming_Reports")
trimmed_multiqc_dir = os.path.join(trimmed_dir, "MultiQC_Reports")

bismark_alignment_reports_dir = os.path.join(bismark_alignments_dir, "Alignment_Reports")
if args.deduped:
    bismark_dedupe_reports_dir = os.path.join(bismark_alignments_dir, "Deduplication_Reports")

bismark_meth_call_reports_dir = os.path.join(bismark_meth_calls_dir, "Methylation_Call_Reports")
bismark_summary_dir = os.path.join(bismark_meth_call_reports_dir, "Methylation_Call_Summary_Reports")
align_multiqc_dir = os.path.join(bismark_meth_calls_dir, "MultiQC_Reports")


# file prefix setting (if present)
output_prefix = str(args.output_prefix)

# filename definitions
processing_zip_file = f"{genelab_dir}/{output_prefix}processing_info{args.assay_suffix}.zip"
readme_file = f"{genelab_dir}/{output_prefix}README{args.assay_suffix}.txt"

raw_multiqc_zip = f"{output_prefix}raw_multiqc{args.assay_suffix}_data.zip"
raw_multiqc_html = f"{output_prefix}raw_multiqc{args.assay_suffix}_report.html"
trimmed_multiqc_zip = f"{output_prefix}trimmed_multiqc{args.assay_suffix}_data.zip"
trimmed_multiqc_html = f"{output_prefix}trimmed_multiqc{args.assay_suffix}_report.html"
align_multiqc_zip = f"{output_prefix}align_and_bismark_multiqc{args.assay_suffix}_data.zip"
align_multiqc_html = f"{output_prefix}align_and_bismark_multiqc{args.assay_suffix}_report.html"
raw_suffix = "_raw.fastq.gz"
raw_R1_suffix = "_R1_raw.fastq.gz"
raw_R2_suffix = "_R2_raw.fastq.gz"
trimmed_suffix = "_trimmed.fastq.gz"
trimmed_R1_suffix = "_R1_trimmed.fastq.gz"
trimmed_R2_suffix = "_R2_trimmed.fastq.gz"

if args.output == "":
    validation_log = f"{args.GLDS_ID}-{output_prefix}methylseq-validation.log"
else:
    validation_log = args.output

# expected directories to check
expected_dirs = [genelab_dir, raw_multiqc_dir,
                trimmed_dir, trimmed_multiqc_dir, trimming_reports_dir,
                bismark_alignments_dir, bismark_alignment_reports_dir,
                bismark_meth_calls_dir, bismark_meth_call_reports_dir, bismark_summary_dir, align_multiqc_dir,
                methylkit_dir]
if args.include_raw_fastq:
    expected_dirs.insert(1, raw_dir)
if args.deduped:
    expected_dirs.insert(7, bismark_dedupe_reports_dir)

# expected file types
bismark_alignment_files_expected_suffixes = [".nucleotide_stats.txt", "_qualimap.zip", "_sorted.bam", "_sorted.bam.bai",
                                            "_report.txt"]

if args.deduped:
    bismark_alignment_files_expected_suffixes.extend(["_sorted.deduplicated.bam", "_sorted.deduplicated.bam.bai", 
                                                      ".deduplication_report.txt"])

bismark_meth_call_files_expected_suffixes = [".bedGraph.gz", ".bismark.cov.gz", ".M-bias.txt", "_splitting_report.txt",
                                            "cytosine_context_summary.txt", "_report.html"]
bismark_meth_call_files_expected_prefixes = ["CHG", "CHH", "CpG"]

bismark_summary_expected_files = [f"bismark_summary_report{args.assay_suffix}.html",
                                f"bismark_summary_report{args.assay_suffix}.txt"]

methylkit_expected_files = [f"differential_methylation_tiles{args.assay_suffix}.csv",
                            f"differential_methylation_bases{args.assay_suffix}.csv",
                            f"contrasts{args.assay_suffix}.csv",
                            f"SampleTable{args.assay_suffix}.csv"]

def main():

    # initializing log file
    setup_log()

    append_message_to_log(f"Summary of checks:")

    if 'RNA' in args.assay_suffix:
        rna = True
    else:
        rna = False

    check_for_file_and_contents(readme_file)
    append_message_to_log(f"    - populated {readme_file} detected")

    append_message_to_log("    - checking directories")
    check_expected_directories(expected=expected_dirs)
    append_message_to_log("\n    - all expected directories found.")

    sample_names = read_samples(args.runsheet)

    (samples_in_raw_fastqc, raw_fastqc_counts) = fastqc_info_from_multiqc_json(raw_multiqc_dir, raw_multiqc_zip)
    (samples_in_trimmed_fastqc, trimmed_fastqc_counts) = fastqc_info_from_multiqc_json(trimmed_multiqc_dir,
                                                                                         trimmed_multiqc_zip)
    (samples_in_align_qc, align_qc_counts) = fastqc_info_from_multiqc_json(align_multiqc_dir,
                                                                                         align_multiqc_zip, 
                                                                                         count_key="total_reads")

    append_message_to_log(f"    - checking samples in multiqc data in '{raw_multiqc_zip}'")
    check_multiqc_samples(sample_names, samples_in_raw_fastqc)
    append_message_to_log(f"    - all expected samples found in raw multiqc")

    append_message_to_log(f"    - checking samples in multiqc data in '{trimmed_multiqc_zip}'")
    check_multiqc_samples(sample_names, samples_in_trimmed_fastqc)
    append_message_to_log(f"    - all expected samples found in trimmed multiqc")

    append_message_to_log(f"    - checking samples in multiqc data in '{align_multiqc_zip}'")
    check_multiqc_samples(sample_names, samples_in_align_qc, per_read=False)
    append_message_to_log(f"    - all expected samples found in alignment and Bismark multiqc")

    append_message_to_log(f"    - checking fastq files")
    check_fastq_files(sample_names)
    if args.include_raw_fastq:
        append_message_to_log(f"    - all expected read files found in '{raw_dir}'")
    append_message_to_log(f"    - all expected read files found in '{trimmed_dir}'")

    # checking for all expected files and folders in bismark alignments directory
    check_dir_for_file_patterns(bismark_alignments_dir, bismark_alignment_files_expected_suffixes, sample_names,
                                rna=rna)
    append_message_to_log(
        f"    - all expected bismark alignment files found under the '{bismark_alignments_dir}' directory")

    # checking for all expected files and folders in bismark meth calls directory
    check_dir_for_file_patterns(bismark_meth_calls_dir, bismark_meth_call_files_expected_suffixes, sample_names,
                                expected_prefixes=bismark_meth_call_files_expected_prefixes, rna=rna)

    # checking for all expected files in bismark summary dir under bismark meth calls directory
    check_dir_for_files(bismark_summary_expected_files, bismark_summary_dir)

    append_message_to_log(f"    - all expected Bismark files found under the '{bismark_meth_calls_dir}' directory")

    # checking for all expected base outputs in methylkit dir
    check_dir_for_files(methylkit_expected_files, methylkit_dir)
    append_message_to_log(f"    - all expected differential methylation files found in the '{methylkit_dir}' directory",
                          one_return=True)

    # checking for processing_info.zip
    check_for_file_and_contents(processing_zip_file)
    append_message_to_log(f"\n    - populated {processing_zip_file} detected")

    report_success()

    get_read_count_stats(raw_fastqc_counts, read_type="Raw")
    get_read_count_stats(trimmed_fastqc_counts, read_type="Trimmed")
    get_read_count_stats(align_qc_counts, read_type="Aligned")

    # add newline at end of log file
    append_message_to_log("")


if __name__ == "__main__":
    main()
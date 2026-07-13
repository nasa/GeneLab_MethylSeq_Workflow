#!/usr/bin/env python

"""
This is a program for generating a README.txt file for GeneLab processed MethylSeq datasets.
"""

import os
import sys
import argparse
import textwrap
import zipfile
import datetime

parser = argparse.ArgumentParser(
    description="This program generates the corresponding README file for GeneLab processed Methylseq datasets.")

required = parser.add_argument_group('required arguments')

required.add_argument("-a", "--assay_suffix", choices=['_GLMethylSeq', '_GLRNAMethylSeq'],
                      help="Specifies what assay suffix to use. Parsed from assay_suffix in ch_meta", action="store", required=True)
required.add_argument("--osd-id", 
                      help='OSDR osd_id ID (e.g. "OSD-69")', action="store", required=True)
parser.add_argument("--output", "--output-name", default="",
                    help='Name of output file')
parser.add_argument("-p", "--output-prefix", action="store", default="",
                    help="Output additional file prefix if there is one")
parser.add_argument("--paired-end", action="store_true", 
                    help="Passed only when data is paired-end, parsed from paired_end in ch_meta")
parser.add_argument("--dedupe-skipped", action="store_true",
                    help="Passed only when deduplication is skipped (for RRBS from ch_meta or user choice)")
parser.add_argument("--include-raw-reads", action="store_true",
                    help="Passed only when raw read data (Merged sequence data) should be included")
parser.add_argument("--name", default="", required=True,
                    help='Name of individual who performed the processing (default: "")')
parser.add_argument("--email", default="", required=True,
                    help='Email address of individual who performed the processing (default: "")')
parser.add_argument("--protocol-ID", default="",
                    help='Protocol document ID followed (default: assay dependent)')
parser.add_argument("--date", action="store", default=datetime.date.today(), type=datetime.date.fromisoformat,
                    help="Date the processed data was generated in ISO format (YYYY-MM-DD). (default: today's date)")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_args()

latest_methylseq_DPPD = "GL-DPPD-7113"


################################################################################

def main():
    with open(output_file, "w") as output:
        write_header(output, args.osd_id, args.name, args.email, args.protocol_ID, args.date)

        write_methylseq_body(output)


################################################################################

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

    print(textwrap.fill(text, width=80, initial_indent="  ",
                        subsequent_indent="  ", break_on_hyphens=False))


def report_failure(message, color="yellow"):
    print("")
    wprint(color_text(message, color))
    print("\nREADME-generation failed.\n")

    sys.exit(1)


def check_for_file_and_contents(file_path):
    """ used by get_processing_zip_contents function """

    if not os.path.exists(file_path):
        report_failure(f"The expected file '{file_path}' does not exist.")
    if not os.path.getsize(file_path) > 0:
        report_failure("The file '{file_path}' is empty.")


def get_processing_zip_contents():
    """ this gets the filenames that are in the processing_info.zip to add them to the readme """

    check_for_file_and_contents(processing_zip_file)

    with zipfile.ZipFile(processing_zip_file) as zip_obj:
        entries = zip_obj.namelist()
        entries.sort()

    return entries


def write_header(output, osd_id, name, email, protocol_id, date=datetime.date.today()):
    if not args.osd_id.startswith("OSD-"):
        print("OSDR osd_id must start with 'OSD-'")
        sys.exit(1)

    header = ["################################################################################\n",
              "{:<77} {:>0}".format(f"## This directory holds processed data for NASA {osd_id}", "##\n"),
              "{:<77} {:>0}".format(f"## https://osdr.nasa.gov/bio/repo/data/studies/{osd_id}/", "##\n"),
              "{:<77} {:>0}".format("##", "##\n"),
              "{:<77} {:>0}".format(f"## Processed by {name} ({email}) on {date}", "##\n"),
              "{:<77} {:>0}".format(f"## Based on {protocol_id}", "##\n"),
              "################################################################################\n\n",
              "Summary of contents:\n\n"]

    output.writelines(header)


# up_and_left: '┌',
# up_and_right: '┐',
# down_and_left: '└',
# down_and_right: '┘',
# vertical: '│',
# horizontal: '─',
# vertical_and_horizontal: '┼',
# down_and_horizontal: '┬',
# up_and_horizontal: '┴',
# top_connection: None,
# bottom_connection: None,
# HorizT='├'
# HLine='─'


def format_string(file_name, file_description, max_offset, level, connect="top", is_terminal=False):
    """
    Return a formatted string for a directory structure with connectors similar to unix "tree" command.

    Args:
        file_name (str): a file or directory name
        file_description (str): a description of the file
        max_offset (int): maximum line length provided to ensure that all formatted strings 
            produced will space the descriptions at the same level.
        level (int): Position in the directory tree
        connect (str, optional): Type of tree connector. "top" indicates an elbow connector that 
            connects to the row above. "top_bottom" indicates a "T" connector that connects to the 
            rows above and below. Defaults to "top".
        is_terminal (bool, optional): Thsi directory/folder is a terminal leaf in the directory 
            tree, do not connect it to upper levels. Defaults to False.

    Returns:
        str: A formatted string with a tree connector at the start, a filename, and a file description
    """
    if connect == "top_bottom":
        start_type = '├──'
    elif connect == "top":
        start_type = '└──'
    else:
        start_type = "───"

    if level == 1:
        s = '   {}'.format(start_type)
    elif level == 2:
        s = '   │   {}'.format(start_type)
    elif level == 3:
        if is_terminal:
            s = '   │       {}'.format(start_type)
        else:
            s = '   │   │   {}'.format(start_type)
    elif level == 4:
        if is_terminal:
            s = '   │           {}'.format(start_type)
        else:
            s = '   │   │       {}'.format(start_type)
    else:
        s = ''
    # NOTE: length of the indent that was actually added needs to be subtracted 
    # from the max padding to get the correct padding value for the current row
    return "{start} {file:<{spacing}}   - {desc:>0}\n".format(start=s, file=file_name, desc=file_description,
                                                              spacing=max_offset - len(s))


def add_level_one(file_name, file_description, max_offset, output):
    output.write(format_string(file_name, file_description, max_offset, 
                               level=1, connect="top_bottom", is_terminal=False))


def add_level_one_last(file_name, file_description, max_offset, output):
    output.write(format_string(file_name, file_description, max_offset, 
                               level=1, connect="top", is_terminal=True))


def add_level_two(file_name, file_description, max_offset, output):
    output.write(format_string(file_name, file_description, max_offset, 
                               level=2, connect="top_bottom", is_terminal=False))


def add_level_two_last(file_name, file_description, max_offset, output):
    output.write(format_string(file_name, file_description, max_offset,
                               level=2, connect="top", is_terminal=True))


def add_level_three(file_name, file_description, max_offset, output):
    output.write(format_string(file_name, file_description, max_offset,
                               level=3, connect="top_bottom", is_terminal=False))


def add_level_three_mid(file_name, file_description, max_offset, output):
    output.write(format_string(file_name, file_description, max_offset,
                               level=3, connect="top_bottom", is_terminal=True))


def add_level_three_mid_last(file_name, file_description, max_offset, output):
    output.write(format_string(file_name, file_description, max_offset,
                               level=3, connect="top", is_terminal=False))


def add_level_three_last(file_name, file_description, max_offset, output):
    output.write(format_string(file_name, file_description, max_offset,
                               level=3, connect="top", is_terminal=True))


def add_level_four(file_name, file_description, max_offset, output):
    output.write(format_string(file_name, file_description, max_offset,
                               level=4, connect="top_bottom", is_terminal=False))


def add_level_four_mid_last(file_name, file_description, max_offset, output):
    output.write(format_string(file_name, file_description, max_offset,
                               level=4, connect="top", is_terminal=False))


def add_spacer(output):
    output.write("   │\n")


def write_methylseq_body(output):
    longest_filename = f"align_and_bismark_multiqc{args.assay_suffix}_data.zip"
    # length of padding is the length of the longest file + the indent_level + some extra
    pad = len(longest_filename) + (4 * 3) + 2

    # this file
    add_level_one(output_file, "this file", pad, output)

    add_spacer(output)

    # Merged Sequence Data File
    # raw reads (usually "Merged sequence data")
    if args.include_raw_reads:
        add_level_one(merged_reads_dir, "initial read fastq files", pad, output)
        if args.paired_end:
            add_level_two("*_R1_raw.fastq.gz", "read1 fastq files", pad, output)
            add_level_two("*_R2_raw.fastq.gz", "read2 fastq files", pad, output)
        else:
            add_level_two("*_raw.fastq.gz", "fastq files", pad, output)
        # Merged Sequence Data/MultiQC Reports
        add_level_two_last(multiqc_dir, "multiQC summary reports of raw FastQC runs", pad, output)
        add_level_three_mid(f"raw_multiqc{args.assay_suffix}.html", "multiQC FastQC summary report", pad, output)
        add_level_three_last(f"raw_multiqc{args.assay_suffix}_data.zip", "multiQC FastQC summary report data", pad, output)

    add_spacer(output)

    # quality-filtered and trimmed reads
    # Trimmed Sequence Data
    add_level_one(filtered_reads_dir, "quality-filtered and trimmed fastq files", pad, output)
    if args.paired_end:
        add_level_two("*_R1_trimmed.fastq.gz", "read1 trimmed fastq files", pad, output)
        add_level_two("*_R2_trimmed.fastq.gz", "read2 trimmed fastq files", pad, output)
    else:
        add_level_two("*_trimmed.fastq.gz", "trimmed fastq files", pad, output)
    # Trimmed Sequence Data/Trimming Reports
    add_level_two(trimming_reports_dir, "per sample trimming reports", pad, output)
    if args.paired_end:
        add_level_three("*_R1_trimming_report.txt", "read1 trimming report", pad, output)
        add_level_three_mid_last("*_R2_trimming_report.txt", "read2 trimming report", pad, output)
    else:
        add_level_three_mid_last("*_trimming_report.txt", "trimming report", pad, output)
    # Trimmed Sequence Data/MultiQC Reports
    add_level_two_last(multiqc_dir, "multiQC summary reports of trimming statistics and filtered and trimmed FastQC runs", pad, output)
    add_level_three_mid(f"trimmed_multiqc{args.assay_suffix}.html", "multiQC trimming summary report", pad, output)
    add_level_three_last(f"trimmed_multiqc{args.assay_suffix}_data.zip", "multiQC trimming summary report data", pad, output)

    add_spacer(output)

    # bismark alignment files
    # Aligned Sequence Data
    add_level_one(bismark_alignments_dir, "bismark alignment files", pad, output)
    add_level_two("*_sorted.bam", "bam files", pad, output)
    add_level_two("*_sorted.bam.bai", "bam index files", pad, output)
    if not args.dedupe_skipped:
        add_level_two("*.deduplicated_sorted.bam", "deduplicated bam files", pad, output)
        add_level_two("*.deduplicated_sorted.bam.bai", "deduplicated bam index files", pad, output)
    # Aligned Sequence Data/Alignment Reports
    if args.dedupe_skipped:
        add_level_two_last(bismark_alignment_report_dir, "bismark alignment reports and qualimap reports", pad, output)
        add_level_three_mid("*.nucleotide_stats.txt", "genome-wide nucleotide statistics", pad, output)
        add_level_three_mid("*_report.txt", "bismark alignment report", pad, output)
        add_level_three_last("*qualimap.zip", "qualimap alignment report", pad, output)
    else:
        add_level_two(bismark_alignment_report_dir, "bismark alignment reports and qualimap reports", pad, output)
        add_level_three("*.nucleotide_stats.txt", "genome-wide nucleotide statistics", pad, output)
        add_level_three("*_report.txt", "bismark alignment report", pad, output)
        add_level_three_mid_last("*qualimap.zip", "qualimap alignment report", pad, output)

    # deduplication files
    if not args.dedupe_skipped:
        # Aligned Sequence Data/Deduplication Reports
        add_level_two_last(bismark_deduplication_report_dir, "bismark deduplication reports", pad, output)
        add_level_three_last("*.deduplication_report.txt", "bismark deduplication report", pad, output)

    add_spacer(output)

    # bismark methylation calls
    # "Methylation Call Data"
    add_level_one(bismark_meth_calls_dir, "methylation-call files from bismark_methylation_extractor", pad, output)
    add_level_two("CHG_context_*.txt.gz", "CHG context methylation calls", pad, output)
    add_level_two("CHH_context_*.txt.gz", "CHH context methylation calls", pad, output)
    add_level_two("CpG_context_*.txt.gz", "CpG context methylation calls", pad, output)
    add_level_two("*.bedGraph.gz", "methylation percentage in bedgraph format", pad, output)
    add_level_two("*.bismark.cov.gz", "read coverage and methylation percentage of methylated and unmethylated reads", pad, output)
    add_level_two("*.CpG_report.txt.gz", "cytosine methylation report for cytosines in CpG context", pad, output)
    # "Methylation Call Data/Methylation Call Reports"
    add_level_two(bismark_meth_call_reports_dir, "methylation-call reports", pad, output)
    add_level_three("*.cytosine_context_summary.txt", "cytosine context summary", pad, output)
    add_level_three("*.M-bias.txt", "methylation bias report", pad, output)
    add_level_three("*_splitting_report.txt", "splitting report", pad, output)
    add_level_three("*_bismark_*_report.html", "bismark invidual sample report", pad, output)
    # bismark summary files
    # Methylation Call Data/Methylation Call Reports/Methylation Call Summary Reports
    add_level_three_mid_last(bismark_summary_dir, "bismark summary reports and html files", pad, output)
    add_level_four(f"bismark_summary_report{args.assay_suffix}.html", "bismark summary report html", pad, output)
    add_level_four_mid_last(f"bismark_summary_report{args.assay_suffix}.txt", "bismark summary report data", pad, output)
    # Methylation Call Data/MultiQC Reports
    add_level_two_last(multiqc_dir, "multiQC summary reports of bismark alignment and methylation calls", pad, output)
    add_level_three_mid(f"align_and_bismark_multiqc{args.assay_suffix}.html", "multiQC bismark and alignment summary report", pad, output)
    add_level_three_last(f"align_and_bismark_multiqc{args.assay_suffix}_data.zip", "multiQC bismark and alignment summary report data", pad, output)

    add_spacer(output)

    # methylkit outputs
    # Differential Methylation Analysis Data
    add_level_one(methylkit_outputs_dir, "MethylKit differential methylation outputs", pad, output)
    add_level_two(f"SampleTable{args.assay_suffix}.csv", "list of samples and conditions", pad, output)
    add_level_two(f"contrasts{args.assay_suffix}.csv", "list of contrasts", pad, output)
    add_level_two(f"differential_methylation_bases{args.assay_suffix}.csv", "differential methylation at each base with sufficient coverage", pad, output)
    add_level_two_last(f"differential_methylation_tiles{args.assay_suffix}.csv", "differential methylation between tiled genomic regions", pad, output)

    add_spacer(output)

    # processing info
    add_level_one_last(processing_zip_file, "zip archive holding info related to processing (workflow files and metadata)", pad, output)

    output.write("\n")


# variable setup #

# universal settings
output_prefix = str(args.output_prefix)

processing_zip_file = f"{output_prefix}processing_info{args.assay_suffix}.zip"

merged_reads_dir = "Merged Sequence Data/"
multiqc_dir = "MultiQC Reports/"
filtered_reads_dir = "Trimmed Sequence Data/"
trimming_reports_dir = "Trimming Reports/"
bismark_alignments_dir = "Aligned Sequence Data/"
bismark_alignment_report_dir = "Alignment Reports/"
bismark_deduplication_dir = "Deduplicated Data/"
bismark_deduplication_report_dir = "Deduplication Reports/"
bismark_meth_calls_dir = "Methylation Call Data/"
bismark_meth_call_reports_dir = "Methylation Call Reports/"
bismark_summary_dir = "Methylation Call Summary Reports/"
methylkit_outputs_dir = "Differential Methylation Analysis Data/"

if args.protocol_ID == "":
    args.protocol_ID = latest_methylseq_DPPD

if args.output == "":
    output_file = f"{output_prefix}README{args.assay_suffix}.txt"
else:
    output_file = args.output

if __name__ == "__main__":
    main()
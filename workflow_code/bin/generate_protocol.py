#!/usr/bin/env python
"""
This script generates a protocol text file for GeneLab Methyl-seq data processing.
It reads software versions from a YAML file generated upstream.
"""

import argparse
import yaml
import os
import sys
from datetime import datetime
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description='Generate protocol file for GeneLab Methyl-Seq pipeline')
    parser.add_argument('--mode', choices=['DNA', 'RNA'], default='DNA',
                        help='Processing mode (DNA methylation or RNA methylation)')
    parser.add_argument('--outdir', required=True,
                        help='Output directory for the protocol file')
    parser.add_argument('--software_table', required=True,
                        help='Path to YAML file containing software versions')
    parser.add_argument('--assay_suffix', default='',
                        help='Suffix for the assay type')
    parser.add_argument('--workflow_version', default='unknown',
                        help='Version of the NF_RCP workflow manifest')
    parser.add_argument('--organism', required=False,
                        help='Organism name for the reference genome')
    parser.add_argument('--dedupe_skipped', action="store_true",
                        help='Passed only when deduplication is skipped (for RRBS or user choice)')
    parser.add_argument('--reference_source', required=False,
                        help='Source of the reference genome')
    parser.add_argument('--reference_version', required=False,
                        help='Version of the reference genome')
    parser.add_argument('--reference_fasta', required=True,
                        help='Path to the reference genome FASTA file')
    parser.add_argument('--reference_gtf', required=True,
                        help='Path to the reference genome GTF file')
    return parser.parse_args()

def read_software_versions(yaml_file):
    try:
        with open(yaml_file, 'r') as f:
            return yaml.safe_load(f)
    except Exception as e:
        sys.stderr.write(f"Error reading software versions file: {e}\n")
        sys.exit(1)

def generate_protocol_content(args, software_versions):
    # Get current date
    current_date = datetime.now().strftime("%Y-%m-%d")
    
    # Create header
    header = f"# GeneLab Methyl-Seq Pipeline Protocol{args.assay_suffix}\n"
    header += f"# Date: {current_date}\n\n"
    
    # Start building the description as a single paragraph
    header += "Data were processed as described in GL-DPPD-7113 "
    header += "(https://github.com/nasa/GeneLab_Data_Processing/blob/master/Methyl-Seq/Pipeline_GL-DPPD-7113_Versions/GL-DPPD-7113.md), " # Verify link
    header += f"using NF_MSCP version {args.workflow_version} " # Verify NF_MSCP
    header += f"(https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MSCP-{args.workflow_version}/Methyl-Seq/Workflow_Documentation/NF_MSCP). " # Verify link
    
    # Add processing description with software versions
    trim_galore_version = software_versions.get('TrimGalore!', 'unknown')
    cutadapt_version = software_versions.get('Cutadapt', 'unknown')
    fastqc_version = software_versions.get('FastQC', 'unknown')
    multiqc_version = software_versions.get('MultiQC', 'unknown')
    
    description = f"In short, raw fastq files were filtered using Trim Galore! (version {trim_galore_version}) powered by Cutadapt (version {cutadapt_version}). "
    description += f"Trimmed fastq file quality was evaluated with FastQC (version {fastqc_version}), and MultiQC (version {multiqc_version}) was used to generate MultiQC reports. "
    
    # Add reference description based on reference source
    reference_description = ""
    
    # Format organism name if provided - replace underscores with spaces and title case
    organism_name = ""
    organism_name_italics = ""  # For use in text (with underscores for Jira)
    if hasattr(args, 'organism') and args.organism:
        organism_name = args.organism.replace('_', ' ').title()
        organism_name_italics = f"_{organism_name}_"  # Surrounded by underscores for Jira italics
    
    # Get reference FASTA and GTF basenames
    ref_fasta_name = ""
    ref_gtf_name = ""
    genome_assembly = ""
    if args.reference_fasta and args.reference_fasta != "null":
        ref_fasta_name = os.path.basename(args.reference_fasta)
        
        # Handle different reference naming conventions
        if args.reference_source and "ensembl" in args.reference_source.lower():
            # For Ensembl: extract after first dot (e.g., Bacillus_subtilis.ASM904v1.dna.toplevel.fa)
            if '.' in ref_fasta_name:
                parts = ref_fasta_name.split('.')
                if len(parts) > 1:
                    genome_assembly = parts[1]
        elif args.reference_source and "ncbi" in args.reference_source.lower():
            # For NCBI: extract only the assembly name (e.g., ASM746v2 from GCF_000007465.2_ASM746v2_genomic.fna.gz)
            if '_genomic' in ref_fasta_name:
                # Get the part right before _genomic
                parts = ref_fasta_name.split('_genomic')[0].split('_')
                # The assembly name is typically the last part before _genomic (after the GCF_accession_)
                if len(parts) > 2:
                    # Skip the GCF part and accession, take the assembly name
                    genome_assembly = parts[-1]
            # Fallback to old method if _genomic is not in the name
            elif '_' in ref_fasta_name:
                parts = ref_fasta_name.split('_')
                if len(parts) > 2:  # Should have at least 3 parts
                    # The assembly name is usually the second part after the GCF_accession
                    genome_assembly = parts[1]
                    # If there are more parts before "genomic", include them in assembly name
                    for i in range(2, len(parts)):
                        if "genomic" in parts[i]:
                            break
                        genome_assembly += "_" + parts[i]
        else:
            # For other sources, try to extract anything that looks like an assembly version
            if '_' in ref_fasta_name and '.' in ref_fasta_name:
                # First try the Ensembl pattern (after first dot)
                parts = ref_fasta_name.split('.')
                if len(parts) > 1:
                    genome_assembly = parts[1]
                
                # If that didn't work, try the NCBI pattern (after first underscore)
                if not genome_assembly:
                    parts = ref_fasta_name.split('_')
                    if len(parts) > 1:
                        genome_assembly = parts[1]
    
    if args.reference_gtf and args.reference_gtf != "null":
        ref_gtf_name = os.path.basename(args.reference_gtf)
    
    # Get software versions
    bismark_version = software_versions.get('Bismark', 'unknown')
    if args.mode == 'RNA':
        hisat2_version = software_versions.get('Hisat2', 'unknown')
    elif args.mode == 'DNA':
        bowtie2_version = software_versions.get('Bowtie2', 'unknown')
    
    # Check if reference source contains "ensembl" or "ncbi"
    is_ensembl = False
    is_ncbi = False
    ref_source_formatted = ""
    if args.reference_source:
        if "ensembl" in args.reference_source.lower():
            is_ensembl = True
            # Format ensembl source nicely (e.g., ensembl_bacteria -> Ensembl Bacteria)
            ref_source_formatted = args.reference_source.replace('_', ' ').title()
        elif "ncbi" in args.reference_source.lower():
            is_ncbi = True
            # Format NCBI source nicely (e.g., ncbi_refseq -> NCBI RefSeq)
            ref_source_formatted = args.reference_source.replace('_', ' ').upper()
    
    # Generate reference description
    reference_description = f"{organism_name_italics} Bismark references were built using Bismark (version {bismark_version}) and "
    if args.mode == 'RNA':
        reference_description += f"Hisat2 (version {hisat2_version})"
    elif args.mode == 'DNA':
        reference_description += f"Bowtie 2 (version {bowtie2_version})"
        
    # Add source and version information
    if is_ensembl and args.reference_version:
        reference_description += f", {ref_source_formatted} release {args.reference_version}"
    elif is_ncbi and args.reference_version:
        reference_description += f", {ref_source_formatted} version {args.reference_version}"
      
    # Add assembly information if available
    if genome_assembly:
        if is_ncbi:
            reference_description += f", NCBI genome assembly {genome_assembly}"
        else:
            reference_description += f", genome assembly {genome_assembly}"
         
    reference_description += f" ({ref_fasta_name}). "
    
    # Add reference description to the protocol
    description += reference_description
    
    # Add Bismark info based on mode and RRBS/non-RRBS data
    description += f"Trimmed reads were aligned to the {organism_name_italics} Bismark reference using Bismark (version {bismark_version}) and "
    if args.mode == 'RNA':
        description += f"Hisat2 (version {hisat2_version})"
    elif args.mode == 'DNA':
        description += f"Bowtie 2 (version {bowtie2_version})" 
    
    samtools_version = software_versions.get('SAMtools', 'unknown')
    qualimap_version = software_versions.get('Qualimap', 'unknown')
    
    description += f", the aligned data were sorted with samtools (version {samtools_version}), and alignment quality was assessed with QualiMap (version {qualimap_version}). "
    
    if args.dedupe_skipped:
        description += f"Bismark (version {bismark_version}) was used to extract methylation calls from aligned (unsorted) data, "
    else:
        description += f"Bismark (version {bismark_version}) was used to deduplicate the aligned (unsorted) data, then methylation calls were extracted from the deduplicated aligned data, "
    
    description += f"and bismark reports were generated for each sample with Bismark (version {bismark_version}). Bismark reports for each sample were combined into a bismark summary report with Bismark (version {bismark_version}). Alignment QualiMap reports and Bismark reports were compiled with MultiQC (version {multiqc_version}). "
    
    dp_tools_version = software_versions.get('dp_tools', 'unknown')
    r_version = software_versions.get('R', 'unknown')
    methylkit_version = software_versions.get('methylKit', 'unknown')
    genomation_version = software_versions.get('genomation', 'unknown')
    gtfToGenePred_version = software_versions.get('gtfToGenePred', 'unknown')
    genePredToBed_version = software_versions.get('genePredToBed', 'unknown')
    
    # Add runsheet generation sentence
    description += f"A runsheet containing sample group information was generated with dp_tools (version {dp_tools_version}) and imported to R (version {r_version}), and methylation call data was imported into R using methylkit (version {methylkit_version}). "
            
    # Add differential methylation analysis sentence
    description += f"Pair-wise differential methylation analysis was performed in R at both the base and tile level using methylkit (version {methylkit_version}). "
    
    # Add gene annotations sentence
    description += f"Gene annotations were assigned to the differential methylation results using genomation (version {genomation_version}) and custom annotation tables generated using UCSC gtfToGenePred (version {gtfToGenePred_version}) and UCSC genePredToBed (version {genePredToBed_version}) from the following gtf annotation file: {ref_gtf_name} "
    
	# Add source and version information
    if is_ensembl and args.reference_version:
        description += f"({ref_source_formatted} release {args.reference_version})."
    elif is_ncbi and args.reference_version:
        description += f"({ref_source_formatted} version {args.reference_version})."
    
    # Add paragraph break after the protocol text
    description += "\n"
    
    # Combine all sections
    content = header + description
    
    return content

def main():
    args = parse_args()
    
    # Read software versions from YAML file
    software_versions = read_software_versions(args.software_table)
    
    # Generate protocol content
    protocol_content = generate_protocol_content(args, software_versions)
    
    # Write to output file
    output_file = os.path.join(args.outdir, f"protocol{args.assay_suffix}.txt")
    try:
        with open(output_file, 'w') as f:
            f.write(protocol_content)
        print(f"Protocol file generated successfully: {output_file}")
    except Exception as e:
        sys.stderr.write(f"Error writing protocol file: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()

# GeneLab Methylation Sequencing Data Processing Workflow

> GeneLab, part of [NASA's Open Science Data Repository (OSDR)](https://www.nasa.gov/osdr), has wrapped each step of the Methylation Sequencing (Methyl-Seq) consensus pipeline ([MethylSeq](https://github.com/nasa/GeneLab_Data_Processing/tree/master/Methyl-Seq)) into a Nextflow workflow with validation and verification of output files built in after each step. This repository contains the Nextflow workflow code (NF_MethylSeq) along with instructions for installation and usage. Exact workflow run info and MethylCP version used to process specific datasets that have been released are available in the \*nextflow\_processing\_info.txt file on the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/), which can be found under 'Files' -> 'GeneLab Processed Methyl-Seq Files' -> 'Supplemental Materials'.

<br>

# NF\_MethylSeq Workflow Information and Usage Instructions

## General workflow info
The current GeneLab Methyl-Seq sequencing data processing pipeline (MethylSeq), [GL-DPPD-7113.md](https://github.com/nasa/GeneLab_Data_Processing/tree/master/Methyl-Seq/Pipeline_GL-DPPD-7113_Versions/GL-DPPD-7113.md), is implemented as a [Nextflow](https://nextflow.io/) DSL2 workflow and can utilize images managed through either [Docker](https://docs.docker.com/get-started/) or [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html), or [conda])https://docs.conda.io/en/latest/) for environments (GeneLab uses Singularity when processing GLDS datasets). This workflow (NF_MethylSeq) is run using the command line interface (CLI) of any unix-based system. While knowledge of creating workflows in Nextflow is not required to run the workflow as is, [the Nextflow documentation](https://nextflow.io/docs/latest/index.html) is a useful resource for users who want to modify and/or extend this workflow.   

## Utilizing the Workflow

1. [Install Nextflow and Singularity](#1-install-nextflow-and-singularity)  
   1a. [Install Nextflow](#1a-install-nextflow)  
   1b. [Install Singularity](#1b-install-singularity)
2. [Download the Workflow Files](#2-download-the-workflow-files)  
3. [Fetch Singularity Images](#3-fetch-singularity-images)  
4. [Run the Workflow](#4-run-the-workflow)  
   4a. [Approach 1: Run the workflow on a GeneLab MethylSeq dataset with automatic retrieval of reference fasta and gtf files](#4a-approach-1-run-the-workflow-on-a-genelab-methylseq-dataset-with-automatic-retrieval-of-reference-fasta-and-gtf-files)  
   4b. [Approach 2: Run the workflow on a GeneLab MethylSeq dataset with custom reference fasta and gtf files](#4b-approach-2-run-the-workflow-on-a-genelab-methylseq-dataset-with-custom-reference-fasta-and-gtf-files)  
   4c. [Approach 3: Run the workflow on a non-GeneLab dataset using a user-created runsheet with automatic retrieval of reference fasta and gtf files](#4c-approach-3-run-the-workflow-on-a-non-genelab-dataset-using-a-user-created-runsheet-with-automatic-retrieval-of-reference-fasta-and-gtf-files)  
   4d. [Approach 4: Run the workflow on a non-GeneLab dataset using a user-created runsheet with custom reference fasta and gtf files](#4d-approach-4-run-the-workflow-on-a-non-genelab-dataset-using-a-user-created-runsheet-with-custom-reference-fasta-and-gtf-files)
5. [Additional Output Files](#5-additional-output-files)  

<br>

---

### 1. Install Nextflow and Singularity 

#### 1a. Install Nextflow

Nextflow can be installed either through [Anaconda](https://anaconda.org/bioconda/nextflow) or as documented on the [Nextflow documentation page](https://www.nextflow.io/docs/latest/getstarted.html).

> Note: If you want to install Anaconda, we recommend installing a Miniforge version appropriate for your system, as documented on the [conda-forge website](https://conda-forge.org/download/), where you can find basic binaries for most systems. More detailed miniforge documentation is available in the [miniforge github repository](https://github.com/conda-forge/miniforge).
> 
> Once conda is installed on your system, you can install the latest version of Nextflow by running the following commands:
> 
> ```bash
> conda install -c bioconda nextflow
> nextflow self-update
> ```

<br>

#### 1b. Install Singularity

Singularity is a container platform that allows usage of containerized software. This enables the GeneLab MSCP workflow to retrieve and use all software required for processing without the need to install the software directly on the user's system.

We recommend installing Singularity on a system wide level as per the associated [documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).

> Note: Singularity is also available through [Anaconda](https://anaconda.org/conda-forge/singularity).

> Note: Alternatively, Docker can be used in place of Singularity. See the [Docker CE installation documentation](https://docs.docker.com/engine/install/).

<br>

---

### 2. Download the Workflow Files

All files required for utilizing the NF_MSCP GeneLab workflow for processing MethylSeq data are in the [workflow_code](workflow_code) directory. To get a 
copy of latest NF_MSCP version on to your system, the code can be downloaded as a zip file from the release page then unzipped after downloading by running the following commands: 

```bash
wget https://github.com/nasa/GeneLab_MethylSeq_Workflow/releases/download/NF_MSCP_1.0.0/NF_MSCP_1.0.0.zip

unzip NF_MSCP_1.0.0.zip
```

<br>

---

### 3. Fetch Singularity Images

Although Nextflow can fetch Singularity images from a url, doing so may cause issues as detailed [here](https://github.com/nextflow-io/nextflow/issues/1210).

To avoid this issue, run the following command to fetch the Singularity images prior to running the NF_MSCP workflow:
> Note: This command should be run in the location containing the `NF_MSCP_2.1.0` directory that was downloaded in [step 2](#2-download-the-workflow-files) above. Depending on your network speed, fetching the images will take ~20 minutes. Approximately 8GB of RAM is needed to download and build the Singularity images.

```bash
bash NF_MSCP_1.0.0/bin/prepull_singularity.sh NF_MSCP_1.0.0/config/by_docker_image.config
```


Once complete, a `singularity` folder containing the Singularity images will be created. Run the following command to export this folder as a Nextflow configuration environment variable to ensure Nextflow can locate the fetched images:

```bash
export NXF_SINGULARITY_CACHEDIR=$(pwd)/singularity
```

<br>

---

### 4. Run the Workflow

While in the location containing the `NF_MSCP_1.0.0` directory that was downloaded in [step 2](#2-download-the-workflow-files), you are now able to run the workflow.

The workflow automatically loads reference files and organism-specific gene annotation files from the [GeneLab annotations table](https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv). For organisms not listed in the table or to use alternative reference files, additional workflow parameters can be specified.

The workflow also parses a dataset's metadata to determine how to process the data. Specifically, it detects whether the library type is DNA or RNA MethylSeq, which determines the aligner used in the alignment step and the assay suffix assigned to output files (`_GLMethylSeq` or `_GLRNAMethylSeq`). It also detects whether the dataset was generated using an RRBS (Reduced Representation Bisulfite Sequencing) protocol, which enables or disables RRBS-specific alignment flags accordingly.

In addition, the workflow uses the library preparation kit identified in the metadata to automatically apply a matching preset from [`library_presets.config`](config/library_presets.config). Each preset defines kit-specific parameters for trimming, aligning, and extracting methylation calls as well as whether the library is non-directional. 

If the library preparation kit is not recognized or no preset match is found, the workflow falls back to default parameters. These presets can also be overridden by specifying the relevant parameters directly at runtime.

 Below are four examples of how to run the NF_MSCP workflow:
> Note: Nextflow commands use both single hyphen arguments (e.g. -help) that denote general nextflow arguments and double hyphen arguments (e.g. --reference_version) that denote workflow specific parameters.  Make sure you are using the proper number of hyphens for each argument.

> Note: To use Docker instead of Singularity, use `-profile docker` in the Nextflow run command. Nextflow will automatically pull images as needed.

> Note: The `-resume` parameter can be used to resume a previously interrupted workflow from where it left off (see [Nextflow documentation](https://www.nextflow.io/docs/latest/getstarted.html#modify-and-resume)) or to restart the workflow from a specific point by changing relevant parameters, which will re-execute that process and all downstream affected processes.

<br>

#### 4a. Approach 1: Run the workflow on a GeneLab MethylSeq dataset with automatic retrieval of reference fasta and gtf files

```bash
nextflow run NF_MSCP_1.0.0/main.nf \ 
   -profile singularity,local \
   --accession OSD-47 
```

<br>

#### 4b. Approach 2: Run the workflow on a GeneLab MethylSeq dataset with custom reference fasta and gtf files

```bash
nextflow run NF_MSCP_1.0.0/main.nf \ 
   -profile singularity,local \
   --accession OSD-47 \
   --ref_version 112 \
   --ref_source ensembl \ 
   --ref_fasta <url/or/path/to/fasta> \ 
   --ref_gtf <url/or/path/to/gtf>
```

> Note: The `--ref_source` and `--ref_version` parameters should match the reference source and version number of the reference fasta and gtf files used. 

> Note: For gene annotations in the differential methylation output table, see the optional `--gene_annotations_file` parameter described in the [Optional Parameters](#optional-parameters) section.

<br>

#### 4c. Approach 3: Run the workflow on a non-GeneLab dataset using a user-created runsheet with automatic retrieval of reference fasta and gtf files

```bash
nextflow run NF_MSCP_1.0.0/main.nf \ 
   -profile singularity,local \
   --runsheet_path </path/to/runsheet> 
```

> Note: Specifications for creating a runsheet manually are described [here](examples/runsheet/README.md).

<br>

#### 4d. Approach 4: Run the workflow on a non-GeneLab dataset using a user-created runsheet with custom reference fasta and gtf files

```bash
nextflow run NF_MSCP_1.0.0/main.nf \ 
   -profile singularity \
   --runsheet_path </path/to/runsheet> \
   --ref_version 112 \
   --ref_source ensembl \ 
   --ref_fasta <url/or/path/to/fasta> \ 
   --ref_gtf <url/or/path/to/gtf> 
```
> Note: This approach should be used for organisms not listed in the [GeneLab annotations table](https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv). 

> Note: The `--ref_source` and `--ref_version` parameters should match the reference source and version number of the reference fasta and gtf files used. 

> Note: For gene annotations in the differential methylation output table, see the optional `--gene_annotations_file` parameter described in the [Optional Parameters](#optional-parameters) section.

<br>


#### Required Parameters For All Approaches:

* `NF_MSCP_1.0.0/main.nf` - Instructs Nextflow to run the NF_MSCP workflow 

* `-profile` - Specifies the configuration profile(s) to load, `singularity` instructs Nextflow to setup and use singularity for all software called in the workflow; use `local` for local execution ([local.config](workflow_code/conf/local.config)) or `slurm` for SLURM cluster execution ([slurm.config](workflow_code/conf/slurm.config))

<br>

**Additional Required Parameters For [Approach 1](#4a-approach-1-run-the-workflow-on-a-genelab-methylseq-dataset-with-automatic-retrieval-of-reference-fasta-and-gtf-files):**

* `--accession` - The OSD or GLDS ID for the dataset to be processed, eg. `GLDS-194` or `OSD-194`

<br>

**Additional Required Parameters For [Approach 2](#4b-approach-2-run-the-workflow-on-a-genelab-methylseq-dataset-with-custom-reference-fasta-and-gtf-files):**

* `--accession` - The OSD or GLDS ID for the dataset to be processed, eg. `GLDS-194` or `OSD-194`

* `--ref_version` - specifies the reference source version to use for the reference genome (Ensembl release `112` is used in this example); only needed when using Ensembl as the reference source

* `--ref_source` - specifies the source of the reference files used (the source indicated in the Approach 2 example is `ensembl`) 

* `--ref_fasta` - specifies the URL or path to a fasta file 

* `--ref_gtf` - specifies the URL or path to a gtf file

<br>

**Additional Required Parameters For [Approach 3](#4c-approach-3-run-the-workflow-on-a-non-genelab-dataset-using-a-user-created-runsheet-with-automatic-retrieval-of-reference-fasta-and-gtf-files):**

* `--runsheet_path` - specifies the path to a local runsheet; if not provided, a runsheet is automatically generated using OSDR metadata (type: string, default: null)

<br>

**Additional Required Parameters For [Approach 4](#4d-approach-4-run-the-workflow-on-a-non-genelab-dataset-using-a-user-created-runsheet-with-custom-reference-fasta-and-gtf-files):**

* `--runsheet_path` - specifies the path to a local runsheet; if not provided, a runsheet is automatically generated using OSDR metadata (type: string, default: null)

* `--ref_version` - specifies the reference source version to use for the reference genome (Ensembl release `112` is used in this example); only needed when using Ensembl as the reference source

* `--ref_source` - specifies the source of the reference files used (the source indicated in the Approach 2 example is `ensembl`) 

* `--ref_fasta` - specifies the URL or path to a fasta file 

* `--ref_gtf` - specifies the URL or path to a gtf file

<br>

#### Optional Parameters:

* `--assay_suffix` - Specifies a string that should be used to label the output file names. If not supplied, the workflow sets it to either `_GLMethylSeq` or `_GLRNAMethylSeq` based on runsheet's metadata (type: string, default: null).

* `--outdir` - Specifies the base directory where the output directory will be created (type: string, default: "results")

* `--reference_store_path` - specifies the directory to store the reference fasta and gtf files (type: string, default: "./References")  

* `--derived_store_path` - specifies the directory to store the tool-specific indices created during processing (type: string, default: "./DerivedReferences")

* `--gene_annotations_file` - Specifies the URL or path to a gene annotation file that adds additional gene annotation columns to the differential methylation output table. This can be:

  - The file listed in the `genelab_annots_link` column of the [GeneLab annotations table](https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)
  - A custom gene annotation file where:
    - For organisms listed in the GeneLab annotations table: gene IDs must be in a column with the same name as column 1 of the GeneLab organism-specific gene annotation file
    - For organisms not listed in the table: gene IDs must be in a column named `gene_id`
  
  Only genes included in the specified annotations file will receive additional annotations in the output.

* `--non_directional` - Enables appropriate processing of non-directional bisulfite libraries at trimming and alignemnt levels (type: boolean, default: false)

* `--skip_dedupe` - Skips deduplication regardless of the library type (type: boolean, default: false) 

* `--force_single_end` - forces the analysis to use single end processing; for paired end datasets, this means only R1 is used; for single end datasets, this should have no effect (type: boolean, default: false)  


<br>

**Additional Optional Parameters:**

All parameters listed above and additional optional arguments for the MSCP workflow, including tool-specific parameters and debug related options that may not be immediately useful for most users, can be viewed by running the following command:

```bash
nextflow run NF_MSCP_1.0.0/main.nf --help
```

See `nextflow run -h` and [Nextflow's CLI run command documentation](https://nextflow.io/docs/latest/cli.html#run) for more options and details common to all nextflow workflows.

<br>

---

### 5. Additional Output Files

> Note: The outputs from the MethylSeq Consensus Pipeline workflow are documented in the [GL-DPPD-7113.md](https://github.com/nasa/GeneLab_Data_Processing/tree/master/Methyl-Seq/Pipeline_GL-DPPD-7113_Versions/GL-DPPD-7113.md) processing protocol.

> Note: All output files are published inside `--outdir`, which is set to `results/` by default.

The additional outputs from the Analysis Staging and Processing are described below:

**Analysis Staging Subworkflow**

   - Metadata/\*_methylSeq_v2_runsheet.csv (table containing metadata required for processing, including the raw reads files location)
   - Metadata/\*-ISA.zip (the ISA archive of the OSD datasets to be processed, downloaded from the OSDR)
   
   
**Version Capturing and Processed-Data Protocol Generation**
   - GeneLab/software_versions_<assay_suffix>.md 
   - GeneLab/processed_data_protocol_<assay_suffix>.txt 


Standard Nextflow resource usage logs are also produced as follows:
> Further details about these logs can also found within [this Nextflow documentation page](https://www.nextflow.io/docs/latest/tracing.html#execution-report).

**Nextflow Resource Usage Logs**

   - nextflow_info/execution_report_{timestamp}.html (an html report that includes metrics about the workflow execution including computational resources and exact workflow process commands)
   - nextflow_info/execution_timeline_{timestamp}.html (an html timeline for all processes executed in the workflow)
   - nextflow_info/execution_trace_{timestamp}.txt (an execution tracing file that contains information about each process executed in the workflow, including: submission time, start time, completion time, cpu and memory used, machine-readable output)
   - nextflow_info/pipeline_dag_{timestamp}.html (a visualization of the workflow process DAG)

<br>

---

### 6. Post-processing


The post-processing workflow is designed to validate processed data generated by the main workflow, package processing information, purge file paths when necessary, and generate an updated assay table that can be hosted on OSDR, in addition to other output files, such as README, md5sums table, and raw data protocol if needed.

Before running the post-processing workflow, set the parameters nested under the `post_processing` block in [`nextflow.config`](workflow_code/nextflow.config) — including the analyst's name and email, protocol ID, OSD and GLDS accession numbers, and whether raw data files should be included or were merged. Once set, run the following command:

```bash
nextflow run post_processing.nf \
   -profile singularity
``` 

For more details on the parameters for the post-processing workflow, run the following command:

```bash
nextflow run main.nf --help post_processing
```

The outputs of the post-processing workflow are described below:
> *Note: All post_processing outputs are saved `--outdir`, which is set to `results/` by default.*

**Individual Output Files** 
 
 - GeneLab/updated_curation_tables/a*.txt (Updated assay table that complies with OSDR's assay table)
 - GeneLab/<GLDS_accession>-methylseq-validation.log (Automated verification and validation log file)
 - GeneLab/README_<assay_suffix>.txt (README file listing and describing the outputs of the workflow)
 - GeneLab/processed_md5sum_<assay_suffix>.tsv (md5sums for the processed data published on OSDR)

 **Purged Output Files**
 > *Note: These files are saved/overwritten directly under their parent directories*

 - *\_multiqc\_<assay_suffix>_data.zip (Zipped MultiQC data folders after purging file paths for raw, trimmed, and align_and_bismark MultiQC data)
 - \*\_bismark\*\_report.txt (Sample-specific Bismark alignment reports after purging file paths)
 
 **Processing Information Archive**

   - GeneLab/processing_info_<assay_suffix>.zip (Archive containing workflow execution metadata)
     - processing_info/samples.txt (single column list of all sample names in the dataset)
     - processing_info/nextflow_log_<assay_suffix>.txt (Nextflow execution logs captured via `nextflow log`)
     - processing_info/nextflow_run_command_<assay_suffix>.txt (Exact command line used to initiate the workflow) 

**Optional Output Files**

   - GeneLab/raw_md5sum_<assay_suffix>.tsv (md5sums for the raw data published on OSDR)
   - GeneLab/raw_data_protocol_<assay_suffix>.txt  (Protocol for raw reads generation and QC)

<br>

---

## Licenses

The software for the Methyl-Seq pipeline and workflow is released under the [NASA Open Source Agreement (NOSA) Version 1.3](License/Methylation_Sequencing_NOSA_License.pdf).


### 3rd Party Software Licenses

Licenses for the 3rd party open source software utilized in the Methyl-Seq pipeline and workflow can be found in the [3rd_Party_Licenses sub-directory](License/3rd_Party_Licenses/). 

<br>

---

## Notices

Copyright © 2024 United States Government as represented by the Administrator of the National Aeronautics and Space Administration.  All Rights Reserved. 

### Disclaimers

No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."

Waiver and Indemnity: RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT. 

The “GeneLab Methylation Sequencing Processing Pipeline and Workflow” software also makes use of 3rd party Open Source software, released under the licenses indicated above. A complete listing of 3rd Party software notices and licenses made use of in "GeneLab Methylation Sequencing Processing Pipeline and Workflow” can be found in the [3rd Party Licenses README.md](License/3rd_Party_Licenses/README.md) file. 

<br>

---
**Developed by:

Michael Lee  
Crystal Han  
Jihan Yehia

**Maintained by:  

Jihan Yehia


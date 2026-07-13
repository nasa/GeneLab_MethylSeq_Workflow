#!/usr/bin/env Rscript

library(optparse)
## other packages loaded below after help menu might be called

############################################
############ handling arguments ############
############################################

parser <- OptionParser(description = "\n  This is really only intented to be called from within the GeneLab MethylSeq workflow.")

parser <- add_option(parser, c("-v", "--verbose"),
    action = "store_true",
    default = FALSE, help = "Print extra output")

parser <- add_option(parser, c("--bismark_methylation_calls_dir"),
    default = "Methylation_Call_Data",
    help = "Directory holding *bismark.cov.gz files")

parser <- add_option(parser, c("--path_to_runsheet"),
    help = "Path to the runsheet")

parser <- add_option(parser, c("--org_name"),
    help = "Organism name (must match species column of *annotations.csv ref table with spaces converted to underscore, ex: Mus_musculus)")

parser <- add_option(parser, c("--ref_bed_path"),
    help = "Path to reference BED file")

parser <- add_option(parser, c("--gene_transcript_map_path"),
    help = "Path to gene-to-transcript.tsv")

parser <- add_option(parser, c("--methylkit_output_dir"),
    default = "MethylKit_Outputs",
    help = "Directory for methylkit output files")

parser <- add_option(parser, c("--limit_samples_to"),
    default = "all",
    help = "Limits the number of samples being processed (won't do real factor comparisons if set)")

parser <- add_option(parser, c("--ref_ensemblVersion"),
    help = "Reference Ensmbl version (must match ensmblVersion column of *annotations.csv ref table)")

parser <- add_option(parser, c("--ref_annotations_tab_link"),
    help = "Link to reference-genome annotations table")

parser <- add_option(parser, c("--methRead_mincov"),
    default = 10, type = "integer",
    help = "Passed to mincov argument of methRead() call")

parser <- add_option(parser, c("--mc_cores"),
    default = 4, type = "integer",
    help = "Passed to mc.cores argument of calculateDiffMeth() call")

parser <- add_option(parser, c("--getMethylDiff_difference"),
    default = 25, type = "integer",
    help = "Passed to difference argument of getMethylDiff() call")

parser <- add_option(parser, c("--getMethylDiff_qvalue"),
    default = 0.01, type = "double",
    help = "Passed to qvalue argument of getMethylDiff() call")

parser <- add_option(parser, c("--tileMethylCounts_mincov"),
    default = 10, type = "integer",
    help = "Passed to cov.bases argument of tileMethylCounts() call")

parser <- add_option(parser, c("--tileMethylCounts_winsize"),
    default = 1000, type = "integer",
    help = "Passed to win.size argument of tileMethylCounts() call")

parser <- add_option(parser, c("--tileMethylCounts_stepsize"),
    default = 1000, type = "integer",
    help = "Passed to step.size argument of tileMethylCounts() call")

parser <- add_option(parser, c("--primary_keytype"),
    help = "The keytype to use for mapping annotations (usually 'ENSEMBL' for most things; 'TAIR' for plants)")

parser <- add_option(parser, c("--test"),
    action = "store_true",
    default = FALSE,
    help = "Provide solely this flag to run with a small test dataset that will be downloaded")

parser <- add_option(parser, c("--file_suffix"),
	help = "File suffix to add to output files: _GLMethylSeq for DNA methylation, _GLRNAMethylSeq for RNA methylation")

parser <- add_option(parser, c("--normalize"),
    action = "store_true",
    default = FALSE,
    help = "Enable normalization. Prior to normalization using a median-based scaling factor, filters regions with coverage exceeding 99.9 percentile to reduce PCR bias effects.")

myargs <- parse_args(parser)


############################################
########### for testing purposes ###########
############################################


if (myargs$test) {
    test_meth_data_link <- "https://figshare.com/ndownloader/files/39246758"
    test_meth_data_tarball <- "MethylSeq-test-meth-call-cov-files.tar"
    test_ref_data_link <- "https://figshare.com/ndownloader/files/38616860"
    test_ref_data_tarball <- "MethylSeq-test-ref-files.tar"
    
    cat("\n  Running test data from:\n\n")
    cat("    - test methylation calls:   ", test_meth_data_link, "\n", sep = "")
    cat("    - test ref files:           ", test_ref_data_link, "\n", sep = "")
    
    cat("\n  NOTICE\n  Any other parameters are ignored when '--test' is provided.\n\n")

    suppressWarnings(suppressMessages(library(curl)))
    
    curl_download(url = test_meth_data_link, destfile = test_meth_data_tarball, quiet = TRUE)
    curl_download(url = test_ref_data_link, destfile = test_ref_data_tarball, quiet = TRUE)

    
    untar(test_meth_data_tarball, restore_times = FALSE)
    untar(test_ref_data_tarball, restore_times = FALSE)
    
    file.remove(test_meth_data_tarball, test_ref_data_tarball)

    myargs$v <- TRUE
    myargs$org_name <- "Mus_musculus"
    myargs$bismark_methylation_calls_dir <- "test-meth-calls"
    myargs$path_to_runsheet <- "test-meth-calls/test-runsheet.csv"
    myargs$ref_bed_path <- "test-ref-files/*.bed"
    myargs$gene_transcript_map_path <- "test-ref-files/*.-gene-to-transcript-map.tsv"
    myargs$methylkit_output_dir <- "test-MethylKit_Outputs"
    myargs$ref_ensemblVersion <- "112"
    myargs$ref_annotations_tab_link <- "https://api.figshare.com/v2/file/download/36597114"
    myargs$methRead_mincov <- 2
    myargs$getMethylDiff_difference <- 1
    myargs$getMethylDiff_qvalue <- 0.5
    myargs$primary_keytype <- "ENSEMBL"
    myargs$tileMethylCounts_mincov <- 10
    myargs$tileMethylCounts_stepsize <- 1000
    myargs$tileMethylCounts_winsize <- 1000

}

if (myargs$v) {
    cat("\n  NOTICE\n  Verbose logging has been specified.\n\n")
}

############################################


############################################
########### checking on arguments ##########
############################################

# checking required arguments were set
required_args <- c("path_to_runsheet" = "--path_to_runsheet",
                   "ref_ensemblVersion" = "--ref_ensemblVersion",
                   "ref_annotations_tab_link" = "--ref_annotations_tab_link",
                   "primary_keytype" = "--primary_keytype")

for (arg in names(required_args)) {
    tryCatch(
        {
            get(arg, myargs)
        },
        error = function(e) {
            error_message <- paste0(
                "\n  The '", as.character(required_args[arg]),
                "' argument must be provided.\nCannot proceed.\n\n"
            )
            stop(error_message, call. = FALSE)
        }
    )
}

# checking primary keytype is what's expected
currently_accepted_keytypes <- c("ENSEMBL", "TAIR")
if (!myargs$primary_keytype %in% currently_accepted_keytypes) {
    error_message <- paste0(
        "\n  The current potential --primary_keytypes are: ", paste(currently_accepted_keytypes, collapse = ", "),
        "\nCannot proceed with: ", myargs$primary_keytype, "\n\n"
    )

    stop(error_message, call. = FALSE)
}

# checking runsheet file exists
if (!file.exists(myargs$path_to_runsheet)) {
    stop("\n  The specified --path_to_runsheet does not seem to point to an actual file.\nCannot proceed.\n\n", call. = FALSE)
}

############################################
############# loading packages #############
############################################

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(methylKit)))
suppressWarnings(suppressMessages(library(genomation)))
suppressWarnings(suppressMessages(library(dplyr)))

# Set timeout for internet operations
options(timeout = 600)

############################################


############################################
############# helper functions #############
############################################

#' Order the coverage files in the same order as the corresponding sample names.
#'
#' @param sample_names A list of sample_names.
#' @param file_suffix The file suffix to look for in the coverage file names, to ensure correct matching of sample names to coverage files.
#' @param paths A list of paths to coverage files.
#' @returns A list of coverage files in the same order as the sample names.
order_input_files <- function(sample_names, file_suffix, paths) {
    ordered_paths <- c()

    for (sample in sample_names) {
        search_pattern <- paste0(sample, file_suffix, "_bismark")
        hits <- paths[grep(search_pattern, paths)]

        # making sure there is exactly one match
        if (length(hits) != 1) {
            stop("\n  There was a problem matching up sample names with their coverage files.\nCannot proceed.\n\n", call. = FALSE)
        }

        ordered_paths <- c(ordered_paths, hits)
    }

    return(ordered_paths)
}

#' Load study metadata.
#'
#' @param runsheet_path A full path to the study runsheet.csv.
#' @returns A data.frame containing the samplenames and corresponding factors from the runsheet.
compare_csv_from_runsheet <- function(runsheet_path) {
    df <- read.csv(runsheet_path)
    factors <- df %>%
        dplyr::select(starts_with("Factor.Value", ignore.case = TRUE)) %>%
        rename_with(~ paste0("factor_", seq_along(.)))
    result <- df %>%
        dplyr::select(sample_id = "Sample.Name") %>%
        bind_cols(factors)
    return(result)
}

############################################


############################################
############# setting things up ############
############################################

### Create output directory
dir.create(myargs$methylkit_output_dir, showWarnings = FALSE)

### Pull in the GeneLab annotation table ###

gene_annotation_link <- myargs$ref_annotations_tab_link

# Download with user-agent
local_file <- "gene_annotations.tsv"
download.file(url = gene_annotation_link,
              destfile = local_file,
              method = "libcurl",
              headers = c("User-Agent" = "Mozilla/5.0"))

functional_annots_tab <- read.table(local_file, sep = "\t", quote = "", header = TRUE)


### Pull in reference annotation data ###

# Read in gene to transcript mapping file
gene_transcript_map <- read.table(myargs$gene_transcript_map_path, sep = "\t",
                                  col.names = c("gene_ID", "feature.name"))

# Read in reference bed file
gene.obj <- readTranscriptFeatures(myargs$ref_bed_path, up.flank = 1000,
                                   down.flank = 1000, remove.unusual = TRUE,
                                   unique.prom = TRUE)


##### Set variable with organism and ensembl version number for methylKit #####

org_and_ensembl_version <- paste(myargs$org_name, myargs$ref_ensemblVersion, sep = "_")


### Load metadata from runsheet csv file and create data frame containing all samples and respective factors ###
study <- compare_csv_from_runsheet(myargs$path_to_runsheet) %>% column_to_rownames(var = "sample_id")

if (myargs$v) {
    cat("\nStudy df:\n")
    print(study)
    cat("\n")
}

##### Format groups and indicate the group that each sample belongs to #####
group <- apply(study, 1, paste, collapse = " & ")
group_names <- paste0("(", group, ")") ## human readable group names
# generate group naming compatible with R models, maintains the default behaviour of make.names but ensures 'X' is never prepended
group <- sub("^BLOCKER_", "", make.names(paste0("BLOCKER_", group)))
names(group) <- group_names
rm(group_names)

# create a lookup table to map back from R naming to human readable group names
group_name_lookup <- data.frame(group = group, group_names = names(group)) %>%
    distinct() %>%
    column_to_rownames(var = "group")

##### Format contrasts table, defining pairwise comparisons for all groups, one-way contrasts #####

# generate matrix of pairwise group combinations for comparison
contrast.names <- combn(levels(factor(names(group))), 2)

# create computationally friendly contrast values
contrasts <- apply(
    contrast.names,
    MARGIN = 2,
    function(col) {
        # limited make.names call for each group (also removes leading parentheses)
        sub("^BLOCKER_", "", make.names(paste0("BLOCKER_", stringr::str_sub(col, 2, -2))))
    }
)

# format contrast combinations for output tables
colnames(contrasts) <- paste(contrast.names[1, ], contrast.names[2, ], sep = "v")

if (myargs$v) {
    cat("\nContrasts table:\n")
    print(contrasts)
    cat("\n")
}

## Output contrasts table
write.csv(
    contrasts,
    file = file.path(myargs$methylkit_output_dir,
                     paste0("contrasts", myargs$file_suffix, ".csv"))
)

# Create sampleTable
sampleTable <- data.frame("sample_id" = rownames(study), "condition" = group)

if (myargs$v) {
    cat("\nSamples and factors df:\n")
    print(sampleTable)
    cat("\n")
}

## Output Sample table
write.csv(
    sampleTable,
    file = file.path(myargs$methylkit_output_dir,
                     paste0("SampleTable", myargs$file_suffix, ".csv")),
    row.names = FALSE
)

# Get list of *.bismark.cov.gz files
bismark_cov_paths <- list.files(
    myargs$bismark_methylation_calls_dir,
    pattern = ".*.bismark.cov.gz",
    full.names = TRUE,
    recursive = TRUE
)

# Make sure files list matches number of samples in runsheet
if (dim(sampleTable)[1] != length(bismark_cov_paths)) {
    error_message <- paste0(
        "\n  The number of '*.bismark.cov.gz' files found in the ",
        myargs$bismark_methylation_calls_dir,
        " directory\n  does not match the number of samples specified in the runsheet.\nCannot proceed.\n\n"
    )

    stop(error_message, call. = FALSE)
}

if (myargs$v) {
    cat("\nBismark coverage files:\n")
    print(bismark_cov_paths)
    cat("\n")
}

# Make sure coverage-file-paths vector is in the same order as the sample names
bismark_cov_paths <- order_input_files(sampleTable$sample_id, myargs$file_suffix, bismark_cov_paths)

# Make a single table with info needed for methylkit
sample_meth_info_df <- full_join(
    data.frame("sample_id" = sampleTable$sample_id, "coverage_file_path" = bismark_cov_paths),
    sampleTable,
    "sample_id")

if (myargs$v) {
    cat("\nPrimary methylkit df:\n")
    print(sample_meth_info_df)
    cat("\n")
}


# ??? subsetting runsheet if specified, not mentioned in Barbara's script ???
if ( myargs$limit_samples_to != "all" ) {
    
    myargs$limit_samples_to <- as.integer(myargs$limit_samples_to)

    runsheet <- runsheet[ 1:myargs$limit_samples_to, ]
    
    # making up factors for testing if things were subset (so that we have multiple factor types even if not really in the subset samples)
    len_type_A <- ceiling(myargs$limit_samples_to / 2)
    len_type_B <- myargs$limit_samples_to - len_type_A
    
    factor_1 <- c(rep("A", len_type_A), rep("B", len_type_B))
    
    factor_df <- data.frame(sample_id = runsheet %>% pull("Sample.Name"), factor_1)

}


############################################


############################################
########## moving onto methylkit ##########
############################################

# Read in methylation counts
# NOTE: treatment vector is composed of integers corresponding to the groups.
#       It is used during the min.per.group filtering in the unite command
meth_obj <- methRead(
    location = as.list(sample_meth_info_df %>% pull(coverage_file_path)),
    sample.id = as.list(sample_meth_info_df %>% pull(sample_id)),
    treatment = as.list(sample_meth_info_df %>% mutate(group_no = as.integer(factor(condition))) %>% pull(group_no)),
    pipeline = "bismarkCoverage",
    assembly = org_and_ensembl_version,
    header = FALSE,
    mincov = myargs$methRead_mincov,
    dbtype = "tabix"
)

### Normalize the data if option was selected ###
if (myargs$normalize) {
    if (myargs$v) {
        cat("\nRunning normalization\n\n")
    }
    # First, filter samples by coverage to account for PCR bias or over-amplification
    meth_obj <- filterByCoverage(meth_obj, lo.count = myargs$methRead_mincov, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
    # Normalize coverage between samples using a scaling factor derived from the median coverage distributions
    meth_obj <- normalizeCoverage(meth_obj, method = "median")
}

### Base analysis ###

### Function for computing differential methylation per base for each contrasts ###
#' Compute differential methylation per base for each constrast.
#' @param i The index of the contrast in the constrasts vector
#' @return A data.frame containing the results of methylKit::calculateDiffMeth with human-readable
#'         column names ready for combining with data from other contrasts
compute_contrast_bases <- function(i) {
    # Get current contrast
    curr_contrasts_vec <- contrasts[, i]

    # Get subset sample info table relevant to current contrast
    curr_sample_info_df <- sample_meth_info_df %>% filter(condition %in% curr_contrasts_vec)

    # Get which samples are relevant to current contrast
    curr_samples_vec <- curr_sample_info_df %>% pull(sample_id)

    # Make binary vector for treatment argument to methRead()
    curr_treatment_vec <- c()
    for (value in curr_sample_info_df$condition) {
        if (value == curr_sample_info_df$condition[1]) {
            curr_treatment_vec <- c(curr_treatment_vec, 1)
        } else {
            curr_treatment_vec <- c(curr_treatment_vec, 0)
        }
    }

    # return methylation statistics for bases
    return(
        as.data.frame(getData(
            calculateDiffMeth(
                methylKit::unite(reorganize(
                    meth_obj,
                    sample.ids = as.list(curr_samples_vec),
                    treatment = curr_treatment_vec
                ),
                min.per.group = 1L), # Keep only bases with coverage in at least one sample per group
                mc.cores = myargs$mc_cores)
            )
        ) %>%
        relocate("meth.diff", .before = "pvalue") %>% # move methyl diff value before stats values
        rename_with(~ paste0(.x, "_", colnames(contrasts)[i]), any_of(c("pvalue", "qvalue", "meth.diff")))
    )
}


# create initial output object at base resolution
#   - compute pair-wise differential methylation for each contrast
#   - add percent methylation values for all samples
#   - add All.mean and All.stdev stats
# NOTE: when retrieving percent methylation, keep all bases (min.per.group = 0L) to
#       include bases that have coverage in some groups (and data for only some
#       contrasts) but no coverage in others
output_bases_df <- lapply(seq_len(dim(contrasts)[2]), compute_contrast_bases) %>%
    purrr::reduce(full_join, by = c("chr", "start", "end", "strand")) %>%
    left_join(as.data.frame(percMethylation(methylKit::unite(meth_obj, min.per.group = 0L), rowids = TRUE)) %>%
        rownames_to_column() %>%
        separate_wider_regex(rowname, c(chr = ".*", "\\.", start = ".*", "\\.", end = ".*")) %>%
        mutate_at(vars(start, end), as.integer) %>%
        rowwise() %>%
        mutate(
            All.mean = mean(c_across(all_of(sample_meth_info_df$sample_id)), na.rm = TRUE),
            All.stdev = sd(c_across(all_of(sample_meth_info_df$sample_id)), na.rm = TRUE)
        ), by = c("chr", "start", "end"))

# Add group means and stdev
for (g in rownames(group_name_lookup)) {
    colnames_to_process <- sample_meth_info_df %>% filter(condition == g) %>% pull(sample_id)
    output_bases_df <- output_bases_df %>%
        rowwise() %>%
        mutate(
            gmean = mean(c_across(all_of(colnames_to_process)), na.rm = TRUE),
            gstdev = sd(c_across(all_of(colnames_to_process)), na.rm = TRUE)
        ) %>%
        rename("gmean" = paste0("Group.Mean_", group_name_lookup[g, ]),
               "gstdev" = paste0("Group.Stdev_", group_name_lookup[g, ]))
}

# Annotate with gene parts (promoter/exon/intron)
suppressWarnings(output_ann_df <- annotateWithGeneParts(as(output_bases_df, "GRanges"), gene.obj))

# Associate the nearest TSS annotations with gene parts
# Add annotations to output bases dataframe
df_list <- list(
        tibble::as_tibble(getAssociationWithTSS(output_ann_df) %>% dplyr::rename(rowid = target.row)),
        tibble::rowid_to_column(as.data.frame(genomation::getMembers(output_ann_df))),
        tibble::rowid_to_column(output_bases_df))

# Create one combined table with all annotations, methylation values, and statistics
# reorder the columns to final expected format
bases_tab_with_features_and_annots <- df_list %>%
        purrr::reduce(full_join, by = "rowid") %>%
        mutate(rowid = NULL) %>%
    left_join(gene_transcript_map) %>%
    left_join(functional_annots_tab, by = c("gene_ID" = myargs$primary_keytype)) %>%
    rename("gene_ID" = myargs$primary_keytype) %>%
    relocate(c("ENSEMBL", "SYMBOL", "GENENAME", "REFSEQ", "ENTREZID", "STRING_id", "GOSLIM_IDS", "feature.name", "chr", "start", "end", "strand")) %>%
    relocate(sample_meth_info_df$sample_id, .after = "intron")
rm(df_list)

# write out differentially methylated bases
write.csv(
    bases_tab_with_features_and_annots,
    row.names = FALSE,
    file = file.path(myargs$methylkit_output_dir,
        paste0("differential_methylation_bases", myargs$file_suffix, ".csv")
    )
)
rm(bases_tab_with_features_and_annots)

### Tile analysis ###

### Function for computing differential methylation per tile for each contrasts ###
#' Compute differential methylation per tile for each constrast.
#' @param i The index of the contrast in the constrasts vector
#' @return A data.frame containing the results of methylKit::calculateDiffMeth with human-readable
#'         column names ready for combining with data from other contrasts
compute_contrast_tiles <- function(i) {
    # Get current contrast
    curr_contrasts_vec <- contrasts[, i]

    # Get subset sample info table relevant to current contrast
    curr_sample_info_df <- sample_meth_info_df %>% filter(condition %in% curr_contrasts_vec)

    # Get which samples are relevant to current contrast
    curr_samples_vec <- curr_sample_info_df %>% pull(sample_id)

    # Make binary vector for treatment argument to methRead()
    curr_treatment_vec <- c()
    for (value in curr_sample_info_df$condition) {
        if (value == curr_sample_info_df$condition[1]) {
            curr_treatment_vec <- c(curr_treatment_vec, 1)
        } else {
            curr_treatment_vec <- c(curr_treatment_vec, 0)
        }
    }

    # return methylation statistics for tiles
    return(
        as.data.frame(getData(
            calculateDiffMeth(
                methylKit::unite(reorganize(
                    meth_tiles_obj,
                    sample.ids = as.list(curr_samples_vec),
                    treatment = curr_treatment_vec
                ),
                min.per.group = 1L), # Keep only tiles with coverage in at least one sample per group
                mc.cores = myargs$mc_cores)
            )
        ) %>%
        relocate("meth.diff", .before = "pvalue") %>% # move methyl diff value before stats values
        rename_with(~ paste0(.x, "_", colnames(contrasts)[i]), any_of(c("pvalue", "qvalue", "meth.diff")))
    )
}

# Group bases into tiles
meth_tiles_obj <- tileMethylCounts(meth_obj,
                                   win.size = myargs$tileMethylCounts_winsize,
                                   step.size = myargs$tileMethylCounts_stepsize,
                                   cov.bases = myargs$tileMethylCounts_mincov)


# create initial output object at tile resolution
#   - compute pair-wise differential methylation for each contrast
#   - add percent methylation values for all samples
#   - add All.mean and All.stdev stats
# NOTE: when retrieving percent methylation, keep all tiles (min.per.group = 0L) to
#       include tiles that have coverage in some groups (and data for only some
#       contrasts) but no coverage in others
output_tiles_df <- lapply(seq_len(dim(contrasts)[2]), compute_contrast_tiles) %>%
    purrr::reduce(full_join, by = c("chr", "start", "end", "strand")) %>%
    left_join(as.data.frame(percMethylation(methylKit::unite(meth_tiles_obj, min.per.group = 0L), rowids = TRUE)) %>%
        rownames_to_column() %>%
        separate_wider_regex(rowname, c(chr = ".*", "\\.", start = ".*", "\\.", end = ".*")) %>%
        mutate_at(vars(start, end), as.integer) %>%
        rowwise() %>%
        mutate(
            All.mean = mean(c_across(all_of(sample_meth_info_df$sample_id)), na.rm = TRUE),
            All.stdev = sd(c_across(all_of(sample_meth_info_df$sample_id)), na.rm = TRUE)
        ), by = c("chr", "start", "end"))

# Add group means and stdev
for (g in rownames(group_name_lookup)) {
    colnames_to_process <- sample_meth_info_df %>% filter(condition == g) %>% pull(sample_id)
    output_tiles_df <- output_tiles_df %>%
        rowwise() %>%
        mutate(
            gmean = mean(c_across(all_of(colnames_to_process)), na.rm = TRUE),
            gstdev = sd(c_across(all_of(colnames_to_process)), na.rm = TRUE)
        ) %>%
        rename("gmean" = paste0("Group.Mean_", group_name_lookup[g, ]),
               "gstdev" = paste0("Group.Stdev_", group_name_lookup[g, ]))
}

# Annotate with gene parts (promoter/exon/intron)
suppressWarnings(output_tiles_ann_df <- annotateWithGeneParts(as(output_tiles_df, "GRanges"), gene.obj))

# Associate the nearest TSS annotations with gene parts
# Add annotations to output bases dataframe
df_list <- list(
        tibble::as_tibble(getAssociationWithTSS(output_tiles_ann_df)) %>% dplyr::rename(rowid = target.row),
        tibble::rowid_to_column(as.data.frame(genomation::getMembers(output_tiles_ann_df))),
        tibble::rowid_to_column(output_tiles_df))

# Create one combined table with all annotations, methylation values, and statistics
# reorder the columns to final expected format
tiles_tab_with_features_and_annots <- df_list %>%
        purrr::reduce(full_join, by = "rowid") %>%
        mutate(rowid = NULL) %>%
    left_join(gene_transcript_map) %>%
    left_join(functional_annots_tab, by = c("gene_ID" = myargs$primary_keytype)) %>%
    rename("gene_ID" = myargs$primary_keytype) %>%
    relocate(c("ENSEMBL", "SYMBOL", "GENENAME", "REFSEQ", "ENTREZID", "STRING_id", "GOSLIM_IDS", "feature.name", "chr", "start", "end", "strand")) %>%
    relocate(sample_meth_info_df$sample_id, .after = "intron")
rm(df_list)

# Write out differentially methylated tiles
write.csv(
    tiles_tab_with_features_and_annots,
    row.names = FALSE,
    file = file.path(myargs$methylkit_output_dir,
        paste0("differential_methylation_tiles", myargs$file_suffix, ".csv")
    )
)
rm(tiles_tab_with_features_and_annots)
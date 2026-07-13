#!/usr/bin/env bash
set -e

# Usage: clean_paths.sh <file>
# Removes cluster-specific paths from a file
# Built for use on N288 cluster only

file="$1"

sed -i \
    -e 's|/global/data/Data_Processing/MethylSeq_Datasets/GLDS_Datasets/||g' \
    -e 's|/global/smf/miniconda38_admin/envs/[^/]*/||g' \
    -e 's|/[^ ]*/main.nf|main.nf|g' \
    -e 's|/[^ ]*/References/||g' \
    -e 's|/[^ ]*/OSD-|OSD-|g' \
    -e 's|"/global/[^ ]*"|"<path-removed-for-security-purposes>"|g' \
    -e 's|/global/[^ ]*|<path-removed-for-security-purposes>|g' \
    "$file"
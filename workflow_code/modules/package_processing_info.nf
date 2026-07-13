process PACKAGE_PROCESSING_INFO { 
    tag "Purging paths of log files and creating zipped processing info"
    label "zip"

    input:
        path(processing_scripts)
        val(assay_suffix)

    output:
        path("processing_info${assay_suffix}.zip"), emit: zip

    script:
    """
    for f in ${processing_scripts}/nextflow*.txt; do
        echo "Purging file paths from \$f"
        clean_paths.sh "\$f"

        # Add assay suffix to the end of nextflow log files before .txt extension, if not already present
        if [[ "\$f" != *"${assay_suffix}.txt" ]]; then
            mv "\$f" "\${f%.txt}${assay_suffix}.txt"
        fi
    done
        
    # Zip 
    zip -r processing_info${assay_suffix}.zip processing_scripts
    """
}
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process PACKAGE_PROCESSING_INFO { 
    //used for run_command and processing_commands,
    // I have it under processing_scripts/ (assigned in launch_MethylSeq.slurm)
    //  check if it can be applied here too, if it does copy clean-paths.sh to bin/

    tag "Purging file paths and zipping processing info"

    input:
        val(files_and_dirs) 
    output:
        path("${params.cleaned_prefix}processing_info${params.assay_suffix}.zip"), emit: zip

    script:
        """
        cat `which clean-paths.sh` > clean-paths.sh
        chmod +x ./clean-paths.sh
        mkdir processing_info/ && \\
        cp -r ${files_and_dirs.join(" ")} processing_info/

        echo "Purging file paths"
        find processing_info/ -type f -exec bash ./clean-paths.sh '{}' ${params.baseDir} \\;
        
        # Purge file paths and then zip
        zip -r ${params.cleaned_prefix}processing_info${params.assay_suffix}.zip processing_info/
        """
}
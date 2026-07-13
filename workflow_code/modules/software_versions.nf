process SOFTWARE_VERSIONS {

    input:
        path(versions_file)
        val(suffix)
    
    output:
        path "software_versions*.md", emit: software_versions_md
        path "software_versions*.yaml", emit: software_versions_yaml

    script:
    def assay_suffix = params.assay_suffix ? params.assay_suffix : "${suffix}"
    """
    software_versions.py ${versions_file} software_versions${assay_suffix}.md --workflow NF_MethylSeq --workflow_version ${workflow.manifest.version} --assay methylseq
    """
}
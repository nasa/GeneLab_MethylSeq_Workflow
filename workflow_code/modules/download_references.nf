process DOWNLOAD_REFERENCES {
  // Download and decompress genome and annotation files
  tag "Organism: ${organism_sci}, Reference Source: ${ref_source}${ref_source.toLowerCase().contains('ensembl') ? ', Reference Version: ' + ref_version : ''}"
  label 'networkBound'
  storeDir "${params.reference_store_path}/${ref_source}/${ref_source.toLowerCase().contains('ensembl') ? ref_version + '/' : ''}${organism_sci}"

  input:
    val(organism_sci)
    val(fasta_url)
    val(gtf_url)
    val(ref_source)
    val(ref_version)
  
  output:
    tuple path("{*.fa,*.fna}"), path("*.gtf"), emit: reference_files

  script:
  """
  # Create temp directories for processing
  mkdir -p temp_fasta temp_gtf

  # Handle fasta file
  if [[ "${fasta_url}" == http* ]]; then
    wget --directory-prefix temp_fasta "${fasta_url}"
  else
    # It's a file path - just copy it directly to main directory
    cp "${fasta_url}" temp_fasta
  fi
  
  if [[ "${gtf_url}" == http* ]]; then
    wget --directory-prefix temp_gtf "${gtf_url}"
  else
    # It's a file path - just copy it directly to main directory
    cp "${gtf_url}" temp_gtf
  fi
    
  # Handle decompression if needed
  if ls temp_fasta/*.gz &> /dev/null; then
    gunzip temp_fasta/*.gz
  fi
  # Move processed files to main directory
  mv temp_fasta/*.fa temp_fasta/*.fna ./ 2>/dev/null || true

  if ls temp_gtf/*.gz &> /dev/null; then
    gunzip temp_gtf/*.gz
  fi
  # Move processed files to main directory
  mv temp_gtf/*.gtf ./ 2>/dev/null || true
  
  # Clean up temp directories
  rm -rf temp_fasta temp_gtf
  """
}
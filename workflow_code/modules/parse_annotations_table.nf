process PARSE_ANNOTATIONS_TABLE {
  // Extracts data from this kind of table: 
  // https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv

  input:
    val(annotations_csv_url_string)
    val(organism_sci)
  
  output:
    val(fasta_url), emit: reference_fasta_url
    val(gtf_url), emit: reference_gtf_url
    val(gene_annotations_url), emit: gene_annotations_url
    val(ref_source), emit: reference_source
    val(ref_version), emit: reference_version
    val(simple_organism_name), emit: simple_organism_name
  
  exec:
    def colorCodes = [
        c_line: "┅" * 70,
        c_back_bright_red: "\u001b[41;1m",
        c_bright_green: "\u001b[32;1m",
        c_blue: "\033[0;34m",
        c_yellow: "\u001b[33;1m",
        c_reset: "\033[0m"
    ]

    def organisms = [:]
    println "${colorCodes.c_yellow}Fetching table from ${annotations_csv_url_string}${colorCodes.c_reset}"
    
    // Check if input is a URL or a local file path
    if (annotations_csv_url_string.startsWith('http://') || annotations_csv_url_string.startsWith('https://')) {
      // For URLs: use toURL() method
      annotations_csv_url_string.toURL().splitEachLine(",") {fields ->
            organisms[fields[1]] = fields
      }
    } else {
      // For local files: use File class
      new File(annotations_csv_url_string).splitEachLine(",") {fields ->
            organisms[fields[1]] = fields
      }
    }
    
    // Extract required fields
    organism_key = organism_sci.capitalize().replace("_"," ")
    
    // Check if the organism exists in the table
    if (organisms.containsKey(organism_key)) {
      simple_organism_name = organisms[organism_key][0]
      fasta_url = organisms[organism_key][5]
      gtf_url = organisms[organism_key][6]
      gene_annotations_url = organisms[organism_key][10]
      
      // Convert figshare ndownloader URL to API endpoint
      if (gene_annotations_url != null && gene_annotations_url.contains('figshare.com/ndownloader/files/')) {
        file_id = (gene_annotations_url =~ /.*\/files\/([a-zA-Z0-9]+).*/)[0][1]
        gene_annotations_url = "https://api.figshare.com/v2/file/download/${file_id}"
      }

      ref_version = organisms[organism_key][3]
      ref_source = organisms[organism_key][4]

      println "${colorCodes.c_blue}Annotation table values parsed for '${organism_key}':${colorCodes.c_bright_green}"
      println "--------------------------------------------------"
      println "- fasta_url: ${fasta_url}"
      println "- gtf_url: ${gtf_url}"
      println "- gene_annotations_url: ${gene_annotations_url}"
      println "- ref_source: ${ref_source}${colorCodes.c_reset}"
      if (ref_source.toLowerCase().contains('ensembl')) {
            println "${colorCodes.c_bright_green}- ref_version: ${ref_version}${colorCodes.c_reset}"
        }
      println "--------------------------------------------------"
    }
    else {
      fasta_url = null
      gtf_url = null
      gene_annotations_url = null
      ref_source = null
      ref_version = null
      simple_organism_name = null
      println "${colorCodes.c_back_bright_red}WARNING: Organism '${organism_key}' not found in annotations table.${colorCodes.c_reset}"
      println "${colorCodes.c_yellow}Returning null values for all outputs.${colorCodes.c_reset}"
        
    }
}

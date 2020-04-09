
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' post_process_salmon
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title post_process_salmon 
#' 
#' @description 
#' Joins individual sample files into one tsv file.  Also converts ucsc names to hgnc_symbol|entrez_ids
#' and outputs the gene level counts.
#' 
#' @param input_file_paths Named character vector of paths to the pipeline output data. Samples will be named by the names of the vector. If not named thye will be named after the quant file.
#' @param output_dir Path to the output folder.
#' @param gene_biotypes The type of biomaRt gene_biotypes that should be output to the gene level output.
#' @param output_entrez_id_matrix T/F
#' @param output_hgnc_matrix T/F
#' @param output_log2_upper_quartile_norm T/F
#' @param output_piped_hugo_entrez_id_matrix T/F
#' @param output_transcript_matrix T/F
#' @param output_upper_quartile_norm T/F
#' @param counts_or_tpm 'counts' or 'tpm'
#' @param ref options 'grch38' with ensemble output or 'hg38' with ucsc output
#' @param thread_num Integer number of threads to run mclapply statements
#' @param this_script_path Path to script that runs this function for documentation urposes
#' @return A list of paths to the files.
#' 
#' @family mart
#' 
#' @export
post_process_salmon = function(
  this_script_path = '',
  input_file_paths,# = system(paste0("ls ", RAW_DATA_DIR, "/pipeline_output/star_salmon/*/*_quant.sf"), intern = TRUE)
  output_dir,# = file.path(base_dir, "post_processing", "star_salmon")
  ref = "grch38", # options 'grch38' with ensemble output or 'hg38' with ucsc output
  gene_biotypes = c('protein_coding', 
                    'IG_C_gene','IG_D_gene', 'IG_J_gene', 'IG_V_gene',
                    'TR_C_gene', 'TR_D_gene', 'TR_J_gene','TR_V_gene'),
  thread_num = 8,
  output_transcript_matrix = TRUE,
  output_hgnc_matrix = TRUE,
  output_entrez_id_matrix = FALSE,
  output_piped_hugo_entrez_id_matrix = FALSE,
  output_upper_quartile_norm = TRUE,
  output_log2_upper_quartile_norm = TRUE,
  counts_or_tpm = "counts"
){
  
  library(magrittr)
  
  dir.create(output_dir, showWarnings = F)
  
  readme_path = file.path(output_dir, "readme.txt")
  if(file.exists(readme_path)){ file.remove(readme_path)}
  
  a = function(...){
    my_output = paste0(...)
    if(!is.null(readme_path)){
      write(my_output, readme_path, append = TRUE)
    }
    cat(paste0(my_output,"\n"))
  }
  
  a(paste0("## Making isoform counts matrix: ", this_script_path))
  a("")
  
  if(length(names(input_file_paths)) == 0){
    names(input_file_paths) = gsub(".quant.sf$", "", basename(input_file_paths))
  }
  
  
  a("Reading in files from input_file_paths:")
  counts_or_tpm = tolower(counts_or_tpm)
  if(counts_or_tpm == "counts"){
    salmon_column = "NumReads"
  } else if (counts_or_tpm == "tpm"){
    salmon_column = "TPM"
  }
  
  read_data = lapply(input_file_paths, function(input_file_path){
    a("  ", input_file_path)
    counts_df = fread(input_file_path, select = c("Name", salmon_column), data.table = F)
    return_list = lapply(counts_df[[salmon_column]], function(x)x) # rbindlist needs a list so we turn this into a list
    names(return_list) = counts_df[["Name"]] # Name the list items so they get assigned to the right column
    return(return_list)
  })
  a("")
  
  # convert to a data table according to the names of the list 
  #   (ie, use.names = T won't assume the genes are in the same order)
  dat = rbindlist(read_data, use.names = TRUE, fill = TRUE)


  # add the sample names
  dat = data.frame(Sample_ID = names(read_data), dat)
  row.names(dat) = NULL
  
  if(output_transcript_matrix){
    
    if (ref == 'grch38'){
      isoform_output_path = file.path(output_dir, paste0("grch38_transcript_", counts_or_tpm,".tsv"))
    } else if (ref == 'hg38') {
      isoform_output_path = file.path(output_dir, paste0("hg38_transcript_", counts_or_tpm,".tsv"))
    }
    
    fwrite(dat, isoform_output_path, sep = "\t")
    
    if(output_upper_quartile_norm || output_log2_upper_quartile_norm){
      norm_dat = binfotron::normalize_rows_by_quartile(data.table(dat))
      fwrite(norm_dat,  file.path(output_dir, paste0("norm_transcript_", counts_or_tpm,".tsv")), sep = "\t")
      if(output_log2_upper_quartile_norm){
        log2_norm_dat = binfotron::log_transform_plus(norm_dat)
        fwrite(log2_norm_dat,  file.path(output_dir, paste0("log2_norm_transcript_", counts_or_tpm,".tsv")), sep = "\t")
      }
    }
  }
  
  
  a(paste0("## Making gene counts matrix: ", this_script_path))
  a("")
  a("Using saved biomaRt.")
  if (ref == 'grch38'){
    BM_results = readRDS(StarSalmon::get_biomart_grch38_path())
    #BM_results = readRDS("/rstudio-common/dbortone/packages/StarSalmon/inst/biomart/grch38/bm_result.rds")
  } else if (ref == 'hg38') {
    BM_results = readRDS(StarSalmon::get_biomart_hg38_path())
    #BM_results = readRDS("/rstudio-common/dbortone/packages/StarSalmon/inst/biomart/hg38/bm_result.rds")
    
  } 
  
  a("")
  
  BM_results$gene_biotype = factor(BM_results$gene_biotype)
  
  # for debugging this it's best to look at the BM_results matrix and make sure
  #   each step is filling in data the right way...
  a(paste0("* Starting with gene_biotypes: ", paste0(unique(BM_results$gene_biotype), collapse = ", ")))
  
  a(paste0("* Restricting gene_biotypes to: ", paste0(gene_biotypes, collapse = ", ")))
  BM_results = BM_results[ BM_results$gene_biotype %in% gene_biotypes, ]
  
  missing_ids = which(is.na(BM_results$entrezgene))
  
  #lots of blanks in symbols
  missing_symbols = which(BM_results$hgnc_symbol == "")
  
  fill_in_ids_using_symbol_indices = missing_ids[missing_ids %ni% missing_symbols]
  fill_in_symbols_using_id_indices = missing_symbols[missing_symbols %ni% missing_ids]
  
  symbols_to_use = BM_results$hgnc_symbol[fill_in_ids_using_symbol_indices]
  
  a("* Used hgnc symbols to fill in missing entrez ids.")

  # mclapply takes 3x longer
  fill_in_ids = unlist(lapply(1:length(symbols_to_use), function(x){ # very slow but mclapply runs into errors
    gene_id = as.character(paste0(unlist(annotate::lookUp(symbols_to_use[x],'org.Hs.eg','SYMBOL2EG')), collapse = ","))
    return(gene_id)
  }))
  
  BM_results$symbols_to_use = NA
  BM_results$symbols_to_use[fill_in_ids_using_symbol_indices] = symbols_to_use # just marking these to make sure we are converting the correct things
  BM_results$fill_in_ids = NA
  BM_results$fill_in_ids[fill_in_ids_using_symbol_indices] = fill_in_ids
  
  
  ids_to_use = as.character(BM_results$entrezgene[fill_in_symbols_using_id_indices])
  
  a("* Used entrez ids to fill in missing hgnc symbols.")

  # mclapply takes 3x longer
  fill_in_symbols = unlist(lapply(1:length(ids_to_use), function(x){ # very slow but mclapply runs into errors
    gene_name = as.character(paste0(unlist(annotate::lookUp(ids_to_use[x],'org.Hs.eg','SYMBOL')), collapse = ","))
    return(gene_name)
  }))
  
  BM_results$fill_in_symbols = NA
  BM_results$fill_in_symbols[fill_in_symbols_using_id_indices] = fill_in_symbols
  BM_results$ids_to_use = NA
  BM_results$ids_to_use[fill_in_symbols_using_id_indices] = ids_to_use
  
  BM_results$fin_symbols = BM_results$hgnc_symbol
  BM_results$fin_symbols[fill_in_symbols_using_id_indices] = fill_in_symbols
  
  BM_results$fin_ids = BM_results$entrezgene
  BM_results$fin_ids[fill_in_ids_using_symbol_indices] = fill_in_ids
  
  # do any of the hgnc or entrez id's convert to multiple?  I've never seen this but it's possible
  if(sum(grepl(",", BM_results$fill_in_ids)) > 0){
    warning("Some hgnc converted to more than one entrez id.  Need to figure out how to resolve this...")
  }
  if(sum(grepl(",", BM_results$fill_in_symbols)) > 0){
    warning("Some entrez ids converted to more than one symbol.  Need to figure out how to resolve this...")
  }
  
  
  BM_results$fin_ids[BM_results$fin_ids == "NA"] = NA
  BM_results$fin_symbols[BM_results$fin_symbols == "NA"] = NA
  BM_results$fin_symbols[BM_results$fin_symbols == ""] = NA
  
  a("* Dropped ensembl data columns with no information on hgnc or entrez.")
  BM_results = BM_results[ ((!is.na(BM_results$fin_ids)) & (!is.na(BM_results$fin_symbols))), ]

  
  
  BM_results = tidyr::unite(BM_results, combined_names, fin_symbols, fin_ids, sep = "|", remove = FALSE)
  # drop dat columns that aren't in BM_results
  dat = dat[, c("Sample_ID", names(dat)[names(dat) %in% BM_results$ucsc])] # went from 190K to 130K probably from loosing gene biotypes
  
  make_matrix_of_cols = c()
  if(output_hgnc_matrix) make_matrix_of_cols = c(make_matrix_of_cols, "fin_symbols")
  if(output_entrez_id_matrix) make_matrix_of_cols = c(make_matrix_of_cols, "fin_ids")
  if(output_piped_hugo_entrez_id_matrix) make_matrix_of_cols = c(make_matrix_of_cols, "combined_names")
  
  # need to transpose so we can have multiple gene names in the same row (can't do that with columns)
  my_samples = as.character(dat$Sample_ID)
  my_genes = names(dat)[-1]
  dat = dat[,-1]
  dat = data.table(Gene = my_genes, t(dat))
  names(dat) = c("Gene", my_samples)
  
  
  for(make_matrix_of_col in make_matrix_of_cols){
    
    if(make_matrix_of_col == "fin_symbols"){
      file_prefix = "hgnc_"
    } else if (make_matrix_of_col == "fin_ids"){
      file_prefix = "entrez_"
    } else if (make_matrix_of_col == "combined_names"){
      file_prefix = "hgnc_entrez_"
    } else {
      stop("Unknown make_matrix_of_col")
    }
    
    combined_name_lut = BM_results[[make_matrix_of_col]]
    names(combined_name_lut) = BM_results[["ucsc"]]
    
    my_dt = dat
    lut_new_names = combined_name_lut[my_dt$Gene] # faster than mclapply and lapply
    my_dt$Gene = lut_new_names
    
    my_dt = my_dt[,-1]
    my_dt = apply(my_dt, 2, function(x){tapply(x, lut_new_names, sum)}) # tapply reorders the items
    my_dt %<>% t %>% as.data.table 
    
    my_dt = data.table(Sample_ID = my_samples, my_dt)
    
    fwrite(my_dt, file.path(output_dir, paste0(file_prefix, counts_or_tpm,".tsv")), sep = "\t")
    
    if(output_upper_quartile_norm || output_log2_upper_quartile_norm){
      norm_dat = binfotron::normalize_rows_by_quartile(data.table(my_dt))
      if(output_upper_quartile_norm) fwrite(norm_dat,  file.path(output_dir, paste0(file_prefix, "norm_",counts_or_tpm,".tsv")), sep = "\t")
      if(output_log2_upper_quartile_norm){
        log2_norm_dat = binfotron::log_transform_plus(norm_dat)
        fwrite(log2_norm_dat,  file.path(output_dir, paste0(file_prefix, "log2_norm_",counts_or_tpm,".tsv")), sep = "\t")
      }
    }
  }
}

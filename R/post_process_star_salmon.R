
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' post_process_star_salmon
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title post_process_star_salmon 
#' 
#' @description 
#' Joins individual sample files into one tsv file.  Also converts ucsc names to hgnc_symbol|entrez_ids
#' and outputs the gene level counts.
#' 
#' @param input_file_paths Character vector of paths to the pipeline output data. 
#' @param output_dir Path to the output folder.
#' @param sample_data_path Path to the sample data which should contatin the sample_id_column and sample_folder_column
#' @param gene_biotypes The type of biomaRt gene_biotypes that should be output to the gene level output.
#' @param sample_folder_column The name of the column that has sample folder names
#' @param sample_id_column The name of the column that has sample ids
#' @param thread_num Integer number of threads to run mclapply statements
#' 
#' @return A path to the rds file.
#' 
#' @family mart
#' 
#' @export
post_process_star_salmon = function(
  # this_script_path = get_script_dir_path(include_file_name = T)
  # init_path = find_file_along_path(this_script_path, "_init.R")
  # source(init_path)
  # base_dir = dirname(init_path)
  
  input_file_paths,# = system(paste0("ls ", RAW_DATA_DIR, "/pipeline_output/star_salmon/*/*_quant.sf"), intern = TRUE)
  output_dir,# = file.path(base_dir, "post_processing", "star_salmon")
  sample_data_path,# = file.path(base_dir, "sample_data", "sample_data.tsv")
  
  gene_biotypes = c('IG_C_gene','IG_D_gene', 'IG_J_gene', 
                    'IG_V_gene','misc_RNA', 'processed_transcript',  'protein_coding',
                    'TR_C_gene', 'TR_D_gene', 'TR_J_gene',  'TR_V_gene',
                    'lincRNA'),
  sample_folder_column = "Sample_Folder",
  sample_id_column = "Sample_ID",
  thread_num = 16
){
  
  library(org.Hs.eg.db)
  library(annotate)
  library(binfotron)

  isoform_output_dir = file.path(output_dir,"ucsc_isoform_counts")
  isoform_output_path = file.path(isoform_output_dir, "ucsc_isoform_counts.tsv")
  dir.create(isoform_output_dir, showWarnings = F)
  
  readme_path = file.path(isoform_output_dir, "readme.txt")
  if(file.exists(readme_path)){ file.remove(readme_path)}
  
  a = function(...){
    my_output = paste0(...)
    if(!is.null(readme_path)){
      write(my_output, readme_path, append = TRUE)
    }
    message(my_output)
  }
  
  a(paste0("Making isoform counts matrix: ", this_script_path) %>% as.header1)
  a("")
  
  # Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # limit_import_to_num_rows = -1 # set to int to limit import for debuging, '-1' to import all data
  # # should_convert_NA_to_zeroes = T
  # should_log2_transform = T
  # # should_drop_rare_genes = T
  # #required_gene_expression = 0.7 # genes will be dropped if expressed in < required_gene_expression of samples
  # upper_quartile_normalize = TRUE
  # 
  # thread_num = 16
  
  
  # Input path ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # input_file_name = "McRee_salmon_quant.txt"
  a("Reading in files from input_file_paths:")
  read_data = mclapply(input_file_paths, function(file_path){
    a("  ", file_path)
    sample_data = fread(file_path, select = c("Name", "NumReads"))%>% as.data.frame
    return_vector = sample_data$NumReads
    names(return_vector) = sample_data$Name
    return(return_vector)
  }, mc.cores = thread_num)
  a("")
  
  # want to name these after the names on the sample sheet
  #   will use the sample folder names to look up the Sample_ID's
  sample_dat = fread(sample_data_path) %>% as.data.frame
  rownames(sample_dat) = sample_dat[,sample_folder_column]
  
  file_folders = lapply(input_file_paths, function(x){
    path_parts = unlist(strsplit(x, "/", fixed = T))
    return(path_parts[length(path_parts )-1])
  }) %>% unlist()
  
  # read_data order is based on input file path order & file_folder order is based on input file folder order
  #  so read data can be named with the sample_Ids in the order of the file folders
  names(read_data) = sample_dat[file_folders, "Sample_ID"]
  rm(sample_dat)
  
  # don't want to assume that all the data have the same column names
  unique_names = names(read_data[[1]])
  invisible( 
    mclapply(read_data, function(a_sample){
      unique_names <<- unique(c(unique_names, names(a_sample)))
    }, mc.cores = thread_num)
  )
  unique_names = sort(unique_names)
  # dat = data.frame(matrix(nrow = length(read_data), ncol = length(unique_names)))
  # names(dat) = c(unique_names)
  
  # make sure the order of columns is the same
  my_rows = mclapply(1:length(read_data), function(sample_index){
    read_data[[sample_index]][unique_names] %>% as.numeric
  }, mc.cores = thread_num)
  
  dat = data.frame(matrix(unlist(my_rows), ncol=length(unique_names), byrow=T, dimnames = list(c(), unique_names)))
  
  dat = data.frame(Sample_ID = names(read_data), dat)
  
  # manually checked that the order of the samples was right at this point by looking at a couple of quant.sf files.
  fwrite(dat, isoform_output_path, sep = "\t")
  
  # isoform output done
  
  
  
  # now do gene level counts
  
  gene_output_dir = file.path(output_dir, "hgnc_entrezid_gene_counts")
  gene_output_path = file.path(gene_output_dir, "hgnc_entrezid_gene_counts.tsv")
  dir.create(gene_output_dir, showWarnings = F)
  
  readme_path = file.path(gene_output_dir, "readme.txt")
  if(file.exists(readme_path)){file.remove(readme_path)}
  
  a = function(...){
    my_output = paste0(...)
    if(!is.null(readme_path)){
      write(my_output, readme_path, append = TRUE)
    }
    message(my_output)
  }
  
  
  a(paste0("Making gene counts matrix: ", this_script_path) %>% as.header1)
  a("")
  a("Converting ucsc ID's to Entrez to fit gene set collections and immune gene signatures...")
  a("Looking up hgnc_symbol, entrezgene and gene_biotype of ucsc id's." %>% as.bullet)
  a("")
  a("Using saved biomaRt from binfotron.")
  BM_results = readRDS(StarSalmon::get_biomart_hsa_ucsc_path())
  a("")
  
  BM_results$gene_biotype = factor(BM_results$gene_biotype)
  
  # for debugging this it's best to look at the BM_results matrix and make sure
  #   each step is filling in data the right way...
  
  a(paste0("Restricting gene_biotypes to: ", paste0(gene_biotypes, collapse = ", ")) %>% as.bullet)
  BM_results = BM_results[ BM_results$gene_biotype %in% gene_biotypes, ]
  
  missing_ids = which(is.na(BM_results$entrezgene))
  
  #lots of blanks in symbols
  missing_symbols = which(BM_results$hgnc_symbol == "")
  
  fill_in_ids_using_symbol_indices = missing_ids[missing_ids %ni% missing_symbols]
  fill_in_symbols_using_id_indices = missing_symbols[missing_symbols %ni% missing_ids]
  
  symbols_to_use = BM_results$hgnc_symbol[fill_in_ids_using_symbol_indices]
  
  a("Used hgnc symbols to fill in missing entrez ids." %>% as.bullet)
  # mclapply takes 3x longer
  fill_in_ids = unlist(lapply(1:length(symbols_to_use), function(x){ # very slow but mclapply runs into errors
    gene_id = as.character(paste0(unlist(lookUp(symbols_to_use[x],'org.Hs.eg','SYMBOL2EG')), collapse = ","))
    # gene_id = as.character(unlist(lookUp(symbols_to_use[x],'org.Hs.eg','SYMBOL2EG'))[1])# some return more than one value.  just taking the first one.
    return(gene_id)
  }))
  
  BM_results$symbols_to_use = NA
  BM_results$symbols_to_use[fill_in_ids_using_symbol_indices] = symbols_to_use # just marking these to make sure we are converting the correct things
  BM_results$fill_in_ids = NA
  BM_results$fill_in_ids[fill_in_ids_using_symbol_indices] = fill_in_ids
  
  
  ids_to_use = BM_results$entrezgene[fill_in_symbols_using_id_indices] %>% as.character
  
  a("Used entrez ids to fill in missing hgnc symbols." %>% as.bullet)
  # mclapply takes 3x longer
  fill_in_symbols = unlist(lapply(1:length(ids_to_use), function(x){ # very slow but mclapply runs into errors
    gene_name = as.character(paste0(unlist(lookUp(ids_to_use[x],'org.Hs.eg','SYMBOL')), collapse = ","))
    #gene_name = as.character(unlist(lookUp(ids_to_use[x],'org.Hs.eg','SYMBOL'))[1])# some return more than one value.  just taking the first one.
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
  
  a("Dropped ensembl data columns with no information on hgnc or entrez." %>% as.bullet)
  BM_results = BM_results[ ((BM_results$fin_ids %>% is_not_na) & (BM_results$fin_symbols %>% is_not_na)), ]
  
  
  BM_results = tidyr::unite(BM_results, combined_names, fin_symbols, fin_ids, sep = "|", remove = FALSE)
  
  BM_results = BM_results[, c("ucsc", "fin_symbols", "fin_ids", "combined_names")]
  
  # Dup_Ensemble could be Dup_UCSC... it's just the original gene names 
  BM_results$Dup_Ensembl = (BM_results$ucsc %in% BM_results$ucsc[duplicated(BM_results$ucsc)])
  BM_results$Dup_Comb = (BM_results$combined_names %in% BM_results$combined_names[duplicated(BM_results$combined_names)])
  
  # we have duplicated ensembl and duplicated combined names at this point.  some are even both!
  # the next step is to get the sorted levels of the combined names
  #  these will be our column names for our unnormalized data columns
  # then go through the names and get the ensembls that match them
  # add these columns and put the result in the unnormalized data colunn
  
  combined_names = unique(BM_results$combined_names)
  combined_names = combined_names[order(combined_names)]
  
  converted_dat = data.frame(matrix(nrow = nrow(dat), ncol = length(combined_names) +1 ))
  names(converted_dat) = c("Sample_ID", combined_names)
  converted_dat$Sample_ID = dat$Sample_ID
  
  a("Made combined names from hgnc|entrez" %>% as.bullet)
  a("Used combined_name to lookup all ensembl id producing that combined_name and summed those data columns to get the combined_name data." %>% as.bullet)
  a("")
  
  for(name in combined_names){
    ensemble_subdat = BM_results[BM_results$combined_names == name,]
    ensembls = ensemble_subdat$ucsc # these are the columns we need to add to get our data column
    
    data_subdat = dat[,ensembls, drop = F]
    return_data = apply(data_subdat, 1, sum, na.rm = F)
    converted_dat[[name]] = return_data
  }
  
  dat = converted_dat
  rm(converted_dat)
  
  dat = dat[order(dat$Sample_ID), ]
  
  # drop columns which have a sum of 0
  col_sums = apply(dat[,2:ncol(dat)], 2, sum)
  dat = dat[, c(TRUE, col_sums > 0)]
  
  fwrite(dat, gene_output_path, sep = "\t")
}
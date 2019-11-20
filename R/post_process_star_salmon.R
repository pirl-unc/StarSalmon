
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
  library(tximport)

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

  a("Reading in files from input_file_paths:")
  read_data = mclapply(input_file_paths, function(input_file_path){
    a("  ", input_file_path)
    counts_df = fread(input_file_path, select = c("Name", "NumReads"), data.table = F)
    return_list = lapply(counts_df[["NumReads"]], function(x)x) # rbindlist needs a list so we turn this into a list
    names(return_list) = counts_df[["Name"]] # Name the list items so they get assigned to the right column
    return(return_list)
  }, mc.cores = thread_num)
  a("")

  # convert to a data table according to the names of the list
  #   (ie, use.names = T won't assume the genes are in the same order)
  dat = rbindlist(read_data, use.names = TRUE, fill = TRUE)

  # now lets make a lookup table (lut) to get the sample names fom the folder names
  sample_dat = fread(sample_data_path, data.table = F)
  rownames(sample_dat) = sample_dat[[sample_folder_column]]
  file_folders = basename(dirname(input_file_paths))
  sample_lut = sample_dat$Sample_ID
  names(sample_lut) = sample_dat[[sample_folder_column]]
  
  # add the sample names
  dat = data.frame(Sample_ID = sample_lut[file_folders], dat)
  row.names(dat) = NULL

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
  
  tx2gene = BM_results[,1:2]

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

  # gene count output done

  

  # set up TPM output

  TPM_output_dir = file.path(output_dir, "TPM_averages")
  TPM_output_path = file.path(TPM_output_dir, "TPM_averages.tsv")
  dir.create(TPM_output_dir, showWarnings = F)

  readme_path = file.path(TPM_output_dir, "readme.txt")
  if(file.exists(readme_path)){ file.remove(readme_path)}

  a = function(...){
    my_output = paste0(...)
    if(!is.null(readme_path)){
      write(my_output, readme_path, append = TRUE)
    }
    message(my_output)
  }

  a(paste0("Use tximport pkg to import transcript-level abundance, estimate counts/transcript lengths, summarize into matrix", this_script_path) %>% as.header1)
  a("")

  temp_files_df = as.data.table(matrix(input_file_paths, nrow=length(input_file_paths), byrow=T), stringsAsFactors = FALSE)
  colnames(temp_files_df) = c("input_file_path")
  temp_files_df[, ID := gsub("_quant.sf", "", basename(temp_files_df$input_file_path))]
  if("Sample_Folder" %in% colnames(sample_dat)) {  #Riaz
    names_df = merge(temp_files_df, sample_dat, by.x = "ID", by.y = "Sample_Folder")
  }  else if ("run_accession" %in% colnames(sample_dat)) {  #HugoLo
    names_df = merge(temp_files_df, sample_dat, by.x = "ID", by.y = "run_accession")
  } else {
    rnaseq_file_path = paste0(RAW_DATA_DIR, "/rnaseq_fastq_paths/", "rnaseq_fastq_paths.tsv")
    SRR_dat = fread(rnaseq_file_path)
    if("run_accession" %in% colnames(SRR_dat)) {  #Prins
      names_df = merge(temp_files_df, SRR_dat, by.x = "ID", by.y = "run_accession")
    } else {
      names_df = merge(temp_files_df, SRR_dat, by.x = "ID", by.y = "Sample_ID", all.x = TRUE)  #Gide
      setnames(names_df, "ID", "Sample_ID")
    }
  }

  files <- input_file_paths
  names(files) = temp_files_df$ID
  txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
  
  # Write to output file
  write.table(data.frame("Sample"=rownames(txi.salmon$abundance),txi.salmon$abundance), TPM_output_path, row.names=FALSE, sep = "\t", quote = FALSE)
  
  # fwrite(txi.salmon$counts, TPM_output_path, sep = "\t")

}

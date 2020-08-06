#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_biomart_hg38_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_biomart_hg38_path 
#' 
#' @description 
#' Returns the path a saved mart
#' 
#' @param none
#' 
#' @return A path to the rds file.
#' 
#' @family mart
#' 
#' @export
get_biomart_hg38_path = function(){
  return(system.file("biomart", "hg38", "hsa_ensembl_ucsc.tsv", package = "StarSalmon"))
}


#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_biomart_grch38_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_biomart_grch38_path 
#' 
#' @description 
#' Returns the path a saved mart
#' 
#' @param none
#' 
#' @return A path to the rds file.
#' 
#' @family mart
#' 
#' @export
get_biomart_grch38_path = function(){
  return(system.file("biomart", "grch38", "bm_result.tsv", package = "StarSalmon"))
}


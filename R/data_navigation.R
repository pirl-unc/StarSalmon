#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' get_biomart_hsa_ucsc_path
#' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @title get_biomart_hsa_ucsc_path 
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
get_biomart_hsa_ucsc_path = function(){
  return(system.file("biomart", "hsa_ensembl_ucsc", "hsa_ensembl_ucsc.rds", package = "StarSalmon"))
}

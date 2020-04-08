Obtained 202004 via:

library("biomaRt")
# biomaRt::listMarts()
#ensembl <- biomaRt::useMart("ensembl") # host=useast??
gene_biotypes = c('IG_C_gene','IG_D_gene', 'IG_J_gene', 
                  'IG_V_gene','misc_RNA', 'processed_transcript',  'protein_coding',
                  'TR_C_gene', 'TR_D_gene', 'TR_J_gene',  'TR_V_gene',
                  'lincRNA')

mart = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl", 
                  host="useast.ensembl.org") # uk, useast, uswest, asia was intermitant, www was hardley ever working
#https://useast.ensembl.org/info/website/archives/index.html
# other hosts failed to find the data !@$!#$
# failed submisions
# failures mid submission

# datasets <- biomaRt::listDatasets(ensembl)
unique_names = my_dt[[1]]
# my_filters = listFilters(mart)
# my_attr = listAttributes(mart)
# ucsc  |  UCSC Stable ID(s) [e.g. ENST00000000233.9] # these aren't ucsc!?!
# had to run this a ton to get it to go all the way through.  was in the process of breaking it up into a loop when it ran through so I saved it.
failed = FALSE
bm_result = tryCatch({
  biomaRt::getBM(
    filters= "ucsc",
    attributes= c("ucsc", "hgnc_symbol", "entrezgene_id", "gene_biotype"),
    values= unique_names,
    mart= mart
  )
}, warning = function(w) {
  message("Got warning")
  failed = TRUE
}, error = function(e) {
  message("Got error")
  failed = TRUE
})

# to match the names of the other hg38 matrix I made
names(bm_result)[3] = "entrezgene"

saveRDS(bm_result, "bm_result.rds")

Post processes output from the star salmon pipeline.

## Example
``` r
post_process_salmon(
  this_script_path = housekeeping::get_script_dir_path(include_file_name = T),
  input_file_paths = <input_file_paths>,
  output_dir = <output_dir>,
  ref = "grch38",
  gene_biotypes = c('protein_coding', 
                    'IG_C_gene','IG_D_gene', 'IG_J_gene', 'IG_V_gene',
                    'TR_C_gene', 'TR_D_gene', 'TR_J_gene','TR_V_gene'),
  thread_num = 8,
  output_transcript_matrix = F,
  output_hgnc_matrix = T,
  output_entrez_id_matrix = F,
  output_piped_hugo_entrez_id_matrix = F,
  output_upper_quartile_norm = F,
  output_log2_upper_quartile_norm = F,
  counts_or_tpm = "counts"
)
```

## Assembling this package
In R:
``` r
housekeeping::assemble_package(package_name = "StarSalmon", my_version = "0.1-03",
  my_dir = "/datastore/alldata/shiny-server/rstudio-common/dbortone/packages/StarSalmon")
```

## Push changes
In bash:
``` bash
cd /datastore/alldata/shiny-server/rstudio-common/dbortone/packages/StarSalmon
my_comment="Minor comment changes.  Added output of all gene_biotypes to readme.txt"
git commit -am "$my_comment"; git push origin master
git tag -a 0.1-03 -m "$my_comment"; git push -u origin --tags
```

## Install
Restart R
In R (local library, packrat library):
``` r
devtools::install_github("Benjamin-Vincent-Lab/StarSalmon")
```

Or for a specific version:
``` r
devtools::install_github("Benjamin-Vincent-Lab/StarSalmon", ref = "0.1-03")
```

## Previous locations
https://sc.unc.edu/dbortone/starsalmon
https://sc.unc.edu/benjamin-vincent-lab/starsalmon
Moved to github so that the package could be accessed without a token.

## Creating BM_results
The following code was used to make the grch38 bm_results
``` r
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
BM_results = tryCatch({
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
```
Oddly here, the 'ucsc' column isn't 'ucsc' but 'ensembl'.  I kept the column with that name as it made switching to the real ucsc table easier.  Additionally I renamed the 'entrezgene_id' to 'entrezgene.' Connecting with the dataabse above was very problematic.  It failed to connect 1 out of 5 times and when it did connnect it didn't finish.  I was giving up on it and was going to write a loop to keep sending smaller batches using the try catch statement when finally the whole thing went through.  For future uses a loop is the way to go.  Also don't expect the column names to stay stable.  They change these on almost a monthly basis.  I'd love to switch to something other than biomaRt, but unfortunately AFAIK there isn't anything else.


## Comments
Using test_code.R, I checked if hgnc or entrez was better for not having one ensemble map to them.  Almost all of the duplicates for the ensembl id were from one ensembl id mapping to multiple entrez. Very few of the hgnc caused multimappings.  
I also checked if using entrez to lookup hgnc and visa versa caused more multimappings.  It wasn't a huge contributor and looking up the genes did find a lot fo new genes:  ~70 for hgnc and ~300 for entrez.
I skiped using tximport.  I don't see the value of this package, since I'd have to make the tx2gene matrix anyway. That's the hard part.
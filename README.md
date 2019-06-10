Functions for writing commands and post processing the Star Salmon pipeline.

Given a sample data matrix that indicates the sample name, sample folder name and input paths
this script can write out pipeline commands to run star salmon on the cluster and post process the
data into ucsc_isoform_counts.tsv and hgnc_entrezid_gene_counts.tsv.  Only ucsc pipeline outputs 
are handled at this time.


## Assembling this package
In R:
``` r
housekeeping::assemble_package(package_name = "StarSalmon", my_version = "0.0-04",
  my_dir = "/datastore/alldata/shiny-server/rstudio-common/dbortone/packages/StarSalmon")
```

## Push changes
In bash:
``` bash
cd /datastore/alldata/shiny-server/rstudio-common/dbortone/packages/StarSalmon
my_comment="Added binfotron to dependencies for using the header and footer formatting."
git commit -am "$my_comment"; git push origin master
git tag -a 0.0-04 -m "$my_comment"; git push -u origin --tags
```

## Install
Restart R
In R (local library, packrat library):
``` r
devtools::install_github("DanteBortone/StarSalmon")
```

Or for a specific version:
``` r
devtools::install_github("DanteBortone/StarSalmon", ref = "0.0-04")
```
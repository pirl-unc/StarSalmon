# this script was run on 20200805 to recover some of the hgnc and entrez id's that were missing from the 
#  biomart conversion from enst and ucsc to hgnc and entrez id
#  Briefly I took the AnnotationDbi ENSTs that had neither the hgnc nor entrez in my biomart results and 
#  added them to my biomart results.  This recovered 1882 hgnc & entrez genes.  I could have gotten more if 
#  I added those without hgnc | entrez but that would have added a lot of duplicates that I don't have time 
#  to deal with at the moent.  Just using AnnotationDbi wasn't an option either as this was missing alot of 
#  genes that biomart had, as in ALL of the TCR and BCR genes.
#  I tried to apply the same to hgnc but was having trouble finding the right db for this and hg38 isn't the 
#  priority right now.  for hgnc I just changed the header. also changed to .tsv output so the repos can 
#  handle changes better in the futre.


# BiocManager::install("AnnotationDbi")
# BiocManager::install('org.Hs.eg.db')
# BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

`%ni%`<- Negate(`%in%`)

grch38_path = "/Users/dantebortone/Desktop/bm_result.rds"
hg38_path = "/Users/dantebortone/Desktop/hsa_ensembl_ucsc.rds"

bm_df = readRDS(grch38_path)
bm_df$hgnc_symbol[bm_df$hgnc_symbol == ""] = NA

# drop those which have no info for either
bm_df = bm_df[!(is.na(bm_df$hgnc_symbol) & is.na(bm_df$entrezgene)),]

# rename this column
names(bm_df)[names(bm_df) == "ucsc"] = 'transcript'
# drop the decimals so we don't miss any that might be off by just the decimal
bm_df$transcript = substr(bm_df$transcript, 1,15) # this didn't make any new duplicates


# hg38_bm_df = readRDS(hg38_path)
# names(hg38_bm_df)[1] = 'transcript'
# split protein and not protein coding, so we can focus on these
# not_prot_bm_df = bm_df[bm_df$gene_biotype != "protein_coding", ]
# prot_bm_df = bm_df[bm_df$gene_biotype == "protein_coding", ]
# rm(bm_df)
# hg38_bm_df = hg38_bm_df[hg38_bm_df$gene_biotype == "protein_coding", ]

# ls("package:org.Hs.eg.db")
# columns(org.Hs.eg.db)
# keytypes(org.Hs.eg.db)
# head(keys(org.Hs.eg.db, keytype="SYMBOL"))
# from our pipeline
sf_file = "/Users/dantebortone/Downloads/IMVigor210-SAM1f66db567eb5-GO29293_ngs_rna_targrna_SAM1f66db567eb5_20160128.quant.sf"
sf_df = fread(sf_file)
enstrans_df = select(org.Hs.eg.db, keys=sf_df$Name, columns = c("SYMBOL","ENTREZID"), keytype="ENSEMBLTRANS")
enstrans_df = enstrans_df[!(is.na(enstrans_df$SYMBOL) & is.na(enstrans_df$ENTREZID)),]
# these are all complete cases now
# grab the genes that don't have a SYMBOL and ID in our biomart results
missing_df = enstrans_df[(enstrans_df$SYMBOL %ni% bm_df$hgnc_symbol) & (enstrans_df$ENTREZID %ni% bm_df$entrezgene), ]
names(missing_df) = c("transcript","hgnc_symbol","entrezgene")
missing_df$gene_biotype = "from AnnotationDbi"


# recovered 1882 entrez ids and hgnc symbols
length(missing_df$hgnc_symbol[missing_df$hgnc_symbol %ni% bm_df$hgnc_symbol])
length(missing_df$entrezgene[missing_df$entrezgene %ni% bm_df$entrezgene])

out_bm_df = rbindlist(list(bm_df, missing_df))
fwrite(out_bm_df, "/Users/dantebortone/Desktop/repaired/bm_result.tsv", sep = "\t")
# saveRDS(out_bm_df, "/Users/dantebortone/Desktop/repaired/bm_result.rds")



# tried to fix hg38 but was having trouble finding the right db for this and hg38 isn't the priority right now
# just changed the header
bm_df = readRDS(hg38_path)
bm_df$hgnc_symbol[bm_df$hgnc_symbol == ""] = NA

# drop those which have no info for either
bm_df = bm_df[!(is.na(bm_df$hgnc_symbol) & is.na(bm_df$entrezgene)),]

# rename this column
names(bm_df)[names(bm_df) == "ucsc"] = 'transcript'
# drop the decimals so we don't miss any that might be off by just the decimal
# bm_df$transcript = substr(bm_df$transcript, 1,15) # this didn't make any new duplicates

# hg38_bm_df = readRDS(hg38_path)
# names(hg38_bm_df)[1] = 'transcript'
# split protein and not protein coding, so we can focus on these
# not_prot_bm_df = bm_df[bm_df$gene_biotype != "protein_coding", ]
# prot_bm_df = bm_df[bm_df$gene_biotype == "protein_coding", ]
# rm(bm_df)
# hg38_bm_df = hg38_bm_df[hg38_bm_df$gene_biotype == "protein_coding", ]

# ls("package:org.Hs.eg.db")
# columns(org.Hs.eg.db)
# keytypes(org.Hs.eg.db)
# head(keys(org.Hs.eg.db, keytype="SYMBOL"))
# hg38 referenced sample from our pipeline
# sf_file = "/Users/dantebortone/Downloads/ipiPD1_Patient10_PRE_quant.csv"
# sf_df = read.csv(sf_file)
# enstrans_df = select(org.Hs.eg.db, keys=as.character(sf_df$Name), columns = c("SYMBOL","ENTREZID"), keytype="UCSCKG")
# enstrans_df = enstrans_df[!(is.na(enstrans_df$SYMBOL) & is.na(enstrans_df$ENTREZID)),]
# # these are all complete cases now
# # grab the genes that don't have a SYMBOL and ID in our biomart results
# missing_df = enstrans_df[(enstrans_df$SYMBOL %ni% bm_df$hgnc_symbol) & (enstrans_df$ENTREZID %ni% bm_df$entrezgene), ]
# names(missing_df) = c("transcript","hgnc_symbol","entrezgene")
# missing_df$gene_biotype = "from AnnotationDbi"
# 
# # saved 1920 hgnc and 3007 entrez ids
# length(missing_df$hgnc_symbol[missing_df$hgnc_symbol %ni% bm_df$hgnc_symbol])
# length(missing_df$entrezgene[missing_df$entrezgene %ni% bm_df$entrezgene])
# 
# out_bm_df = rbindlist(list(bm_df, missing_df))
fwrite(bm_df, "/Users/dantebortone/Desktop/repaired/hsa_ensembl_ucsc.tsv", sep = "\t")

# saveRDS(bm_df, "/Users/dantebortone/Desktop/repaired/hsa_ensembl_ucsc.rds")








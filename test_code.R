# using this code I wanted to see if hgnc or entrez was better for not having one ensemble map to them.
# also checked if ensemble had better one to one with either hgnc or entrez.  
# I also wanted to see if using entrez id's to lookup missing hgnc and visa vesa is generally helpful
#   I don't want to do it if that means I make an entrez that was mapping to 1 hgnc now map to 2
#  all of the measurements got closer to one whch means for the most part I'm adding information and one-to-one mappings

# next I looked at the ucsc
# ucsc was a tiny bit worse but not by any meaningful amount
# ucsc did find more genes about 1000 more genes

# next I check how bad the combined names made things
# it really wasn't a big deal
# still I will avoid using the combined set because of the impact it will have on gene signatures

# library(data.table)
# library(org.Hs.eg.db)
# `%ni%`<- Negate(`%in%`)
this_script_path = housekeeping::get_script_dir_path(include_file_name = T)
ref = "grch38" # options 'grch38' with ensemble output or 'hg38' with ucsc output
gene_biotypes = c('protein_coding',
                  'IG_C_gene','IG_D_gene', 'IG_J_gene', 'IG_V_gene',
                  'TR_C_gene', 'TR_D_gene', 'TR_J_gene','TR_V_gene')
thread_num = 8


input_file_paths = c("/datastore/nextgenout5/share/labs/Vincent_Lab/members/dbortone/raft/analyses/rna_quant_test/outputs/Hugo_Pt1/SRR3184279/salmon-1.1.0=82305/Homo_sapiens.GRCh38.cdna.all.fa=a20ef/Hugo_Pt1_SRR3184279.Aligned.sortedByCoord.out.bam=226e9/Hugo_Pt1_SRR3184279.quant.sf", "/datastore/nextgenout5/share/labs/Vincent_Lab/members/dbortone/raft/analyses/rna_quant_test/outputs/Hugo_Pt2/SRR3184280/salmon-1.1.0=cf105/Homo_sapiens.GRCh38.cdna.all.fa=a20ef/Hugo_Pt2_SRR3184280.Aligned.sortedByCoord.out.bam=2b56f/Hugo_Pt2_SRR3184280.quant.sf", "/datastore/nextgenout5/share/labs/Vincent_Lab/members/dbortone/raft/analyses/rna_quant_test/outputs/Hugo_Pt4/SRR3184281/salmon-1.1.0=fa50c/Homo_sapiens.GRCh38.cdna.all.fa=a20ef/Hugo_Pt4_SRR3184281.Aligned.sortedByCoord.out.bam=9c653/Hugo_Pt4_SRR3184281.quant.sf")
output_dir = "/datastore/nextgenout5/share/labs/Vincent_Lab/members/dbortone"
sample_names = sapply(input_file_paths, function(input_file_path){
  path_segments = strsplit(input_file_path, split = "/")[[1]]
  sample_name = path_segments[which(path_segments == "outputs") + 1]
  return(sample_name)
  }, USE.NAMES = FALSE)

names(input_file_paths) = sample_names


# BM_results = readRDS("/rstudio-common/dbortone/packages/StarSalmon/inst/biomart/grch38/bm_result.rds")
# gene_biotypes = c('protein_coding',
#                   'IG_C_gene','IG_D_gene', 'IG_J_gene', 'IG_V_gene',
#                   'TR_C_gene', 'TR_D_gene', 'TR_J_gene','TR_V_gene')
# 
# output_transcript_matrix = FALSE
# output_hgnc_matrix = TRUE
# output_entrez_id_matrix = TRUE
# output_piped_hugo_entrez_id_matrix = TRUE
# output_upper_quartile_norm = TRUE
# output_log2_upper_quartile_norm = TRUE
# counts_or_tpm = "counts"


StarSalmon::post_process_salmon(
  this_script_path = housekeeping::get_script_dir_path(include_file_name = T),
  input_file_paths = input_file_paths,
  output_dir = output_dir,
  ref = "grch38",
  gene_biotypes = c('protein_coding', 
                    'IG_C_gene','IG_D_gene', 'IG_J_gene', 'IG_V_gene',
                    'TR_C_gene', 'TR_D_gene', 'TR_J_gene','TR_V_gene'),
  thread_num = 8,
  output_transcript_matrix = F,
  output_hgnc_matrix = TRUE,
  output_entrez_id_matrix = F,
  output_piped_hugo_entrez_id_matrix = F,
  output_upper_quartile_norm = F,
  output_log2_upper_quartile_norm = F,
  counts_or_tpm = "counts"
)


# now looking at duplicate ENST's only of complete BM_resutls

duplicate_bm = BM_results[BM_results$ucsc %in% unique(BM_results$ucsc[duplicated(BM_results$ucsc)]),]
duplicate_bm = tidyr::unite(duplicate_bm, ucsc_hgnc, ucsc, fin_symbols, sep = "|", remove = FALSE)
duplicate_bm = tidyr::unite(duplicate_bm, ucsc_entrez, ucsc, fin_ids, sep = "|", remove = FALSE)
duplicate_bm = tidyr::unite(duplicate_bm, ucsc_both, ucsc, combined_names, sep = "|", remove = FALSE)

sum(is.na(duplicate_bm$ucsc_hgnc))
sum(is.na(duplicate_bm$ucsc_entrez))
sum(is.na(duplicate_bm$ucsc_both))


unique_ensembl = unique(duplicate_bm$ucsc) # grch 700
length(unique_ensembl)
unique_hgnc = unique(duplicate_bm$ucsc_hgnc) # grch 711
length(unique_hgnc)
unique_entrez = unique(duplicate_bm$ucsc_entrez) # grch 1523
length(unique_entrez)
unique_both = unique(duplicate_bm$ucsc_both) # grch 1541
length(unique_both)

# which results in less multimapping hgnc or entrez
unique_ucsc_counts_per_hgnc = tapply(BM_results$ucsc, BM_results$fin_symbols, function(x){length(unique(x))})
mean(unique_ucsc_counts_per_hgnc) # grch 7.397086
unique_hgnc_counts_per_ucsc = tapply(BM_results$fin_symbols, BM_results$ucsc, function(x){length(unique(x))})
mean(unique_hgnc_counts_per_ucsc) # grch 1.000081
unique_ucsc_counts_per_entrez = tapply(BM_results$ucsc, BM_results$fin_ids, function(x){length(unique(x))})
mean(unique_ucsc_counts_per_entrez) # 7.412676
unique_entrez_counts_per_ucsc = tapply(BM_results$fin_ids, BM_results$ucsc, function(x){length(unique(x))})
mean(unique_entrez_counts_per_ucsc) # 1.006072

unique_ucsc_counts_per_hgnc = tapply(duplicate_bm$ucsc, duplicate_bm$fin_symbols, function(x){length(unique(x))})
mean(unique_ucsc_counts_per_hgnc) # grch 6.836538
unique_hgnc_counts_per_ucsc = tapply(duplicate_bm$fin_symbols, duplicate_bm$ucsc, function(x){length(unique(x))})
mean(unique_hgnc_counts_per_ucsc) # grch 1.015714
unique_ucsc_counts_per_entrez = tapply(duplicate_bm$ucsc, duplicate_bm$fin_ids, function(x){length(unique(x))})
mean(unique_ucsc_counts_per_entrez) # 7.502463
unique_entrez_counts_per_ucsc = tapply(duplicate_bm$fin_ids, duplicate_bm$ucsc, function(x){length(unique(x))})
mean(unique_entrez_counts_per_ucsc) # 2.175714

# hgnc results in less multimapping



# does making combined names make multimaping worse?
unique_ucsc_counts_per_combined = tapply(BM_results$ucsc, BM_results$combined_names, function(x){length(unique(x))})
mean(unique_ucsc_counts_per_combined) # grch 7.391556
unique_combined_counts_per_ucsc = tapply(BM_results$combined_names, BM_results$ucsc, function(x){length(unique(x))})
mean(unique_combined_counts_per_ucsc) # grch 1.006204

unique_ucsc_counts_per_combined = tapply(duplicate_bm$ucsc, duplicate_bm$combined_names, function(x){length(unique(x))})
mean(unique_ucsc_counts_per_combined) # grch 6.758772
unique_combined_counts_per_ucsc = tapply(duplicate_bm$combined_names, duplicate_bm$ucsc, function(x){length(unique(x))})
mean(unique_combined_counts_per_ucsc) # grch 2.201429
# yes, making combined names make multimaping worse



# did looking up the ids make matters worse
hgnc_original_bm = BM_results[BM_results$hgnc_symbol != "", ]

unique_ucsc_counts_per_hgnc = tapply(hgnc_original_bm$ucsc, hgnc_original_bm$hgnc_symbol, function(x){length(unique(x))})
mean(unique_ucsc_counts_per_hgnc) # grch 7.39711 compared to 7.397086
unique_hgnc_counts_per_ucsc = tapply(hgnc_original_bm$hgnc_symbol, hgnc_original_bm$ucsc, function(x){length(unique(x))})
mean(unique_hgnc_counts_per_ucsc) # grch 1.000081 compared to 1.000081

entrez_original_bm = BM_results[!is.na(BM_results$entrezgene), ]
unique_ucsc_counts_per_entrez = tapply(entrez_original_bm$ucsc, entrez_original_bm$entrezgene, function(x){length(unique(x))})
mean(unique_ucsc_counts_per_entrez) # 7.594286 compared to 7.412676
unique_entrez_counts_per_ucsc = tapply(entrez_original_bm$entrezgene, entrez_original_bm$ucsc, function(x){length(unique(x))})
mean(unique_entrez_counts_per_ucsc) # 1.006101 compared to 1.006072
# not much and more id's were found


length(unique(BM_results$hgnc_symbol)) # 18273
length(unique(BM_results$fin_symbols)) # 18326


length(unique(BM_results$entrezgene)) # 17851
length(unique(BM_results$fin_ids)) #18397



post_process_salmon(
  this_script_path = housekeeping::get_script_dir_path(include_file_name = T),
  input_file_paths = input_file_paths,# = system(paste0("ls ", RAW_DATA_DIR, "/pipeline_output/star_salmon/*/*_quant.sf"), intern = TRUE)
  output_dir = output_dir,# = file.path(base_dir, "post_processing", "star_salmon")
  ref = "grch38", # options 'grch38' with ensemble output or 'hg38' with ucsc output
  gene_biotypes = c('protein_coding',
                    'IG_C_gene','IG_D_gene', 'IG_J_gene', 'IG_V_gene',
                    'TR_C_gene', 'TR_D_gene', 'TR_J_gene','TR_V_gene'),
  thread_num = 8,
  output_transcript_matrix = TRUE,
  output_hgnc_matrix = TRUE,
  output_entrez_id_matrix = TRUE,
  output_piped_hugo_entrez_id_matrix = TRUE,
  output_upper_quartile_norm = TRUE,
  output_log2_upper_quartile_norm = TRUE,
  counts_or_tpm = "counts"
)

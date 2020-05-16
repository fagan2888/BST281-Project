library(ggplot2)
library(ComplexHeatmap)

gencode <- read.delim('gencode.v22.genes.txt') %>% select(gene_id, gene_name)

# gene_list = c('IGF2', 'NOTUM', 'ACAN', 'ALDH3A1', 'CCDC146', 'CES1', 'CHGB',
#               'COL22A1', 'CYP4A11', 'FABP4')

# gene_list = c('DDX43', 'PROM1', 'SLC22A12', 'ASTN1', 'CCDC146', 'CDHR1', 'CTNNA2', 'CYP4A11',
#   'CYP4A22','GLYAT', 'GRIK5', 'HOXA13', 'HOXD11', 'PCP4', 'SLC5A8', 'TRIM50', 'UNC5D')

gene_list = c('IGF2', 'NOTUM', 'ACAN','ALDH3A1', 'CCDC146', 'CES1', 'CHGB', 
'COL22A1', 'CYP4A11', 'FABP4', 'FDCSP', 'GPR64', 'HOXA11', 'L1CAM','NTN1','PROM2','SHISA3','SLC38A11', 'TMEM100', 'TNFSF4','TNNT3','TREX2','WNK2')

gencode_genes = gencode %>% filter(gencode$gene_name %in% gene_list)

df_tot = NULL

for (dds_idx in ls(dds_env)){
    counts_norm <- counts(dds_env[[dds_idx]], normalize=TRUE)
    counts_norm <- as.data.frame(counts_norm)
    counts_norm <- counts_norm[gencode_genes$gene_id,]
    
    counts_norm <- dplyr::left_join(counts_norm %>% mutate(gene_id = rownames(counts_norm)),
                     gencode_genes,
                     by = 'gene_id') %>%
      select(-c('gene_id')) %>%
      tibble::column_to_rownames(var="gene_name")
    
    cancer_type <- unlist(strsplit(dds_idx, 'ddsMat'))
    metadata_fnam <- paste0("data/metadata/metadata.files_", cancer_type, ".txt")
    metadata <- read.delim(metadata_fnam, header=T) %>%
      select('cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id',
              'cases.0.demographic.gender') %>%
      mutate('sample'=gsub('-', '.', cases.0.samples.0.portions.0.analytes.0.aliquots.0.submitter_id))
    

    meta_male <- metadata %>% filter(cases.0.demographic.gender == 'male')
    meta_female <- metadata %>% filter(cases.0.demographic.gender == 'female')
    # # 
    # counts_norm_male <- counts_norm %>%
    #   select(meta_male$sample) %>%
    #   rowMeans() %>%
    #   as.data.frame() %>%
    #   rename(male='.') %>%
    #   t() %>%
    #   as.data.frame() %>%
    #   tibble::rownames_to_column('gender')
    # 
    counts_norm_male <- counts_norm %>%
      select(meta_male$sample) %>%
      rowMeans() %>%
      as.data.frame() %>%
      rename(!!cancer_type:='.') %>%
      t() %>%
      as.data.frame() %>%
      rename_all(paste0, ' male')
# 
#     counts_norm_female <- counts_norm %>%
#       select(meta_female$sample) %>%
#       rowMeans() %>%
#       as.data.frame() %>%
#       rename(female='.') %>%
#       t() %>%
#       as.data.frame() %>%
#       tibble::rownames_to_column('gender')

    counts_norm_female <- counts_norm %>%
      select(meta_female$sample) %>%
      rowMeans() %>%
      as.data.frame() %>%
      rename(!!cancer_type:='.') %>%
      t() %>%
      as.data.frame() %>%
      rename_all(paste0, ' female')
    
    df <- cbind(counts_norm_male, counts_norm_female)
    
    df_tot <- rbind(df_tot, df)
    
    }

df_tot <- df_tot[,order(desc(colnames(df_tot)))]

df_tot <- log(df_tot + 1)

gene_num = length(gene_list)


Heatmap(as.matrix(df_tot), name = "log(count)", row_names_side = "left", row_dend_side = "right", 
        column_names_side = "top", cluster_columns = FALSE, column_split = rep(1:gene_num, each=2),
        col = viridis(10, direction=-1))

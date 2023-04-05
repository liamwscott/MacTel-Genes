library(tidyverse) 
library(BiocManager)
library(biomaRt)
library(org.Hs.eg.db)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

########### L O A D  S O U R C E  F I L E S ################

published_genes = read_delim("source_files/mactel_genes.txt")
gwas3_genes = read_delim("source_files/GWAS3_genes_fuma.txt")


################ Get Significant Genes only #################

gwas3_sig = gwas3_genes %>% filter(minGwasP < 5e-08)

############### Get Ensemble IDs ###########################

published_genes$ensembl_id = NA
gwas3_sig$source = "GWAS3"

gwas3_sig = gwas3_sig %>% dplyr::select(symbol, source, ensg)
colnames(gwas3_sig) = colnames(published_genes)
mactel_genes = rbind(gwas3_sig, published_genes)
mactel_genes = mactel_genes[!duplicated(mactel_genes$gene),]

rows_to_replace = mactel_genes[is.na(mactel_genes$ensembl_id),]


rows_to_replace$ensembl_id <- mapIds(org.Hs.eg.db, keys = rows_to_replace$gene, 
                      keytype = "SYMBOL", column="ENSEMBL") 

mactel_genes[is.na(mactel_genes$ensembl_id),] = rows_to_replace

###################### Add Locus Information #########################  
  
locus = as.data.frame(getBM(attributes = c( 'chromosome_name', 'start_position',
                                            'end_position', 'ensembl_gene_id', 
                                            'phenotype_description'),
                            filters = 'ensembl_gene_id', 
                            values = mactel_genes$ensembl_id, 
                            mart = ensembl))

locus = (locus %>%
            group_by(chromosome_name, start_position, end_position, ensembl_gene_id) %>%
            summarize(phenotype = str_c(phenotype_description, collapse = "; ")))

mactel_gene_table = merge(mactel_genes, locus, by.x = "ensembl_id",
                          by.y = "ensembl_gene_id", all.x = T)

mactel_gene_table_cavalier = mactel_gene_table
mactel_gene_table_cavalier$list_id
mactel_gene_table_cavalier$list_name 

####################### Write Files #################################

write_delim(mactel_gene_table, file = "mactel_genes.txt"). # full list



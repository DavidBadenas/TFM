# Cargar las librerías necesarias
library(biomaRt)
library(dplyr)
library(openxlsx)

# Leer el archivo CSV con los SNPs significativos para próstata
SNPs_prostate <- read.csv(file = "Significant_SNPs_prostate.csv", sep = ";")

# Limpiar los identificadores de SNPs, extrayendo la primera parte antes de ":"
SNPs_prostate$SNP <- sapply(strsplit(SNPs_prostate$SNP, ":"), `[`, 1)

# Conectarse a Ensembl utilizando el biomart para SNPs humanos
snpmart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")

# Opcional: listar los atributos y filtros disponibles en el biomart (para exploración)
listAttributes(snpmart)
listFilters(snpmart)

# Obtener los genes asociados a los SNPs
snpsGenes <- getBM(attributes = c('refsnp_id','ensembl_gene_stable_id'), 
                   filters = 'snp_filter', 
                   values = SNPs_prostate$SNP, 
                   mart = snpmart)

# Conectarse a Ensembl utilizando el biomart para genes humanos
genemart <- useEnsembl(biomart = "genes", 
                       dataset = "hsapiens_gene_ensembl")

# Obtener los símbolos de los genes asociados
genesSymbol <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = snpsGenes$ensembl_gene_stable_id,
                     mart = genemart)

# Combinar la información de SNPs con los símbolos de los genes
snpsSymbols <- merge(snpsGenes, genesSymbol, 
                     by.x = "ensembl_gene_stable_id", 
                     by.y = "ensembl_gene_id")

# Renombrar la columna 'refsnp_id' a 'SNP'
snpsSymbols <- snpsSymbols %>%
  rename("SNP" = "refsnp_id")

# Combinar los datos originales de los SNPs con la información de los genes
result_prostate <- left_join(SNPs_prostate, snpsSymbols, by = "SNP")

# Guardar el resultado en un archivo Excel
write.xlsx(result_prostate, file = "SNPs_prostate_genes.xlsx")

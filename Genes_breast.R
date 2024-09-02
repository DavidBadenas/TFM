# Cargar las librerías necesarias
library(biomaRt)
library(dplyr)
library(openxlsx)

# Leer los datos de los SNPs significativos para el cáncer de mama
SNPs_breast <- read.csv(file = "Significant_SNPs_breast.csv")

# Extraer el primer componente del identificador de SNP (antes de ':')
SNPs_breast$SNP <- sapply(strsplit(SNPs_breast$SNP, ":"), `[`, 1)

# Conectar a la base de datos Ensembl para SNPs humanos
snpmart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")

# (Opcional) Listar atributos y filtros disponibles en el mart de SNPs
listAttributes(snpmart)
listFilters(snpmart)

# Obtener los identificadores de genes asociados a los SNPs
snpsGenes <- getBM(attributes = c('refsnp_id', 'ensembl_gene_stable_id'), 
                   filters = 'snp_filter', 
                   values = SNPs_breast$SNP, 
                   mart = snpmart)

# Conectar a la base de datos Ensembl para genes humanos
genemart <- useEnsembl(biomart = "genes", 
                       dataset = "hsapiens_gene_ensembl")

# Obtener los símbolos de los genes asociados a los identificadores de Ensembl
genesSymbol <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = snpsGenes$ensembl_gene_stable_id,
                     mart = genemart)

# Unir la información de los SNPs con los símbolos de los genes
snpsSymbols <- merge(snpsGenes, genesSymbol, 
                     by.x = "ensembl_gene_stable_id", 
                     by.y = "ensembl_gene_id")

# Renombrar la columna 'refsnp_id' a 'SNP' para la unión final
snpsSymbols <- snpsSymbols %>%
  rename("SNP" = "refsnp_id")

# Unir los datos originales de SNPs con los símbolos de los genes
result_breast <- left_join(SNPs_breast, snpsSymbols, by = "SNP")

# Guardar el resultado en un archivo Excel
write.xlsx(result_breast, file = "SNPs_breast_genes.xlsx")

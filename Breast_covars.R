# Cargar las librerías necesarias
library(readxl)
library(dplyr)

# Leer los datos
covar <- read.table("covar_mds.txt", header = TRUE)
mama <- read_excel("imputed_data_breast.xlsx")

# Renombrar la columna 'Sample ID' a 'FID'
mama <- mama %>%
  rename(FID = `Sample ID`)

# Unir los datos por la columna 'FID'
breast_covars <- inner_join(covar, mama, by = "FID")

# Seleccionar las columnas de interés
breast_covars <- breast_covars %>%
  select(FID, IID, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, Site, axillary_surgery, 
         age_at_radiotherapy_start_yrs, chemo, radio_breast_ct_volume_cm3, radio_boost, 
         Breast.cancer.phenotype, Breast.cancer.stage)

# Verificar las categorías en la columna 'Site'
table(breast_covars$Site)

# Mapear los nombres de los sitios a números
site_map <- c("Barcelona" = 1, "Gent" = 2, "Leicester" = 3, "Leuven" = 4, 
              "Mannheim" = 5, "Milan" = 6, "Montpellier" = 7, 
              "Mount Sinai" = 8, "Santiago" = 9)
breast_covars$Site <- site_map[breast_covars$Site]

# Verificar los cambios en la columna 'Site'
table(breast_covars$Site)

# Verificar las categorías en la columna 'Breast.cancer.phenotype'
table(breast_covars$Breast.cancer.phenotype)

# Mapear los fenotipos del cáncer de mama a números
phenotype_map <- c("DCIS" = 1, "HER2 positive" = 2, "Luminal" = 3, 
                   "Luminal B HER2 positive" = 4, "Triple negative" = 5)
breast_covars$Breast.cancer.phenotype <- phenotype_map[breast_covars$Breast.cancer.phenotype]

# Verificar los cambios en la columna 'Breast.cancer.phenotype'
table(breast_covars$Breast.cancer.phenotype)

# Verificar las categorías en la columna 'Breast.cancer.stage'
table(breast_covars$Breast.cancer.stage)

# Mapear las etapas del cáncer de mama a números
stage_map <- c("Stage 0" = 1, "Stage I" = 2, "Stage II" = 3, "Stage III" = 4)
breast_covars$Breast.cancer.stage <- stage_map[breast_covars$Breast.cancer.stage]

# Verificar los cambios en la columna 'Breast.cancer.stage'
table(breast_covars$Breast.cancer.stage)

# Mostrar el dataframe final
breast_covars

# Guardar el dataframe en un archivo de texto
write.table(breast_covars, file = "breast_covars.txt", sep = ",", row.names = FALSE, quote = FALSE)

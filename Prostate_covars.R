# Cargar las librerías necesarias
library(readxl)
library(dplyr)

# Leer los datos
covar <- read.table("covar_mds.txt", header = TRUE)
prostata <- read_excel("imputed_data_prostate 2.xlsx")

# Renombrar la columna 'Sample ID' a 'FID'
prostata <- prostata %>%
  rename(FID = `Sample ID`)

# Unir los datos por la columna 'FID'
prostate_covars <- inner_join(covar, prostata, by = "FID")

# Seleccionar las columnas de interés
prostate_covars <- prostate_covars %>%
  select(FID, IID, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, Site, pre_radio_turp, 
         radio_pelvic, ace_inhibitor, Prostate.Cancer.Stage, other_collagen_vascular_disease, 
         previous_abdominal_surgery, psa_prediagnostic_biopsy_ng_ml, hormone_therapy, 
         radio_dose_per_fraction_Gy, radio_rectum_v60_pc, radio_rectum_v75_pc)

# Verificar las categorías en la columna 'Site'
table(prostate_covars$Site)

# Mapear los nombres de los sitios a números
site_map <- c(
  "Candiolo" = 1, "Freiburg" = 2, "Gent" = 3, "Karlsruhe ST" = 4, "Karlsruhe VN" = 5, 
  "Leicester" = 6, "Leuven" = 7, "Ludwigshafen" = 8, "Maastricht" = 9, "Manchester" = 10, 
  "Mannheim" = 11, "Milan" = 12, "Montpellier" = 13, "Mount Sinai" = 14, "Nimes" = 15, 
  "Santiago" = 16
)
prostate_covars$Site <- site_map[prostate_covars$Site]

# Verificar los cambios en la columna 'Site'
table(prostate_covars$Site)

# Verificar las categorías en la columna 'Prostate.Cancer.Stage'
table(prostate_covars$Prostate.Cancer.Stage)

# Mapear las etapas del cáncer de próstata a números
stage_map <- c("Stage I" = 1, "Stage II" = 2, "Stage III" = 3, "Stage IV" = 4)
prostate_covars$Prostate.Cancer.Stage <- stage_map[prostate_covars$Prostate.Cancer.Stage]

# Verificar los cambios en la columna 'Prostate.Cancer.Stage'
table(prostate_covars$Prostate.Cancer.Stage)

# Mostrar el dataframe final
prostate_covars

# Guardar el dataframe en un archivo de texto
write.table(prostate_covars, file = "prostate_covars.txt", sep = ",", row.names = FALSE, quote = FALSE)

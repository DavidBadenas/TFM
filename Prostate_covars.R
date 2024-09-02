# Cargar las librerías necesarias
library(dplyr)
library(readxl)

# Leer los datos
covar <- read.table("covar_mds.txt", header = TRUE)
prostate <- read_excel("imputed_data_prostate 2.xlsx")

# Renombrar la columna 'Sample ID' a 'FID'
prostate <- prostate %>%
  rename(FID = `Sample ID`)

# Unir los datos por la columna 'FID'
prostate_covars <- inner_join(covar, prostate, by = "FID")

# Seleccionar las columnas de interés
prostate_covars <- prostate_covars %>%
  select(FID, IID, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, Site, pre_radio_turp, radical_prostatectomy, radio_pelvic, ace_inhibitor, Prostate.Cancer.Stage, other_collagen_vascular_disease, previous_abdominal_surgery, psa_prediagnostic_biopsy_ng_ml, hormone_therapy, radio_dose_per_fraction_Gy, radio_rectum_v60_pc, radio_rectum_v75_pc)

# Revisar la estructura del dataframe
str(prostate_covars)

# Convertir columnas de tipo character a factor
prostate_covars <- prostate_covars %>%
  mutate_if(is.character, as.factor)

# Revisar nuevamente la estructura del dataframe
str(prostate_covars)

# Guardar el dataframe resultante en un archivo de texto
write.table(prostate_covars, file = "prostate_covars.txt", sep = "\t", row.names = FALSE, quote = FALSE)

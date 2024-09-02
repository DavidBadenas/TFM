# Cargar las librerías necesarias
library(dplyr)
library(readxl)
library(glmnet)

# MAMA
# Leer los datos de mama
covar_breast <- read_excel("imputed_data_breast.xlsx")

# Mostrar la tabla de frecuencia para 'Recurrence' y 'chemo'
table(covar_breast$Recurrence, covar_breast$chemo)

# Eliminar columnas innecesarias
covar_breast <- covar_breast[, !names(covar_breast) %in% c("Subject.Id", "Sample ID", "Site")]

# Convertir 'Recurrence' en factor y preparar las matrices para el modelo
covar_breast$Recurrence <- factor(covar_breast$Recurrence)
X <- as.matrix(covar_breast[, -which(names(covar_breast) == "Recurrence")])
Y <- covar_breast$Recurrence

# Ajustar el modelo LASSO
mlasso <- glmnet(X, Y, standardize = TRUE, alpha = 1, family = "binomial")

# Graficar el modelo LASSO
plot(mlasso)

# Validación cruzada para seleccionar el mejor lambda
set.seed(1234)
cv.lasso <- cv.glmnet(X, Y, standardize = TRUE, family = "binomial")

# Graficar la validación cruzada
plot(cv.lasso)

# Obtener el valor de lambda que minimiza el error
lambda_min <- cv.lasso$lambda.min
print(lambda_min)

# Obtener los coeficientes asociados a lambda_min
coef_min <- coef(mlasso, s = lambda_min)
print(coef_min)

########################################################

# PROSTATA
# Leer los datos de próstata
covar_prostate <- read_excel("imputed_data_prostate.xlsx")

# Eliminar columnas innecesarias
covar_prostate <- covar_prostate[, !names(covar_prostate) %in% c("Subject.Id", "Sample ID", "Site")]

# Convertir 'Recurrence' en factor y preparar las matrices para el modelo
covar_prostate$Recurrence <- factor(covar_prostate$Recurrence)
X <- as.matrix(covar_prostate[, -which(names(covar_prostate) == "Recurrence")])
Y <- covar_prostate$Recurrence

# Ajustar el modelo LASSO
mlasso <- glmnet(X, Y, standardize = TRUE, alpha = 1, family = "binomial")

# Graficar el modelo LASSO
plot(mlasso)

# Validación cruzada para seleccionar el mejor lambda
set.seed(1234)
cv.lasso <- cv.glmnet(X, Y, standardize = TRUE, family = "binomial")

# Graficar la validación cruzada
plot(cv.lasso)

# Obtener el valor de lambda que minimiza el error
lambda_min <- cv.lasso$lambda.min
print(lambda_min)

# Obtener los coeficientes asociados a lambda_min
coef_min <- coef(mlasso, s = lambda_min)
print(coef_min)

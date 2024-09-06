## FINAL MASTER PROJECT


1º file_conversion.sh --> Script para transformar los datos imputados de formato txt a formato vcf.
2º eliminar_pacientes.sh --> Script para seleccionar los pacientes con fenotipos conocidos.
2º seleccionar_snps_IC<0.3.R --> Script para seleccionar solo aquellos SNPs con un imputation score (IC) inferior a 0.3.
3º eli_snps.sh --> Script para eliminar los SNPs seleccionados en el paso anterios (IC menos a 0.3).
4º QC_breast.sh --> Script que realiza el control de calidad en los pacientes del estudio de cáncer de mama.
5º QC_prostate.sh --> Script que realiza el control de calidad en los pacientes del estudio de cáncer de próstata.
6º LASSO.R --> Script para seleccionar aquellas variables que se van a usar en el GWAS como covariables.
7º Breast_covars.R --> Script para unir las covariables seleccionas por LASSO y las dictadas por oncólogos especialistas y codificar las variables categóricas, en el caso del estudio de pacientes con cáncer de mama.
8º Prostate_covars.R --> Script para unir las covariables seleccionas por LASSO y las dictadas por oncólogos especialistas y codificar las variables categóricas, en el caso del estudio de pacientes con cáncer de próstata.
9º GWAS_breast.sh --> Script para la realización del GWAS con un modelos de regresión logística con las covariables obtenidas en el paso 7, para los pacientes del estudio de mama.
10º GWAS_prostate.sh -->Script para la realización del GWAS con un modelos de regresión logística con las covariables obtenidas en el paso 8, para los pacientes del estudio de prostata.
11º Genes_breast.R --> Script para mirar que genes están asociados con las mutaciones con mayor significancia en el estudio de cáncer de mama.
12º Genes_prostate.R --> Script para mirar que genes están asociados con las mutaciones con mayor significancia en el estudio de cáncer de prostata.
13º PRS_breast.sh --> Script que realiza el PRS en los pacientes del estudio de mama.
14º PRS_prostate.sh --> Scriprt que realiza el PRS en los pacientes del estudio de prostata.

El resto de los script se utilizan el los scripts del Quality control y del GWAS.

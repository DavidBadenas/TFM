# TFM

1. file_conversion.sh --> Script para convertir los datos imputados de formato txt a formato vcf
2. eliminar_pacientes.sh --> Script de R para eliminar aquellos pacientes de los que no sepamos su fenotipo (Cambiar breast por prostate).
3. seleccionar_snps_IC<0.3.R  --> Script de R para seleccionar todos aquellos SNPs con un imputation score inferior a 0.3.
4. eli_snps.sh --> Script para eliminar los SNPs obtenidos en el paso anterior de los datos tanto de mama como de prostata.
5. QC_breast --> Script del control de calidad para los pacientes de cancer de mama. (El script es el mismo que para prostata pero cambiando el nombre de las variables).
6. --> Script que realiza un LASSO para la detección de las covariables, añade aquellas dictadas por oncólogos especialistas y los 10 PCAs obtenidas en el paso de la estratificación poblacional del control de calidad.
7. GWAS_breast --> 
8.  

#!/bin/bash

#SBATCH --job-name= GWAS # Job name
#SBATCH -p highmem
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=davidbadenas@vhio.net     # Where to send mail
#SBATCH --ntasks=1
#SBATCH --mem= 120G                     # Job memory request
#SBATCH --cpus-per-task=2
#SBATCH --output=%x_%j.log   # Standard output and error log
#SBATCH --error=%x_%j.err

singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/GWAS_breast_con_la_mitad/ plink19.sif /plink --bfile breast_GWAS --covar covar_prostate.txt --logistic --hide-covar --pheno tabla_recurrencias.txt --out logistic_results
awk '!/'NA'/' logistic_results_prostate.assoc.logistic > logistic_results_prostate_2.assoc.logistic


singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/GWAS_prostate_todo/ r-immunedeconv-qqman-optparse.sif Rscript Manhattan_plot.R

singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/GWAS_prostate_todo/ r-immunedeconv-qqman-optparse.sif Rscript QQ_plot.R

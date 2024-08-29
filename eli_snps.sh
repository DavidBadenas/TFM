#!/bin/bash
#SBATCH --job-name=eli_snp  # Job name
#SBATCH -p highmem
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=davidbadenas@vhio.net     # Where to send mail
#SBATCH --ntasks=1
#SBATCH --mem=120G                     # Job memory request
#SBATCH --cpus-per-task=1
#SBATCH --output=%x_%j.log   # Standard output and error log
#SBATCH --error=%x_%j.err

singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/previo_al_QC/ bcftools.sif bcftools view -T ^eli_snps.txt imputed_breast.vcf.gz -o imputed_breast_alta_calidad.vcf.gz

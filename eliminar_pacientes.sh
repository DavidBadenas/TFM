#!/bin/bash
#SBATCH --job-name=eliminar  # Job name
#SBATCH -p highmem
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=davidbadenas@vhio.net     # Where to send mail
#SBATCH --ntasks=1
#SBATCH --mem=120G                     # Job memory request
#SBATCH --cpus-per-task=1
#SBATCH --output=%x_%j.log   # Standard output and error log
#SBATCH --error=%x_%j.err

singularity run -H $PWD:/home/ bcftools.sif bcftools view --samples-file pacientes_breast.txt imputed.vcf.gz -o imputed_breast.vcf.gz
singularity run -H $PWD:/home/ bcftools.sif bcftools view --samples-file pacientes_prostate.txt imputed.vcf.gz -o imputed_prostate.vcf.gz
echo "Finish"

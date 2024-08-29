#!/bin/bash

#SBATCH --job-name=test_impute2vcf_chr2  # Job name
#SBATCH -p short
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=albamas@vhio.net     # Where to send mail
#SBATCH --ntasks=1
#SBATCH --mem=5G                     # Job memory request
#SBATCH --cpus-per-task=2
#SBATCH --output=%x_%j.log   # Standard output and error log
#SBATCH --error=%x_%j.err


#Test with chr2

chr="2"

#Convert imputation dosages to genotypes
grep -v "misc snp pos" requite_ids.dose.oncoarray_requite_chr"$chr"_unphased.txt | cut -d " " -f6- | awk '{FS=" "; OFS="\t"}{for (i=1; i <=NF;++i) if ($i < 0.5) {$i = "0/0"} else if ($i >= 0.5 && $i < 1.5) {$i = "0/1"} else {$i = "1/1"}; print $0}' > genotypes_chr"$chr".txt

#Add header with sample names to genotypes table
head -1 requite_ids.dose.oncoarray_requite_chr"$chr"_unphased.txt | cut -d " " -f6- | sed 's/ /\t/g' > sample_names_header.txt
cat sample_names_header.txt genotypes_chr"$chr".txt > genotypes_chr"$chr"_header.txt

#Create table with first 9 columns of vcf file
grep -v "misc snp pos" requite_ids.dose.oncoarray_requite_chr"$chr"_unphased.txt | awk -v chr="$chr" 'BEGIN{FS=OFS=" "}{CHROM=chr; QUAL="."; FILTER=".";INFO="."; FORMAT="GT"; print CHROM,$3,$2,$4,$5,QUAL,FILTER,INFO,FORMAT}' > first_columns_chr"$chr".txt

#Add header to first columns
echo "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT" > header_vcf.txt
cat header_vcf.txt first_columns_chr"$chr".txt | sed 's/ /\t/g' > first_columns_chr"$chr"_header.txt


#Create vcf file
echo -e "##fileformat=VCFv4.2" > chr"$chr".vcf
paste first_columns_chr"$chr"_header.txt genotypes_chr"$chr"_header.txt >> chr"$chr".vcf


rm sample_names_header.txt genotypes_chr"$chr".txt genotypes_chr"$chr"_header.txt first_columns_chr"$chr".txt first_columns_chr"$chr"_header.txt header_vcf.txt
mv chr"$chr".vcf results/

#!/bin/bash

#SBATCH --job-name= QC # Job name
#SBATCH -p highmem
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=davidbadenas@vhio.net     # Where to send mail
#SBATCH --ntasks=1
#SBATCH --mem=120G                     # Job memory request
#SBATCH --cpus-per-task=2
#SBATCH --output=%x_%j.log   # Standard output and error log
#SBATCH --error=%x_%j.err

cut -d ',' -f 2 prostate_GWAS_recurrence_updated.csv > pacientes_prostate.txt

singularity run -H $PWD:/home/ plink19.sif /plink --vcf imputed_prostate.vcf.gz --make-bed --out imputed_prostate
singularity run -H $PWD:/home/ plink19.sif /plink --bfile imputed_prostate --exclude snp_ids.txt --make-bed --out imputed_prostate_new

awk '{$5=1; print}' imputed_prostate_new.fam > archivo_sin_snps_hombres.fam
rm imputed_prostate.*
mv archivo_sin_snps_hombres.fam imputed_prostate.fam
mv imputed_prostate_new.bed imputed_prostate.bed
mv imputed_prostate_new.bim imputed_prostate.bim

# STEP 1
# Investigate missingness per individual and per SNP and make histograms.
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate --missing 
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ r-immunedeconv-qqman-optparse.sif Rscript hist_miss.R

# Delete SNPs with missingness >0.2.
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate --geno 0.2 --make-bed --out imputed_prostate_2
# Delete individuals with missingness >0.02.
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_2 --mind 0.2 --make-bed --out imputed_prostate_3
# Delete SNPs with missingness >0.02.
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_3 --geno 0.02 --make-bed --out imputed_prostate_4
# Delete individuals with missingness >0.02.
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_4 --mind 0.02 --make-bed --out imputed_prostate_5


# STEP 2

singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_5 --check-sex 
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ r-immunedeconv-qqman-optparse.sif Rscript gender_check.R
grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_5 --remove sex_discrepancy.txt --make-bed --out imputed_prostate_6 

# STEP 3
# Generate a bfile with autosomal SNPs only and delete SNPs with a low minor allele frequency (MAF).
# Select autosomal SNPs only
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' imputed_prostate_6.bim > snp_1_22.txt
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_6 --extract snp_1_22.txt --make-bed --out imputed_prostate_7

# Generate a plot of the MAF distribution.
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_7 --freq --out MAF_check
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ r-immunedeconv-qqman-optparse.sif Rscript MAF_check.R

# Remove SNPs with a low MAF frequency. MIRAAR
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_7 --maf 0.05 --make-bed --out imputed_prostate_8


# STEP 4
# Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE).

singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_8 --hardy
awk '{ if ($9 <0.00001) print $0 }' plink.hwe>plinkzoomhwe.hwe
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ r-immunedeconv-qqman-optparse.sif Rscript hwe.R
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_8 --hwe 1e-10 --hwe-all --make-bed --out imputed_prostate_9


# STEP 5
# Generate a plot of the distribution of the heterozygosity rate of your subjects.

singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_9 --exclude inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_9 --extract indepSNP.prune.in --het --out R_check
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ r-immunedeconv-qqman-optparse.sif Rscript check_heterozygosity_rate.R
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ r-immunedeconv-qqman-optparse.sif Rscript heterozygosity_outliers_list.R

sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt

singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_9 --remove het_fail_ind.txt --make-bed --out imputed_prostate_10


# STEP 6

singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_10 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2
awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome 
# singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ r-immunedeconv-qqman-optparse.sif Rscript Relatedness.R
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_10 --remove pacientes_UN.txt --make-bed --out imputed_prostate_11   #pacientes_UN hacerlo a mano

# STEP 7

singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --vcf ALL.2of4intersection.20100804.genotypes.vcf.gz --make-bed --out ALL.2of4intersection.20100804.genotypes
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile ALL.2of4intersection.20100804.genotypes --set-missing-var-ids @:#[b37]\$1,\$2 --make-bed --out ALL.2of4intersection.20100804.genotypes_no_missing_IDs

singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile ALL.2of4intersection.20100804.genotypes_no_missing_IDs --geno 0.2 --allow-no-sex --make-bed --out 1kG_MDS
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile 1kG_MDS --mind 0.2 --allow-no-sex --make-bed --out 1kG_MDS2
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile 1kG_MDS2 --geno 0.02 --allow-no-sex --make-bed --out 1kG_MDS3
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile 1kG_MDS3 --mind 0.02 --allow-no-sex --make-bed --out 1kG_MDS4
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile 1kG_MDS4 --maf 0.05 --allow-no-sex --make-bed --out 1kG_MDS5


awk '{print$2}' imputed_prostate_11.bim > HapMap_SNPs.txt
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile 1kG_MDS5 --extract HapMap_SNPs.txt --make-bed --out 1kG_MDS6
awk '{print$2}' 1kG_MDS6.bim > 1kG_MDS6_SNPs.txt
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_11 --extract 1kG_MDS6_SNPs.txt --recode --make-bed --out imputed_prostate_MDS

awk '{print$2,$4}' imputed_prostate_MDS.map > buildhapmap.txt
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile 1kG_MDS6 --update-map buildhapmap.txt --make-bed --out 1kG_MDS7
awk '{print$2,$5}' 1kG_MDS7.bim > 1kg_ref-list.txt
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_MDS --reference-allele 1kg_ref-list.txt --make-bed --out imputed_prostate-adj

awk '{print$2,$5,$6}' 1kG_MDS7.bim > 1kGMDS7_tmp
awk '{print$2,$5,$6}' imputed_prostate-adj.bim > imputed_prostate-adj_tmp
sort 1kGMDS7_tmp imputed_prostate-adj_tmp |uniq -u > all_differences.txt

awk '{print$1}' all_differences.txt | sort -u > flip_list.txt
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate-adj --flip flip_list.txt --reference-allele 1kg_ref-list.txt --make-bed --out corrected_imputed_prostate

awk '{print$2,$5,$6}' corrected_imputed_prostate.bim > corrected_imputed_prostate_tmp
sort 1kGMDS7_tmp corrected_imputed_prostate_tmp |uniq -u  > uncorresponding_SNPs.txt
awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exlusion.txt

singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile corrected_imputed_prostate --exclude SNPs_for_exlusion.txt --make-bed --out imputed_prostate_MDS2
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile 1kG_MDS7 --exclude SNPs_for_exlusion.txt --make-bed --out 1kG_MDS8
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_MDS2 --bmerge 1kG_MDS8.bed 1kG_MDS8.bim 1kG_MDS8.fam --allow-no-sex --make-bed --out MDS_merge2
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile MDS_merge2 --extract indepSNP.prune.in --genome --out MDS_merge2
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile MDS_merge2 --read-genome MDS_merge2.genome --cluster --mds-plot 10 --out MDS_merge2

awk '{print$1,$1,$2}' 20100804.ALL.panel > race_1kG.txt
sed 's/JPT/ASN/g' race_1kG.txt>race_1kG2.txt
sed 's/ASW/AFR/g' race_1kG2.txt>race_1kG3.txt
sed 's/CEU/EUR/g' race_1kG3.txt>race_1kG4.txt
sed 's/CHB/ASN/g' race_1kG4.txt>race_1kG5.txt
sed 's/CHD/ASN/g' race_1kG5.txt>race_1kG6.txt
sed 's/YRI/AFR/g' race_1kG6.txt>race_1kG7.txt
sed 's/LWK/AFR/g' race_1kG7.txt>race_1kG8.txt
sed 's/TSI/EUR/g' race_1kG8.txt>race_1kG9.txt
sed 's/MXL/AMR/g' race_1kG9.txt>race_1kG10.txt
sed 's/GBR/EUR/g' race_1kG10.txt>race_1kG11.txt
sed 's/FIN/EUR/g' race_1kG11.txt>race_1kG12.txt
sed 's/CHS/ASN/g' race_1kG12.txt>race_1kG13.txt
sed 's/PUR/AMR/g' race_1kG13.txt>race_1kG14.txt

awk '{print$1,$2,"OWN"}' imputed_prostate_MDS.fam>racefile_own.txt
cat race_1kG14.txt racefile_own.txt | sed -e '1i\FID IID race' > racefile.txt
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ r-immunedeconv-qqman-optparse.sif Rscript MDS_merged.R
awk '{ if ($4 < -0.005 && $5 <0.03) print $1,$2 }' MDS_merge2.mds > EUR_MDS_merge2
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_11 --keep EUR_MDS_merge2 --make-bed --out imputed_prostate_12


singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_12 --extract indepSNP.prune.in --genome --out imputed_prostate_12
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_12 --read-genome imputed_prostate_12.genome --cluster --mds-plot 10 --out imputed_prostate_12_mds
awk '{print$1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' imputed_prostate_12_mds.mds > covar_mds.txt


# STEP 8

awk '{print$2}' imputed_prostate_12.bim > HapMap_SNPs_sub.txt
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile 1kG_MDS5 --extract HapMap_SNPs_sub.txt --make-bed --out 1kG_MDS6_sub
awk '{print$2}' 1kG_MDS6_sub.bim > 1kG_MDS6_sub_SNPs.txt
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_12 --extract 1kG_MDS6_sub_SNPs.txt --recode --make-bed --out imputed_prostate_sub_MDS

awk '{print$2,$4}' imputed_prostate_sub_MDS.map > buildhapmap_sub.txt
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile 1kG_MDS6_sub --update-map buildhapmap_sub.txt --make-bed --out 1kG_MDS7_sub
awk '{print$2,$5}' 1kG_MDS7_sub.bim > 1kg_ref-list_sub.txt
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_sub_MDS --reference-allele 1kg_ref-list_sub.txt --make-bed --out imputed_prostate_sub-adj

awk '{print$2,$5,$6}' 1kG_MDS7_sub.bim > 1kGMDS7_sub_tmp
awk '{print$2,$5,$6}' imputed_prostate_sub-adj.bim > imputed_prostate_sub-adj_tmp
sort 1kGMDS7_sub_tmp imputed_prostate_sub-adj_tmp |uniq -u > all_differences_sub.txt

awk '{print$1}' all_differences_sub.txt | sort -u > flip_list_sub.txt
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_sub-adj --flip flip_list_sub.txt --reference-allele 1kg_ref-list_sub.txt --make-bed --out corrected_imputed_prostate_sub

awk '{print$2,$5,$6}' corrected_imputed_prostate_sub.bim > corrected_imputed_prostate_tmp_sub
sort 1kGMDS7_sub_tmp corrected_imputed_prostate_tmp_sub |uniq -u  > uncorresponding_SNPs_sub.txt
awk '{print$1}' uncorresponding_SNPs_sub.txt | sort -u > SNPs_for_exlusion_sub.txt

singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile corrected_imputed_prostate_sub --exclude SNPs_for_exlusion_sub.txt --make-bed --out imputed_prostate_MDS2_sub
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile 1kG_MDS7_sub --exclude SNPs_for_exlusion_sub.txt --make-bed --out 1kG_MDS8_sub
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_MDS2_sub --bmerge 1kG_MDS8_sub.bed 1kG_MDS8_sub.bim 1kG_MDS8_sub.fam --allow-no-sex --make-bed --out MDS_merge2_sub
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile MDS_merge2_sub --extract indepSNP.prune.in --genome --out MDS_merge2_sub
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile MDS_merge2_sub --read-genome MDS_merge2_sub.genome --cluster --mds-plot 10 --out MDS_merge2_sub

awk '{print$1,$1,$2}' 20100804.ALL.panel > race_1kG.txt
sed 's/JPT/ASN/g' race_1kG.txt>race_1kG2_sub.txt
sed 's/ASW/AFR/g' race_1kG2_sub.txt>race_1kG3_sub.txt
sed 's/CHB/ASN/g' race_1kG3_sub.txt>race_1kG4_sub.txt
sed 's/CHD/ASN/g' race_1kG4_sub.txt>race_1kG5_sub.txt
sed 's/YRI/AFR/g' race_1kG5_sub.txt>race_1kG6_sub.txt
sed 's/LWK/AFR/g' race_1kG6_sub.txt>race_1kG7_sub.txt
sed 's/MXL/AMR/g' race_1kG7_sub.txt>race_1kG8_sub.txt
sed 's/CHS/ASN/g' race_1kG8_sub.txt>race_1kG9_sub.txt
sed 's/PUR/AMR/g' race_1kG9_sub.txt>race_1kG10_sub.txt

awk '{print$1,$2,"OWN"}' imputed_prostate_MDS_sub.fam>racefile_own_sub.txt
cat race_1kG10_sub.txt racefile_own_sub.txt | sed -e '1i\FID IID race' > racefile_sub.txt
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ r-immunedeconv-qqman-optparse.sif Rscript MDS_merged_sub.R
awk '{ if ($4 < -0.005 && $5 <0.03) print $1,$2 }' MDS_merge2_sub.mds > EUR_MDS_merge2_sub
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_12 --keep EUR_MDS_merge2 --make-bed --out imputed_prostate_13


singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_13 --extract indepSNP.prune.in --genome --out imputed_prostate_13
singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/QC/ plink19.sif /plink --bfile imputed_prostate_13 --read-genome imputed_prostate_13.genome --cluster --mds-plot 10 --out imputed_prostate_13_mds
awk '{print$1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' imputed_prostate_13_mds.mds > covar_mds_sub.txt

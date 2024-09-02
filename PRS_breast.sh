dat <- read.table("logistic_results_breast.assoc_2.logistic", header=T)
dat$BETA <- log(dat$OR)
write.table(dat, "Breast.Transformed", quote=F, row.names=F)
q()

singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/GWAS_breast_con_la_mitad/ plink19.sif /plink \
    --bfile imputed_breast_13 \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump Breast.Transformed \
    --clump-snp-field SNP \
    --clump-field P \
    --out breast
	
	
awk 'NR!=1{print $3}' breast.clumped >  breast.valid.snp
awk '{print $2,$9}' logistic_results_breast.assoc_2.logistic > SNP.pvalue

echo "0.000001 0 0.000001" > range_list
echo "0.00000215 0 0.00000215" >> range_list
echo "0.00000464 0 0.00000464" >> range_list
echo "0.00001 0 0.00001" >> range_list
echo "0.0000215 0 0.0000215" >> range_list 
echo "0.0000464 0 0.0000464" >> range_list
echo "0.0001 0 0.0001" >> range_list


singularity run -H $PWD:/mnt/bioinfnas/bioinformatics/test/david/GWAS_breast_con_la_mitad/ plink19.sif /plink \
    --bfile imputed_breast_13 \
    --score Breast.Transformed 2 4 10 header \
    --q-score-range range_list SNP.pvalue \
    --extract breast.valid.snp \
    --out breast
	

	
	
En R	
p.threshold <- c(0.000001, 0.00000215, 0.00000464, 0.00001, 0.0000215, 0.0000464, 0.0001)
phenotype <- read.table("breast_phenotype.txt", header=T)
covariate <- read.table("EUR.cov", header=T)
pheno <- merge(phenotype, covariate, by=c("FID", "IID"))
null.r2 <- summary(lm(newly_diagnosed_recurrence_metastasis_cancer~., data=pheno[,!colnames(pheno)%in%c("FID","IID")]))$r.squared
prs.result <- NULL


for(i in p.threshold){
    pheno.prs <- merge(pheno, 
                        read.table(paste0("breast.",i,".profile"), header=T)[,c("FID","IID", "SCORE")],
                        by=c("FID", "IID"))
    model <- summary(lm(newly_diagnosed_recurrence_metastasis_cancer~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")]))
    model.r2 <- model$r.squared
    prs.r2 <- model.r2-null.r2
    prs.coef <- model$coeff["SCORE",]
    prs.result <- rbind(prs.result, 
        data.frame(Threshold=i, R2=prs.r2, 
                    P=as.numeric(prs.coef[4]), 
                    BETA=as.numeric(prs.coef[1]),
                    SE=as.numeric(prs.coef[2])))
}
print(prs.result[which.max(prs.result$R2),])
q() 
	
	
VISUALIZATION

library(ggplot2)
prs.result$print.p <- round(prs.result$P, digits = 3)
prs.result$print.p[!is.na(prs.result$print.p) &
                    prs.result$print.p == 0] <-
    format(prs.result$P[!is.na(prs.result$print.p) &
                            prs.result$print.p == 0], digits = 2)
prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)
# Initialize ggplot, requiring the threshold as the x-axis (use factor so that it is uniformly distributed)
ggplot(data = prs.result, aes(x = factor(Threshold), y = R2)) +
    # Specify that we want to print p-value on top of the bars
    geom_text(
        aes(label = paste(print.p)),
        vjust = -1.5,
        hjust = 0,
        angle = 45,
        cex = 4,
        parse = T
    )  +
    # Specify the range of the plot, *1.25 to provide enough space for the p-values
    scale_y_continuous(limits = c(0, max(prs.result$R2) * 1.25)) +
    # Specify the axis labels
    xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
    ylab(expression(paste("PRS model fit:  ", R ^ 2))) +
    # Draw a bar plot
    geom_bar(aes(fill = -log10(P)), stat = "identity") +
    # Specify the colors
    scale_fill_gradient2(
        low = "dodgerblue",
        high = "firebrick",
        mid = "dodgerblue",
        midpoint = 1e-4,
        name = bquote(atop(-log[10] ~ model, italic(P) - value),)
    ) +
    # Some beautification of the plot
    theme_classic() + theme(
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size =
                                        18),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust =
                                    1)
    )
# save the plot
ggsave("breast.height.bar.png", height = 7, width = 7)
q() # exit R

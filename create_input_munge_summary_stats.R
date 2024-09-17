#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
df_file <- args[1]
pheno <- as.character(args[2])

df <- fread(df_file,header=F)

#CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P
#REGENIE: reference allele (allele 0), alternative allele (allele 1)
#REGENIE:  estimated effect sizes (for allele 1 on the original scale)

colnames(df) <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N","TEST", "BETA", "SE", "CHISQ", "LOG10P")
dfclean <- df %>% select(CHROM,ID,GENPOS,ALLELE1,ALLELE0,A1FREQ,LOG10P,BETA,SE)

dfclean$pval <- 10^(-dfclean$LOG10P)

#select column in order:
dfclean <- dfclean %>% select(ID,CHROM,GENPOS,ALLELE1,ALLELE0,BETA,SE,A1FREQ,pval)

#rename columns :
colnames(dfclean) <- c("snpid", "b37chr", "bp", "a1", "a2", "LOG_ODDS", "se", "eaf", "pval")
#save gwas all chr results:

#write.table(dfclean,paste0("output/",pheno,"_betase_input_mungestat"),quote=FALSE,sep="\t",row.names=FALSE)


#make sure eaf column is numeric:
dfclean$eaf <- as.numeric(dfclean$eaf)

#create MAF columns and filtered out for MAF >= 0.01:
dfclean <- dfclean %>% mutate(MAF = ifelse(dfclean$eaf <= 0.5, dfclean$eaf, 1 - dfclean$eaf))
dfclean_MAF001 <- dfclean %>% filter(MAF >= 0.01)

#save gaws results for variants with MAF >= 0.01:
write.table(dfclean_MAF001,paste0("output/maf001_",pheno,"_betase_input_mungestat"),quote=FALSE,sep="\t",row.names=FALSE)
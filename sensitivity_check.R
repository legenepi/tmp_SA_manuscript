#!/usr/bin/env Rscript


suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

args = commandArgs(trailingOnly=TRUE)

#sumstat1_file = "/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/Disc_VS_Meta_sumstat_6cohorts.txt"
#sumstat2_file = "output/maf001_broad_pheno_nocomob_for_sensitivity"
#output_file = "output/sensitivity_nocomob_suggestive_discovery.txt"

sumstat1_file = args[1]
sumstat2_file = args[2]
output_file = args[3]

#input:
sumstat1 <- fread(sumstat1_file,header=T)
sumstat2 <- fread(sumstat2_file,header=T)

#create odds_ratio ci_lower ci_upper for no-comorbidity GWAS:
sumstat2$odds_ratio <- exp(sumstat2$LOG_ODDS)
sumstat2$ci_lower <- exp(sumstat2$LOG_ODDS - qnorm(0.975)*sumstat2$se)
sumstat2$ci_upper <- exp(sumstat2$LOG_ODDS + qnorm(0.975)*sumstat2$se)
#create unique column:
sumstat2$'OR_95_CI' <- paste0(format(round(sumstat2$odds_ratio,2),nsmall=2)," [",format(round(sumstat2$ci_lower,2),nsmall=2),"-",format(round(sumstat2$ci_upper,2),nsmall=2),"]")
#direction of effect:
sumstat2 <- sumstat2 %>% mutate(direction_of_effect=ifelse(sumstat2$LOG_ODDS > 0, "+","-"))

#check:
sumstat2 <- sumstat2 %>% rename(SNP_UKB=snpid, chr=b37chr, pos_b37=bp, EA=a1, NEA=a2)
colnames(sumstat2)[6:15] <- paste(colnames(sumstat2)[6:15], "no_comob", sep = "_")
check <- left_join(sumstat1, sumstat2, by = c("SNP_UKB","chr","pos_b37","EA","NEA"))
check$concordance_nocomob <- ifelse((check$odds_ratio > 0 & check$odds_ratio_no_comob > 0) | (check$odds_ratio < 0 & check$odds_ratio_no_comob < 0), "yes", "no")

write.table(check,output_file,sep="\t",quote=F,row.names=F)

png("output/sensitivity_nocomob_check_oddsratio.png")
plot(check$odds_ratio, check$odds_ratio_no_comob)
dev.off()
png("output/sensitivity_nocomob_check_pval.png")
plot(check$pval, check$pval_no_comob, xlim= c(0,0.0003), ylim=c(0,0.0003))
dev.off()
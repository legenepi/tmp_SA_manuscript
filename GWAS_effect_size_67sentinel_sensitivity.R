#!/usr/bin/env Rscript


TO BE FINISHED

#Rationale: compare effect size of the 67 sentinel variants between the discovery GWAS and the sensitivity analysis
#allergy-free.

#On bash:
#ungzipped the summary stats for noallergy GWAS:
#tar -xvzf REGENIE_assoc_31082023.tar.gz \
#home/n/nnp5/PhD/PhD_project//REGENIE_assoc/output/maf001_noallergy_pheno_betase_input_mungestat
#cp in Rdrive:
#cp home/n/nnp5/PhD/PhD_project/REGENIE_assoc/output/maf001_noallergy_pheno_betase_input_mungestat \
#/rfs/TobinGroup/nnp5/data/
#filter for sentinel variants:
#dos2unix 67sentinelvars_rsid.txt
#remove the duplicate row for rs762279794:
grep -w -F -f 67sentinelvars_rsid.txt /rfs/TobinGroup/nnp5/data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat | \
    grep -v "CAA" \
    > /scratch/gen1/nnp5/tmp_manuscript/67sentvars_gwas
grep -w -F -f 67sentinelvars_rsid.txt /rfs/TobinGroup/nnp5/data/maf001_noallergy_pheno_betase_input_mungestat | \
    grep -v "CAA" \
    > /scratch/gen1/nnp5/tmp_manuscript/67sentvars_noallergy

#In R:
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
library(readxl)
library(cowplot) # it allows you to save figures in .png file
library(smplot2)
library(ggpubr)

#input file:
gwas <- fread("/scratch/gen1/nnp5/tmp_manuscript/67sentvars_gwas") %>% select(V1, V6)
colnames(gwas) <- c("snpid", "logodds_gwas")
noallergy <- fread("/scratch/gen1/nnp5/tmp_manuscript/67sentvars_noallergy") %>% select(V1, V6)
colnames(noallergy) <- c("snpid", "logodds_noallergy")
df_plot <- gwas %>% left_join(noallergy, by = "snpid")



png("output/67sentinels_DiscVSSensitallergy_effectsize_comparison.png",units="in", width=10, height=10, res=800)
ggplot(data = df_plot, aes(x = logodds_noallergy, y = logodds_gwas)) +
  sm_statCorr(color = "black", corr_method = "pearson", text_size = 3, size = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey", size = 0.25, linetype = "dashed") +
  geom_smooth(method="lm") +
  geom_point(data = df_plot, aes(x = logodds_noallergy,color="#D55E00",size=0.15)) +
  theme_minimal() + geom_vline(xintercept = 0, linetype="dashed", color = "grey", size = 0.25) +
  geom_hline(yintercept = 0, linetype="dashed", color="grey", size = 0.5) +
  ylim(-0.6, +0.4) + xlim(-0.6, + 0.4) +
  xlab("Sensitivity no allergy") + ylab("Discovery")
dev.off()

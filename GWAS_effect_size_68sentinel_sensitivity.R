#!/usr/bin/env Rscript
#Rationale: compare effect size of the 68 sentinel variants between the discovery GWAS and the sensitivity analysis
#allergy-free and comorbidities-free. - update 17/07/2025 wit the sugg sentinel on chr1 that was missed.

#On bash:
#ungzipped the summary stats for noallergy GWAS:
#tar -xvzf REGENIE_assoc_31082023.tar.gz \
#home/n/nnp5/PhD/PhD_project//REGENIE_assoc/output/maf001_noallergy_pheno_betase_input_mungestat
#cp in Rdrive:
#cp home/n/nnp5/PhD/PhD_project/REGENIE_assoc/output/maf001_noallergy_pheno_betase_input_mungestat \
#/rfs/TobinGroup/nnp5/data/
#R:\TobinGroup\nnp5\PhD_DocMiscell_Bridge\sentinel_vars_severeasthma.xlsx copy and paste vars into 68sentvars_gwas
#filter for sentinel variants:
#dos2unix /rfs/TobinGroup/nnp5/data/sugg_68_sentinels.txt
#remove the duplicate row for rs762279794:
#grep -w -F -f /rfs/TobinGroup/nnp5/data/sugg_68_sentinels.txt /rfs/TobinGroup/nnp5/data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat  | \
#    grep -v "CAA" \
#    > /scratch/gen1/nnp5/tmp_manuscript/68sentvars_gwas

#grep -w -F -f /rfs/TobinGroup/nnp5/data/sugg_68_sentinels.txt /rfs/TobinGroup/nnp5/data/maf001_noallergy_pheno_betase_input_mungestat  | \
#    grep -v "CAA" \
#    > /scratch/gen1/nnp5/tmp_manuscript/68sentvars_noallergy

grep -w -F -f /rfs/TobinGroup/nnp5/data/sugg_68_sentinels.txt /rfs/TobinGroup/nnp5/data/maf001_broad_pheno_nocomob_betase_input_mungestat | \
    grep -v "CAA" \
    > /scratch/gen1/nnp5/tmp_manuscript/68sentvars_nocomorb


#In R:
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
library(readxl)
library(cowplot) # it allows you to save figures in .png file
library(smplot2)
library(ggpubr)
library("ggrepel")

#input file:
gwas <- fread("/scratch/gen1/nnp5/tmp_manuscript/68sentvars_gwas") %>% select(V1, V6)
colnames(gwas) <- c("snpid", "logodds_gwas")

##no-allergy:
noallergy <- fread("/scratch/gen1/nnp5/tmp_manuscript/68sentvars_noallergy") %>% select(V1, V6)
colnames(noallergy) <- c("snpid", "logodds_noallergy")
df_plot <- gwas %>% left_join(noallergy, by = "snpid")
df_plot$logOR_diff <- abs(df_plot$logodds_gwas - df_plot$logodds_noallergy)


png("output/68sentinels_DiscVSSensitallergy_effectsize_comparison.png",units="in", width=10, height=10, res=800)
ggplot(data = df_plot, aes(x = logodds_noallergy, y = logodds_gwas)) +
  sm_statCorr(color = "black", corr_method = "pearson", text_size = 3, size = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey", size = 0.25, linetype = "dashed") +
  geom_smooth(method="lm") +
  geom_point(data = df_plot, aes(x = logodds_noallergy,color="#D55E00",size=0.15)) +
  theme_minimal() + geom_vline(xintercept = 0, linetype="dashed", color = "grey", size = 0.25) +
  geom_hline(yintercept = 0, linetype="dashed", color="grey", size = 0.5) +
  ylim(-0.6, +0.4) + xlim(-0.6, + 0.4) +
  xlab("Sensitivity no allergy") + ylab("Discovery") +
  geom_text_repel(
    data = subset(df_plot, logOR_diff >= 0.1 ),
    aes(label = snpid),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))
dev.off()


##no-comorbidities:
nocomorb <- fread("/scratch/gen1/nnp5/tmp_manuscript/68sentvars_nocomorb") %>% select(V1, V6)
colnames(nocomorb) <- c("snpid", "logodds_nocomorb")
df_plot <- gwas %>% left_join(nocomorb, by = "snpid")
df_plot$logOR_diff <- abs(df_plot$logodds_gwas - df_plot$logodds_nocomorb)


png("output/68sentinels_DiscVSSensitcomob_effectsize_comparison.png",units="in", width=10, height=10, res=800)
ggplot(data = df_plot, aes(x = logodds_nocomorb, y = logodds_gwas)) +
  sm_statCorr(color = "black", corr_method = "pearson", text_size = 3, size = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey", size = 0.25, linetype = "dashed") +
  geom_smooth(method="lm") +
  geom_point(data = df_plot, aes(x = logodds_nocomorb,color="#D55E00",size=0.15)) +
  theme_minimal() + geom_vline(xintercept = 0, linetype="dashed", color = "grey", size = 0.25) +
  geom_hline(yintercept = 0, linetype="dashed", color="grey", size = 0.5) +
  ylim(-0.6, +0.4) + xlim(-0.6, + 0.4) +
  xlab("Sensitivity no comorbidities") + ylab("Discovery") +
  geom_text_repel(
    data = subset(df_plot, logOR_diff >= 0.1 ),
    aes(label = snpid),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))
dev.off()


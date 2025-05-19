#!/usr/bin/env Rscript

#Rationale: compare effect size of the 21 replicated variants between the discovery GWAS and the GBMI.

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
library(readxl)
library(cowplot) # it allows you to save figures in .png file
library(smplot2)
library(ggpubr)


#replicated:
repl <- fread("ReplVars.txt")
repl <- repl %>% rename(rsid = 'SNP_UKB', CHR = 'chr', BP = 'pos_b38')

#gbmi:
eur_gbmi <- fread("Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz")
eur_gbmi <- eur_gbmi %>% rename(CHR = '#CHR', BP =   POS)
#15 replicated variants are present in the GBMI, the others need to find through proxy:
repl_gbmi <- repl %>% left_join(eur_gbmi, by = c('CHR','BP')) %>% rename(rsid = 'rsid.x')
#flip beta sign:
#chek ALT == a1 if not flip effect size:
repl_gbmi <- repl_gbmi %>%
             mutate(eur_gbmi_beta_flipped = ifelse(repl_gbmi$EA == repl_gbmi$ALT, inv_var_meta_beta, -inv_var_meta_beta))

#two variants not present, find proxy:
##rs201499805 on chromosome 10: use rs144536148 which is the proxy used from SNP_GUU
proxy <- eur_gbmi %>% filter(CHR == 10, BP == 9001864, rsid == "rs144536148")
#17:38073838_CCG_C chr 17, pos38 is 39917585:
proxy2 <- eur_gbmi %>% filter(CHR == 17, BP > 39917585-200, BP < 39917585+200)
proxy2 <- proxy2[3,]
proxy_all <- rbind(proxy, proxy2) %>% rename(proxy = rsid)
proxy_all$rsid <- c("rs201499805", "17:38073838_CCG_C")

#add manually beta and se for the proxies:
repl_gbmi$eur_gbmi_beta_flipped[c(12,21)] <- proxy_all$inv_var_meta_beta
repl_gbmi$inv_var_meta_sebeta[c(12,21)] <- proxy_all$inv_var_meta_sebeta

#upload severe GWAS summary stats to have beta (logOR):
gwas <- fread("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat")
gwas <- gwas %>% rename(rsid = snpid)
repl_gbmi <- repl_gbmi %>% left_join(gwas, by = 'rsid')

df_plot <- repl_gbmi %>% select(rsid, LOG_ODDS, se, eur_gbmi_beta_flipped, inv_var_meta_sebeta)

png("output/21replvars_DiscVSeurgbmi_effectsize_comparison.png",units="in", width=10, height=10, res=800)
ggplot(data = df_plot, aes(x = eur_gbmi_beta_flipped, y = LOG_ODDS)) +
  sm_statCorr(color = "black", corr_method = "pearson", text_size = 3, size = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey", size = 0.25, linetype = "dashed") +
  geom_smooth(method="lm") +
  geom_point(data = df_plot, aes(x = eur_gbmi_beta_flipped,color="#D55E00",size=0.4)) +
  theme_minimal() + geom_vline(xintercept = 0, linetype="dashed", color = "grey", size = 0.25) +
  geom_hline(yintercept = 0, linetype="dashed", color="grey", size = 0.5) +
  ylim(-0.25, +0.25) + xlim(-0.25, + 0.25) + xlab("GBMI") + ylab("Discovery")
dev.off()

#Welch's t-test: Bonferroni corrected pvalue: 0.05/21
df_plot <- df_plot %>% mutate(eur_z_stat=(eur_gbmi_beta_flipped-LOG_ODDS)/sqrt(inv_var_meta_sebeta^2+se^2),eur_p=2*pnorm(-abs(eur_z_stat)))
df_plot %>% filter(eur_p < (0.05/nrow(df_plot))) #9 variants at Bonferroni p-value
df_plot %>% filter(eur_p < (0.05) #19 variants at nominal p-value
#2 variants that difference in effect size is not significant.
fwrite(df_plot,"output/21replvars_DiscVSeurgbmi_effectsize_welchttest.txt", quote = F)

#plot beta comparison by rsid:
df_plot_gwas <- df_plot %>% select(rsid, LOG_ODDS, se)
colnames(df_plot_gwas) <- c("snpid","beta","se")
df_plot_gwas$study <- as.factor("Discovery")
df_plot_gbmi <- df_plot %>% select(rsid, eur_gbmi_beta_flipped, inv_var_meta_sebeta)
colnames(df_plot_gbmi) <- c("snpid","beta","se")
df_plot_gbmi$study <- as.factor("GBMI")
slope_df <- rbind(df_plot_gwas, df_plot_gbmi)
png(file="output/21replvars_DiscVSeurgbmi_effectsize_pointrange.png", res=600, width=4800, height=4800, pointsize=10)
ggplot(slope_df, aes(colour = study)) +
geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
geom_pointrange(aes(x = snpid, y = beta, ymin = beta-se, ymax = beta+se ), lwd = 1, position = position_dodge(width = 1/2), fatten = .5, size = 1) +
coord_flip() +
theme_bw()
dev.off()
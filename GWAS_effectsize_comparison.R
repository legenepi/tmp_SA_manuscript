#!/usr/bin/env Rscript

#Rationale: compare the 128 EUR-GBMI index variants' effect size against our discovery GWAS effect size.
#179 GBMI: 46 not in EUR, 3 on chr23 which we did not analysed, so a total of 128 in EUR to be compared.

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
library(readxl)
library(cowplot) # it allows you to save figures in .png file
library(smplot2)
library(ggpubr)

##Comparison with sentinel variants as found by GBMI to avoid winner's course from our discovery analysis.
#179 variants from GBMI - remove sentinel on chr23 as we did not analysed it: from TableS2 of GBMI paper PMID 36778051.
#awk '{print $1}' GBMI_sentinel_vars.txt | awk -F ":" '{print $1":"$2"-"$2+1}' |  grep -v "chr23:" | tail -n +2 \
#> GBMI_sentinel_vars_b38_liftover_input
#online liftover:
#awk -F '-' '{print $1}' GBMI_sentinel_vars_b37_liftover_output.bed > GBMI_sentinel_vars_b37_key
#gbmi_loci <- read_excel("GBMI_sentinel_vars.xlsx",sheet = "Sheet1") %>% filter(!grepl("chr23:",SNPID))
#b37_key <- fread("GBMI_sentinel_vars_b37_key",header=F)
#colnames(b37_key) <- "b37"
#gbmi_loci_b37_key <- cbind(gbmi_loci,b37_key)

#upload severe GWAS summary stats:
#gwas <- fread("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/maf001_broad_pheno_1_5_ratio_betase_input_mungestat")
#gwas$b37 <- paste0("chr",gwas$b37chr,":",gwas$bp)
#gbmi_gwas <- left_join(gbmi_loci_b37_key, gwas, by = "b37")
#gbmi_gwas <- gbmi_gwas %>% rename(gwas_beta = LOG_ODDS, gwas_sebeta = se)

#chek ALT == a1 if not flip effect size:
#gbmi_gwas <- gbmi_gwas %>%
#             mutate(gbmi_beta_flipped = ifelse(gbmi_gwas$a1 == gbmi_gwas$ALT, gbmi_beta, -gbmi_beta))

#using effect size from GBMI EUR-only:
#eur_gbmi <- fread("Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz")
#eur_gbmi <- eur_gbmi %>% rename(CHR = '#CHR', BP =   POS)
#eur_gbmi_gwas <- left_join(gbmi_gwas, eur_gbmi, by = c("CHR", "BP", "REF", "ALT"))
#eur_gbmi_gwas <- eur_gbmi_gwas %>% rename(eur_gbmi_beta = "inv_var_meta_beta", eur_gbmi_sebeta = "inv_var_meta_sebeta")

#chek ALT == a1 if not flip effect size:
#eur_gbmi_gwas <- eur_gbmi_gwas %>%
#             mutate(eur_gbmi_beta_flipped = ifelse(eur_gbmi_gwas$a1 == eur_gbmi_gwas$ALT, eur_gbmi_beta, -eur_gbmi_beta))

#Table S8 for the 49 (3 from chr23 that I did not analysed) gwas non-significant in EUR-GBMI:
#eur_notsig <- read_excel("TableS8_EUR_GBMI_notsign.xlsx", sheet = "Sheet1")
#eur_notsig_df <- semi_join(eur_gbmi_gwas, eur_notsig, by = c("CHR", "BP", "REF", "ALT"))
#they are 46, it's right.

#EUR-GBMI significant only and save the file because it has the 128 variants significnat in EUR-GBMI and all columns for plots:
#Table S8 for the 49 (3 from chr23 that I did not analysed) gwas non-significant in EUR-GBMI:
#eur_notsig <- read_excel("TableS8_EUR_GBMI_notsign.xlsx", sheet = "Sheet1")
#eur_gbmi_gwas_onlysig <- anti_join(eur_gbmi_gwas, eur_notsig, by = c("CHR", "BP", "REF", "ALT"))
#Check which variants have opposite effect size:
#eur_gbmi_gwas_onlysig <- eur_gbmi_gwas_onlysig %>% mutate(beta_concordance = ifelse(((eur_gbmi_gwas_onlysig$eur_gbmi_beta_flipped > 0 & eur_gbmi_gwas_onlysig$gwas_beta > 0) | (eur_gbmi_gwas_onlysig$eur_gbmi_beta_flipped < 0 & eur_gbmi_gwas_onlysig$gwas_beta < 0)),"+","-"))
#fwrite(eur_gbmi_gwas_onlysig,"EUR_GBMI_significant.txt", quote = F)
eur_gbmi_gwas_onlysig <- fread("EUR_GBMI_significant.txt")

#Welch's t-test: Bonferroni corrected pvalue: 0.05/174 = 0.0002873563
eur_gbmi_gwas_onlysig <- eur_gbmi_gwas_onlysig %>% mutate(eur_z_stat=(eur_gbmi_beta_flipped-gwas_beta)/sqrt(eur_gbmi_sebeta^2+gwas_sebeta^2),eur_p=2*pnorm(-abs(eur_z_stat)))
eur_sig <- eur_gbmi_gwas_onlysig %>% filter(eur_p < (0.05/nrow(eur_gbmi_gwas_onlysig)))

#Add replicated variant with colorcoded way in the first plot:
#awk -F ';' '{print $1}' /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/meta_analysis_bonferroni_Nsuggestive_replicated_6cohorts > /home/n/nnp5/PhD/PhD_project/tmp_manuscript/rsid_replicated
replicated <- fread("/home/n/nnp5/PhD/PhD_project/tmp_manuscript/rsid_replicated")
colnames(replicated) <- "snpid"
replicated_gwas_gbmi <- eur_gbmi_gwas_onlysig %>% filter(snpid %in% replicated$snpid)
print(replicated_gwas_gbmi$snpid)
fwrite(eur_gbmi_gwas_onlysig,"EUR_GBMI_welchttest.txt", quote = F)

png("output/Eur_gbmi_gwas_effectsize_comparison_onlysign.png",units="in", width=10, height=10, res=800)
ggplot(data = eur_gbmi_gwas_onlysig, aes(x = gwas_beta, y = eur_gbmi_beta_flipped)) +
  sm_statCorr(color = "black", corr_method = "pearson", text_size = 3, size = 0.5) +
  #geom_errorbar(aes(ymin = eur_gbmi_beta_flipped - eur_gbmi_sebeta,ymax = eur_gbmi_beta_flipped + eur_gbmi_sebeta), colour="lightblue", alpha=0.5, size=0.3) +
  #geom_errorbarh(aes(xmin = gwas_beta - gwas_sebeta,xmax = gwas_beta + gwas_sebeta), colour="lightblue", alpha=0.5, size=0.3) +
  geom_point(shape = 21, fill = "#0072B2", color = "#0072B2", size = 0.4) + ylab("EUR-GBMI beta") + xlab("GWAS beta") +
  theme_minimal() + geom_vline(xintercept = 0, linetype="dashed", color = "grey", size = 0.25) +
  geom_hline(yintercept = 0, linetype="dashed", color="grey", size = 0.25) +
  ylim(-0.3, +0.3) + xlim(-0.3, + 0.3) +
  geom_abline(intercept = 0, slope = 1, color = "grey", size = 0.25, linetype = "dashed") +
  geom_smooth(method="lm") +
  geom_point(data=eur_sig,aes(x=gwas_beta,y=eur_gbmi_beta),color="#D55E00",size=0.4) +
  geom_point(data=replicated_gwas_gbmi ,aes(x=gwas_beta,y=eur_gbmi_beta),color="#CC79A7",size=0.4)
dev.off()

#palette.colors(palette = "Okabe-Ito")
#        black        orange       skyblue   bluishgreen        yellow
#    "#000000"     "#E69F00"     "#56B4E9"     "#009E73"     "#F0E442"
#         blue    vermillion reddishpurple          gray
#    "#0072B2"     "#D55E00"     "#CC79A7"     "#999999"

#Slope value with EUR-GBMI, discovery GWAS - for the 174-46=128 significant variants EUR-GBMI:
res <- lm(eur_gbmi_gwas_onlysig$gwas_beta ~ eur_gbmi_gwas_onlysig$eur_gbmi_beta_flipped)
summary(res)
paste0("Slope of EUR-GBMI ~ GWAS replicated: ", round(res$coefficients[[2]],2))
#This means that when EUR-GBMI is equal 1, GWAS is equal to 1*1.34, so 0.34% higher in GWAS.

#Slope plot:
eur_gbmi_slope <- eur_gbmi_gwas_onlysig %>% select(snpid, eur_gbmi_beta_flipped, eur_gbmi_sebeta) %>% rename(beta = eur_gbmi_beta_flipped, se = eur_gbmi_sebeta)
eur_gbmi_slope$study <- as.factor("EUR_GBMI")
gwas_slope <- eur_gbmi_gwas_onlysig %>% select(snpid, gwas_beta, gwas_sebeta) %>% rename(beta = gwas_beta, se = gwas_sebeta)
gwas_slope$study  <- as.factor("GWAS")
slope_df <- rbind(eur_gbmi_slope, gwas_slope)
slope_plot <- ggplot(data = slope_df, aes(x = study, y = beta, fill = study)) +
  sm_slope(labels = c("EUR_GBMI", "GWAS"), group = snpid, line.params = list(linewidth = 0.2), point.params = list(size = 1.5)) +
  scale_fill_manual(values = sm_color("blue", "orange")) +
  geom_hline(yintercept = 0, linetype="dashed", color="grey", size = 0.25)
ggsave("output/Eur_gbmi_gwas_effectsize_slope_onlysig.pdf", plot = last_plot(), width = 20, height = 20, units = "cm")
save_plot("output/Eur_gbmi_gwas_effectsize_slope_onlysign.png", slope_plot, base_asp = 1.4)

#plot beta comparison by rsid:
beta_comparison_repl_pointrange <- ggplot(slope_df, aes(colour = study)) +
geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
geom_pointrange(aes(x = snpid, y = beta, ymin = beta-se, ymax = beta+se ), lwd = 1, position = position_dodge(width = 1/2)) +
coord_flip() +
theme_bw()
ggsave("output/Eur_gbmi_gwas_effectsize_pointrange_onlysign.pdf", plot = last_plot(), width = 20, height = 50, units = "cm")
save_plot("output/Eur_gbmi_gwas_effectsize_pointrange_onlysign.png", beta_comparison_repl_pointrange, base_asp = 1.4)
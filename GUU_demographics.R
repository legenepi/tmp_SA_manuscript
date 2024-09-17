#!/usr/bin/env Rscript

#Rationale: demographic of GASP/U-BIOPRED/UKBiobank case-control cohort used in the replication meta-analysis

NB: U-BIOPRED only adult, so it will be GASP with maybe some childhood onset asthma

library(tidyverse)
library(dplyr)
library(data.table)
library(hablar)

#load input
##demographic table:
#cp /home/n/nnp5/PhD/PhD_project/Post_GWAS/input/ubio_gasp_gwas/pheno_cov_ubio_gasp_ukb.txt \
#    /home/n/nnp5/PhD/PhD_project/tmp_manuscript/
demo <- fread("/home/n/nnp5/PhD/PhD_project/Post_GWAS/input/ubio_gasp_gwas/pheno_cov_ubio_gasp_ukb.txt")
cases <- demo %>% filter(pheno == 1)
controls <- demo %>% filter(pheno == 0)

#Split participants into their cohorts to give IID to Ian Sayers and Ian Adcock for demographic traits in GASP and
U-BIOPRED respectively:

gasp_cases_IID <- demo %>% filter(severe_cohort == "GASP_SEVEREasthma_batch2" | severe_cohort == "GASP_SEVEREasthma") %>%
                  select(FID, IID)
fwrite(gasp_cases_IID,"gasp_cases_IID.txt",sep="\t",quote=F)

ubiopred_cases <- demo %>% filter(severe_cohort == "UBIOPRED_adult_SEVEREasthma_case" | severe_cohort == "UBIOPRED_adult_SEVEREasthma_case_batch2") %>%
                  select(FID, IID)
fwrite(ubiopred_cases, "ubiopred_cases_IID.txt",sep="\t",quote=F)

ubiopred_controls_IID <- demo %>% filter(severe_cohort == "UBIOPRED_adult_asthma_control") %>%
                  select(FID, IID)
fwrite(ubiopred_controls_IID, "ubiopred_controls_IID.txt",sep="\t",quote=F)

#Genetic sex:
# is female, 1 is male
table(cases$sex,exclude=NULL)
prop.table(table(cases$sex,exclude=NULL))

table(controls$sex,exclude=NULL)
prop.table(table(controls$sex,exclude=NULL))

#Age:
print(sum(is.na(cases$age)))
print(mean(cases$age,na.rm=TRUE))
print(sd(cases$age,na.rm=TRUE))

print(sum(is.na(controls$age)))
print(mean(controls$age,na.rm=TRUE))
print(sd(controls$age,na.rm=TRUE))

#UK-Biobank controls data:
ukbb_IID <- demo %>% filter(severe_cohort == "UKBIOBANK_adult_controls") %>%
                  select(IID)
ukbb_IID$IID <- as.factor(ukbb_IID$IID)
eur_ukbb <- read.table("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/demo_EUR_pheno_cov_broadasthma.txt",header=T,sep=" ")
eur_ukbb$IID <- as.factor(eur_ukbb$IID)
eur_ukbb_in_meta <- eur_ukbb %>% inner_join(ukbb_IID, by = c("IID")) %>%
                    select(IID,eid,BMI,smoking_status,category_onset)
#check for cateory onset: it is NAs for all of them because they are controls:
table(eur_ukbb_in_meta$category_onset,exclude=NULL)

##lung funciton for UK-Biobank:
#Lung function from Kath file:
bridge_app648_8389 <- fread("/data/gen1/UKBiobank/application_648/mapping_to_app8389.txt",header=T)
bridge_app648_8389$app8389 <- as.character(bridge_app648_8389$app8389)
bridge_app648_8389$app648 <- as.character(bridge_app648_8389$app648)
#awk '{print $1, $52}' /data/gen1/UKBiobank_500K/severe_asthma/data/ukbiobank_master_app56607.sample \
#    > /data/gen1/UKBiobank_500K/severe_asthma/data/bridge_app648_56607
bridge_app648_56607 <- fread("/data/gen1/UKBiobank_500K/severe_asthma/data/bridge_app648_56607",sep=" ",header=T)
colnames(bridge_app648_56607) <- c("app648","app56607")
bridge_app648_56607$app56607 <- as.character(bridge_app648_56607$app56607)
bridge_app648_56607$app648 <- as.character(bridge_app648_56607$app648)
bridge_648_8389_56607 <- inner_join(bridge_app648_8389,bridge_app648_56607,by="app648")

perc_pred_fev1_kath <- fread("/data/gen1/UKBiobank_500K/severe_asthma/data/percent_pred_fev1.txt",header=T)
fev1_fvc_kath <- fread("/data/gen1/UKBiobank_500K/severe_asthma/data/ff.txt",header=T)
LF_kath <- inner_join(perc_pred_fev1_kath,fev1_fvc_kath,by="ID_1") %>% rename(app648=ID_1)
LF_kath$app648 <- as.character(LF_kath$app648)
LF_kath_app56607 <- left_join(LF_kath, bridge_648_8389_56607, by="app648") %>% select(app56607,fev1_perc_pred,ff.best) %>% rename(eid=app56607)
LF_kath_app56607$eid <- as.factor(LF_kath_app56607$eid)
eur_ukbb_in_meta$eid <- as.factor(eur_ukbb_in_meta$eid)
eur_ukbb_in_meta <- eur_ukbb_in_meta %>% left_join(LF_kath_app56607,by="eid")

#eosinophil:
eos <- fread("Eosinophill_count_30150.csv")
colnames(eos) <- c("eid","eos1","eos2","eso3")
eos$eid <- as.factor(eos$eid)
demo_eos <- left_join(eur_ukbb_in_meta,eos,by="eid") %>% select(eos1,eos2,eso3)
demo_eos <- demo_eos %>%
  rowwise() %>%
  mutate(min_eos = min_(c_across()),
         max_eos = max_(c_across()))


#neutrophil:
neu <- fread("Neutrophill_count_30140.csv")
colnames(neu) <- c("eid","neu1","neu2","neu3")
neu$eid <- as.factor(neu$eid)
demo_neu <- left_join(eur_ukbb_in_meta,neu,by="eid") %>% select(neu1,neu2,neu3)
demo_neu <- demo_neu %>%
  rowwise() %>%
  mutate(min_neu = min_(c_across()),
         max_neu = max_(c_across()))

##add eosinophils and neutrophils to the main table:
eur_ukbb_in_meta_demo <- cbind(eur_ukbb_in_meta,demo_eos$max_eos,demo_neu$max_neu) %>%
                         rename(eos_count = 'demo_eos$max_eos', neu_count = 'demo_neu$max_neu')
fwrite(eur_ukbb_in_meta_demo,"ukbb_ctrl_in_GUU_demotraits.txt",sep="\t",quote=F)
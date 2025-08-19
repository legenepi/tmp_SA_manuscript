#!/usr/bin/env Rscript

#Rationale: create update version of smoking status, age at onset asthma, prednisolone use.
#For smoking status, co-author raised that the number of unknown is odd
#For age at onset, it might be that I used two Data-Field when only one was fine.
#For prednisolone use: strange the difference in unknown. - I think prednisolone use it's fine as it is. It comes from prescription data and not all individuals have primary care records.

library(tidyverse)
library(dplyr)
library(data.table)
library(hablar)

#Smoking status: use the smoking status from the smoking behaviour paper.
#/data/gen1/UKBiobank/Smoking/clean_phenotypes_jan17/ukb648_smoking_behaviour_pheno_covar_270617
smk_status_648 <- fread("/data/gen1/UKBiobank/Smoking/clean_phenotypes_jan17/ukb648_smoking_behaviour_pheno_covar_270617") %>%
                        select(app_id.648,SI) %>% rename(app648 = "app_id.648")
smk_status_648$app648 <- as.character(smk_status_648$app648)

#bridge file, app648 to app56607:
bridge_app648_56607 <- fread("/data/gen1/UKBiobank_500K/severe_asthma/data/bridge_app648_56607",sep=" ",header=T)
colnames(bridge_app648_56607) <- c("app648","app56607")
bridge_app648_56607$app56607 <- as.character(bridge_app648_56607$app56607)
bridge_app648_56607$app648 <- as.character(bridge_app648_56607$app648)
smk_status_56607 <- smk_status_648 %>% left_join(bridge_app648_56607, by="app648") %>% select(-app648)

#Retrieve all-comer asthma European NOT cases (64218 individuals):
##European demographic table:
demo_eur <- fread("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/demo_EUR_pheno_cov_broadasthma.txt",header=T,sep=" ")
##all-comer asthma IDs:
allasthma_ID <- read.table("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/Eid_intersection_asthma_diagnosis_ATLEAST_1_evidence.txt")
##all asthma european demographics table:
demo <- demo_eur %>% filter(IID %in% allasthma_ID$V1)
demo$eid <- as.character(demo$eid)
##withdrawn up to feb22:
#copy from
withdr <- fread("Eid_withdrawn_participants_upFeb2022.txt")
withdr_IID <- withdr$V1
demo <- demo %>% filter(!IID %in% withdr_IID)

#filter out cases to have the asthma not cases:
asthma_notcases_eur <- demo %>% filter(is.na(cases_broad_EUR)) %>% rename(app56607 = "FID")
asthma_notcases_eur$app56607 <- as.character(asthma_notcases_eur$app56607)
#smoking status as SI from smoking behaviour paper:
asthma_notcases_eur <- asthma_notcases_eur %>% left_join(smk_status_56607, by = "app56607")
table(asthma_notcases_eur$SI,exclude=NULL)
print(prop.table(table(asthma_notcases_eur$SI,exclude=NULL)))

#cases:
cases_eur <- demo %>% filter(!is.na(cases_broad_EUR)) %>% rename(app56607 = "FID")
cases_eur$app56607 <- as.character(cases_eur$app56607)
cases_eur <- cases_eur %>% left_join(smk_status_56607, by = "app56607")
table(cases_eur$SI,exclude=NULL)
print(prop.table(table(cases_eur$SI,exclude=NULL)))

#controls:
controls <- demo_eur %>% filter(broad_pheno_1_5_ratio == 0) %>% rename(app56607 = "FID")
controls$app56607 <- as.character(controls$app56607)
controls <- controls %>% left_join(smk_status_56607, by = "app56607")
table(controls$SI,exclude=NULL)
print(prop.table(table(controls$SI,exclude=NULL)))

#smoking status statistics:
nocases_smk <- asthma_notcases_eur %>% select(SI)
nocases_smk$tmp_col <- as.factor(0)
cases_smk <- cases_eur %>% select(SI)
cases_smk$tmp_col <- as.factor(1)
controls_smk <- controls %>% select(SI)
controls_smk$tmp_col <- as.factor(0)

print("Smoking status - cases VS asthma not cases")
asthma_smk <- rbind(nocases_smk, cases_smk)
chisq.test(table(asthma_smk$SI,asthma_smk$tmp_col))

print("Smoking status - cases VS controls")
cases_controls_smk <- rbind(controls_smk, cases_smk)
chisq.test(table(cases_controls_smk$SI,cases_controls_smk$tmp_col))

#Age at onset asthma:


#!/usr/bin/env Rscript

#Rationale: check lung function, BMI of cases with hesin/death with comorbidities and cases with hesin/death without comorbidities.


library(tidyverse)
library(dplyr)
library(data.table)

all <- read.table("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/demo_EUR_pheno_cov_broadasthma.txt", header=T)

all_asthma <- read.table("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/Eid_intersection_asthma_diagnosis_ATLEAST_1_evidence.txt")

# create variable == 1 if individual is in all comer asthma group, else 0
all$allcomer <- ifelse(all$IID %in% all_asthma$V1, 1, 0)

comorbid <- read.table("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/eid_emphchronCOPD_union.txt")
# create variable ==1 if individual has comorbidities, else 0
all$comorbid <- ifelse(all$IID %in% comorbid$V1, 1, 0)

all$eid <- as.factor(all$eid)

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

all_LF_kath <- left_join(all,LF_kath_app56607,by="eid")

#summarising command and comparision between different case/comorbidity status:

print("FEV1")
all_LF_kath %>% group_by(broad_pheno_1_5_ratio,pheno_1_5_ratio,comorbid) %>% dplyr::summarize(fev1 = mean(fev1_perc_pred, na.rm=TRUE))

print("FEV1/FVC")
all_LF_kath %>% group_by(cases_broad,pheno_1_5_ratio,comorbid) %>% dplyr::summarize(ratio_fev1_fvc = mean(ff.best, na.rm=TRUE))

print("BMI")
all_LF_kath %>% group_by(cases_broad,pheno_1_5_ratio,comorbid) %>% dplyr::summarize(BMI = mean(BMI, na.rm=TRUE))

print("sex")
all_LF_kath %>% group_by(cases_broad,pheno_1_5_ratio,comorbid) %>% dplyr::summarize(sex = mean(genetic_sex, na.rm=TRUE))


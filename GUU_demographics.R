#!/usr/bin/env Rscript

#Rationale: demographic of GASP/U-BIOPRED/UKBiobank case-control cohort used in the replication meta-analysis

NB: U-BIOPRED only adult, so it will be GASP with maybe some childhood onset asthma

library(tidyverse)
library(dplyr)
library(data.table)
library(hablar)
library(readxl)

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
#check for category onset: it should be NAs for all of them because they are controls:
table(eur_ukbb_in_meta$category_onset,exclude=NULL)

##lung function for UK-Biobank:
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

#Upload the UKBB data:
ukbb <- fread("ukbb_ctrl_in_GUU_demotraits.txt")
##Need to create smoking status as done for the discovery cohort (with cigarette pack per year and smoking status):
eur_ukbb <- read.table("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/demo_EUR_pheno_cov_broadasthma.txt",header=T,sep=" ") %>%
            select(IID, cigarette_pack_years)
ukbb_demo <- ukbb %>% left_join(eur_ukbb, by = "IID")
ukbb_demo <- ukbb_demo %>% mutate(pack_per_year_threshold = case_when(cigarette_pack_years < 5 ~ "less_than_5",
                                                                          cigarette_pack_years >= 5 ~ "equal_or_more_than_5"))

ukbb_demo <- ukbb_demo %>% mutate(ubiopred_smk = ifelse((ukbb_demo$pack_per_year_threshold == "less_than_5" & ukbb_demo$smoking_status == 1) | ukbb_demo$smoking_status == 0 , "Never_smoker",
                                                      ifelse(ukbb_demo$smoking_status == 2, "Ever_smoker",
                                                      ifelse(ukbb_demo$pack_per_year_threshold == "equal_or_more_than_5" & ukbb_demo$smoking_status == 1, "Ever_smoker", NA))))

ukbb_demo <- ukbb_demo %>% rename(GWAS_ID = IID, Sm_Status = ubiopred_smk, age_onset = category_onset,
                        FEV1_Percent_Predicted = fev1_perc_pred, FEV1_FVC_Ratio = ff.best,
                        Eosinophils = eos_count, Neutrophils = neu_count) %>%
                        select(-eid, -smoking_status, -cigarette_pack_years, -pack_per_year_threshold)
ukbb_demo$pheno <- as.factor(0)

#Upload and clean demographic data for U-BIOPRED:
#Saved the file in tab format (Excel --> export as txt)
ubiopred <- fread("UBIOPRED_data.txt")
#filter for columns:
col2keep <- c("Patient", "cohort", "Body_Mass_Index_(kg/m2)","eosinophils_(x10^3/uL)", "Smoking_Status", "FEV1/FVC_Ratio_Predicted", "FEV1_Predicted_(L)", "Onset_OR_First_Diagnosis_Age_(years)", "neutrophils_(x10^3/uL)")
ubiopred <- ubiopred %>% select(all_of(col2keep))

#bridge file: use the column “Baseline_visit_kitID” to bridge with my case/controls IDs, and “Patient” and “cohort”
#to bridge with UBIOPRED collabs data:
bridge_file <- read_excel("/rfs/TobinGroup/GWAtraits/FEV/AIRPROM/UBIOPRED/phenotype_data/UBIOPRED_pheno.xlsx") %>%
               select("Baseline_visit_kitID","Patient","cohort") %>% rename(IID = Baseline_visit_kitID)
bridge_file$IID <- as.character(bridge_file$IID)

#cases and controls:
ubiopred_cases <- fread("ubiopred_cases_IID.txt")
ubiopred_cases$pheno <- as.factor(1)
ubiopred_controls <- fread("ubiopred_controls_IID.txt")
ubiopred_controls$pheno <- as.factor(0)

#ubiopred and bridge file:
ubiopred_bridge_file <- ubiopred %>% left_join(bridge_file, by = c("Patient","cohort"))

#cases and bridge file:
ubiopred_cases_bridge_file <- ubiopred_cases %>% left_join(bridge_file, by = "IID")
#controls and bridge file:
ubiopred_controls$IID <- as.character(ubiopred_controls$IID)
ubiopred_controls_bridge_file <- ubiopred_controls %>% left_join(bridge_file, by = "IID")

#info I need: BMI, Smoking status, Age onset asthma, FEV1 predicted, FEV1/FVC, Eosinophil count, Neutrophil count
ubiopred_cases_demo <- ubiopred_cases_bridge_file %>% left_join(ubiopred_bridge_file, by = c("IID", "Patient", "cohort"))
ubiopred_controls_demo <- ubiopred_controls_bridge_file %>% left_join(ubiopred_bridge_file, by = c("IID", "Patient", "cohort"))

#merge ubiopred case and controls:
ubiopred_demo <- rbind(ubiopred_cases_demo, ubiopred_controls_demo)
#Age onset asthma: need to create adult and childhood variable: if age on set < 18, childhood, if >= adult.
ubiopred_demo <- ubiopred_demo %>% mutate(age_onset = ifelse(ubiopred_demo$'Onset_OR_First_Diagnosis_Age_(years)' < 18, "onset_early", "onset_adult"))
ubiopred_demo <- ubiopred_demo %>% rename(GWAS_ID = "IID", BMI = "Body_Mass_Index_(kg/m2)", Sm_Status = Smoking_Status,
                        FEV1_Percent_Predicted = "FEV1_Predicted_(L)", FEV1_FVC_Ratio = "FEV1/FVC_Ratio_Predicted",
                        Eosinophils = "eosinophils_(x10^3/uL)" , Neutrophils = "neutrophils_(x10^3/uL)")
ubiopred_demo <- ubiopred_demo %>% select(GWAS_ID, pheno, BMI, Sm_Status, FEV1_Percent_Predicted,
                                          FEV1_FVC_Ratio, age_onset, Eosinophils, Neutrophils)
#change smoking status into ever and never:
ubiopred_demo <- ubiopred_demo %>% mutate(Sm_Status = ifelse(ubiopred_demo$Sm_Status == "current_smoker", "Ever_smoker",
 ifelse(ubiopred_demo$Sm_Status == "ex_smoker", "Ever_smoker",
 ifelse(ubiopred_demo$Sm_Status == "non_smoker", "Never_smoker", "NA"))))
ubiopred_demo <- ubiopred_demo %>% mutate(FEV1_FVC_Ratio = ubiopred_demo$FEV1_FVC_Ratio/100)

#GASP demographics:
gasp <- read_excel("GASP_demographics.xlsx", sheet = "gasp")
#info I need: BMI, Smoking status, Age onset asthma, FEV1 predicted, FEV1/FVC, Eosinophil count, Neutrophil count
gasp_IID <- read_table("gasp_cases_IID.txt") %>% rename(GWAS_ID = IID)
gasp_demo <- gasp_IID %>% left_join(gasp, by = "GWAS_ID")
gasp_demo$pheno <- as.factor(1)
gasp_demo <- gasp_demo %>% mutate(age_onset = ifelse(gasp_demo$'Age_of_Onset' < 18, "onset_early", "onset_adult"))
#change smoking status into ever and never:
gasp_demo <- gasp_demo %>% mutate(Sm_Status = ifelse(gasp_demo$Sm_Status == "Current smoker", "Ever_smoker",
 ifelse(gasp_demo$Sm_Status == 'Ex-smoker', "Ever_smoker",
 ifelse(gasp_demo$Sm_Status == "Never", "Never_smoker", "NA"))))
gasp_demo <- gasp_demo %>% select(GWAS_ID, pheno, BMI, Sm_Status, FEV1_Percent_Predicted,
                                  FEV1_FVC_Ratio, age_onset, Eosinophils, Neutrophils)

#Merge all the three datasets:
all_demo <- rbind(ukbb_demo, gasp_demo, ubiopred_demo)
write.table(all_demo, "all_GUU_demographics.txt", sep = "\t", quote = F)
write.table(all_demo, "/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/all_GUU_demographics.txt", sep = "\t", quote = F)
all_demo$Eosinophils <- as.numeric(all_demo$Eosinophils)
cases_all_demo <- all_demo %>% filter(pheno == 1)
controls_all_demo <- all_demo %>% filter(pheno == 0)

#BMI:
print(mean(cases_all_demo$BMI, na.rm=TRUE))
print(sd(cases_all_demo$BMI, na.rm=TRUE))
print(mean(controls_all_demo$BMI, na.rm=TRUE))
print(sd(controls_all_demo$BMI, na.rm=TRUE))

#Smoking status:
print(table(cases_all_demo$Sm_Status, useNA = "always"))
print(table(controls_all_demo$Sm_Status, useNA = "always"))

#Category onset asthma:
print(table(cases_all_demo$age_onset, useNA = "always"))
print(table(controls_all_demo$age_onset, useNA = "always"))

#FEV1 Predicted:
print(mean(cases_all_demo$FEV1_Percent_Predicted, na.rm=TRUE))
print(sd(cases_all_demo$FEV1_Percent_Predicted, na.rm=TRUE))
print(mean(controls_all_demo$FEV1_Percent_Predicted, na.rm=TRUE))
print(sd(controls_all_demo$FEV1_Percent_Predicted, na.rm=TRUE))

#FEV1_FVC_Ratio
print(mean(cases_all_demo$FEV1_FVC_Ratio, na.rm=TRUE))
print(sd(cases_all_demo$FEV1_FVC_Ratio, na.rm=TRUE))
print(mean(controls_all_demo$FEV1_FVC_Ratio, na.rm=TRUE))
print(sd(controls_all_demo$FEV1_FVC_Ratio, na.rm=TRUE))

#Eosinophil count:
print(mean(cases_all_demo$Eosinophils, na.rm=TRUE))
print(sd(cases_all_demo$Eosinophils, na.rm=TRUE))
print(mean(controls_all_demo$Eosinophils, na.rm=TRUE))
print(sd(controls_all_demo$Eosinophils, na.rm=TRUE))

#Neutrophil count:
print(mean(cases_all_demo$Neutrophils, na.rm=TRUE))
print(sd(cases_all_demo$Neutrophils, na.rm=TRUE))
print(mean(controls_all_demo$Neutrophils, na.rm=TRUE))
print(sd(controls_all_demo$Neutrophils, na.rm=TRUE))
#!/usr/bin/env Rscript


#run and save it as : Rscript Descriptive_difference_test.R > descriptive_difference_test_report

library(tidyverse)
library(dplyr)
library(data.table)
library(hablar)
library(VennDiagram)
library(RColorBrewer)

#load input
demo <- read.table("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/demo_EUR_pheno_cov_broadasthma.txt",header=T,sep=" ")
demo$eid <- as.factor(demo$eid)
demo <- demo %>% rename(clustered_ancestry = clustered.ethnicity)
demo$clustered_ancestry  <- as.factor(demo$clustered_ancestry)
demo$smoking_status <- as.factor(demo$smoking_status)
demo$broad_pheno_1_5_ratio <- as.factor(demo$broad_pheno_1_5_ratio)

#all-comer asthma, all genetic ancestry:
#cp Eid_intersection_asthma_diagnosis_ATLEAST_1_evidence.txt /home/n/nnp5/PhD/PhD_project/tmp_manuscript/
allcomer_allanc <- fread("Eid_intersection_asthma_diagnosis_ATLEAST_1_evidence.txt") %>% rename(eid = V1)
allcomer_allanc$allcomer_asthma_allanc <- as.factor(1)
allcomer_allanc$eid <- as.factor(allcomer_allanc$eid)
demo <- demo %>% left_join(allcomer_allanc, by = "eid")

#check of ECB/COPD comorbidities:
#cp eid_emphchronCOPD_union.txt /home/n/nnp5/PhD/PhD_project/tmp_manuscript/
copd_ecb <- fread("eid_emphchronCOPD_union.txt") %>% rename(eid = V1)
copd_ecb$comob <- as.factor(1)
demo_comob <- demo %>% left_join(copd_ecb, by = "eid")

#hesin:
#app56607:
hes <- fread("QC_hesin_diag_asthma.txt") %>% rename(eid = app56607_ids, level_hesin = level) %>% select(eid, level_hesin) %>% unique()
hes$eid <- as.factor(hes$eid)
hesin_1 <- hes %>% filter(level_hesin == 1)
#death
#grep "J45\|J46" /rfs/TobinGroup/data/UKBiobank/application_56607/hes/death_cause.txt \
#    > eid_asthma_death
death <- fread("eid_asthma_death")
colnames(death) <- c("eid","ins_index","arr_index","level_death","cause_icd10")
death <- death %>% select(eid, level_death) %>% unique()
death$eid <- as.factor(death$eid)
death_1 <- death %>% filter(level_death == 1)
demo_comob_hesin <- demo_comob %>% left_join(hesin_1, by = "eid")
demo <- demo_comob_hesin %>% left_join(death_1, by = "eid")

#Total of 71,899 individuals with asthma either with or wihtout comorbidities, european-like genetic ancestry:
eur_all_comer <- demo %>% filter(allcomer_asthma_allanc == 1)

#smoking as per ubiopred:
#Split participants (both case and controls) into
#non-smokers (according to Shaw et al. 2015: previous/non current smoker with < 5 pack per year history) and
#smokers (according to Shaw et al. 2015: current smoker or ex-smokers with >= 5 pack year history)

eur_all_comer <- eur_all_comer %>% mutate(pack_per_year_threshold = case_when(cigarette_pack_years < 5 ~ "less_than_5",
                                                                          cigarette_pack_years >= 5 ~ "equal_or_more_than_5"))

eur_all_comer <- eur_all_comer %>% mutate(ubiopred_smk = ifelse((eur_all_comer$pack_per_year_threshold == "less_than_5" & eur_all_comer$smoking_status == 1) | eur_all_comer$smoking_status == 0 , "non_smoker",
                                                      ifelse(eur_all_comer$smoking_status == 2, "smoker",
                                                      ifelse(eur_all_comer$pack_per_year_threshold == "equal_or_more_than_5" & eur_all_comer$smoking_status == 1, "smoker", NA))))

##Hay fever-rhinitis, eczema/atopic dermatitis:
allergy <- fread("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/eid_union_hayfev_rhinitis_eczema_derma_ATLEAST_1_evidence.txt",header=F)
allergy <- unique(allergy)
colnames(allergy)[1] <- "eid"
allergy$eid <- as.character(allergy$eid)
allergy$allergy <- as.factor("1")
eur_all_comer <- left_join(eur_all_comer,allergy,by="eid")
#eur_all_comer_allergy <- left_join(eur_all_comer,allergy,by="eid")

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
#eur_all_comer_LF_kath <- left_join(eur_all_comer,LF_kath_app56607,by="eid")
eur_all_comer <- left_join(eur_all_comer,LF_kath_app56607,by="eid")

#create column with the different categories I want to compare:
#Check demographics in five different group:
#eur_asthma_comob
#eur_asthma_no_comob
#cases_prescript_only
#cases_hesin_death_comob
#cases_hesin_death_no_comob
#different sub-groups:
eur_all_comer_comob <- eur_all_comer %>% filter(comob == 1)
eur_all_comer_nocomob <- eur_all_comer %>% filter(is.na(comob))
eur_all_comer_comob_cases <- eur_all_comer_comob %>% filter(broad_pheno_1_5_ratio == 1)
eur_all_comer_nocomob_cases <- eur_all_comer_nocomob %>% filter(broad_pheno_1_5_ratio == 1)
eur_all_comer_nocomob_cases_scripts <- eur_all_comer_nocomob %>% filter(broad_pheno_1_5_ratio == 1 & (is.na(level_hesin) | is.na(level_death)))

##eaosinophils and neutrophils need a little of more thinking to create the dataset.
#eosinophil:
eos <- fread("Eosinophill_count_30150.csv")
colnames(eos) <- c("eid","eos1","eos2","eso3")
eos$eid <- as.factor(eos$eid)
eur_all_comer_eos <- left_join(eur_all_comer,eos,by="eid") %>% select(broad_pheno_1_5_ratio,eos1,eos2,eso3)
eur_all_comer_eos <- eur_all_comer_eos %>%
  rowwise() %>%
  mutate(min_eos = min_(c_across(-broad_pheno_1_5_ratio)),
         max_eos = max_(c_across(-broad_pheno_1_5_ratio)))
#neutrophil:
neu <- fread("Neutrophill_count_30140.csv")
colnames(neu) <- c("eid","neu1","neu2","eso3")
neu$eid <- as.factor(neu$eid)
eur_all_comer_neu <- left_join(eur_all_comer,neu,by="eid") %>% select(broad_pheno_1_5_ratio,neu1,neu2,eso3)
eur_all_comer_neu <- eur_all_comer_neu %>%
  rowwise() %>%
  mutate(min_neu = min_(c_across(-broad_pheno_1_5_ratio)),
         max_neu = max_(c_across(-broad_pheno_1_5_ratio)))




#Total of 71,899 individuals with asthma either with or wihtout comorbidities, european-like genetic ancestry:
eur_all_comer <- eur_all_comer_comob %>% filter(allcomer_asthma_allanc == 1)
print("Demographics analysis in five different groups of European-like genetic ancestry with diagnosis for asthma:")


#age :
print("Age all-comer asthma")
mean(eur_all_comer$age_at_recruitment)
sd(eur_all_comer$age_at_recruitment)
mean(eur_all_comer_comob$age_at_recruitment)
mean(eur_all_comer_nocomob$age_at_recruitment)
mean(eur_all_comer_comob_cases$age_at_recruitment)
mean(eur_all_comer_nocomob_cases$age_at_recruitment)

wilcox.test(eur_all_comer$age_at_recruitment~eur_all_comer$broad_pheno_1_5_ratio)

#BMI
print("BMI")
mean(eur_all_comer_comob$BMI, na.rm=TRUE)
mean(eur_all_comer_nocomob$BMI, na.rm=TRUE)
mean(eur_all_comer_comob_cases$BMI, na.rm=TRUE)
mean(eur_all_comer_nocomob_cases$BMI, na.rm=TRUE)
wilcox.test(eur_all_comer$BMI~eur_all_comer$broad_pheno_1_5_ratio)

#FEV1 % predicted
print("FEV1 % predicted")
mean(eur_all_comer_comob$fev1_perc_pred, na.rm=TRUE)
mean(eur_all_comer_nocomob$fev1_perc_pred, na.rm=TRUE)
mean(eur_all_comer_comob_cases$fev1_perc_pred, na.rm=TRUE)
mean(eur_all_comer_nocomob_cases$fev1_perc_pred, na.rm=TRUE)
wilcox.test(eur_all_comer_LF_kath$fev1_perc_pred~eur_all_comer_LF_kath$broad_pheno_1_5_ratio)


#ratio_FEV1_FVC
print("ratio FEV1 FVC")
mean(eur_all_comer_comob$ff.best, na.rm=TRUE)
mean(eur_all_comer_nocomob$ff.best, na.rm=TRUE)
mean(eur_all_comer_comob_cases$ff.best, na.rm=TRUE)
mean(eur_all_comer_nocomob_cases$ff.best, na.rm=TRUE)
wilcox.test(eur_all_comer_LF_kath$ff.best~eur_all_comer_LF_kath$broad_pheno_1_5_ratio)

#eosinophil:
print("eosinophil")
wilcox.test(eur_all_comer_eos$max_eos~eur_all_comer_eos$broad_pheno_1_5_ratio)

#neutrophil:
print("neutrophil")
wilcox.test(eur_all_comer_neu$max_neu~eur_all_comer_eos$broad_pheno_1_5_ratio)

print("Test categorical variable with chi-square test:")
print("NB: Category onset, Prednisolone use, hospitalisation cannot be test as specific to cases only")

#sex:
print("sex")
chisq.test(table(eur_all_comer$genetic_sex,eur_all_comer$broad_pheno_1_5_ratio))

#Smoking status (ubiopred):
print("smokig status")
chisq.test(table(eur_all_comer$ubiopred_smk,eur_all_comer$broad_pheno_1_5_ratio))

#Hay fever-rhinitis, eczema/atopic dermatitis:
print("allergic condition")
chisq.test(table(eur_all_comer_allergy$allergy,eur_all_comer_allergy$broad_pheno_1_5_ratio))

print("Test continuous variable with Wilcox Mann-Whitney test:")


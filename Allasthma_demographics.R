#!/usr/bin/env Rscript

#Rationale: demographic of all-asthma europeans - to be compared with the severe asthma phenotype for the manuscript

library(tidyverse)
library(dplyr)
library(data.table)
library(hablar)

#load input
##demographic table:
demo_eur <- read.table("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/demo_EUR_pheno_cov_broadasthma.txt",header=T,sep=" ")
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

#hesin:
#app56607:
hes <- fread("QC_hesin_diag_asthma.txt") %>% rename(eid = app56607_ids, level_hesin = level) %>% select(eid, level_hesin) %>% unique()
hes$eid <- as.factor(hes$eid)
hesin_1 <- hes %>% filter(level_hesin == 1)
demo_hes <- demo %>% left_join(hes, by = "eid")

#death
#grep "J45\|J46" /rfs/TobinGroup/data/UKBiobank/application_56607/hes/death_cause.txt \
#    > eid_asthma_death
death <- fread("eid_asthma_death")
colnames(death) <- c("eid","ins_index","arr_index","level_death","cause_icd10")
death <- death %>% select(eid, level_death) %>% unique()
death$eid <- as.factor(death$eid)
death_1 <- death %>% filter(level_death == 1)
demo_hes_death <- demo_hes %>% left_join(death, by = "eid")


#function outliers:
outliers <- function(x) {

  Q1 <- quantile(x, probs=.25)
  Q3 <- quantile(x, probs=.75)
  iqr = Q3-Q1

 upper_limit = Q3 + (iqr*1.5)
 lower_limit = Q1 - (iqr*1.5)

 x > upper_limit | x < lower_limit
}

remove_outliers <- function(df, cols = names(df)) {
  for (col in cols) {
    df <- df[!outliers(df[[col]]),]
  }
  df
}

#eosinophil:
eos <- fread("Eosinophill_count_30150.csv")
colnames(eos) <- c("eid","eos1","eos2","eso3")
eos$eid <- as.factor(eos$eid)
demo_eos <- left_join(demo,eos,by="eid") %>% select(eos1,eos2,eso3)
demo_eos <- demo_eos %>%
  rowwise() %>%
  mutate(min_eos = min_(c_across()),
         max_eos = max_(c_across()))


#neutrophil:
neu <- fread("Neutrophill_count_30140.csv")
colnames(neu) <- c("eid","neu1","neu2","neu3")
neu$eid <- as.factor(neu$eid)
demo_neu <- left_join(demo,neu,by="eid") %>% select(broad_pheno_1_5_ratio,neu1,neu2,neu3)
demo_neu <- demo_neu %>%
  rowwise() %>%
  mutate(min_neu = min_(c_across()),
         max_neu = max_(c_across()))



#smoking as per ubiopred:
#Split participants (both case and controls) into
#non-smokers (according to Shaw et al. 2015: previous/non current smoker with < 5 pack per year history) and
#smokers (according to Shaw et al. 2015: current smoker or ex-smokers with >= 5 pack year history)
demo <- demo %>% mutate(pack_per_year_threshold = case_when(cigarette_pack_years < 5 ~ "less_than_5",
                                                                          cigarette_pack_years >= 5 ~ "equal_or_more_than_5"))
demo <- demo %>% mutate(ubiopred_smk = ifelse((demo$pack_per_year_threshold == "less_than_5" & demo$smoking_status == 1) | demo$smoking_status == 0 , "non_smoker",
                                                      ifelse(demo$smoking_status == 2, "smoker",
                                                      ifelse(demo$pack_per_year_threshold == "equal_or_more_than_5" & demo$smoking_status == 1, "smoker", NA))))


#Prednisolone use:
#from UKBiobank_asthmaMeds_stratification_27062023.tar.gz, copy data/all_UKBB_with_prednisolone_gp_scripts_edit
#and data/all_UKBB_with_any_gp_scripts_edit into tmp_manuscript/
eid_pred <- fread("all_UKBB_with_prednisolone_gp_scripts_edit")
eid_scripts <- fread("all_UKBB_with_any_gp_scripts_edit")
eid_pred$pred_use <- as.factor(1)
eid_scripts$scripts <- as.factor(1)
eid_pred$V1 <- as.factor(eid_pred$V1)
eid_scripts$V1 <- as.factor(eid_scripts$V1)
eid_script_pred <- left_join(eid_scripts,eid_pred,by="V1")
eid_script_pred <- eid_script_pred %>% mutate(pred_use = ifelse(is.na(eid_script_pred$pred_use), 0, 1))
colnames(eid_script_pred)[1] <- "eid"
eid_script_pred$eid <- as.factor(eid_script_pred$eid)
eid_script_pred <- eid_script_pred %>% select(eid,pred_use)
demo <- left_join(demo,eid_script_pred,by="eid")

##Hay fever-rhinitis, eczema/atopic dermatitis:
allergy <- fread("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/eid_union_hayfev_rhinitis_eczema_derma_ATLEAST_1_evidence.txt",header=F)
allergy <- unique(allergy)
colnames(allergy)[1] <- "eid"
allergy$eid <- as.character(allergy$eid)
allergy$allergy <- as.factor("1")
demo <- left_join(demo,allergy,by="eid")

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
demo <- demo %>% left_join(LF_kath_app56607,by="eid")

# ECB/COPD comorbidities:
#cp eid_emphchronCOPD_union.txt /home/n/nnp5/PhD/PhD_project/tmp_manuscript/
copd_ecb <- fread("eid_emphchronCOPD_union.txt") %>% rename(eid = V1)
copd_ecb$comob <- as.factor(1)
demo <- demo %>% left_join(copd_ecb, by = "eid")

#Save descriptive tables for all-asthma Europeans:
fwrite(demo,"all_asthma_demographics_20102024",sep="\t",quote=F)
fwrite(demo_hes_death,"all_asthma_demographics_hesdeath_20102024",sep="\t",quote=F)
fwrite(demo_eos,"all_asthma_demographics_eos_20102024",sep="\t",quote=F)
fwrite(demo_neu,"all_asthma_demographics_neu_20102024",sep="\t",quote=F)

#Create descriptive function - calculate summary stats for demographic info:
descriptive <- function(demo,demo_hes_death,demo_eos,demo_neu){

#sex:
print("Sex : count and percentage")
print(table(demo$genetic_sex,exclude=NULL))
print(prop.table(table(demo$genetic_sex,exclude=NULL)))

#age : mean and SD:
print("Age: mean and SD")
print(sum(is.na(demo$age_at_recruitment)))
print(mean(demo$age_at_recruitment,na.rm=TRUE))
print(sd(demo$age_at_recruitment,na.rm=TRUE))

#BMI : mean and SD:
print("BMI: mean and SD")
print(sum(is.na(demo$BMI)))
bmi <- demo %>% filter(!is.na(BMI))
bmi <- remove_outliers(bmi, 'BMI')
print(mean(bmi$BMI,na.rm=TRUE))
print(sd(bmi$BMI,na.rm=TRUE))

##blood cell count
#eosinophil:
print("Eosinophils: mean and SD")
print(sum(is.na(demo_eos$max_eos)))
df <- demo_eos %>% filter(!is.na(max_eos))
df <- remove_outliers(df, 'max_eos')
print(mean(df$max_eos))
print(sd(df$max_eos))

#neutrophils:
print("Neutrophils: mean and SD")
print(sum(is.na(demo_neu$max_neu)))
df <- demo_neu %>% filter(!is.na(max_neu))
df <- remove_outliers(df, 'max_neu')
print(mean(df$max_neu))
print(sd(df$max_neu))

#Smoking status (ubiopred):
print("Smoking status (ubiopred): count and percentage")
print("Cases")
print(table(demo$ubiopred_smk,exclude=NULL))
print(prop.table(table(demo$ubiopred_smk,exclude=NULL)))

#fev1_perc_pred:
print("fev1_perc_pred (From Kath data): mean and SD")
print(sum(is.na(demo$fev1_perc_pred)))
df <- demo %>% filter(!is.na(fev1_perc_pred))
df <- remove_outliers(df, 'fev1_perc_pred')
print(mean(df$fev1_perc_pred))
print(sd(df$fev1_perc_pred))

#ff.best:
print("ff.best (From Kath data): mean and SD")
print(sum(is.na(demo$ff.best)))
df <- demo %>% filter(!is.na(ff.best))
df <- remove_outliers(df, 'ff.best')
print(mean(df$ff.best))
print(sd(df$ff.best))

#Hospitalisation:
if (demo_hes_death == "NA"){
print("controls does not need further descriptive")
}
else {
print("Hospitalisation: count and percentage")
print(table(demo_hes_death$level_hesin,exclude=NULL))
print(prop.table(table(demo_hes_death$level_hesin,exclude=NULL)))

#Death:
print("Death: count and percentage")
print(table(demo_hes_death$level_death,exclude=NULL))
print(prop.table(table(demo_hes_death$level_death,exclude=NULL)))

#Category onset:
print("Category onset: count and percentage")
print("Cases")
print(table(demo$category_onset,exclude = NULL))
print(prop.table(table(demo$category_onset,exclude = NULL)))

#Prednisolone use:
print("Prednisolone use: count and percentage")
print("Cases")
print(table(demo$pred_use,exclude = NULL))
print(prop.table(table(demo$pred_use,exclude = NULL)))

#Hay fever-rhinitis, eczema/atopic dermatitis:
print("Hay fever-rhinitis, eczema/atopic dermatitis:")
print(table(demo$allergy,exclude=NULL))
print(prop.table(table(demo$allergy,exclude=NULL)))

#comorbidities:
print("Comorbidities:")
print(table(demo$comob,exclude=NULL))
print(prop.table(table(demo$comob,exclude=NULL)))
}
}

#Descriptive all-comer asthma:
descriptive(demo,demo_hes_death,demo_eos,demo_neu)

#Descriptive all-comer asthma not cases:
cases <- demo %>% filter(cases_broad_EUR == 1)
cases_IID <- cases$IID
demo_notcases <- demo %>% filter(! IID %in% cases_IID)
demo_notcases_hes_death <- demo_hes_death %>% filter(! IID %in% cases_IID)
demo_notcases_eos <- left_join(demo_notcases,eos,by="eid") %>% select(eos1,eos2,eso3)
demo_notcases_eos <- demo_notcases_eos %>%
  rowwise() %>%
  mutate(min_eos = min_(c_across()),
         max_eos = max_(c_across()))
demo_notcases_neu <- left_join(demo_notcases,neu,by="eid") %>% select(neu1,neu2,neu3)
demo_notcases_neu <- demo_notcases_neu %>%
  rowwise() %>%
  mutate(min_neu = min_(c_across()),
         max_neu = max_(c_across()))

#Save descriptive tables for all-asthma Europeans:
fwrite(demo_notcases,"notcases_asthma_demographics_20102024",sep="\t",quote=F)
fwrite(demo_notcases_hes_death,"notcases_asthma_demographics_hesdeath_20102024",sep="\t",quote=F)
fwrite(demo_notcases_eos,"notcases_asthma_demographics_eos_20102024",sep="\t",quote=F)
fwrite(demo_notcases_neu,"notcases_asthma_demographics_neu_20102024",sep="\t",quote=F)
descriptive(demo_notcases,demo_notcases_hes_death,demo_notcases_eos,demo_notcases_neu)


#Descriptive asthma cases:
cases <- demo %>% filter(cases_broad_EUR == 1)
cases_hes_death <- demo_hes_death %>% filter(IID %in% cases_IID)
cases_eos <- left_join(cases,eos,by="eid") %>% select(eos1,eos2,eso3)
cases_eos <- cases_eos %>%
  rowwise() %>%
  mutate(min_eos = min_(c_across()),
         max_eos = max_(c_across()))
cases_neu <- left_join(cases,neu,by="eid") %>% select(neu1,neu2,neu3)
cases_neu <- cases_neu %>%
  rowwise() %>%
  mutate(min_neu = min_(c_across()),
         max_neu = max_(c_across()))

#Save descriptive tables for all-asthma Europeans:
fwrite(cases,"cases_asthma_demographics_20102024",sep="\t",quote=F)
fwrite(cases_hes_death,"cases_asthma_demographics_hesdeath_20102024",sep="\t",quote=F)
fwrite(cases_eos,"cases_asthma_demographics_eos_20102024",sep="\t",quote=F)
fwrite(cases_neu,"cases_asthma_demographics_neu_20102024",sep="\t",quote=F)
descriptive(cases,cases_hes_death,cases_eos,cases_neu)

#Descriptive controls:
df_tmp <- read.table("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/demo_EUR_pheno_cov_broadasthma.txt",header=T,sep=" ")
controls <- df_tmp %>% filter(broad_pheno_1_5_ratio == 0)
controls$eid <- as.factor(controls$eid)
controls <- controls %>% left_join(LF_kath_app56607,by="eid")
controls_eos <- left_join(controls,eos,by="eid") %>% select(eos1,eos2,eso3)
controls_eos <- controls_eos %>%
  rowwise() %>%
  mutate(min_eos = min_(c_across()),
         max_eos = max_(c_across()))
controls_neu <- left_join(controls,neu,by="eid") %>% select(neu1,neu2,neu3)
controls_neu <- controls_neu %>%
  rowwise() %>%
  mutate(min_neu = min_(c_across()),
         max_neu = max_(c_across()))

fwrite(controls_eos,"controls_asthma_demographics_eos_23102024",sep="\t",quote=F)
fwrite(controls_neu,"controls_asthma_demographics_neu_23102024",sep="\t",quote=F)

controls <- controls %>% mutate(pack_per_year_threshold = case_when(cigarette_pack_years < 5 ~ "less_than_5",
                                                                          cigarette_pack_years >= 5 ~ "equal_or_more_than_5"))
controls$pack_per_year_threshold <- as.factor(controls$pack_per_year_threshold)
controls <- controls %>% mutate(ubiopred_smk = ifelse((controls$pack_per_year_threshold == "less_than_5" & controls$smoking_status == 1) | controls$smoking_status == 0 , "non_smoker",
                                                      ifelse(controls$smoking_status == 2, "smoker",
                                                      ifelse(controls$pack_per_year_threshold == "equal_or_more_than_5" & controls$smoking_status == 1, "smoker", NA))))


fwrite(controls,"controls_asthma_demographics_23102024",sep="\t",quote=F)
descriptive(controls,"NA",controls_eos,controls_neu)


#Statistical difference with wilcox.test() between asthma-non cases and cases:
demo_notcases$NC_C <- as.factor(0)
cases$NC_C <- as.factor(1)
df <- rbind(demo_notcases,cases)

##create a column in df: asthma-non cases and cases:
print("genetic sex")
wilcox.test(df$genetic_sex~df$NC_C,)$p.value
print("age at recruitment")
wilcox.test(df$age_at_recruitment~df$NC_C)$p.value
print("BMI")
wilcox.test(df$BMI~df$NC_C)$p.value
print("FEV1 % predicted")
wilcox.test(df$fev1_perc_pred~df$NC_C)$p.value
print("FEV1/FVC")
wilcox.test(df$ff.best~df$NC_C)$p.value
print("Smoking status")
chisq.test(table(df$ubiopred_smk,df$NC_C))
print("Category onset")
chisq.test(table(df$category_onset,df$NC_C))
print("prednisolone use")
chisq.test(table(df$pred_use,df$NC_C))
print("allergy")
chisq.test(table(df$allergy,df$NC_C))
print("comorbidities")
chisq.test(table(df$comob,df$NC_C))

print("eosinophils")
demo_notcases_eos$NC_C <- as.factor(0)
cases_eos$NC_C <- as.factor(1)
df_eos <- rbind(demo_notcases_eos,cases_eos)
wilcox.test(df_eos$max_eos~df_eos$NC_C)$p.value

print("neutrophils")
demo_notcases_neu$NC_C <- as.factor(0)
cases_neu$NC_C <- as.factor(1)
df_neu <- rbind(demo_notcases_neu,cases_neu)
wilcox.test(df_neu$max_neu~df_neu$NC_C)$p.value

##statistical difference between cases/control:
cases <- fread("cases_asthma_demographics_20102024")
cases_eos <- fread("cases_asthma_demographics_eos_20102024")
cases_neu <- fread("cases_asthma_demographics_neu_20102024")

cases_tmp <- cases %>% select(broad_pheno_1_5_ratio,genetic_sex,age_at_recruitment,BMI,fev1_perc_pred,ff.best) #ubiopred_smk, to add when trobleshooting it
controls_tmp <- controls %>% select(broad_pheno_1_5_ratio,genetic_sex,age_at_recruitment,BMI,fev1_perc_pred,ff.best) #ubiopred_smk, to add when trobleshooting it
df <- rbind(cases_tmp,controls_tmp)
print("genetic sex")
wilcox.test(df$genetic_sex~df$broad_pheno_1_5_ratio,)$p.value
print("age at recruitment")
wilcox.test(df$age_at_recruitment~df$broad_pheno_1_5_ratio)$p.value
print("BMI")
wilcox.test(df$BMI~df$broad_pheno_1_5_ratio)$p.value
print("FEV1 % predicted")
wilcox.test(df$fev1_perc_pred~df$broad_pheno_1_5_ratio)$p.value
print("FEV1/FVC")
wilcox.test(df$ff.best~df$broad_pheno_1_5_ratio)$p.value
print("Smoking status")
chisq.test(table(df$ubiopred_smk,df$broad_pheno_1_5_ratio))
print("eosinophils")
cases_eos$broad_pheno_1_5_ratio <- as.factor(0)
controls_eos$broad_pheno_1_5_ratio <- as.factor(1)
df_eos <- rbind(cases_eos,controls_eos)
wilcox.test(df_eos$max_eos~df_eos$broad_pheno_1_5_ratio)$p.value
print("neutrophils")
cases_neu$broad_pheno_1_5_ratio <- as.factor(0)
controls_neu$broad_pheno_1_5_ratio <- as.factor(1)
df_neu <- rbind(cases_neu,controls_neu)
wilcox.test(df_neu$max_neu~df_neu$broad_pheno_1_5_ratio)$p.value





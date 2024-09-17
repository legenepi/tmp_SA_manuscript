#!/usr/bin/env Rscript

#Rationale: Create phenotype-covariate file for the sensitivity analysis without cases with comorbidities

library(tidyverse)
library(dplyr)
library(data.table)

all <- read.table("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/demo_EUR_pheno_cov_broadasthma.txt", header=T)
comorbid <- read.table("/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/eid_emphchronCOPD_union.txt")
# create variable ==1 if individual has comorbidities, else 0
all$comorbid <- ifelse(all$IID %in% comorbid$V1, 1, 0)
all <- all %>% mutate(broad_pheno_nocomob = ifelse(broad_pheno_1_5_ratio == 1 & comorbid == 1, NA,
broad_pheno_1_5_ratio))
summary(as.factor(all$broad_pheno_nocomob))

write.table(all,"/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/demo_EUR_pheno_cov_broadasthma_nocomob.txt",
row.names = FALSE, col.names = TRUE ,quote=FALSE, sep=" ", na = "NA")
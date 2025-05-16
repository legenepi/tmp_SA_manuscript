#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(readxl)

args = commandArgs(trailingOnly=TRUE)
locus <- as.character(args[2])
credsetnumber <- args[3]
rsid <- args[4]

gwas <- fread(args[1])
credset <- read_excel("credsetSNPs_forlocuszoom.xlsx") %>%
           rename(snpid = `SNP ID`, Credible_set = `Credible set`, Replicated_locus = `Replicated locus`) %>%
           filter(Replicated_locus == locus, Credible_set == credsetnumber)
gwas_credset <- gwas %>% left_join(credset, by = "snpid") %>%
                mutate(ANNOT = ifelse(is.na(ANNOT), 0, ANNOT))
table(gwas_credset$ANNOT)
gwas_credset <- gwas_credset[,c(1:10,17)]
fwrite(gwas_credset,paste0("/scratch/gen1/nnp5/tmp_manuscript/header_",locus,"_",rsid,"_annot"), quote=F, sep="\t")
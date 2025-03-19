#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

#susie
susie_PIP <- fread(args[1]) %>% rename(PIP=V2, vars.variable=V1)
susie_credset_snp <- as.data.frame(fread(args[2])) %>% rename(vars.variable=V1) %>% mutate(vars.variable=paste0("V",vars.variable))
susie_sumstat <- fread(args[3])
susie_df <- cbind(susie_sumstat, susie_PIP)
susie_df_credset <- inner_join(susie_df,susie_credset_snp,by="vars.variable")
susie_df_credset <- susie_df_credset %>% rename(vars.variable_prob=PIP)
susie_df_credset$locus <- as.character(args[4])
susie_df_credset <- susie_df_credset %>% select(locus, rsid, chromosome, position, allele1, allele2, vars.variable_prob,vars.variable)
susie_df_credset <- susie_df_credset %>% rename(PIP=vars.variable_prob)
fwrite(susie_df_credset,args[5],row.names=F,quote=F,sep="\t")

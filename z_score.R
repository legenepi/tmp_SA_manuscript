#!/usr/bin/env Rscript

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

args = commandArgs(trailingOnly=TRUE)

gwas <- fread(args[1])
gwas$z_score <- gwas$beta / gwas$se
write.table(gwas$z_score,args[2],row.names=FALSE,quote=FALSE,col.names=F)

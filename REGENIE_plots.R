#!/usr/bin/env Rscript
library(qqman) #qq plot and manhattan plot
library(ggplot2)
library(data.table)
library(grid)
library(tidyverse)
#call the plot_functions.R for QQ-plot,Manhattan and lambda:
source("src/plot_functions.R")

args = commandArgs(TRUE)


#inputfile:
input_file = args[1] #input mungestat file
pheno = args[2]
ldsc_intercept = as.numeric(args[3])
#input
meta <- fread(input_file, header=T, fill=T)

title_plot <- paste0("GWAS_",pheno)


#QQ plot with qqman package: one genome-wide for each test-statistic p-val:
plot.qqplot(pval_vec = meta$pval, title= title_plot)


#Manhattan plot:
# require these columns: rs,chr,ps,pval
#df <- meta %>% select(snpid,b37chr,bp,pval)
#colnames(df) <- c("rs","chr","ps","pval")
#df$chr <- as.numeric(df$chr)
#df$ps <- as.numeric(df$ps)
##magnify for p<10-15 association:
#plot.Manha(df, title = title_plot)
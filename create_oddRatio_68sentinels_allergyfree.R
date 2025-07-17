#!/usr/bin/env Rscript

#Rationale: re-create the OR and (%% CI to add teh rs on chr1 that was missed.

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

library(data.table)
library(tidyverse)
sumstat <- fread("/scratch/gen1/nnp5/tmp_manuscript/68sentvars_noallergy")
sumstat$odds_ratio <- exp(sumstat$V6)
sumstat$ci_lower <- exp(sumstat$V6 - qnorm(0.975)*sumstat$V7)
sumstat$ci_upper <- exp(sumstat$V6 + qnorm(0.975)*sumstat$V7)
sumstat$'OR_95_CI' <- paste0(format(round(sumstat$odds_ratio,2),nsmall=2)," [",format(round(sumstat$ci_lower,2),nsmall=2),"-",format(round(sumstat$ci_upper,2),nsmall=2),"]")

fwrite(sumstat,"output/68sentinel_allergyfree_sensitivity_OR", sep="\t", quote =F)

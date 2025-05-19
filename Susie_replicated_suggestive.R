#remotes::install_github("stephenslab/susieR")

args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(data.table)
library(susieR)

print("Reading data...")

data <- as.matrix(fread(args[1], header=FALSE))

print("Data read.")

print("Creating correlation matrix...")

data.cor = cor(data)

print("Correlation matrix complete.")

sumstat <- read.table(args[2])

print("Running Susie...")

#coverage default is 0.95:
#fitted_rss <- susie_rss(z = sumstat$V1, R = data.cor, L = 10, n=46086)
#lower coverage parameter:
#3_rs35570272_32547662_33547662: coverage: 50 for credset
#5_rs1837253_rs3806932_109901872_110905675: yes two different credsets with coverage 0.9
#12_rs73107993_47695873_48695873: 55 for credset
#15_rs11071559_60569988_61569988: no credset with coverage 0.80 (min)
#16_rs3024619_26864806_27864806: no credset with coverage 0.80 (min)
#17_17:38073838_CCG_C_37573838_38573838: no credset with coverage 0.80 (min)
print("coverage 0.60")
fitted_rss <- susie_rss(z = sumstat$V1, R = data.cor, L = 10, n=46086, coverage=0.95)


summary(fitted_rss)$cs

print(fitted_rss$pip)

write.table(summary(fitted_rss)$cs,paste0(args[3],"credset"),row.names=TRUE,quote=FALSE,col.names=TRUE,sep="\t")
write.table(fitted_rss$pip,args[3],row.names=TRUE,quote=FALSE,col.names=FALSE,sep="\t")

jpeg(args[4])

susie_plot(fitted_rss, y="PIP")

dev.off()
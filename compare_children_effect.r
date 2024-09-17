suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
library(ggpubr)

# children
signal1 <- fread("/data/gen1/LF_HRC_transethnic/effects_children/all_signal_effect_in_children.txt")
fev1_signal <- signal1[trait=="FEV1",]
fvc_signal <- signal1[trait=="FVC",]
pef_signal <- signal1[trait=="PEF",]
ff_signal <- signal1[trait=="FF",]

fev1 <- fread("/data/gen1/LF_HRC_transethnic/untransformed_original/UKB_EUR/fev1_res_EUR_all.bgen.stats.gz")
fev1 <- mutate(fev1,SNPid=paste(CHR,BP,pmin(toupper(ALLELE1),toupper(ALLELE0)),pmax(toupper(ALLELE1),toupper(ALLELE0)),sep="_"))
fev1_use <- fev1[SNPid%in%fev1_signal$MarkerName,]

fvc <- fread("/data/gen1/LF_HRC_transethnic/untransformed_original/UKB_EUR/fvc_res_EUR_all.bgen.stats.gz")
fvc <- mutate(fvc,SNPid=paste(CHR,BP,pmin(toupper(ALLELE1),toupper(ALLELE0)),pmax(toupper(ALLELE1),toupper(ALLELE0)),sep="_"))
fvc_use <- fvc[SNPid%in%fvc_signal$MarkerName,]

pef <- fread("/data/gen1/LF_HRC_transethnic/untransformed_original/UKB_EUR/pef_res_EUR_all.bgen.stats.gz")
pef <- mutate(pef,SNPid=paste(CHR,BP,pmin(toupper(ALLELE1),toupper(ALLELE0)),pmax(toupper(ALLELE1),toupper(ALLELE0)),sep="_"))
pef_use <- pef[SNPid%in%pef_signal$MarkerName,]

ratio <- fread("/data/gen1/LF_HRC_transethnic/untransformed_original/UKB_EUR/ratio_res_EUR_all.bgen.stats.gz")
ratio <- mutate(ratio,SNPid=paste(CHR,BP,pmin(toupper(ALLELE1),toupper(ALLELE0)),pmax(toupper(ALLELE1),toupper(ALLELE0)),sep="_"))
ratio_use <- ratio[SNPid%in%ff_signal$MarkerName,]

# adults - UKB
signal2 <- rbind(fev1_use,fvc_use,pef_use,ratio_use)
write.table(signal2,"/scratch/gen1/jc824/transethnic/effects_children/signal_effect_in_adults_UKB.txt",col.names=T, row.names=F, sep="\t", quote=F)
signal2 <- fread("/scratch/gen1/jc824/transethnic/effects_children/signal_effect_in_adults_UKB.txt")
signal2_sig <- signal2[P_BOLT_LMM<5e-8,]
# adults - nonUKB
signal3 <- fread("/scratch/gen1/jc824/transethnic/effects_children/signal_effect_in_adults_nonUKB.txt")
setnames(signal3,"P-value","P")
signal3_sig <- signal3[P<5e-8,] 
# discovery: signal2 UKB
SNP_children <- signal1[MarkerName%in%signal2_sig$SNPid,] %>%
                select(SNPid="MarkerName",,trait="trait",EA_children="Allele1",Beta_children="Effect",se_children="StdErr")
SNP_adults <- signal3[MarkerName%in%signal2_sig$SNPid,] %>%
              select(SNPid="MarkerName",EA_adults="Allele1",Beta_adults="Effect",se_adults="StdErr")
merged <- merge(SNP_children,SNP_adults,by="SNPid") %>%
          mutate(ea=ifelse(EA_children==EA_adults,TRUE,FALSE),z_stat=(Beta_children-Beta_adults)/sqrt(se_children^2+se_adults^2),p=2*pnorm(-abs(z_stat)))
### observed 90 (out of 701) signals with nominal evidence (p<0.05) of age-dependent effect compared to that expected by chance (binomial test p = 9.226e-16)
merged0.05 <- merged[p<0.05,]
BinomTest <- binom.test(nrow(merged0.05),nrow(merged),0.05)

sig <- merged[merged$p<0.05/756,]
write.table(sig,"/scratch/gen1/jc824/transethnic/effects_children/signals_Bonforroni_children_vs_adults_DiscoveryUKB.txt",col.names=T,row.names=F,sep="\t",quote=F)

png("/scratch/gen1/jc824/transethnic/effects_children/compare_children_effect_DiscoveryUKB.png",units="in", width=10, height=10, res=800)
ggplot(merged, aes(x=Beta_children, y=Beta_adults)) + geom_point(alpha=0.3) +
geom_point(data=sig,aes(x=Beta_children,y=Beta_adults),color='red',size=3) +
geom_smooth(method="lm") +
stat_cor(method="pearson") +
facet_wrap(~trait,nrow=2,scale="free") +
xlab("Effect size in children") +
ylab("Effect size in adults")
dev.off()

# discovery: signal3 nonUKB
SNP_children <- signal1[MarkerName%in%signal3_sig$MarkerName,] %>%
                select(SNPid="MarkerName",,trait="trait",EA_children="Allele1",Beta_children="Effect",se_children="StdErr")
SNP_adults <- signal2[SNPid%in%signal3_sig$MarkerName,] %>%
              select(SNPid="SNPid",EA_adults="ALLELE1",Beta_adults0="BETA",se_adults="SE")
merged <- merge(SNP_children,SNP_adults,by="SNPid") %>%
          mutate(ea=ifelse(toupper(EA_children)==toupper(EA_adults),TRUE,FALSE),Beta_adults=ifelse(ea==TRUE,Beta_adults0,-1*Beta_adults0),z_stat=(Beta_children-Beta_adults)/sqrt(se_children^2+se_adults^2),p=2*pnorm(-abs(z_stat)))
### observed 19 (out of 55) signals with nominal evidence (p<0.05) of age-dependent effect compared to that expected by chance (binomial test p = 9.32e-12)
merged0.05 <- merged[p<0.05,]
BinomTest <- binom.test(nrow(merged0.05),nrow(merged),0.05)

sig <- merged[merged$p<0.05/756,]
write.table(sig,"/scratch/gen1/jc824/transethnic/effects_children/signals_Bonforroni_children_vs_adults_DiscoveryNonUKB.txt",col.names=T,row.names=F,sep="\t",quote=F)

png("/scratch/gen1/jc824/transethnic/effects_children/compare_children_effect_DiscoveryNonUKB.png",units="in", width=10, height=10, res=800)
ggplot(merged, aes(x=Beta_children, y=Beta_adults)) + geom_point(alpha=0.3) +
geom_point(data=sig,aes(x=Beta_children,y=Beta_adults),color='red',size=3) +
geom_smooth(method="lm") +
stat_cor(method="pearson") +
facet_wrap(~trait,nrow=2,scale="free") +
xlab("Effect size in children") +
ylab("Effect size in adults")
dev.off()

######################################################################################
adults_fev1 <- fread("/scratch/gen1/nrgs1/transethnic/R2/FEV1_EUR_nonUKB_untransformed_meta.1tbl") %>%
               mutate(trait="FEV1")
adults_fvc <- fread("/scratch/gen1/nrgs1/transethnic/R2/FVC_EUR_nonUKB_untransformed_meta.1tbl") %>%
               mutate(trait="FVC")
adults_ff <- fread("/scratch/gen1/nrgs1/transethnic/R2/RATIO_EUR_nonUKB_untransformed_meta_noNEO.1tbl") %>%
               mutate(trait="FEV1/FVC")
adults_pef <- fread("/scratch/gen1/nrgs1/transethnic/R2/PEF_EUR_nonUKB_untransformed_meta.1tbl") %>%
               mutate(trait="PEF")


children_fev1 <- fread("/scratch/gen1/jc824/transethnic/effects_children/children_EUR_meta/FEV1_EUR_untransformed_meta.1tbl") 
children_fev1 <- children_fev1[MarkerName%in%adults_fev1$MarkerName,]
children_fvc <- fread("/scratch/gen1/jc824/transethnic/effects_children/children_EUR_meta/FVC_EUR_untransformed_meta.1tbl") 
children_fvc <- children_fvc[MarkerName%in%adults_fvc$MarkerName,]
children_ff <- fread("/scratch/gen1/jc824/transethnic/effects_children/children_EUR_meta/FF_EUR_untransformed_meta.1tbl") 
children_ff <- children_ff[MarkerName%in%adults_ff$MarkerName,]
children_pef <- fread("/scratch/gen1/jc824/transethnic/effects_children/children_EUR_meta/PEF_EUR_untransformed_meta.1tbl") 
children_pef <- children_pef[MarkerName%in%adults_pef$MarkerName,]

adults <- rbind(adults_fev1,adults_fvc,adults_ff,adults_pef) %>%
          mutate(ID=paste0(MarkerName,"_",trait)) %>%
          select(ID="ID",variant="MarkerName",trait="trait",EA_adults="Allele1",Beta_adults0="Effect",se_adults="StdErr")

children_fev1$trait <- "FEV1"
children_fvc$trait <- "FVC"
children_ff$trait <- "FEV1/FVC"
children_pef$trait <- "PEF"          
children <- rbind(children_fev1,children_fvc,children_ff,children_pef) %>%
            mutate(ID=paste0(MarkerName,"_",trait)) %>%
            select(ID="ID",EA_children="Allele1",Beta_children="Effect",se_children="StdErr")

merged <- merge(adults,children,by="ID") %>%
          mutate(ea=ifelse(toupper(EA_children)==toupper(EA_adults),TRUE,FALSE),Beta_adults=ifelse(ea==TRUE,Beta_adults0,-1*Beta_adults0),z_stat=(Beta_children-Beta_adults)/sqrt(se_children^2+se_adults^2),p=2*pnorm(-abs(z_stat)))
write.table(merged,"/scratch/gen1/jc824/transethnic/effects_children/allSNPs_children_vs_adults_DiscoveryUKB.txt",col.names=T,row.names=F,sep="\t",quote=F)
### observed 113 (out of 1077) signals with nominal evidence (p<0.05) of age-dependent effect compared to that expected by chance (binomial test p-value = 2.56e-13)
merged0.05 <- merged[p<0.05,]
BinomTest <- binom.test(nrow(merged0.05),nrow(merged),0.05,alternative="greater")

sig <- merged[merged$p<0.05/nrow(merged),]
write.table(sig,"/scratch/gen1/jc824/transethnic/effects_children/signals_Bonforroni_children_vs_adults_DiscoveryUKB.txt",col.names=T,row.names=F,sep="\t",quote=F)

png("/scratch/gen1/jc824/transethnic/effects_children/compare_children_effect_DiscoveryUKB.png",units="in", width=10, height=10, res=800)
ggplot(merged, aes(x=Beta_children, y=Beta_adults)) + geom_point(alpha=0.3) +
geom_point(data=sig,aes(x=Beta_children,y=Beta_adults),color='red',size=3) +
geom_smooth(method="lm") +
stat_cor(method="pearson") +
facet_wrap(~trait,nrow=2,scale="free") +
xlab("Effect size in children") +
ylab("Effect size in adults")
dev.off()
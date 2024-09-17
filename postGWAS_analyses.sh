#!/bin/bash

#SBATCH --output=/scratch/gen1/nnp5/sensitivity_comorb_tmp_data/logerror/%x-%j.out
#SBATCH --time=15:0:0
#SBATCH --mem=200gb
#SBATCH --account=gen1
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4

PATH_OUT="/home/n/nnp5/PhD/PhD_project/tmp_manuscript/output"
PHENO="broad_pheno_nocomob"

mkdir ${PATH_OUT}/allchr
GWAS="/home/n/nnp5/PhD/PhD_project/tmp_manuscript/output/allchr"

#merge all assoc file:
zcat ${PATH_OUT}/${PHENO}.1.regenie.step2_${PHENO}.regenie.gz | tail -n +2 | \
     awk -F " " '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' \
    > ${GWAS}/${PHENO}_allchr.assoc.txt

for i in {2..22}
    do zcat ${PATH_OUT}/${PHENO}.${i}.regenie.step2_${PHENO}.regenie.gz | \
    tail -n +2 | awk -F " " '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' \
    >> ${GWAS}/${PHENO}_allchr.assoc.txt
done


gzip ${GWAS}/${PHENO}_allchr.assoc.txt
cp ${GWAS}/${PHENO}_allchr.assoc.txt.gz /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/

#create input ldsc file both for all vars and maf >= 0.01.
module load R
chmod o+x src/create_input_munge_summary_stats.R
dos2unix src/create_input_munge_summary_stats.R
Rscript src/create_input_munge_summary_stats.R \
    ${GWAS}/${PHENO}_allchr.assoc.txt.gz \
    ${PHENO}

cp ${PATH_OUT}/maf001_broad_pheno_nocomob_betase_input_mungestat /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/

#LDSC interactively.
#The results are the same for set with all vars and fitlered maf 0.01., because LDSC uses only vars > 0.01. So run
#on the filtered set
cd /home/n/nnp5/software/ldsc
#to create the env:
#conda env create --file environment.yml
#source activate ldsc

tot_n=44598
PHENO="maf001_broad_pheno_nocomob"

awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' ${PATH_OUT}/${PHENO}_betase_input_mungestat \
    > ${PATH_OUT}/${PHENO}_betase_input_mungestat_clean

#interactively:
conda activate ldsc
/home/n/nnp5/software/ldsc/munge_sumstats.py \
--sumstats ${PATH_OUT}/${PHENO}_betase_input_mungestat_clean \
--N ${tot_n}   \
--out ${PATH_OUT}/${PHENO}_allchr_step2_regenie \
--merge-alleles /data/gen1/UKBiobank/Smoking/Meta_Analysis_AWI_UGR_UKBAFR/files/w_hm3.snplist \
--chunksize 500000

/home/n/nnp5/software/ldsc/ldsc.py \
--h2 ${PATH_OUT}/${PHENO}_allchr_step2_regenie.sumstats.gz \
--ref-ld-chr /data/gen1/reference/ldsc/eur_w_ld_chr/ \
--w-ld-chr /data/gen1/reference/ldsc/eur_w_ld_chr/ \
--out ${PATH_OUT}/${PHENO}_allchr_step2_regenie_h2

conda deactivate
cd /home/n/nnp5/PhD/PhD_project/tmp_manuscript/
cp ${PATH_OUT}/${PHENO}_allchr_step2_regenie_h2.log /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/

#Plots:
ldsc_intercept='1.02'

#run Manhattan, qqplot, lambda for vars with maf >= 0.01:
PHENO="maf001_broad_pheno_nocomob"
chmod o+x src/plot_functions.R
dos2unix src/plot_functions.R
chmod o+x src/REGENIE_plots.R
dos2unix src/REGENIE_plots.R
Rscript src/REGENIE_plots.R ${PATH_OUT}/${PHENO}_betase_input_mungestat ${PHENO} ${ldsc_intercept}


#Miami plot:
to find a new library as the 'hudson' do not work on R 4.3.1

#look at suggestive variants from discovery GWAS to see if they have similar effect size and same direction of effect:
PHENO="maf001_broad_pheno_nocomob"
awk '{print $1}' ${PATH_OUT}/maf001_pheno_1_5_ratio_sentinel_variants.txt | \
    grep -w -F -f - ${PATH_OUT}/${PHENO}_betase_input_mungestat > ${PATH_OUT}/${PHENO}_BTS45_sentinel_vars

#create OddsRatio and direction of effect using Create_oddsratio.R


#Compare broad_pheno suggestive variants with broad_pheno_nocomob: direction of effect
awk '{print $1}' /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/UKBiobank_severeasthma_sentinel_suggestive_to_replicate.txt | \
    grep -F -f - -w /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/maf001_broad_pheno_nocomob_betase_input_mungestat \
    > ${PATH_OUT}/maf001_broad_pheno_nocomob_for_sensitivity

chmod o+x src/sensitivity_check.R
dos2unix src/sensitivity_check.R
Rscript src/sensitivity_check.R /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/Disc_VS_Meta_sumstat_6cohorts.txt \
    ${PATH_OUT}/maf001_broad_pheno_nocomob_for_sensitivity \
    ${PATH_OUT}/sensitivity_nocomob_suggestive_discovery.txt

cp ${PATH_OUT}/sensitivity_nocomob_suggestive_discovery.txt /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/
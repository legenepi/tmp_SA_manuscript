#!/bin/bash

#SBATCH --job-name=regenie1
#SBATCH --output=/scratch/gen1/nnp5/sensitivity_comorb_tmp_data/logerror/%x-%j.out
#SBATCH --time=15:0:0
#SBATCH --mem=200gb
#SBATCH --account=gen1
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4


PATH_DATA="/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data"
work_dir="/home/n/nnp5/PhD/PhD_project/tmp_manuscript"
geno_dir="/data/ukb/genotyped"
scratch_dir="/scratch/gen1/nnp5/sensitivity_comorb_tmp_data"

OUT_DIR="/home/n/nnp5/PhD/PhD_project/tmp_manuscript/output"
pheno="broad_pheno_nocomob"

/home/n/nnp5/software/regenie_v3/regenie_v3.2.1.gz_x86_64_Centos7_mkl \
  --step 1 \
  --bed ${scratch_dir}/ukb_cal_allchr_v2 \
  --extract ${scratch_dir}/ukb_cal_allchr_eur_qc.snplist \
  --keep ${scratch_dir}/ukb_cal_allchr_eur_qc.id \
  --phenoFile ${PATH_DATA}/demo_EUR_pheno_cov_broadasthma_nocomob.txt \
  --phenoCol ${pheno} \
  --covarFile ${PATH_DATA}/demo_EUR_pheno_cov_broadasthma_nocomob.txt \
  --covarColList age_at_recruitment,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,genetic_sex \
  --bt \
  --bsize 1000 \
  --loocv \
  --threads 4 \
  --gz \
  --out ${OUT_DIR}/${pheno}.regenie.step1
#!/bin/bash

#SBATCH --output=/scratch/gen1/nnp5/sensitivity_comorb_tmp_data/logerror/%x-%j.out
#SBATCH --time=15:0:0
#SBATCH --mem=200gb
#SBATCH --account=gen1
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-22

#run as: sbatch --array=1-22 src/run_regenie_step2.sh
chr=$SLURM_ARRAY_TASK_ID
PATH_DATA="/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data"
work_dir="/home/n/nnp5/PhD/PhD_project/tmp_manuscript"
geno_dir="/data/ukb/genotyped"
scratch_dir="/scratch/gen1/nnp5/sensitivity_comorb_tmp_data"
sample_DIR="/data/gen1/UKBiobank_500K/severe_asthma/data"
OUT_DIR="/home/n/nnp5/PhD/PhD_project/tmp_manuscript/output"
pheno="broad_pheno_nocomob"

/home/n/nnp5/software/regenie_v3/regenie_v3.2.1.gz_x86_64_Centos7_mkl \
  --step 2 \
  --bgen /data/ukb/imputed_v3/ukb_imp_chr${chr}_v3.bgen \
  --ref-first \
  --sample ${sample_DIR}/ukbiobank_app56607_for_regenie.sample \
  --keep ${scratch_dir}/ukb_cal_allchr_eur_qc.id \
  --phenoFile ${PATH_DATA}/demo_EUR_pheno_cov_broadasthma_nocomob.txt \
  --phenoCol ${pheno} \
  --covarFile ${PATH_DATA}/demo_EUR_pheno_cov_broadasthma_nocomob.txt \
  --covarColList age_at_recruitment,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,genetic_sex \
  --bt \
  --gz \
  --threads 4 \
  --minMAC 10 \
  --minINFO 0.3 \
  --firth --approx --pThresh 0.01 \
  --pred ${OUT_DIR}/${pheno}.regenie.step1_pred.list \
  --bsize 1000 \
  --out ${OUT_DIR}/${pheno}.${chr}.regenie.step2
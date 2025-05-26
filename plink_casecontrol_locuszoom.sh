#!/bin/bash

#Rationale: create plink file for case/control cohort (as a job array) in plink2 for LocusZoom plots

#SBATCH --job-name=bgenix
#SBATCH --output=/scratch/gen1/nnp5/tmp_manuscript/%x-%j.out
#SBATCH --time=72:0:0
#SBATCH --mem=100gb
#SBATCH --account=gen1
#SBATCH --export=NONE

#run as: sbatch --array=2,5,6,8,9,10,11,12,15 src/plink_casecontrol_locuszoom.sh
#run as: sbatch --array=3,16,17 src/plink_casecontrol_locuszoom.sh
#dos2unix src/plink_casecontrol_locuszoom.sh

chr=$SLURM_ARRAY_TASK_ID

##A.Create plink file for case/control cohort (as a job array) in plink2
module load plink2
plink2 \
  --bgen /data/ukb/nobackup/imputed_v3/ukb_imp_chr${chr}_v3.bgen ref-first \
  --keep /home/n/nnp5/PhD/PhD_project/Post_GWAS/input/broadasthma_individuals \
  --make-bed \
  --out /scratch/gen1/nnp5/tmp_manuscript/broad_pheno_plink_file_v3_chr${chr} \
  --sample /data/gen1/UKBiobank_500K/severe_asthma/data/ukbiobank_app56607_for_regenie.sample

##B.Calculate ld in plink v1.9b:
#rsid="rs12470864"
#rsid='rs6761047'
#rsid='rs10455025'
#rsid='rs1837253'
#rsid='rs1986009'
#rsid='rs2070729'
#rsid='rs148639908'
#rsid='rs7824394'
#rsid='rs992969'
#rsid='rs12413578'
#rsid='rs1444789'
#rsid='rs10160518'
#rsid='rs705705'
#rsid='rs3024971'
#rsid='rs17293632'

#module unload plink2
#module load plink
#plink \
#    --bfile /scratch/gen1/nnp5/tmp_manuscript/broad_pheno_plink_file_v3_chr${chr} \
#    --r2 \
#    --ld-snp ${rsid} \
#    --ld-window-kb 1000 \
#    --ld-window 99999 \
#    --ld-window-r2 0 \
#    --out /scratch/gen1/nnp5/tmp_manuscript/ld_chr${chr}_${rsid}
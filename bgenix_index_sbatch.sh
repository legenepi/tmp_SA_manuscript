#!/bin/bash

#SBATCH --job-name=bgenix
#SBATCH --output=/scratch/gen1/nnp5/Fine_mapping/%x-%j.out
#SBATCH --time=72:0:0
#SBATCH --mem=50gb
#SBATCH --account=gen1
#SBATCH --export=NONE

#run as: sbatch --array=3,5,12,15,16,17 src/bgenix_index_sbatch.sh
#dos2unix src/bgenix_index_sbatch.sh

chr=$SLURM_ARRAY_TASK_ID

cd /scratch/gen1/nnp5/Fine_mapping/tmp_data/

module load qctool
qctool -g /data/ukb/imputed_v3/ukb_imp_chr${chr}_v3.bgen \
    -s /data/gen1/UKBiobank_500K/severe_asthma/data/ukbiobank_app56607_for_regenie.sample \
    -incl-samples /home/n/nnp5/PhD/PhD_project/Post_GWAS/input/broadasthma_individuals \
    -og /scratch/gen1/nnp5/Fine_mapping/tmp_data/sevasthma_chr${chr}_v3.bgen

/data/gen1/bin/bgenix -g /scratch/gen1/nnp5/Fine_mapping/tmp_data/sevasthma_chr${chr}_v3.bgen -index
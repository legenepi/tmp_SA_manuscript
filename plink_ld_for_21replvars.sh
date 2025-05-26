#!/bin/bash

#Rationale: create in-sample ld measures with respect to the 21 replicated variants for LocusZoom Plots

#SBATCH --job-name=ld_measure_21replvrs
#SBATCH --output=/scratch/gen1/nnp5/tmp_manuscript/%x-%j.out
#SBATCH --time=72:0:0
#SBATCH --mem=100gb
#SBATCH --account=gen1
#SBATCH --export=NONE

#run as: sbatch src/plink_ld_for_21replvars.sh
#dos2unix src/plink_ld_for_21replvars.sh
#dos2unix src/21replvars_rsid_chr.txt

#crate file with rsid and chr
##B.Calculate ld in plink v1.9b:
module unload plink2
module load plink

file="/home/n/nnp5/PhD/PhD_project/tmp_manuscript/21replvars_rsid_chr.txt"
while IFS= read -r line; do
  rsid=$(echo "$line" | awk -F "\t" '{print $1}')
  chr=$(echo "$line" | awk -F "\t" '{print $2}')
  echo $rsid
  echo $chr
  plink \
    --bfile /scratch/gen1/nnp5/tmp_manuscript/broad_pheno_plink_file_v3_chr${chr} \
    --r2 \
    --ld-snp ${rsid} \
    --ld-window-kb 1000 \
    --ld-window 99999 \
    --ld-window-r2 0 \
    --out /scratch/gen1/nnp5/tmp_manuscript/ld_chr${chr}_${rsid}
done < "$file"
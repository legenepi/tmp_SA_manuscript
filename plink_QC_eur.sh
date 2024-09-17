#!/bin/bash

#SBATCH --job-name=plink_QC_eur
#SBATCH --output=/scratch/gen1/nnp5/sensitivity_comorb_tmp_data/logerror/%x-%j.out
#SBATCH --time=1:0:0
#SBATCH --mem=500gb
#SBATCH --account=gen1
#SBATCH --export=NONE


PATH_DATA="/data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data"
work_dir="/home/n/nnp5/PhD/PhD_project/tmp_manuscript"
geno_dir="/data/ukb/genotyped"
scratch_dir="/scratch/gen1/nnp5/sensitivity_comorb_tmp_data"

#mkdir /scratch/gen1/nnp5/sensitivity_comorb_tmp_data

#to create list of eur ids
#awk {'print $1, $2'} ${PATH_DATA}/demo_EUR_pheno_cov_broadasthma.txt | tail -n +2 > ${scratch_dir}/ukb_eur_ids

#cp from /rfs to my home the .fam
#cp /rfs/TobinGroup/data/UKBiobank/application_56607/ukb56607_cal_chr1_v2_s488239.fam ${work_dir}/

#rm ${work_dir}/list_plink_files
#for chr in {2..22}
#do echo "${geno_dir}/ukb_cal_chr${chr}_v2.bed \
#     ${geno_dir}/ukb_cal_chr${chr}_v2.bim \
#     ${work_dir}/ukb56607_cal_chr1_v2_s488239.fam" >> ${work_dir}/list_plink_files
#done

#module load gcc/12.3.0-yxgv2bl
#module load plink

#plink --bed ${geno_dir}/ukb_cal_chr1_v2.bed \
#	--bim ${geno_dir}/ukb_cal_chr1_v2.bim \
#	--fam ${work_dir}/ukb56607_cal_chr1_v2_s488239.fam \
#	--merge-list ${work_dir}/list_plink_files \
#	--make-bed --out ${scratch_dir}/ukb_cal_allchr_v2

#&&

awk '{print $1, $1}' /home/n/nnp5/UKBiobank_datafields/data/Eid_withdrawn_participants_upFeb2022.txt \
    > ${scratch_dir}/ukb_withdrawns_ids &&

module load plink2
plink2 --bfile ${scratch_dir}/ukb_cal_allchr_v2 \
	--maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
	--keep ${scratch_dir}/ukb_eur_ids \
	--remove ${scratch_dir}/ukb_withdrawns_ids \
	--mind 0.1 \
	--write-snplist --write-samples --no-id-header \
	--out ${scratch_dir}/ukb_cal_allchr_eur_qc &&

#cp ${scratch_dir}/ukb_cal_allchr_eur_qc.id ${work_dir}


#Updated fine-mapping analysis using Susie only for 6 signals
#!/bin/bash

#PBS -N SuSie_replicated_suggestive
#PBS -j oe
#PBS -o SuSie_replicated_suggestive
#PBS -l walltime=20:0:0
#PBS -l vmem=50gb
#PBS -l nodes=1:ppn=4
#PBS -d .
#PBS -W umask=022

##create bgen and index file:
##data.bgen and data.bgen.bgi (chr of interest are: 3,5,12,15,16,17):
#sbatch --array=3,5,12,15,16,17 src/bgenix_index_sbatch.sh

PATH_finemapping="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma"
PATH_OUT="/home/n/nnp5/PhD/PhD_project/Fine_mapping_severe_asthma/output/susie_replicated_suggestive_Jan2025"
PATH_ASSOC="/home/n/nnp5/PhD/PhD_project/REGENIE_assoc/output/allchr"
##set working directory:
cd ${PATH_finemapping}

##load required tools:
module load gcc
module load R
module load plink2

#take variants to be re-ran:
grep "rs73107993\|rs11071559\|rs3024619\|17:38073838_CCG_C\|rs35570272\|rs1837253_rs3806932" \
    ${PATH_finemapping}/input/fine_mapping_regions_replicated_suggestive_input > \
    /home/n/nnp5/PhD/PhD_project/tmp_manuscript/fine_mapping_regions_replicated_suggestive_input_Jan2025

##May 2025:
#rs2188962_rs848 and rs1444789_rs201499805  loci to be re-run:
grep "rs2188962_rs848\|rs201499805_rs1444789" \
    ${PATH_finemapping}/input/fine_mapping_regions_replicated_suggestive_input > \
    /home/n/nnp5/PhD/PhD_project/tmp_manuscript/fine_mapping_regions_replicated_suggestive_input_May2025

##output folder:
mkdir ${PATH_finemapping}/output/susie_replicated_suggestive_Jan2025
mkdir ${PATH_OUT}

##Analysis:
for line in {1..6}
do
#line=6
##Input data:
SNP=$(awk -v row="$line" ' NR == row {print $1 }' /home/n/nnp5/PhD/PhD_project/tmp_manuscript/fine_mapping_regions_replicated_suggestive_input_Jan2025)
chr=$(awk -v row="$line" ' NR == row {print $2 }' /home/n/nnp5/PhD/PhD_project/tmp_manuscript/fine_mapping_regions_replicated_suggestive_input_Jan2025)
start=$(awk -v row="$line" 'NR == row {print $4}' /home/n/nnp5/PhD/PhD_project/tmp_manuscript/fine_mapping_regions_replicated_suggestive_input_Jan2025)
end=$(awk -v row="$line" 'NR == row {print $5}' /home/n/nnp5/PhD/PhD_project/tmp_manuscript/fine_mapping_regions_replicated_suggestive_input_Jan2025)

for line in {1..2}
do
line=1
##Input data:
SNP=$(awk -v row="$line" ' NR == row {print $1 }' /home/n/nnp5/PhD/PhD_project/tmp_manuscript/fine_mapping_regions_replicated_suggestive_input_May2025)
chr=$(awk -v row="$line" ' NR == row {print $2 }' /home/n/nnp5/PhD/PhD_project/tmp_manuscript/fine_mapping_regions_replicated_suggestive_input_May2025)
start=$(awk -v row="$line" 'NR == row {print $4}' /home/n/nnp5/PhD/PhD_project/tmp_manuscript/fine_mapping_regions_replicated_suggestive_input_May2025)
end=$(awk -v row="$line" 'NR == row {print $5}' /home/n/nnp5/PhD/PhD_project/tmp_manuscript/fine_mapping_regions_replicated_suggestive_input_May2025)


##Creating region bgen
#if chr has double digit, I do not need the '0' in the chromosome name, so I need two different string for the
#-incl-range argument:
#mkdir /scratch/gen1/nnp5/Fine_mapping #if it does not exist already
#mkdir /scratch/gen1/nnp5/Fine_mapping/tmp_data #it it does nto exist already
cd /scratch/gen1/nnp5/Fine_mapping/tmp_data/
if [[ ${chr} -lt 10 ]]
then
/data/gen1/bin/bgenix -g /scratch/gen1/nnp5/Fine_mapping/tmp_data/sevasthma_chr${chr}_v3.bgen \
    -incl-range 0${chr}:${start}-${end} \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.bgen
fi

if [[ ${chr} -gt 9 ]]
then
  /data/gen1/bin/bgenix -g /scratch/gen1/nnp5/Fine_mapping/tmp_data/sevasthma_chr${chr}_v3.bgen \
    -incl-range ${chr}:${start}-${end} \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.bgen
fi

cd ${PATH_finemapping}

##Exclude multi-allelic variants and find the common SNP IDs for the genotyped matrix and the zscore input files:
#use the file for each regions created by FINEMAP.sh:
grep -v -w -F -f /data/gen1/UKBiobank_500K/imputed/multiallelic.snps \
    ${PATH_finemapping}/input/ldstore_chr${chr}_${SNP}.z \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt

awk 'NR > 1 {print $1}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_snps.txt

##zscore:
#dos2unix /home/n/nnp5/PhD/PhD_project/tmp_manuscript/src/z_score.R
Rscript src/z_score.R /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_z_scores.txt


##Format region data for input to R
plink2 \
    --bgen /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.bgen ref-first \
    --sample ${PATH_finemapping}/input/ldstore.sample \
    --export A \
    --extract /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_snps.txt \
    --out /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}

cut -f7- /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.raw \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols.raw

awk 'NR>1 {print}' /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols.raw \
    > /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols_no_header.raw

Rscript /home/n/nnp5/PhD/PhD_project/tmp_manuscript/src/Susie_replicated_suggestive.R \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}.cols_no_header.raw \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_z_scores.txt \
    ${PATH_OUT}/susie_replsugg_${chr}_${SNP}_${start}_${end}.txt \
    ${PATH_OUT}/susie_replsugg_${chr}_${SNP}_${start}_${end}.jpeg

awk -F "\t" '{print $6}' ${PATH_OUT}/susie_replsugg_${chr}_${SNP}_${start}_${end}.txtcredset | \
    tr , '\n' | tail -n +2 > ${PATH_OUT}/susie_replsugg_${chr}_${SNP}_${start}_${end}.credset.indx

Rscript /home/n/nnp5/PhD/PhD_project/tmp_manuscript/src/credset_susie.R \
    ${PATH_OUT}/susie_replsugg_${chr}_${SNP}_${start}_${end}.txt\
    ${PATH_OUT}/susie_replsugg_${chr}_${SNP}_${start}_${end}.credset.indx \
    /scratch/gen1/nnp5/Fine_mapping/tmp_data/${SNP}_no_ma_GWAS_sumstats.txt \
    ${chr}_${SNP}_${start}_${end} \
    ${PATH_OUT}/susie_replsugg_credset.${SNP}.$chr.$start.$end
done

#Check how many credible sets were identified - if more than one, noted down:

#Merge credset into a unique file:
head -n 1 ${PATH_OUT}/susie_replsugg_credset.rs35570272.3.32547662.33547662 \
    > ${PATH_finemapping}/output/susie_replicated_suggestive_Jan2025/susie_replsugg_all_credset_additional_March2025.txt && \
    tail -n +2 -q ${PATH_OUT}/susie_replsugg_credset.*[0-9] \
    >> ${PATH_finemapping}/output/susie_replsugg_all_credset_additional_March2025.txt

#Create src/report/Finemapping_repl_sugg_lowcov_March2025.xlsx with the ${PATH_finemapping}/output/susie_replsugg_all_credset_additional_March2025.txt
#Copy src/report/Finemapping_repl_sugg_lowcov_March2025.xlsx into /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/

##Remove intermediate files:
#/scratch/gen1/nnp5/Fine_mapping/tmp_data/*_no_ma_GWAS_sumstats.txt
#/scratch/gen1/nnp5/Fine_mapping/tmp_data/*.cols_no_header.raw
#/scratch/gen1/nnp5/Fine_mapping/tmp_data/*.cols.raw
#/scratch/gen1/nnp5/Fine_mapping/tmp_data/*.raw
#/scratch/gen1/nnp5/Fine_mapping/tmp_data/*_no_ma_snps.txt
#/scratch/gen1/nnp5/Fine_mapping/tmp_data/*_no_ma_z_scores.txt
#${PATH_OUT}/susie_replsugg_*.credset.
#${PATH_OUT}/susie_replsugg_*.txtcredset
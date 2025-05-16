#!/bin/env/bash

#Rationale: Make LocusZoom Plot for fine-mapped loci and credible sets for the manuscript
#NB: if multiple credsets, do multiple LocusZoom Plot and highlight the different credset in each.

module load R
module load plink
#download LocusZoom:
wget https://statgen.sph.umich.edu/locuszoom/download/locuszoom_1.4_srconly.tgz --no-check-certificate


PATH_OUT="/home/n/nnp5/PhD/PhD_project/tmp_manuscript/output"
PHENO="maf001_broad_pheno_1_5_ratio"

mkdir ${PATH_OUT}/Locuszoom_builtin
mkdir /scratch/gen1/nnp5/tmp_manuscript

#input: change i and rsid for each credible set:
2_rs12470864_102426362_103426362	rs12470864	2	102309902	IL1RL1
2_rs6761047_242192858_243192858	rs6761047	2	241753443	D2HGDH
5_rs1837253_rs3806932_109901872_110905675	rs10455025	5	111069301	TSLP
5_rs1837253_rs3806932_109901872_110905675	rs1837253	5	111066174	TSLP
5_rs2188962_rs848_131270805_132496500	rs1986009	5	132552294	IL5
5_rs2188962_rs848_131270805_132496500	rs2070729	5	132484229	C5orf56
6_rs148639908_90463614_91463614	rs148639908	6	90253895	BACH2
8_rs7824394_80792599_81792599	rs7824394	8	80380364	AC034114.2

9_rs992969_5709697_6709697	rs992969	9	6209697	IL33 # to be done

10_rs201499805_rs1444789_8542744_9564361	rs12413578	10	9007290	AC044784.1
10_rs201499805_rs1444789_8542744_9564361	rs1444789	10	9022398	AC044784.1
11_rs10160518_75796671_76796671	rs10160518	11	76585627	AP001189.5
12_rs705705_55935504_56935504	rs705705	12	56041720	AC034102.4
12_rs3024971_56993727_57993727	rs3024971	12	57099944	STAT6
15_rs17293632_66942596_67942596	rs17293632	15	67150258	SMAD3

####for webtool implementation with LD with 1000GP-EUR:
head -n 1 /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/${PHENO}_betase_input_mungestat \
    >  /scratch/gen1/nnp5/tmp_manuscript/header
#chr2:
awk '$2==2 {print $0}' /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/${PHENO}_betase_input_mungestat \
    > /scratch/gen1/nnp5/tmp_manuscript/chr2_input_mungestat
cat /scratch/gen1/nnp5/tmp_manuscript/header /scratch/gen1/nnp5/tmp_manuscript/chr2_input_mungestat \
    > /scratch/gen1/nnp5/tmp_manuscript/header_chr2_input_mungestat
####end of webtool implementation

####built-in implementation with LD with case/control cohort and credset highlighted:
##i, rsid, locus to change for each credset:
i=9
rsid="rs992969"
locus="9_rs992969_5709697_6709697"

#LD file with respect to highest PIP variant for each credset:
#change the rsid in the .sh file:
sbatch --array=2,5,6,8,9,10,11,12,15 src/plink_casecontrol_locuszoom.sh

#the User-supplied LD should have columns:
#snp1	Any SNP in your plotting region.
#snp2	Should always be the reference SNP in the region.
#dprime	D' between snp2 (reference SNP) and snp1.
#rsquare	r2 between snp2 (reference SNP) and snp1.
#The dprime column can be all missing if it is not known. Rsquare must be present, and must be valid data.
#The file should be whitespace delimited, and the header (column names shown above) must exist.
echo "snp1 snp2 dprime rsquare" > /scratch/gen1/nnp5/tmp_manuscript/header_ld
awk '{print $6, $3, "NA", $7}' /scratch/gen1/nnp5/tmp_manuscript/ld_chr${i}_${rsid}.ld | \
    tail -n +2 | \
    cat /scratch/gen1/nnp5/tmp_manuscript/header_ld - \
    > /scratch/gen1/nnp5/tmp_manuscript/ld_chr${i}_${rsid}_locuszoom

#GWAS summary stats file (--metal option):
#header file, for all:
head -n 1 /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/${PHENO}_betase_input_mungestat \
    >  /scratch/gen1/nnp5/tmp_manuscript/header
#for each chr:
awk -v x=$i '$2==x {print $0}' /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/${PHENO}_betase_input_mungestat \
    > /scratch/gen1/nnp5/tmp_manuscript/chr${i}_input_mungestat
cat /scratch/gen1/nnp5/tmp_manuscript/header /scratch/gen1/nnp5/tmp_manuscript/chr${i}_input_mungestat \
    > /scratch/gen1/nnp5/tmp_manuscript/header_chr${i}_input_mungestat

#CREDSET ANNOTATION: add the ANNOT column to the metal file:
dos2unix src/credsetannot_metal.R
chmod u+x src/credsetannot_metal.R
#if credset number 1:
Rscript src/credsetannot_metal.R \
    "/scratch/gen1/nnp5/tmp_manuscript/header_chr${i}_input_mungestat" \
    ${locus} \
    1 \
    ${rsid}
#if credset number 2:
Rscript src/credsetannot_metal.R \
    "/scratch/gen1/nnp5/tmp_manuscript/header_chr${i}_input_mungestat" \
    ${locus} \
    2 \
    ${rsid}


#LOCUS ZOOM PLOT:
#7.3: -log10(5e-8) 5.3: -log10(5e-6)
#I need python 2.7:
#conda create --name env_name python=2.7
conda activate env_name
/data/gen1/reference/locuszoom-standalone/bin/locuszoom \
    --metal /scratch/gen1/nnp5/tmp_manuscript/header_${locus}_${rsid}_annot \
    --refsnp ${rsid} --flank 1Mb \
    --plotonly \
    signifLine="7.3,5.3" \
    --ld /scratch/gen1/nnp5/tmp_manuscript/ld_chr${i}_${rsid}_locuszoom \
    --delim tab --pvalcol pval --markercol snpid \
    --build hg19 \
    --verbose \
    showAnnot=TRUE \
    annotCol='ANNOT' \
    annotPch='1,21'

#Nick's script is here: /scratch/gen1/nrgs1/adonyo/region_plots/region_plot.sh

#For locuszoom plot of 21 replicated variants:
for i in {2,3,5,6,8,9,10,11,12,15,16,17}
do
  awk -v x=$i '$2==x {print $0}' /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/${PHENO}_betase_input_mungestat \
    > /scratch/gen1/nnp5/tmp_manuscript/chr${i}_input_mungestat
  cat /scratch/gen1/nnp5/tmp_manuscript/header /scratch/gen1/nnp5/tmp_manuscript/chr${i}_input_mungestat \
    > /home/n/nnp5/PhD/PhD_project/tmp_manuscript/header_chr${i}_input_mungestat
done
#And then I used the webtool:
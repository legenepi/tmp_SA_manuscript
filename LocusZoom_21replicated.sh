#!/bin/env/bash

#Rationale: Make LocusZoom Plot for the 21 replicated variants


module load R
module load plink

mkdir /home/n/nnp5/PhD/PhD_project/tmp_manuscript/output/21replvars_locuszoom
PHENO="maf001_broad_pheno_1_5_ratio"


#input: change i and rsid for each credible set:
#rs12470864	2 - done
#rs6761047	2 - done
#rs35570272	 - done
#rs1837253	5 - done - weak LD
#rs3806932	5 - done
#rs2188962	5 - done
#rs848	5 - done
#rs9271365	6 - done
#rs148639908	6 - done
#rs7824394	8 - done
#rs992969	9 - done
#rs201499805	10 - done
#rs1444789	10 - done
#rs10160518	11 - done
#rs73107993	12 - done
#rs705705	12 - done
#rs3024971	12 - done
#rs11071559	15 - done
#rs17293632	15 - done
#rs3024619	16
#17:38073838_CCG_C	17 - do the LD manually


#bgen files:
#sbatch --array=2,5,6,8,9,10,11,12,15 src/plink_casecontrol_locuszoom.sh
#sbatch --array=3,16,17 src/plink_casecontrol_locuszoom.sh

####built-in implementation with LD with case/control cohort and credset highlighted:
##i, rsid, locus to change for each credset:
i=17
rsid="17:38073838_CCG_C"

#LD with respect to the 21 replicated variant:
sbatch src/plink_ld_for_21replvars.sh
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

#GWAS summary stats file (input file for --metal option):
#header file, for all:
head -n 1 /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/${PHENO}_betase_input_mungestat \
    >  /scratch/gen1/nnp5/tmp_manuscript/header
#for each chr:
awk -v x=$i '$2==x {print $0}' /data/gen1/UKBiobank_500K/severe_asthma/Noemi_PhD/data/${PHENO}_betase_input_mungestat \
    > /scratch/gen1/nnp5/tmp_manuscript/chr${i}_input_mungestat
cat /scratch/gen1/nnp5/tmp_manuscript/header /scratch/gen1/nnp5/tmp_manuscript/chr${i}_input_mungestat \
    > /scratch/gen1/nnp5/tmp_manuscript/header_chr${i}_input_mungestat



#LOCUS ZOOM PLOT:
#7.3: -log10(5e-8) 5.3: -log10(5e-6)
#I need python 2.7:
#conda create --name env_name python=2.7
conda activate env_name

/data/gen1/reference/locuszoom-standalone/bin/locuszoom \
    --metal /scratch/gen1/nnp5/tmp_manuscript/header_chr${i}_input_mungestat \
    --refsnp ${rsid} --flank 1Mb \
    --plotonly \
    signifLine="7.3,5.3" \
    --ld /scratch/gen1/nnp5/tmp_manuscript/ld_chr${i}_${rsid}_locuszoom \
    --delim tab --pvalcol pval --markercol snpid \
    --build hg19 \
    --verbose \
    --prefix /home/n/nnp5/PhD/PhD_project/tmp_manuscript/output/21replvars_locuszoom/
#!/bin/sh

#module load python/intel/27
conda init
conda activate py27

export PATH=${PATH}:/data/gen1/reference/locuszoom-standalone/bin/


METAL=$1
PREFIX=$2
LD=$3
#DB=/data/gen1/LF_HRC_transethnic/locuszoom/UKB_HRC.db
DB=/data/gen1/AIRPROM/assoc/region_plots/UK10K+1000G_Phase3.db
CHR=$4
START=$5
END=$6
SNP=$7
SUGGESTIVE=$8
HISTART=$9
HIEND=${10}
TITLE=${11}
#ARGS="--metal $METAL --chr $CHR --start $START --end $END --refsnp $SNP --prefix $PREFIX --build hg19 --ld $LD --db $DB --plotonly --no-date"
ARGS="--metal $METAL --chr $CHR --start $START --end $END --refsnp $SNP --prefix $PREFIX --build hg19 --ld $LD --plotonly --no-date"
#[ $# -gt 10 ] && ARGS="$ARGS --denote-markers-file ${11}"

locuszoom $ARGS \
    ymax=4 \
    ldcuts="0,.2,.5,.8,1" \
    ldColors="white,gray50,yellow,orange,red,blue" \
    recombColor=lightskyblue \
    signifLine="${SUGGESTIVE},8.3" \
    signifLineColor="green,red" \
    rfrows=8 \
    geneFontSize=.7 \
    showPartialGenes=TRUE \
    refsnpTextSize=0.7 \
    title="${TITLE}" \
    titleCex=0.8 \
    titleFontFace="italic"
#    legend='none' \
#    width=7 height=7 \
#    showAnnot=FALSE \
#    hiStart=$HISTART \
#    hiEnd=$HIEND \
#    hiColor=green \
#    hiAlpha=0.4 \
#    annotCol='ANNOT' \
#    annotPch='1,21' \

#PBS -S /bin/bash
#PBS -q batch
#PBS -N csv
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=480:00:00
#PBS -l mem=80gb

#above is the sapelo2 que parameters
#this changes the output to current working directory
cd $PBS_O_WORKDIR

#this loads the required modules on the cluster
module load SAMtools/1.6-foss-2016b
module load FASTX-Toolkit/0.0.14-foss-2016b

#Set reference genome, input file path, and output file path
GENOME=/work/cemlab/reference_genomes/97103_v2.fa
INPUT=/home/lfa81121/primerizer/testsnps.csv
OUTPUT=/home/lfa81121/primerizer/potentialprimers.csv

#Stores the current value for $IFS then changes the internal field seperator to a comma in order to read .csv columns.
OLDIFS=$IFS
IFS=,

#Checks if $INPUT is a valid file.
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }

#Creates the output file and writes the headers into the first row.
echo CHROM,POS,REF,ALT,AD_REFLOW,AD_ALTLOW,DPLOW,GQLOW,PLLOW,SNPINDEXLOW,AD_REFHIGH,AD_ALTHIGH,DPHIGH,GQHIGH,PLHIGH,SNPINDEXHIGH,REF_FRQ,DELTASNP,NSNPS,TRICUBEDELTASNP,G,GPRIME,PVALUE,NEGLOG10PVAL,QVALUE,MINDP,TRICUBEDP,CL95,CL99,POSLESS250,POSPLUS250,SNP1,SNP501,SNP27,SNP24,RSNP27,RSNP24 > $OUTPUT

#Reads each line of the .csv and calculates target positions in geneome. Target positions are extracted from $GENOME and printed into $OUTPUT along with all columns from $INPUT.
# Note that we must change $IFS back to its original value in order to output to the .csv, and then return it to "," again to restart the loop.
while read CHROM POS REF ALT AD_REFLOW AD_ALTLOW DPLOW GQLOW PLLOW SNPINDEXLOW AD_REFHIGH AD_ALTHIGH DPHIGH GQHIGH PLHIGH SNPINDEXHIGH REF_FRQ DELTASNP NSNPS TRICUBEDELTASNP G GPRIME PVALUE NEGLOG10PVAL QVALUE MINDP TRICUBEDP CL95 CL99
  do
    IFS=$OLDIFS
    SNP1="$(samtools faidx $GENOME $CHROM":"$POS"-"$POS | awk '{if(NR>1)print}')"
    POSPLUS250=$((POS+250))
    POSLESS250=$((POS-250))
    SNP501="$(samtools faidx $GENOME $CHROM":"$POSLESS250"-"$POSPLUS250 | awk '{if(NR>1)print}')"
    POSLESS23=$((POS-23))
    POSLESS26=$((POS-26))
    POSPLUS23=$((POS+23))
    POSPLUS26=$((POS+26))
    SNP27="$(samtools faidx $GENOME $CHROM":"$POSLESS26"-"$POS | awk '{if(NR>1)print}')"
    SNP24="$(samtools faidx $GENOME $CHROM":"$POSLESS23"-"$POS | awk '{if(NR>1)print}')"
    RSNP27="$(samtools faidx $GENOME $CHROM":"$POS"-"$POSPLUS26 | awk '{if(NR>1)print}' | tr 'gtca' 'cagt' | rev)"
    RSNP24="$(samtools faidx $GENOME $CHROM":"$POS"-"$POSPLUS23 | awk '{if(NR>1)print}' | tr 'gtca' 'cagt' | rev)"
    LINE=$CHROM","$POS","$REF","$ALT","$AD_REFLOW","$AD_ALTLOW","$DPLOW","$GQLOW","$PLLOW","$SNPINDEXLOW","$AD_REFHIGH","$AD_ALTHIGH","$DPHIGH","$GQHIGH","$PLHIGH","$SNPINDEXHIGH","$REF_FRQ","$DELTASNP","$NSNPS","$TRICUBEDELTASNP","$G","$GPRIME","$PVALUE","$NEGLOG10PVAL","$QVALUE","$MINDP","$TRICUBEDP","$CL95","$CL99","$POSLESS250","$POSPLUS250","$SNP1","$SNP501","$SNP27","$SNP24","$RSNP27","$RSNP24
    echo $LINE >> $OUTPUT
    IFS=,
  done < $INPUT

#Returns $IFS to its original value.
IFS=$OLDIFS

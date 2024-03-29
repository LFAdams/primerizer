#PBS -S /bin/bash
#PBS -q batch
#PBS -N primerizer
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=480:00:00
#PBS -l mem=80gb

#above is the sapelo2 que parameters
#below changes the working directory of the compute node to current directory of
#the login node used to submit the script
cd $PBS_O_WORKDIR

#this loads the required modules on the cluster
module load SAMtools/1.6-foss-2016b

#Set reference genome, input file path, and output file path. Change to your directories.
GENOME=/work/cemlab/reference_genomes/97103_v2.fa
INPUT=/home/lfa81121/primerizer/testsnps.csv
OUTPUT=/home/lfa81121/primerizer/potentialprimers.csv

#Converts line endings to unix format. Might need to use mac2unix instead.
dos2unix $INPUT

#Stores the current value for $IFS then changes the internal field seperator to
#a comma in order to read .csv columns.
OLDIFS=$IFS
IFS=,

#Checks if $INPUT is a valid file.
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }

#Creates the output file and writes the headers into the first row.
#This should be changed to correspond to the columns of your .csv file, make
#sure to include "POSLESS250,POSPLUS250,SNP1,SNP501,SNP27,SNP24,RSNP27,RSNP24"
#at the end, as this sets the labels for the columns that will be created. ex: echo X,Y,Z
echo ,CHROM,POS,REF,ALT,AD_REFLOW,AD_ALTLOW,DPLOW,GQLOW,PLLOW,SNPINDEXLOW,AD_REFHIGH,AD_ALTHIGH,DPHIGH,GQHIGH,PLHIGH,SNPINDEXHIGH,REF_FRQ,DELTASNP,NSNPS,TRICUBEDELTASNP,G,GPRIME,PVALUE,NEGLOG10PVAL,QVALUE,MINDP,TRICUBEDP,CL95,CL99,POSLESS250,POSPLUS250,SNP1,SNP501,SNP27,SNP24,RSNP27,RSNP24 > $OUTPUT

#Reads each line of the .csv file, $INPUT, and calculates target positions in geneome. Target
#positions are extracted from $GENOME and printed into $OUTPUT along with all
#columns from $INPUT.
#You must copy the same column headers from line 36 above into lines 48 and 64.
#Line 48 should have the column headers seperated by spaces rather than commas as above.
#this reads them into bash as variables. ex: while read X Y Z
#The column headers are then called as variables on line 64 preceded by $ to call the variable,
#and seperated by a comma enclosed in quotes. ex: LINE=$X","$Y","$Z
#Note that we must change $IFS back to its original value in order to output to
#the .csv, and then return it to "," again to restart the loop.
while read NUM CHROM POS REF ALT AD_REFLOW AD_ALTLOW DPLOW GQLOW PLLOW SNPINDEXLOW AD_REFHIGH AD_ALTHIGH DPHIGH GQHIGH PLHIGH SNPINDEXHIGH REF_FRQ DELTASNP NSNPS TRICUBEDELTASNP G GPRIME PVALUE NEGLOG10PVAL QVALUE MINDP TRICUBEDP CL95 CL99
  do
    IFS=$OLDIFS
    SNP1="$(samtools faidx $GENOME $CHROM":"$POS"-"$POS | awk '{if(NR>1)print}')"
    POSPLUS250=$((POS+250))
    POSLESS250=$((POS-250))
    SNP501="$(samtools faidx $GENOME $CHROM":"$POSLESS250"-"$POSPLUS250 | awk '{if(NR>1)print}')"
    SNP501NOWHITE="$(echo -e "${SNP501}" | tr -d '[:space:]')"
    POSLESS23=$((POS-23))
    POSLESS26=$((POS-26))
    POSPLUS23=$((POS+23))
    POSPLUS26=$((POS+26))
    SNP27="$(samtools faidx $GENOME $CHROM":"$POSLESS26"-"$POS | awk '{if(NR>1)print}')"
    SNP24="$(samtools faidx $GENOME $CHROM":"$POSLESS23"-"$POS | awk '{if(NR>1)print}')"
    RSNP27="$(samtools faidx $GENOME $CHROM":"$POS"-"$POSPLUS26 | awk '{if(NR>1)print}' | tr 'gtcaGTCA' 'cagtCAGT' | rev)"
    RSNP24="$(samtools faidx $GENOME $CHROM":"$POS"-"$POSPLUS23 | awk '{if(NR>1)print}' | tr 'gtcaGTCA' 'cagtCAGT' | rev)"
    LINE=$NUM","$CHROM","$POS","$REF","$ALT","$AD_REFLOW","$AD_ALTLOW","$DPLOW","$GQLOW","$PLLOW","$SNPINDEXLOW","$AD_REFHIGH","$AD_ALTHIGH","$DPHIGH","$GQHIGH","$PLHIGH","$SNPINDEXHIGH","$REF_FRQ","$DELTASNP","$NSNPS","$TRICUBEDELTASNP","$G","$GPRIME","$PVALUE","$NEGLOG10PVAL","$QVALUE","$MINDP","$TRICUBEDP","$CL95","$CL99","$POSLESS250","$POSPLUS250","$SNP1","$SNP501NOWHITE","$SNP27","$SNP24","$RSNP27","$RSNP24
    echo $LINE >> $OUTPUT
    IFS=,
  done < <(tail -n +2 $INPUT)

#Returns $IFS to its original value.
IFS=$OLDIFS

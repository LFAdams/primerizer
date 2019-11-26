#PBS -S /bin/bash
#PBS -N 189225
#PBS -q batch
#PBS -l nodes=1:ppn=12
#PBS -l walltime=48:00:00
#PBS -l mem=32gb

##before running this script be sure to cat together all forward/reverse reads
##for each bulk into a single forward read and single reverse read file per
##bulk, and download your reference genome.

BASEDIR=/scratch/lfa81121/189225
BULKHIGH1=/home/lfa81121/SRR8751599_1.fastq
BULKHIGH2=/home/lfa81121/SRR8751599_2.fastq
REFGENOME=/work/cemlab/reference_genomes/97103_v2.fa

mkdir $BASEDIR
cd $BASEDIR

module load SAMtools/1.6-foss-2016b

time  samtools faidx $REFGENOME

module load BWA/0.7.15-foss-2016b


time bwa mem -t 8 $REFGENOME $BULKHIGH1 $BULKHIGH2 > Hbulk_aligned.sam

time  samtools view -@ 8 -S -b Hbulk_aligned.sam > H_aligned.bam


time  samtools sort -@ 8 -m 3G H_aligned.bam -o H_sorted.bam


time samtools index H_sorted.bam


module load picard/2.16.0-Java-1.8.0_144

time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar AddOrReplaceReadGroups INPUT=H_sorted.bam OUTPUT=H_sorted.rg.bam RGID=high_bulk RGSM=high_bulk RGLB=high_bulk RGPL=ILLUMINA RGPU=ignore


time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates INPUT= H_sorted.rg.bam OUTPUT=H_sorted_mkdupl_rg.bam METRICS_FILE=H_sorted_mkduplMetrics.txt


time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar BuildBamIndex INPUT=H_sorted_mkdupl_rg.bam


time java -Xmx20g -classpath "/usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144" -jar  /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar CreateSequenceDictionary REFERENCE=$REFGENOME OUTPUT="$REFGENOME".dict

module load GATK/3.8-1-Java-1.8.0_144

java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REFGENOME -I H_sorted_mkdupl_rg.bam -o H_intervals.list


java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $REFGENOME -I H_sorted_mkdupl_rg.bam -targetIntervals H_intervals.list -o H_realigned.bam


java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller  -nct 8 -R $REFGENOME -I H_realigned.bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o H_realigned_reads_raw_variants_gvcf.vcf

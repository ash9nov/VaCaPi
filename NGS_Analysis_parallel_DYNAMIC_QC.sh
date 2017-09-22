#!/bin/bash

# file transfer from miSEQ directory
echo "Please enter the experiment directory name:"
read dir1
echo "Enter the password for file transfer"
read PASS
	#if we need to hide the visibility of password
	# stty -echo; read PASS; stty echo;
#writing the log file.........
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>/data/log_record/$dir1.log 2>&1
#.............................
echo "You have selected $dir1 for analysis"

echo "Analysis starteded at:"; date +"%T" ; date +'%d/%m/%Y';

echo "copying the raw experiment data to the data directory for analysis"
cp -r /mnt/miseq/Medisinsk_Genetikk/Data/$dir1  /data/

echo "total number of samples to be processed:"
smpl=$(find /data/$dir1/Data/Intensities/BaseCalls/ -name "*.fastq.gz" | grep -v _R2_001.fastq.gz |grep -v Undetermind*|sort | sed 's/_R1_001.fastq.gz//' | wc -l)
echo $smpl

a=$((smpl/8))
b=$((smpl%8))

echo " changing the working directory to Reference directory for processing"
cd /data/Data/hg19/
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#============================================||-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Short Read Alignment (DATA pre-processing): ||-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#=========================================== ||-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BWA alignment:..............................
echo "Alignment of the short reads FASTQ files with the Reference genome (hg19):"

i=1
while [ "$i" -le "$a" ]; 
do
	c=$((i*8))
	find /data/$dir1/Data/Intensities/BaseCalls/ -name "*.fastq.gz" | grep -v _R2_001.fastq.gz |grep -v Undetermined*|sort | head -n $c| tail -n 8| sed 's/_R1_001.fastq.gz//' | parallel bwa mem -M ucsc.hg19.fasta {}_R1_001.fastq.gz {}_R2_001.fastq.gz '>' '{}'.sam 
	i=$(($i+1))
done
##
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*.fastq.gz" | grep -v _R2_001.fastq.gz |grep -v Undetermined*|sort | tail -n $b| sed 's/_R1_001.fastq.gz//' | parallel bwa mem -M ucsc.hg19.fasta {}_R1_001.fastq.gz {}_R2_001.fastq.gz '>' '{}'.sam 

#.............................................
# removing the fastq.gz file for making space.
rm /data/$dir1/Data/Intensities/BaseCalls/*.fastq.gz
#.............................................
#.........................................................................................................................................................................................
# Picard Data Preprocessing:..................
#Convert to BAM, sort and mark duplicates
echo "Converting aligned SAM files to BAM files:"
#a: Convert to BAM and sorting

i=1
while [ "$i" -le "$a" ]; 
do
	c=$((i*8))
	find /data/$dir1/Data/Intensities/BaseCalls/ -name "*.sam"| sort | head -n $c| tail -n 8| sed 's/_L001.sam//' | parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar SortSam INPUT= {}_L001.sam OUTPUT= {}_sorted.bam SORT_ORDER=coordinate
	i=$(($i+1))
done
##
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*.sam"| sort | tail -n $b| sed 's/_L001.sam//' | parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar SortSam INPUT= {}_L001.sam OUTPUT= {}_sorted.bam SORT_ORDER=coordinate

#.............................................
# removing the SAM file for making space.
rm /data/$dir1/Data/Intensities/BaseCalls/*.sam
#.............................................

echo "Marking the Duplicats in BAM files:"	
	#b: Marking the Duplicate
i=1
while [ "$i" -le "$a" ]; 
do
	c=$((i*8))
	find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted.bam"| sort | head -n $c| tail -n 8| sed 's/_sorted.bam//' | parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar MarkDuplicates INPUT= {}_sorted.bam OUTPUT= {}_sorted_dedup.bam METRICS_FILE= {}_metrics.txt
	i=$(($i+1))
done
##
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted.bam"| sort | tail -n $b| sed 's/_sorted.bam//' | parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar MarkDuplicates INPUT= {}_sorted.bam OUTPUT= {}_sorted_dedup.bam METRICS_FILE= {}_metrics.txt

#.............................................
# removing the sorted.bam file for making space.
rm /data/$dir1/Data/Intensities/BaseCalls/*_sorted.bam
# removing the metrics.txt file for making space.
rm /data/$dir1/Data/Intensities/BaseCalls/_metrics.txt
#.............................................
#adding header informating to the bam file
echo "adding header informating to the bam file"

i=1
while [ "$i" -le "$a" ]; 
do
	c=$((i*8))
	find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup.bam"| sort |head -n $c| tail -n 8| sed 's/_sorted_dedup.bam//'| parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar AddOrReplaceReadGroups I= {}_sorted_dedup.bam O= {}_sorted_dedup_RG.bam SORT_ORDER=coordinate RGID= {} RGLB=bar RGPL=illumina RGPU=illumina_miSEQ RGSM= {} CREATE_INDEX=True
	i=$(($i+1))
done
##
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup.bam"| sort | tail -n $b| sed 's/_sorted_dedup.bam//'| parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar AddOrReplaceReadGroups I= {}_sorted_dedup.bam O= {}_sorted_dedup_RG.bam SORT_ORDER=coordinate RGID= {} RGLB=bar RGPL=illumina RGPU=illumina_miSEQ RGSM= {} CREATE_INDEX=True

#.............................................
# removing the sorted_deduplicated file for making space.
rm /data/$dir1/Data/Intensities/BaseCalls/*_sorted_dedup.bam
#.............................................

#.........................................................................................................................................................................................
#Indel Realignment:.....................
echo "INDEL Realignment using GATK:"

	#target creator:
echo "running GATK TARGET CREATOR"
	
i=1
while [ "$i" -le "$a" ]; 
do
	c=$((i*8))
	find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG.bam"| sort |head -n $c| tail -n 8| sed 's/_sorted_dedup_RG.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ucsc.hg19.fasta -I {}_sorted_dedup_RG.bam -known /data/Data/GATK_resources/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -o {}_realigner.intervals
	i=$(($i+1))
done
##
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG.bam"| sort | tail -n $b| sed 's/_sorted_dedup_RG.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ucsc.hg19.fasta -I {}_sorted_dedup_RG.bam -known /data/Data/GATK_resources/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -o {}_realigner.intervals
	
# Performing actual realignment
echo "performing actual ReAlignment"

i=1
while [ "$i" -le "$a" ]; 
do
	c=$((i*8))
	find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG.bam"| sort |head -n $c| tail -n 8| sed 's/_sorted_dedup_RG.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I {}_sorted_dedup_RG.bam -known /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -targetIntervals {}_realigner.intervals -o {}_sorted_dedup_RG_IndReAl.bam
	i=$(($i+1))
done
##
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG.bam"| sort | tail -n $b| sed 's/_sorted_dedup_RG.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I {}_sorted_dedup_RG.bam -known /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -targetIntervals {}_realigner.intervals -o {}_sorted_dedup_RG_IndReAl.bam

#.............................................
# removing the sorted_dedup_RG file for making space.
rm /data/$dir1/Data/Intensities/BaseCalls/*_sorted_dedup_RG.bam
# removing the sorted_dedup_RG file for making space.
rm /data/$dir1/Data/Intensities/BaseCalls/*_realigner.intervals
#.............................................

#Base Quality Score Recalibration:.............
echo "Base Quality Score Recalibration step:"
	#a: BaseRecalibrator:
echo "running GATK BaseRecalibrator"

i=1
while [ "$i" -le "$a" ]; 
do
	c=$((i*8))
	find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG_IndReAl.bam"| sort |head -n $c| tail -n 8| sed 's/_sorted_dedup_RG_IndReAl.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl.bam -knownSites /data/Data/GATK_resources/dbsnp_138.hg19.vcf -knownSites /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -o {}_recal.table
	i=$(($i+1))
done
##
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG_IndReAl.bam"| sort | tail -n $b| sed 's/_sorted_dedup_RG_IndReAl.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl.bam -knownSites /data/Data/GATK_resources/dbsnp_138.hg19.vcf -knownSites /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -o {}_recal.table

	#b: PrintReads: 
echo "running GATK PrintReads"

i=1
while [ "$i" -le "$a" ]; 
do
	c=$((i*8))
	find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG_IndReAl.bam"| sort |head -n $c| tail -n 8| sed 's/_sorted_dedup_RG_IndReAl.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T PrintReads -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl.bam -BQSR {}_recal.table -o {}_sorted_dedup_RG_IndReAl_Baserecal.bam
	i=$(($i+1))
done
##
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG_IndReAl.bam"| sort | tail -n $b| sed 's/_sorted_dedup_RG_IndReAl.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T PrintReads -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl.bam -BQSR {}_recal.table -o {}_sorted_dedup_RG_IndReAl_Baserecal.bam

#.............................................
# removing the sorted_dedup_RG_IndReal file for making space.
rm /data/$dir1/Data/Intensities/BaseCalls/*_sorted_dedup_RG_IndReAl.bam
# removing the recal.table file for making space.
rm /data/$dir1/Data/Intensities/BaseCalls/*_recal.table
#.............................................
	# the FINAL BAM file is : "sorted_dedup_RG_IndReAl_Baserecal.bam" , that will be stored in database
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#====================================||------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Quality Control (Coverage repoprt)  ||
#====================================||------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ANALYSIS CALCULATE THE COVERAGE REPORT:
# Tool used in GATK DepthOfCoverage

# generating the .list file of all the "SORTED_DEDUPLICATED_INDEL-REALIGNED_BASE-RECALIBRATED" BAM files.
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG_IndReAl_Baserecal.bam"| sort > /data/$dir1/Data/Intensities/BaseCalls/FINAL_BAMs.list

# making the Quality control directory
mkdir /data/$dir1/Data/Intensities/BaseCalls/coverage_report

java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T DepthOfCoverage -R ucsc.hg19.fasta -I /data/$dir1/Data/Intensities/BaseCalls/FINAL_BAMs.list -o /data/$dir1/Data/Intensities/BaseCalls/coverage_report/$dir1 -L /data/Data/NexteraRapidCapture-71370-targeted-regions_v4.bed

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#==================||-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Variant Calling   ||
#==================||--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Variant calling:..........................
echo "Running GATK HaplotypeCaller:"	
#HaplotypeCaller:

i=1
while [ "$i" -le "$a" ]; 
do
	c=$((i*8))
	find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG_IndReAl_Baserecal.bam"| sort |head -n $c| tail -n 8| sed 's/_sorted_dedup_RG_IndReAl_Baserecal.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl_Baserecal.bam -o {}_raw_SNP_INDEL.g.vcf -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000 --dbsnp /data/Data/GATK_resources/dbsnp_138.hg19.vcf -L /data/Data/NexteraRapidCapture-71370-targeted-regions_v4.bed
	i=$(($i+1))
done
##
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG_IndReAl_Baserecal.bam"| sort | tail -n $b| sed 's/_sorted_dedup_RG_IndReAl_Baserecal.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl_Baserecal.bam -o {}_raw_SNP_INDEL.g.vcf -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000 --dbsnp /data/Data/GATK_resources/dbsnp_138.hg19.vcf -L /data/Data/NexteraRapidCapture-71370-targeted-regions_v4.bed

#Genotyping:.............................
echo "Running GATK GenotypeGVCFs"
#GenotypeGVCFs

i=1
while [ "$i" -le "$a" ]; 
do
	c=$((i*8))
	find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_raw_SNP_INDEL.g.vcf"| sort |head -n $c| tail -n 8| sed 's/_raw_SNP_INDEL.g.vcf//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ucsc.hg19.fasta -V {}_raw_SNP_INDEL.g.vcf -o {}_SNP_INDEL_genotyped.vcf
	i=$(($i+1))
done
##
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_raw_SNP_INDEL.g.vcf"| sort | tail -n $b| sed 's/_raw_SNP_INDEL.g.vcf//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ucsc.hg19.fasta -V {}_raw_SNP_INDEL.g.vcf -o {}_SNP_INDEL_genotyped.vcf

#.............................................
# removing the raw_SNP_INDEL.g.vcf file for making space.
rm /data/$dir1/Data/Intensities/BaseCalls/*_raw_SNP_INDEL.g.vcf
#.............................................
	#SNP_INDEL_genotyped.vcf is the final VCF file
#.............................................

# variant Quality Score Recalibration:
echo "VQSR step will be SKIPPED!! here, as we are running the gene panels"

#.........................................................................................................................................................................................

#variantFiltration:
echo "GATK variantFiltration:"
for i in /data/$dir1/Data/Intensities/BaseCalls/*_SNP_INDEL_genotyped.vcf ;
	do 
		i2=${i%.*}_filtered.vcf;
		java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T VariantFiltration -R ucsc.hg19.fasta --variant $i -o $i2 --clusterWindowSize 10 --filterExpression "MQ0 >= 4 && ((MQ0/(1.0 * DP)) > 0.1)" --filterName "HARD_TO_VALIDATE" --filterExpression "DP < 5 " --filterName "LowCoverage" --filterExpression "QUAL < 30.0 " --filterName "VeryLowQual" --filterExpression "QUAL > 30.0 && QUAL < 50.0 " --filterName "LowQual" --filterExpression "QD < 1.5 " --filterName "LowQD" --filterExpression "SB > -10.0 " --filterName "StrandBias" ;
	done

#.........................................................................................................................................................................................
cd /data/Data/annovar/

#=====================||--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#ANNOVAR annotation:  ||------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#=====================||---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# ANNOVAR annotation of variants:
echo "ANNOVAR annotation of variants:"

for i in /data/$dir1/Data/Intensities/BaseCalls/*_genotyped_filtered.vcf ;
	do 
		i2=${i%.*}_ANNOVAR;
		perl table_annovar.pl $i humandb/ -buildver hg19 -out $i2 -remove -protocol refGene,phastConsElements46way,genomicSuperDups,1000g2015aug_eur,snp138,ljb26_all,esp6500si_all,exac03,clinvar_20150629,cosmic70 -operation g,r,r,f,f,f,f,f,f,f -nastring . -vcfinput 
	done

#.............................................
# removing the raw_SNP_INDEL.g.vcf file for making space.
rm /data/$dir1/Data/Intensities/BaseCalls/*_SNP_INDEL_genotyped_filtered.vcf
#.............................................

#.........................................................................................................................................................................................
echo "VARIANT CALLING done!!!"
echo "FILTER TAGGING done!!!"
echo "ANNOTATION of variants done!!!"


#=====================||-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#File rearrangement:  ||----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#=====================||---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#File rearrangement:
echo "File Re-Arrangement (in to Sub-Directories) respective to their analysis:"

#....................................
mkdir /data/$dir1/Data/Intensities/BaseCalls/6_GATK_BQSR
mv /data/$dir1/Data/Intensities/BaseCalls/*_sorted_dedup_RG_IndReAl_Baserecal.* /data/$dir1/Data/Intensities/BaseCalls/6_GATK_BQSR

#....................................
mkdir /data/$dir1/Data/Intensities/BaseCalls/8_GATK_GenotypeGVCFs
mv /data/$dir1/Data/Intensities/BaseCalls/*_SNP_INDEL_genotyped.vcf* /data/$dir1/Data/Intensities/BaseCalls/8_GATK_GenotypeGVCFs

#....................................
mkdir /data/$dir1/Data/Intensities/BaseCalls/10_ANNOVAR_annotation
mv /data/$dir1/Data/Intensities/BaseCalls/*_filtered_ANNOVAR* /data/$dir1/Data/Intensities/BaseCalls/10_ANNOVAR_annotation


echo "$PASS\n"| sudo -S mkdir /mnt/miseq/Medisinsk_Genetikk/Resultater/$dir1
echo "$PASS\n"| sudo -S cp -r /data/$dir1/Data/Intensities/BaseCalls/6_GATK_BQSR /mnt/miseq/Medisinsk_Genetikk/Resultater/$dir1
echo "$PASS\n"| sudo -S cp -r /data/$dir1/Data/Intensities/BaseCalls/8_GATK_GenotypeGVCFs /mnt/miseq/Medisinsk_Genetikk/Resultater/$dir1
echo "$PASS\n"| sudo -S cp -r /data/$dir1/Data/Intensities/BaseCalls/10_ANNOVAR_annotation /mnt/miseq/Medisinsk_Genetikk/Resultater/$dir1
echo "$PASS\n"| sudo -S cp -r /data/$dir1/Data/Intensities/BaseCalls/coverage_report /mnt/miseq/Medisinsk_Genetikk/Resultater/$dir1
echo "Analysis finished at:"; date +"%T" ; date +'%d/%m/%Y';
echo "$PASS\n"| sudo -S cp /data/log_record/$dir1.log /mnt/miseq/Medisinsk_Genetikk/Resultater/$dir1

#removing the experiment directory from server
cd /data/
rm -r $dir1 

df -lah
#.........................................................................................................................................................................................
#.........................................................................................................................................................................................
#.........................................................................................................................................................................................


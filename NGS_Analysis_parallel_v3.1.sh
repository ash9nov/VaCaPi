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
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*.fastq.gz" | grep -v _R2_001.fastq.gz |sort | sed 's/_R1_001.fastq.gz//' | wc -l
echo "List of samples to be processed:"
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*.fastq.gz" | grep -v _R2_001.fastq.gz |sort | head -8| sed 's/_R1_001.fastq.gz//'

echo " changing the working directory to Reference directory for processing"
cd /data/Data/hg19/
#.........................................................................................................................................................................................
# BWA alignment:..............................
echo "Alignment of the short reads FASTQ files with the Reference genome (hg19):"
find /data/$dir1/Data/Intensities/BaseCalls/  -name "*.fastq.gz" | grep -v _R2_001.fastq.gz |sort | head -n 8| tail -n 8|sed 's/_R1_001.fastq.gz//' | parallel bwa mem -M ucsc.hg19.fasta {}_R1_001.fastq.gz {}_R2_001.fastq.gz '>' '{}'.sam 
find /data/$dir1/Data/Intensities/BaseCalls/  -name "*.fastq.gz" | grep -v _R2_001.fastq.gz |sort | head -n 16| tail -n 8|sed 's/_R1_001.fastq.gz//' | parallel bwa mem -M ucsc.hg19.fasta {}_R1_001.fastq.gz {}_R2_001.fastq.gz '>' '{}'.sam 
find /data/$dir1/Data/Intensities/BaseCalls/  -name "*.fastq.gz" | grep -v _R2_001.fastq.gz |sort | head -n 24| tail -n 8|sed 's/_R1_001.fastq.gz//' | parallel bwa mem -M ucsc.hg19.fasta {}_R1_001.fastq.gz {}_R2_001.fastq.gz '>' '{}'.sam 
#.........................................................................................................................................................................................
# Picard Data Preprocessing:..................
#Convert to BAM, sort and mark duplicates
echo "Converting aligned SAM files to BAM files:"
	#a: Convert to BAM and sorting
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*.sam"| sort | head -n 8| tail -n 8| sed 's/_L001.sam//' | parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar SortSam INPUT= {}_L001.sam OUTPUT= {}_sorted.bam SORT_ORDER=coordinate
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*.sam"| sort | head -n 16| tail -n 8| sed 's/_L001.sam//' | parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar SortSam INPUT= {}_L001.sam OUTPUT= {}_sorted.bam SORT_ORDER=coordinate
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*.sam"| sort | head -n 24| tail -n 8| sed 's/_L001.sam//' | parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar SortSam INPUT= {}_L001.sam OUTPUT= {}_sorted.bam SORT_ORDER=coordinate

echo "Marking the Duplicats in BAM files:"	
	#b: Marking the Duplicate
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted.bam"| sort | head -n 8| tail -n 8| sed 's/_sorted.bam//' | parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar MarkDuplicates INPUT= {}_sorted.bam OUTPUT= {}_sorted_dedup.bam METRICS_FILE= {}_metrics.txt
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted.bam"| sort | head -n 16| tail -n 8| sed 's/_sorted.bam//' | parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar MarkDuplicates INPUT= {}_sorted.bam OUTPUT= {}_sorted_dedup.bam METRICS_FILE= {}_metrics.txt
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted.bam"| sort | head -n 24| tail -n 8| sed 's/_sorted.bam//' | parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar MarkDuplicates INPUT= {}_sorted.bam OUTPUT= {}_sorted_dedup.bam METRICS_FILE= {}_metrics.txt

#adding header informating to the bam file
echo "adding header informating to the bam file"
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup.bam"| sort |head -n 8| tail -n 8| sed 's/_sorted_dedup.bam//'| parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar AddOrReplaceReadGroups I= {}_sorted_dedup.bam O= {}_sorted_dedup_RG.bam SORT_ORDER=coordinate RGID= {} RGLB=bar RGPL=illumina RGPU=illumina_miSEQ RGSM= {} CREATE_INDEX=True
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup.bam"| sort |head -n 16| tail -n 8| sed 's/_sorted_dedup.bam//'| parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar AddOrReplaceReadGroups I= {}_sorted_dedup.bam O= {}_sorted_dedup_RG.bam SORT_ORDER=coordinate RGID= {} RGLB=bar RGPL=illumina RGPU=illumina_miSEQ RGSM= {} CREATE_INDEX=True
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup.bam"| sort |head -n 24| tail -n 8| sed 's/_sorted_dedup.bam//'| parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar AddOrReplaceReadGroups I= {}_sorted_dedup.bam O= {}_sorted_dedup_RG.bam SORT_ORDER=coordinate RGID= {} RGLB=bar RGPL=illumina RGPU=illumina_miSEQ RGSM= {} CREATE_INDEX=True
#.........................................................................................................................................................................................
#Indel Realignment:.....................
echo "INDEL Realignment using GATK:"

	#target creator:
echo "running GATK TARGET CREATOR"
	find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG.bam"| sort |head -n 8| tail -n 8| sed 's/_sorted_dedup_RG.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ucsc.hg19.fasta -I {}_sorted_dedup_RG.bam -known /data/Data/GATK_resources/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -o {}_realigner.intervals
	find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG.bam"| sort |head -n 16| tail -n 8| sed 's/_sorted_dedup_RG.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ucsc.hg19.fasta -I {}_sorted_dedup_RG.bam -known /data/Data/GATK_resources/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -o {}_realigner.intervals
	find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG.bam"| sort |head -n 24| tail -n 8| sed 's/_sorted_dedup_RG.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ucsc.hg19.fasta -I {}_sorted_dedup_RG.bam -known /data/Data/GATK_resources/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -o {}_realigner.intervals
	# Performing actual alignment
echo "performing actual Alignment"	
	find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG.bam"| sort |head -n 8| tail -n 8| sed 's/_sorted_dedup_RG.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I {}_sorted_dedup_RG.bam -known /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -targetIntervals {}_realigner.intervals -o {}_sorted_dedup_RG_IndReAl.bam
	find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG.bam"| sort |head -n 16| tail -n 8| sed 's/_sorted_dedup_RG.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I {}_sorted_dedup_RG.bam -known /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -targetIntervals {}_realigner.intervals -o {}_sorted_dedup_RG_IndReAl.bam
	find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG.bam"| sort |head -n 24| tail -n 8| sed 's/_sorted_dedup_RG.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I {}_sorted_dedup_RG.bam -known /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -targetIntervals {}_realigner.intervals -o {}_sorted_dedup_RG_IndReAl.bam

#Base Quality Score Recalibration:.............
echo "Base Quality Score Recalibration step:"
	#a: BaseRecalibrator:
echo "running GATK BaseRecalibrator"
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG_IndReAl.bam"| sort |head -n 8| tail -n 8| sed 's/_sorted_dedup_RG_IndReAl.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl.bam -knownSites /data/Data/GATK_resources/dbsnp_138.hg19.vcf -knownSites /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -o {}_recal.table
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG_IndReAl.bam"| sort |head -n 16| tail -n 8| sed 's/_sorted_dedup_RG_IndReAl.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl.bam -knownSites /data/Data/GATK_resources/dbsnp_138.hg19.vcf -knownSites /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -o {}_recal.table
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG_IndReAl.bam"| sort |head -n 24| tail -n 8| sed 's/_sorted_dedup_RG_IndReAl.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl.bam -knownSites /data/Data/GATK_resources/dbsnp_138.hg19.vcf -knownSites /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -o {}_recal.table

	#b: PrintReads: 
echo "running GATK PrintReads"
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG_IndReAl.bam"| sort |head -n 8| tail -n 8| sed 's/_sorted_dedup_RG_IndReAl.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T PrintReads -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl.bam -BQSR {}_recal.table -o {}_sorted_dedup_RG_IndReAl_Baserecal.bam
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG_IndReAl.bam"| sort |head -n 16| tail -n 8| sed 's/_sorted_dedup_RG_IndReAl.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T PrintReads -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl.bam -BQSR {}_recal.table -o {}_sorted_dedup_RG_IndReAl_Baserecal.bam
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG_IndReAl.bam"| sort |head -n 24| tail -n 8| sed 's/_sorted_dedup_RG_IndReAl.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T PrintReads -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl.bam -BQSR {}_recal.table -o {}_sorted_dedup_RG_IndReAl_Baserecal.bam

#Variant calling:..........................
echo "Running GATK HaplotypeCaller:"	
#HaplotypeCaller:
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG_IndReAl_Baserecal.bam"| sort |head -n 8| tail -n 8| sed 's/_sorted_dedup_RG_IndReAl_Baserecal.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl_Baserecal.bam -o {}_raw_SNP_INDEL.g.vcf -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000 --dbsnp /data/Data/GATK_resources/dbsnp_138.hg19.vcf -L /data/Data/NexteraRapidCapture-71370-targeted-regions_v4.bed
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG_IndReAl_Baserecal.bam"| sort |head -n 16| tail -n 8| sed 's/_sorted_dedup_RG_IndReAl_Baserecal.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl_Baserecal.bam -o {}_raw_SNP_INDEL.g.vcf -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000 --dbsnp /data/Data/GATK_resources/dbsnp_138.hg19.vcf -L /data/Data/NexteraRapidCapture-71370-targeted-regions_v4.bed
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_sorted_dedup_RG_IndReAl_Baserecal.bam"| sort |head -n 24| tail -n 8| sed 's/_sorted_dedup_RG_IndReAl_Baserecal.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl_Baserecal.bam -o {}_raw_SNP_INDEL.g.vcf -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000 --dbsnp /data/Data/GATK_resources/dbsnp_138.hg19.vcf -L /data/Data/NexteraRapidCapture-71370-targeted-regions_v4.bed

#Genotyping:.............................
echo "Running GATK GenotypeGVCFs"
#GenotypeGVCFs
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_raw_SNP_INDEL.g.vcf"| sort |head -n 8| tail -n 8| sed 's/_raw_SNP_INDEL.g.vcf//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ucsc.hg19.fasta -V {}_raw_SNP_INDEL.g.vcf -o {}_SNP_INDEL_genotyped.vcf
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_raw_SNP_INDEL.g.vcf"| sort |head -n 16| tail -n 8| sed 's/_raw_SNP_INDEL.g.vcf//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ucsc.hg19.fasta -V {}_raw_SNP_INDEL.g.vcf -o {}_SNP_INDEL_genotyped.vcf
find /data/$dir1/Data/Intensities/BaseCalls/ -name "*_raw_SNP_INDEL.g.vcf"| sort |head -n 24| tail -n 8| sed 's/_raw_SNP_INDEL.g.vcf//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ucsc.hg19.fasta -V {}_raw_SNP_INDEL.g.vcf -o {}_SNP_INDEL_genotyped.vcf

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

# ANNOVAR annotation of variants:
echo "ANNOVAR annotation of variants:"

for i in /data/$dir1/Data/Intensities/BaseCalls/*_genotyped_filtered.vcf ;
	do 
		i2=${i%.*}_ANNOVAR;
		perl table_annovar.pl $i humandb/ -buildver hg19 -out $i2 -remove -protocol refGene,phastConsElements46way,genomicSuperDups,1000g2015aug_eur,snp138,ljb26_all,esp6500si_all,exac03,clinvar_20150629,cosmic70 -operation g,r,r,f,f,f,f,f,f,f -nastring . -vcfinput 
	done

#.........................................................................................................................................................................................
echo "VARIANT CALLING done!!!"
echo "FILTER TAGGING done!!!"
echo "ANNOTATION of variants done!!!"
#.........................................................................................................................................................................................

#File rearrangement:
echo "File Re-Arrangement (in to Sub-Directories) respective to their analysis:"
mkdir /data/$dir1/Data/Intensities/BaseCalls/1_BWA_mapped_SAM
mv /data/$dir1/Data/Intensities/BaseCalls/*.sam /data/$dir1/Data/Intensities/BaseCalls/1_BWA_mapped_SAM
#....................................
mkdir /data/$dir1/Data/Intensities/BaseCalls/2_Picard_sorted
mv /data/$dir1/Data/Intensities/BaseCalls/*_sorted.bam* /data/$dir1/Data/Intensities/BaseCalls/2_Picard_sorted
#....................................
mkdir /data/$dir1/Data/Intensities/BaseCalls/3_Picard_DeDuplicated
mv /data/$dir1/Data/Intensities/BaseCalls/*_sorted_dedup.bam* /data/$dir1/Data/Intensities/BaseCalls/3_Picard_DeDuplicated
mv /data/$dir1/Data/Intensities/BaseCalls/*_metrics.txt* /data/$dir1/Data/Intensities/BaseCalls/3_Picard_DeDuplicated
#....................................
mkdir /data/$dir1/Data/Intensities/BaseCalls/4_Picard_withReadGroup
mv /data/$dir1/Data/Intensities/BaseCalls/*_sorted_dedup_RG.* /data/$dir1/Data/Intensities/BaseCalls/4_Picard_withReadGroup
#....................................
mkdir /data/$dir1/Data/Intensities/BaseCalls/5_GATK_Indel_Realigned
mv /data/$dir1/Data/Intensities/BaseCalls/*_realigner.intervals /data/$dir1/Data/Intensities/BaseCalls/5_GATK_Indel_Realigned
mv /data/$dir1/Data/Intensities/BaseCalls/*_sorted_dedup_RG_IndReAl.* /data/$dir1/Data/Intensities/BaseCalls/5_GATK_Indel_Realigned
#....................................
mkdir /data/$dir1/Data/Intensities/BaseCalls/6_GATK_BQSR
mv /data/$dir1/Data/Intensities/BaseCalls/*_recal.table /data/$dir1/Data/Intensities/BaseCalls/6_GATK_BQSR
mv /data/$dir1/Data/Intensities/BaseCalls/*_sorted_dedup_RG_IndReAl_Baserecal.* /data/$dir1/Data/Intensities/BaseCalls/6_GATK_BQSR
#....................................
mkdir /data/$dir1/Data/Intensities/BaseCalls/7_GATK_HaplotypeCaller
mv /data/$dir1/Data/Intensities/BaseCalls/*_raw_SNP_INDEL.g.vcf* /data/$dir1/Data/Intensities/BaseCalls/7_GATK_HaplotypeCaller
#....................................
mkdir /data/$dir1/Data/Intensities/BaseCalls/8_GATK_GenotypeGVCFs
mv /data/$dir1/Data/Intensities/BaseCalls/*_SNP_INDEL_genotyped.vcf* /data/$dir1/Data/Intensities/BaseCalls/8_GATK_GenotypeGVCFs

mkdir /data/$dir1/Data/Intensities/BaseCalls/9_GATK_VariantFiltration
mv /data/$dir1/Data/Intensities/BaseCalls/*_SNP_INDEL_genotyped_filtered.vcf* /data/$dir1/Data/Intensities/BaseCalls/9_GATK_VariantFiltration

mkdir /data/$dir1/Data/Intensities/BaseCalls/10_ANNOVAR_annotation
mv /data/$dir1/Data/Intensities/BaseCalls/*_filtered_ANNOVAR* /data/$dir1/Data/Intensities/BaseCalls/10_ANNOVAR_annotation


echo "$PASS\n"| sudo -S mkdir /mnt/miseq/Medisinsk_Genetikk/Resultater/$dir1
echo "$PASS\n"| sudo -S cp -r /data/$dir1/Data/Intensities/BaseCalls/6_GATK_BQSR /mnt/miseq/Medisinsk_Genetikk/Resultater/$dir1
echo "$PASS\n"| sudo -S cp -r /data/$dir1/Data/Intensities/BaseCalls/8_GATK_GenotypeGVCFs /mnt/miseq/Medisinsk_Genetikk/Resultater/$dir1
echo "$PASS\n"| sudo -S cp -r /data/$dir1/Data/Intensities/BaseCalls/10_ANNOVAR_annotation /mnt/miseq/Medisinsk_Genetikk/Resultater/$dir1

echo "Analysis finished at:"; date +"%T" ; date +'%d/%m/%Y';

echo "$PASS\n"| sudo -S cp /data/log_record/$dir1.log /mnt/miseq/Medisinsk_Genetikk/Resultater/$dir1

#.........................................................................................................................................................................................


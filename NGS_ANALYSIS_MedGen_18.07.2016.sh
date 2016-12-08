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
rm /data/$dir1/Data/Intensities/BaseCalls/*_metrics.txt
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
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
# Running the GATK-DepthOfCoverage command. It will generate 7 files
java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T DepthOfCoverage -R ucsc.hg19.fasta -I /data/$dir1/Data/Intensities/BaseCalls/FINAL_BAMs.list -o /data/$dir1/Data/Intensities/BaseCalls/coverage_report/$dir1 -L /data/Data/NexteraRapidCapture-71370-targeted-regions_v4.bed
# NOTE! Out of these 7 files, we will use only 2 files "$dir1.sample_summary" && "$dir1.sample_interval_summary" for our need.
# Removing the unnecessery 5 files:
#rm /data/$dir1/Data/Intensities/BaseCalls/coverage_report/$dir1
rm /data/$dir1/Data/Intensities/BaseCalls/coverage_report/$dir1.sample_cumulative_coverage_counts
rm /data/$dir1/Data/Intensities/BaseCalls/coverage_report/$dir1.sample_cumulative_coverage_proportions
rm /data/$dir1/Data/Intensities/BaseCalls/coverage_report/$dir1.sample_interval_statistics
rm /data/$dir1/Data/Intensities/BaseCalls/coverage_report/$dir1.sample_statistics

#Customization of coverage report data..
#change of directory (from hg19 to coverage_report)
cd /data/$dir1/Data/Intensities/BaseCalls/coverage_report/
# changing the long names
sed -i "s/\/data\/$dir1\/Data\/Intensities\/BaseCalls\///g" *
#Generating the sample_low_coverage_region file (<30) from the sample_interval_summary file (superset) ............
cut -f 1 "$dir1.sample_interval_summary" > "$dir1.sample_interval_mean_coverage"
cut -f 1 "$dir1.sample_interval_summary" > "$dir1.temp1"
i=1
while [ "$i" -le "$smpl" ] ; 
do
        jj=$(((6*$i)-1))
        cut -f $jj "$dir1.sample_interval_summary" > "$dir1.temp"
        paste -d "\t" "$dir1.sample_interval_mean_coverage" "$dir1.temp" > "$dir1.op"
        mv "$dir1.op" "$dir1.sample_interval_mean_coverage"
        read -r FIRSTLINE < $dir1.temp
        paste -d "\t" "$dir1.temp1" "$dir1.temp" > "$FIRSTLINE.txt"
        awk 'FNR==NR { a[$1]=$2; next } $1 in a { print $1"\t"a[$1]"\t"$2"\t"$3 }' "$FIRSTLINE.txt" /data/Data/Target_ex_n_gn_mdfd_strtPlus1_merged_regions_for_CovRept_2120.bed > "$FIRSTLINE.annotated.txt"
        awk -v x=29 -F "\t" '$2<=x {print}' "$FIRSTLINE.annotated.txt" > "$FIRSTLINE.low_coverage_region.txt"
        rm "$FIRSTLINE.txt"
        rm "$FIRSTLINE.annotated.txt"
        i=$(($i+1));
done
rm $dir1.temp*
sort "$dir1.sample_interval_mean_coverage" > bb1
sort /data/Data/Target_ex_n_gn_mdfd_strtPlus1_merged_regions_for_CovRept_2120.bed > bb2
join -j 1 bb1 bb2 > "$dir1.sample_interval_mean_coverage_annotated.txt"
rm bb1
rm bb2
rm $dir1.sample_interval_mean_coverage
rm $dir1.sample_interval_summary
# Generating the Low_coverage_gap_region files for each sample
cut -f 1 "$dir1" > "$dir1.sample_coverage"
cut -f 1 "$dir1" > "$dir1.temp1"
i=1
while [ "$i" -le "$smpl" ] ; 
do
        jj=$(($i+3))
        cut -f $jj "$dir1" > "$dir1.temp"
        paste -d "\t" "$dir1.sample_coverage" "$dir1.temp" > "$dir1.op"
        mv "$dir1.op" "$dir1.sample_coverage"
        read -r FIRSTLINE < $dir1.temp
        paste -d "\t" "$dir1.temp1" "$dir1.temp" > "$dir1.FIRSTLINE.txt"
        cut -f 2 "$dir1.FIRSTLINE.txt" >"$dir1.FIRSTLINE_coverage"
        paste -d "\t" /data/Data/Target_Nucleotide_annotated.csv "$dir1.FIRSTLINE_coverage" >"$dir1.FIRSTLINE.txt"
        awk -v x=29 -F "\t" '$5<=x {print}' "$dir1.FIRSTLINE.txt" > "$FIRSTLINE.gap_region.txt"
        rm "$dir1.FIRSTLINE.txt"
        i=$(($i+1));
	Rscript ~/my_tools/NGS_pipeline/R_code_for_gap_region_generation.r $FIRSTLINE.gap_region.txt $FIRSTLINE.gap_region_final.csv
rm "$FIRSTLINE.gap_region.txt"
done
rm "$dir1.FIRSTLINE_coverage"
rm "$dir1.sample_coverage"
rm $dir1.temp*
#---------------------------------------------
#Running Samtools-flagstat to generate the mapping statistics
for i in /data/$dir1/Data/Intensities/BaseCalls/*_sorted_dedup_RG_IndReAl_Baserecal.bam;
	do
		i2=${i%.*}_mapping_statistics.txt;
                samtools flagstat $i > $i2
	done
#adding sample name on top of each sample mapping_statistics file
for i in /data/$dir1/Data/Intensities/BaseCalls/*_mapping_statistics.txt;
	do
		echo "$i\n$(cat $i)" > $i # it adds the sample name on to of FLAGSTAT output
                sed -i "s/\/data\/$dir1\/Data\/Intensities\/BaseCalls\///g" $i
                sed -i "s/_sorted_dedup_RG_IndReAl_Baserecal_mapping_statistics.txt//g" $i
	done
#generating tab-seperated table of all-sample mapping summary
cp /data/Data/flag_property "$dir1.mapping_summary.txt"
for i in /data/$dir1/Data/Intensities/BaseCalls/*_mapping_statistics.txt;
	do
                cut -d "+" -f 1 "$i" > "$dir1.temp"
                paste -d "\t" "$dir1.mapping_summary.txt" "$dir1.temp" > "$dir1.map"
                mv "$dir1.map" "$dir1.mapping_summary.txt"
        done
rm $dir1.temp
rm *mapping_statistics.txt
#---------------------------------------------
# Generating PER-RUN total mapping percentage
TOTAL=$(grep 'total' $dir1.mapping_summary.txt | sed 's/\t/+/g' | bc)
echo "total reads= $TOTAL"
MAPPED=$(grep 'mapped' $dir1.mapping_summary.txt | grep -v '_mapped'| sed 's/\t/+/g' | bc)  # sum up the entries in colum
echo "mapped reads = $MAPPED"
Ratio=$(awk "BEGIN {printf \"%.8f\",${MAPPED}/${TOTAL}}")
echo "ratio= $Ratio"
percentage_stat=$(echo "scale=6; $Ratio*100" | bc)
echo "%percentage mapping for the present run= $percentage_stat"
echo "Total reads= $TOTAL" > $dir1.mapped_percentage.txt
echo "Total Mapped reads= $MAPPED" >> $dir1.mapped_percentage.txt
echo "Percentage of mapping for the run $dir1 is= $percentage_stat %" >> $dir1.mapped_percentage.txt
#---------------------------
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#=======================================||---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Gender Control (Copy Number Variation) ||
#=======================================||---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Creating directory for CNV analysis 
mkdir /data/$dir1/Data/Intensities/BaseCalls/cnv_analysis
cd /data/$dir1/Data/Intensities/BaseCalls/cnv_analysis
# step1: Generating Sequencing-accessible regions
# cnvkit.py access /data/Data/hg19/ucsc.hg19.fasta -s 10000 -o /data/$dir1/Data/Intensities/BaseCalls/cnv_analysis/access-10kb.mn10.bed
#running all steps in single "Batch" command
cnvkit.py  batch /data/$dir1/Data/Intensities/BaseCalls/*IndReAl_Baserecal.bam --normal --targets /data/Data/Target_with_exon_and_gene_info.bed --fasta /data/Data/hg19/ucsc.hg19.fasta --split --access /data/Data/access-10kb.mn10.bed --output-reference /data/$dir1/Data/Intensities/BaseCalls/cnv_analysis/flat_reference.cnn -d /data/$dir1/Data/Intensities/BaseCalls/cnv_analysis
# Identifying the gender of pasients
cnvkit.py gender /data/$dir1/Data/Intensities/BaseCalls/cnv_analysis/*cns > /data/$dir1/Data/Intensities/BaseCalls/cnv_analysis/p_gender.txt
#changing the long names
sed -i "s/\/data\/$dir1\/Data\/Intensities\/BaseCalls\/cnv_analysis\///g" /data/$dir1/Data/Intensities/BaseCalls/cnv_analysis/p_gender.txt
sed -i "s/_sorted_dedup_RG_IndReAl_Baserecal.cns//g" /data/$dir1/Data/Intensities/BaseCalls/cnv_analysis/p_gender.txt
cut -f 1-2 /data/$dir1/Data/Intensities/BaseCalls/cnv_analysis/p_gender.txt > /data/$dir1/Data/Intensities/BaseCalls/cnv_analysis/pasient_gender.txt
rm /data/$dir1/Data/Intensities/BaseCalls/cnv_analysis/p_gender.txt
#--------------------------------------------------------------------------
#Change of working directory back to "hg19"
cd /data/Data/hg19/
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


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

#Renaming the fields in INFO column for organized view in FILTUS.
for i in /data/$dir1/Data/Intensities/BaseCalls/*hg19_multianno.vcf ;
	do
		i2=${i%.*}_RENAMED.vcf;
		sed -e 's/;Gene.refGene/;AAA_01_Gene.refGene/g; s/;Func.refGene/;AAA_02_Func.refGene/g; s/;GeneDetail.refGene/;AAA_03_GeneDetail.refGene/g; s/;AAChange.refGene/;AAA_04_AAChange.refGene/g; s/;ExonicFunc.refGene/;AAA_05_ExonicFunc.refGene/g; s/;snp138/;AAA_06_snp138/g; s/;1000g2015aug_eur/;AAA_07_1000g2015aug_eur/g; s/;ExAC_AFR/;AAA_08_ExAC_AFR/g; s/;ExAC_ALL/;AAA_09_ExAC_ALL/g; s/;ExAC_AMR/;AAA_10_ExAC_AMR/g; s/;ExAC_EAS/;AAA_11_ExAC_EAS/g; s/;ExAC_FIN/;AAA_12_ExAC_FIN/g; s/;ExAC_NFE/;AAA_13_ExAC_NFE/g; s/;ExAC_OTH/;AAA_14_ExAC_OTH/g; s/;ExAC_SAS/;AAA_15_ExAC_SAS/g; s/;esp6500si_all/;AAA_16_esp6500si_all/g; s/;clinvar_20150629/;AAA_17_clinvar_20150629/g; s/;cosmic70/;AAA_18_cosmic70/g; s/;Polyphen2_HDIV_pred/;AAA_19_Polyphen2_HDIV_pred/g; s/;Polyphen2_HDIV_score/;AAA_20_Polyphen2_HDIV_score/g; s/;Polyphen2_HVAR_pred/;AAA_21_Polyphen2_HVAR_pred/g; s/;Polyphen2_HVAR_score/;AAA_22_Polyphen2_HVAR_score/g; s/;SIFT_pred/;AAA_23_SIFT_pred/g; s/;SIFT_score/;AAA_24_SIFT_score/g; s/;MutationAssessor_pred/;AAA_25_MutationAssessor_pred/g; s/;MutationAssessor_score/;AAA_26_MutationAssessor_score/g; s/;MutationTaster_pred/;AAA_27_MutationTaster_pred/g; s/;MutationTaster_score/;AAA_28_MutationTaster_score/g; s/;LR_pred/;AAA_29_LR_pred/g; s/;LR_score/;AAA_30_LR_score/g; s/;FATHMM_pred/;AAA_31_FATHMM_pred/g; s/;FATHMM_score/;AAA_32_FATHMM_score/g; s/;GERP++_RS/;AAA_33_GERP++_RS/g; s/;genomicSuperDups/;AAA_34_genomicSuperDups/g; s/;phastConsElements46way/;AAA_35_phastConsElements46way/g; s/;phyloP100way_vertebrate/;AAA_36_phyloP100way_vertebrate/g; s/;phyloP46way_placental/;AAA_37_phyloP46way_placental/g' < $i > $i2
	done

#removing the original VCF file (generated from ANNOVAR)
rm /data/$dir1/Data/Intensities/BaseCalls/*_SNP_INDEL_genotyped_filtered_ANNOVAR.hg19_multianno.vcf
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
echo "$PASS\n"| sudo -S cp /data/$dir1/Data/Intensities/BaseCalls/cnv_analysis/pasient_gender.txt /mnt/miseq/Medisinsk_Genetikk/Resultater/$dir1
echo "Analysis finished at:"; date +"%T" ; date +'%d/%m/%Y';
echo "$PASS\n"| sudo -S cp /data/log_record/$dir1.log /mnt/miseq/Medisinsk_Genetikk/Resultater/$dir1

#removing the experiment directory from server
cd /data/
rm -r $dir1

df -lah
#.........................................................................................................................................................................................
#.........................................................................................................................................................................................
#.........................................................................................................................................................................................


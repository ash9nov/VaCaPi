#!/bin/bash
#extracting the reads with MQ<20 from a BAM file
	#samtools view -h sample.bam| awk '($0 ~ /^$@/ || int($5)<20)' | samtoools view -Sbo sample.low_MQ.bam -
cd /data/low_MQ_analysis_31_10_2016/
for i in /data/low_MQ_analysis_31_10_2016/BAM/*.bam;
        do
                i2=${i%.*}_low_MQ.bam;
                samtools view -h $i| awk '($0 ~ /^$@/ || int($5)<20)' | samtools view -Sbo $i2 -
        done
mkdir ./low_MQ_BAM
mv ./BAM/*_low_MQ.bam ./low_MQ_BAM

#extracting the regions with mapping quality below 20 from bam to bed file
	# bamToBed -i sample.low_MQ.bam > low_MQ_sample.bed
for i in /data/low_MQ_analysis_31_10_2016/low_MQ_BAM/*_low_MQ.bam;
        do
                i2=${i%.*}.bed;
		bamToBed -i $i > $i2
	done
mkdir ./low_MQ_BED
mv ./low_MQ_BAM/*_low_MQ.bed ./low_MQ_BED

#annotation of the extracted bed file
	# bedtools intersect -a low_MQ_sample.bed -wa -wb -b /data/Data/Target_with_exon_and_gene_info.bed > test2.bed
for i in /data/low_MQ_analysis_31_10_2016/low_MQ_BED/*_low_MQ.bed
	do
		i2=${i%.*}_annotated.bed;
		bedtools intersect -a $i -wa -wb -b /data/Data/Target_with_exon_and_gene_info.bed| cut -f 1-6,10,11 > $i2
	done
mkdir ./low_MQ_BED_annotated
mv ./low_MQ_BED/*annotated.bed ./low_MQ_BED_annotated

# counting the number of reads in each gene for each sample
	# cut -f 8 ./D1361-16_S3_sorted_dedup_RG_IndReAl_Baserecal_low_MQ_annotated.bed |sort |uniq -c |sort -nr
for i in /data/low_MQ_analysis_31_10_2016/low_MQ_BED_annotated/*annotated.bed
        do
                i2=${i%.*}_read_frequency_per_gene.csv;
		cut -f 8 $i |sort |uniq -c |sort -nr > $i2
	done
mkdir ./read_frequency_per_gene
mv ./low_MQ_BED_annotated/*_read_frequency_per_gene.csv ./read_frequency_per_gene
# counting the number of reads in each gene for ALL samples together
cut -f 8 ./low_MQ_BED_annotated/* |sort |uniq -c |sort -nr > ./read_frequency_per_gene/ALL_sample_together_read_frequency_per_gene.csv;

#counting the number of reads per exon per gene for top 15 genes
#selecting the top 15 genes:
a=$(awk '{print $2}' ./read_frequency_per_gene/ALL_sample_together_read_frequency_per_gene.csv| head -15)

for i in /data/low_MQ_analysis_31_10_2016/low_MQ_BED_annotated/*annotated.bed
	do
		i2=${i%.*}_top_15_gene.csv;
		grep "$a" $i | grep -vE '(HNF1A|SDHAF2)' > $i2
	done
mkdir ./top_15_gene_with_low_MQ
mv ./low_MQ_BED_annotated/*_top_15_gene.csv ./top_15_gene_with_low_MQ

#counting the number of reads per exon for all samples
for i in /data/low_MQ_analysis_31_10_2016/low_MQ_BED_annotated/*annotated.bed
	do
		i2=${i%.*}_exon_count.csv;
		for (( c=1; c<=100; c++ ));
                         do
				grep "_Exon_${c}_" $i |wc -l
			done
	done
























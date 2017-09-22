#NGS Analsysis Steps:

.......................................................................................
#Step1:
#Getting the refrence genome:
# get the reference genome form 
url: ftp://ftp.broadinstitute.org/bundle/2.8/hg19/
#file_name ucsc.hg19.fasta.gz 
URL: ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta.gz
# .dict file unnecessery to download as it has different directory path, better to generate own .dict file
# .fai index file is also available and can be downloaded as well with fasta file.
# Also get the md5 checksum file.
#check the downlaoded file md5 key (make sure the address of file in md5 key file is corrected to the same directory of the file.)
> md5sum key_file_name.md5
.......................................................................................
#Step2:
#Quality check of FASTQ files:
> fastqc sample.fastq
#!Note: the report files generated are stored in the folder with the sample name with fastqc extension

.......................................................................................
#Step3: Preparing a reference for use with BWA and GATK
		# https://www.broadinstitute.org/gatk/guide/article?id=2798
#a: Creating the BWA index from reference genome
	> bwa index -a bwtsw ucsc.hg19.fasta
	#!Note here bwtsw option is for longe genome, so for human & -p is to replicate teh prefix name index files
	#This creates a collection of files used by BWA to perform the alignment.

#b: Generate the fasta file index
	> samtools faidx hg19.fa 
	# This creates a file called reference.fa.fai, with one record per line for each of the contigs in the FASTA reference file. Each record is composed of the contig name, size, location, basesPerLine and bytesPerLine.

#c: Generate the sequence dictionary
	> java -jar ~/my_tools/picard-tools-1.140/picard.jar CreateSequenceDictionary \
    	REFERENCE=hg19.fa \ 
    	OUTPUT=hg19.dict
    # This creates a file called hg19.dict formatted like a SAM header, describing the contents of reference FASTA file.

.......................................................................................
#Step4:
#Alignment of the short reads with the Reference genome:

#Alignment of paired-end files:
> bwa mem -M -t 8 ref.fa read1.fq read2.fq > aln.sam
# The -M flag causes BWA to mark shorter split hits as secondary (essential for Picard compatibility).
# -t 8 threads are the most optimal (as checked with 1,4,8 & 16 threads)

	#running the program parallely
	> find ./  -name "*.fastq.gz" | grep -v _R2_001.fastq.gz |sort | sed 's/_R1_001.fastq.gz//' | parallel bwa mem -M ref.fa {}_R1_001.fastq.gz {}_R2_001.fastq.gz '>' '{}'.sam 
		#Note! don't write anyaddtional output name before {} in the output file

.......................................................................................
#Step5: Convert to BAM, sort and mark duplicates
#a: Convert to BAM and sorting
> java -jar ~/my_tools/picard-tools-1.140/picard.jar  \
	SortSam  \ 
    INPUT=aligned_reads.sam  \ 
    OUTPUT=sorted_reads.bam  \ 
    SORT_ORDER=coordinate
    #Running the program parallely
    > find ./ -name "*.sam"| sort | sed 's/_L001.sam//' | parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar SortSam INPUT= {}_L001.sam OUTPUT= {}_sorted.bam SORT_ORDER=coordinate

#b: Marking the Duplicate
> java -jar ~/my_tools/picard-tools-1.140/picard.jar \
	MarkDuplicates  \ 
    INPUT=sorted_reads.bam  \ 
    OUTPUT=dedup_reads.bam  \
    METRICS_FILE=metrics.txt
	#Running the program parallely
    > find ./ -name "*_sorted.bam"| sort | sed 's/_sorted.bam//' | parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar MarkDuplicates INPUT= {}_sorted.bam OUTPUT= {}_sorted_dedup.bam METRICS_FILE= {}_metrics.txt

.......................................................................................
# NOTE!This is not part of GATK varient calling step
	#calculation of Depth of Coverage:
		> samtools depth -b /data/QC_analysis/NexteraRapidCapture-71370-targeted-regions_v4.bed sorted_dedup.bam > 1_48_Depth_of_Coverage.txt 

.......................................................................................
#Step6: Indexing the BAM file 

> java -jar picard.jar BuildBamIndex  \ 
    INPUT=dedup_reads.bam 
#Index files are required to visualize the sequence in the genome browsers (e.g. IGV)
			..........................
			## to view the bam files
			>samtools view -h read.bam
			..........................
.......................................................................................
#Step7: Indel Realignment: (It only makes a difference for indels anyhow)
#Note! Realignment around indels helps improve the accuracy of several  of the  downstream processing steps
#a: Realignment target creator: Identify what regions need to be realigned
# preprocessing step to find intervals that may need realignment.
		#URL: https://docs.google.com/file/d/0B2dK2q40HDWeLTFzNndsNDBuVms/preview
#pre-requsit: .dist file and fasta indax file (create using SAMtools)

> java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator \
	   -R ucsc.hg19.fasta \
	   -I original.bam \
	   -known indels.vcf \
	   -o realigner.intervals
	   #input BAM file not necessery if processing only at known indels
	   #using a list of known indels will both speed up processign and improve accuracy, but
	   # required.
	   # http://gatkforums.broadinstitute.org/discussion/1874/effects-of-dbsnp-in-the-step-of-indel-realignment
	   # http://gatkforums.broadinstitute.org/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it
	   #Running the program parallel		
	   > find /data/Blodbak_1_48_analysis/withRG_sorted_dedup/ -name "*_sorted_dedup_RG.bam"| sort | sed 's/_sorted_dedup_RG.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ucsc.hg19.fasta -I {}_sorted_dedup_RG.bam -known /data/Data/GATK_resources/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -o {}_realigner.intervals
	   # 12 sample at a time

	   		#TROUBLE SHOOTING: 
	   		#IF we get Errors about missing read group (RG) information
	   				
	   				#Then add the HEADER read group to the BAM file.
	   				> java -jar picard.jar AddOrReplaceReadGroups \
    		  				I= reads_without_RG.bam \
    						O=  reads_with_RG.bam \
   				    		SORT_ORDER=coordinate \  #default as input file
    						RGID=foo \  #default is 1, and can not be 'null'
    						RGLB=bar \
    						RGPL=illumina \
    						RGPU=illumina_miSEQ #can not be 'null'
    						RGSM=Sample1 \
    						CREATE_INDEX=True
    				#or
    				#> java -jar ~/my_tools/picard-tools-1.140/picard.jar AddOrReplaceReadGroups I=dedup_sorted_aln_t8.bam O=RG_dedup_sorted_aln_t8.bam SORT_ORDER=coordinate RGID= ash9nov RGLB=bar RGPL=illumina RGPU=illumina_miSEQ RGSM=Sample_SS_0284 CREATE_INDEX=True

    				#Running the program parallel
    				> find ./ -name "*_sorted_dedup.bam"| sort | sed 's/_sorted_dedup.bam//'| parallel java -jar ~/my_tools/picard-tools-1.140/picard.jar AddOrReplaceReadGroups I= {}_sorted_dedup.bam O= {}_sorted_dedup_RG.bam SORT_ORDER=coordinate RGID= {} RGLB=bar RGPL=illumina RGPU=illumina_miSEQ RGSM= {} CREATE_INDEX=True


#b: Indel Realigner: perform the actual alignment.
> java -jar GenomeAnalysisTK.jar -T IndelRealigner \
	  -R ucsc.hg19.fasta
	  -I original.bam \
	  -known indels.vcf \
	  -targetIntervals realigner.intervals \
	  -o realigned.bam
#or
#> java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I RG_dedup_sorted_aln_t8.bam -targetIntervals realigner.intervals -o IndelRealigned_RG_dedup_sorted_aln_t8.bam
#Running the progeam parallel
> find /data/Blodbak_1_48_analysis/withRG_sorted_dedup/ -name "*_sorted_dedup_RG.bam"| sort | sed 's/_sorted_dedup_RG.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R ucsc.hg19.fasta -I {}_sorted_dedup_RG.bam -known /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -targetIntervals {}_realigner.intervals -o {}_sorted_dedup_RG_IndReAl.bam
# 16 sample at a time
	  # Must use the same input file(s) used in RealignerTargetCreator step
	  # processing options
	  	 #- Only at known indels: much faster, accurate for ~90-95% of indels
	  	 #- At indels seen in the original BAM alignments: the recommended mode!!!!!
	  	 #- Using full Smith-Waterman realignment: most accurate, but heavy computational
	  	 #  cost and not really necessary with the new techs
.......................................................................................
#Step:8: Base Quality Score Recalibration:
		# URL: https://docs.google.com/file/d/0B2dK2q40HDWeZk1rMXpTYmZzTXc/preview

#a: BaseRecalibrator: Mode the error modes and recalibrate qualities. producing the recalibration table
> java -jar GenomeAnalysisTK.jar -T BaseRecalibrator \
		-R ucsc.hg19.fasta \
		-I realigned.bam \
		-knownSites dbsnp137.vcf \ # its just example.
		-knownSites gold.standard.indels.vcf \ # its just an example.
		-o recal.table
		#Running the progeam parallel
		> find /data/Blodbak_1_48_analysis/withRG_sorted_dedup/ -name "*_sorted_dedup_RG_IndReAl.bam"| sort | sed 's/_sorted_dedup_RG_IndReAl.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl.bam -knownSites /data/Data/GATK_resources/dbsnp_138.hg19.vcf -knownSites /data/Data/GATK_resources/1000G_phase1.indels.hg19.sites.vcf -o {}_recal.table
	   	# 8 sample at a time
		
		#NOTE! make sure the knownSites and reference genome have the contigs in same order.
		#web resources for knownSites:
		# http://gatkforums.broadinstitute.org/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it
		# ftp://ftp.broadinstitute.org/bundle/2.8/hg19/
		# possible error can be reordering of the contigs in reference in same order to the known site contigs.
		# solution: https://www.broadinstitute.org/gatk/guide/article?id=1328
#b: PrintReads: Write the recalibrated data to file.
> java -jar GenomeAnalysisTK.jar -T PrintReads \
		-R ucsc.hg19.fasta \
		-I realigned.bam \
		-BQSR recal.table \   #Original qualities can be retaind with OQ tag.
		-o recal.bam
		# Creates a new bam file using the input table generated previously which has exquisitely accurate base substitution, insertion, and deletion quality scoures
		#Running the progeam parallel
		> find /data/Blodbak_1_48_analysis/withRG_sorted_dedup/ -name "*_sorted_dedup_RG_IndReAl.bam"| sort | sed 's/_sorted_dedup_RG_IndReAl.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T PrintReads -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl.bam -BQSR {}_recal.table -o {}_sorted_dedup_RG_IndReAl_Baserecal.bam

#c: To plot the before/after plot, we need to recalibrate the resultant recalibrated bam file.
	# so we will repeter the step #a once more with the recalibrated bam file.
> java -jar GenomeAnalysisTK.jar -T BaseRecalibrator \
		-R ucsc.hg19.fasta \
		-I realigned.bam \
		-knownSites dbsnp137.vcf \ # its just example.
		-knownSites gold.standard.indels.vcf \ # its just an example.
		-BQSR recal.table \
		-o after_recal.table

#d: AnalyzeCovariates: 
	#Note! It requires the R ​gsalib package)

	#Plot a single recalibration table
>java -jar GenomeAnalysisTK.jar \
      -T AnalyzeCovariates \
      -R ucsc.hg19.fastasta \
      -BQSR recal.table \
      -plots BQSR.pdf
	
	# Plot before (first pass) and after (second pass) recalibration tables to compare them 
> java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates \
		-R ucsc.hg19.fasta \
		-before recal.table \
		-after after_recal.table \
		-plots recal_plots.pdf   # there is option to keep the internediate.csv file used for plotting, if we want to play with the plot data.
		# Error: mostly due to Rscript. Make sure for right version of R and right packages example: ggplots, ggpots2, gsalib.
# Note! Post-recalibrated quality score should fit the empirically-derived quality score very well; no obvious systematic biases should remain..
.......................................................................................
					
#Step:9:					### variant_CALLING ###

# URL: https://docs.google.com/file/d/0B2dK2q40HDWeQUFYUFRmM1hhRUE/preview
# URL: https://docs.google.com/file/d/0B2dK2q40HDWeYzVTUGs0bjM3M1E/preview
# As per GATK people, Joint variant discovery (MULTI-SAMPLE variant DISCOVERY) is more informative
#variant calling of individual samples will miss important informations.


#UnifiedGenotyper: Calls SNPs & Indels on per-locus basis

#HaplotypeCaller: The program determines which regions of the genome it needs to operate on, based on the presence of significant evidence for variation.
				# so it Calls SNPs & Indels simeltenously via re-assembly of haplotypes in an active region.
#URL: https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php

	#- Calls SNPs and indels simultaneously.
	#- Performs local re-assembly to Identify haplotypes.

#Single-sample all-sites calling on DNAseq (for `-ERC GVCF` cohort analysis workflow)
> java -jar GenomeAnalysisTK.jar -T HaplotypeCaller \
		-R ucsc.hg19.fasta \
		-I sample1.bam  # analysis ready samples
		-o output.raw.snps.indels.g.vcf \  #extension must be in .vcf, NOT in .gvcf
		-ERC GVCF \  # --emitRefConfidence GVCF\ , also available in BP_RESOLUTION
		--variant_index_type LINEAR \   #variant index arguments are related to file compression
		--variant_index_parameter 128000 \
		--dbsnp dbsnp137.vcf \ #optional
		-L /data/Data/NexteraRapidCapture-71370-targeted-regions_v4.bed #targets.intervals_list \ #optional, URL:  http://gatkforums.broadinstitute.org/discussion/4133/when-should-i-use-l-to-pass-in-a-list-of-intervals#latest
		#-l INFO #for turning off the warning
#running paralle program on all samples
> find /data/Blodbak_1_48_analysis/base_recalibration/ -name "*_sorted_dedup_RG_IndReAl_Baserecal.bam"| sort | sed 's/_sorted_dedup_RG_IndReAl_Baserecal.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl_Baserecal.bam -o {}_raw_SNP_INDEL.g.vcf -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000 --dbsnp /data/Data/GATK_resources/dbsnp_138.hg19.vcf 

# running the haplotypeCaller using BED file
> find /data/Blodbak_1_48_analysis/base_recalibration/ -name "*_sorted_dedup_RG_IndReAl_Baserecal.bam"| sort | sed 's/_sorted_dedup_RG_IndReAl_Baserecal.bam//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ucsc.hg19.fasta -I {}_sorted_dedup_RG_IndReAl_Baserecal.bam -o {}_raw_SNP_INDEL.g.vcf -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000 --dbsnp /data/Data/GATK_resources/dbsnp_138.hg19.vcf -L /data/Data/NexteraRapidCapture-71370-targeted-regions_v4.bed


#GenotypeGVCFs: performs joint genotyping on all samples together:
> java -jar GenomeAnalysisTK.jar -T GenotypeGVCFs \ 
		-R ucsc.hg19.fasta \
		-V sample1.g.vcf \
		-V sample2.g.vcf \
		-V sampleN.g.vcf \
		-o output.vcf
		# if >200 samples, combine in batches first using CombineGVCFs
		> java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ucsc.hg19.fasta -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak1_S1_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak2_S2_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak3_S3_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak4_S4_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak5_S5_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak6_S6_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak7_S7_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak8_S8_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak9_S9_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak10_S10_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak11_S11_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak12_S12_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak13_S13_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak14_S14_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak15_S15_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak16_S16_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak17_S17_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak18_S18_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak19_S19_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak20_S20_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak21_S21_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak22_S22_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak23_S23_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbak24_S24_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-25_S1_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-26_S2_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-27_S3_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-28_S4_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-29_S5_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-30_S6_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-31_S7_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-32_S8_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-33_S9_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-34_S10_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-35_S11_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-36_S12_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-37_S13_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-38_S14_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-39_S15_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-40_S16_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-41_S17_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-42_S18_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-43_S19_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-44_S20_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-45_S21_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-46_S22_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-47_S23_raw_SNP_INDEL_L.g.vcf -V /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/Blodbank-48_S24_raw_SNP_INDEL_L.g.vcf -o /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/1_to_48_genotype.vcf

		#running parallel for individual sample genotyping
		> find /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/ -name "*_raw_SNP_INDEL_L.g.vcf"| sort | sed 's/_raw_SNP_INDEL_L.g.vcf//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ucsc.hg19.fasta -V {}_raw_SNP_INDEL_L.g.vcf -o {}_SNP_INDEL_genotyped.vcf
......................................................................................

#Step10: variant Quality Score Recalibration: 
			#Assigning accurate confidance scores to each putative mutation call
			#Buiding a model of what true genetic variation looks like will allow us to rank-order variants based on their likelihood of being real.

	#a: VariantRecalibrator:
			#Note! SNPs and Indels must be recalibrated separately!
	> java -jar GenomeAnalysisTK.jar -T VariantRecalibrator \
			-R ucsc.hg19.fasta \
			-input raw_SNPs_Indels.vcf \  #its -input not -I
			-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /data/Data/GATK_resources/hapmap_3.3.hg19.sites.vcf \
   			-resource:omni,known=false,training=true,truth=false,prior=12.0 /data/Data/GATK_resources/1000G_omni2.5.hg19.sites.vcf \
   			-resource:1000G,known=false,training=true,truth=false,prior=10.0 /data/Data/GATK_resources/1000G_phase1.snps.high_confidence.hg19.sites.vcf \ #issue with the file, contig not in same order.
			-resource:dbsnp,known=true,training=false,truth=false,prior=6.0 /data/Data/GATK_resources/dbsnp_138.hg19.vcf \
			-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR --maxGaussians 4  # -an InbreedingCoeff \ #The InbreedingCoeff statistic is a population-level calculation that is only available with 10 or more samples. If you have fewer samples you will need to omit that particular filter statement.
			-mode SNP \   #like wise for Indels
			-recalFile raw_SNPs.recal \  #likewise raw_INDELs.recal
			-tranchesFile raw_SNPs.tranches \ #likewise raw_INDELs.tranches
			-rscriptFile recal.plots.R
			#running parallel for individual sample genotyping
			> find /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/ -name "*_SNP_INDEL_genotyped.vcf" | sort | sed 's/_SNP_INDEL_genotyped.vcf//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T VariantRecalibrator -R ucsc.hg19.fasta -input {}_SNP_INDEL_genotyped.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /data/Data/GATK_resources/hapmap_3.3.hg19.sites.vcf -resource:omni,known=false,training=true,truth=false,prior=12.0 /data/Data/GATK_resources/1000G_omni2.5.hg19.sites.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 /data/Data/GATK_resources/1000G_phase1.snps.high_confidence.hg19.sites.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 /data/Data/GATK_resources/dbsnp_138.hg19.vcf -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR --maxGaussians 4 -mode SNP -recalFile {}_raw_SNPs.recal -tranchesFile {}_raw_SNPs.tranches -rscriptFile {}_recal.plots.R
			#Indelspecific
			--maxGaussians 4 \   #set it to "2" if face error message "NaN LOD value assigned. Clustering with this few variants and these annotations is unsafe. 
   			-resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.b37.sites.vcf \
   			-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.b37.vcf \
   			-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \ #The InbreedingCoeff statistic is a population-level calculation that is only available with 10 or more samples. If you have fewer samples you will need to omit that particular filter statement.
   			-mode INDEL \
# http://gatkforums.broadinstitute.org/discussion/1259/what-vqsr-training-sets-arguments-should-i-use-for-my-specific-project
	#b: ApplyRecalibration:
	>java -jar  GenomeAnalysisTK.jar -T ApplyRecalibration \
			-R ucsc.hg19.fasta \
			-input raw_SNPs_Indels.vcf \   #its -input not -I
			-mode SNP \   # -mode INDEL 
			-recalFile raw_SNPs.recal \
			-tranchesFile raw_SNPs.tranches \
			-o recal_SNPs.vcf \
			-ts_filter_level 99.5    # for indel 99.0
			#running parallel for individual sample genotyping
			> find /data/Blodbak_1_48_analysis/Haplotypecaller_with_BED_option/ -name "*_raw_SNPs.recal" | sort | sed 's/_raw_SNPs.recal//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T ApplyRecalibration -R ucsc.hg19.fasta -input {}_SNP_INDEL_genotyped.vcf -mode SNP -recalFile {}_raw_SNPs.recal -tranchesFile {}_raw_SNPs.tranches -o {]_SNPs_genotyped_recalibrated.vcf -ts_filter_level 99.5

......................................................................................

#Step11: variantFiltration:

#!NOTE this commands needs the .idx files of each sample run the command in same directory

java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T variantFiltration \
-R ucsc.hg19.fasta \
--variant snp.recalibrated.vcf \
-o snp.recalibrated.filtered.vcf \
--clusterWindowSize 10 \
--filterExpression "MQ0 >= 4 && ((MQ0/(1.0 * DP)) > 0.1)" \
--filterName "HARD_TO_VALIDATE" \
--filterExpression "DP < 5 " \
--filterName "LowCoverage" \
--filterExpression "QUAL < 30.0 " \
--filterName "VeryLowQual" \
--filterExpression "QUAL > 30.0 && QUAL < 50.0 " \
--filterName "LowQual" \
--filterExpression "QD < 1.5 " \
--filterName "LowQD" \
--filterExpression "SB > -10.0 " \
--filterName "StrandBias"

#Running the code parallel
> find /data/Data/hg19/ -name "*_SNP_INDEL_genotyped.vcf" | sort | sed 's/_SNP_INDEL_genotyped.vcf//' | parallel java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T VariantFiltration -R ucsc.hg19.fasta --variant {}_SNP_INDEL_genotyped.vcf -o {}_SNP_INDEL_genotyped.filtered.vcf --clusterWindowSize 10 --filterExpression "MQ0 >= 4 && ((MQ0/(1.0 * DP)) > 0.1)" --filterName "HARD_TO_VALIDATE" --filterExpression "DP < 5 " --filterName "LowCoverage" --filterExpression "QUAL < 30.0 " --filterName "VeryLowQual" --filterExpression "QUAL > 30.0 && QUAL < 50.0 " --filterName "LowQual" --filterExpression "QD < 1.5 " --filterName "LowQD" --filterExpression "SB > -10.0 " --filterName "StrandBias"
#Note! its not working parallaly.

#Running the code sequencely
> 

......................................................................................

#Step12: Annotation of variants

# Using ANNOVAR:

#for single sample:
perl table_annovar.pl $i humandb/ -buildver hg19 -out $i2 -remove -protocol refGene,phastConsElements46way,genomicSuperDups,1000g2015aug_eur,snp138,ljb26_all,esp6500si_all,exac03,clinvar_20150629,cosmic70 -operation g,r,r,f,f,f,f,f,f,f -nastring . -vcfinput






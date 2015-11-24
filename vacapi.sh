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
	> bwa index -a bwtsw hg19.fa
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

.......................................................................................
#Step5: Convert to BAM, sort and mark duplicates
#a: Convert to BAM and sorting
> java -jar ~/my_tools/picard-tools-1.140/picard.jar  \
	SortSam  \ 
    INPUT=aligned_reads.sam  \ 
    OUTPUT=sorted_reads.bam  \ 
    SORT_ORDER=coordinate
#b: Marking the Duplicate
> java -jar ~/my_tools/picard-tools-1.140/picard.jar \
	MarkDuplicates  \ 
    INPUT=sorted_reads.bam  \ 
    OUTPUT=dedup_reads.bam  \
    METRICS_FILE=metrics.txt

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
	   -R hg19.fa \
	   -I original.bam \
	   -known indels.vcf \
	   -o realigner.intervals
	   #input BAM file not necessery if processing only at known indels
	   #using a list of known indels will both speed up processign and improve accuracy, but
	   # required.
	   # http://gatkforums.broadinstitute.org/discussion/1874/effects-of-dbsnp-in-the-step-of-indel-realignment
	   # http://gatkforums.broadinstitute.org/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it
	   		
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

#b: Indel Realigner: perform the actual alignment.
> java -jar GenomeAnalysisTK.jar -T IndelRealigner \
	  -R hg19.fa
	  -I original.bam \
	  -known indels.vcf \
	  -targetIntervals realigner.intervals \
	  -o realigned.bam
#or
#> java -jar ~/my_tools/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R hg19.fa -I RG_dedup_sorted_aln_t8.bam -targetIntervals realigner.intervals -o IndelRealigned_RG_dedup_sorted_aln_t8.bam
	  # Must use the same input file(s) used in RealignerTargetCreator step
	  # processing options
	  	 #- Only at known indels: much faster, accurate for ~90-95% of indels
	  	 #- At indels seen in the original BAM alignments: the recommended mode!!!!!
	  	 #- Using full Smith-Waterman realignment: most accurate, but heavu computational
	  	 #  cost and not really necessary with the new techs
.......................................................................................
#Step:8: Base Quality Score Recalibration:
		# URL: https://docs.google.com/file/d/0B2dK2q40HDWeZk1rMXpTYmZzTXc/preview

#a: BaseRecalibrator: Mode the error modes and recalibrate qualities. producing the recalibration table
> java -jar GenomeAnalysisTK.jar -T BaseRecalibrator \
		-R hg19.fa \
		-I realigned.bam \
		-knownSites dbsnp137.vcf \ # its just example.
		-knownSites gold.standard.indels.vcf \ # its just an example.
		-o recal.table
		#NOTE! make sure the knownSites and reference genome have the contigs in same order.
		#web resources for knownSites:
		# http://gatkforums.broadinstitute.org/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it
		# ftp://ftp.broadinstitute.org/bundle/2.8/hg19/
		# possible error can be reordering of the contigs in reference in same order to the known site contigs.
		# solution: https://www.broadinstitute.org/gatk/guide/article?id=1328
#b: PrintReads: Write the recalibrated data to file.
> java -jar GenomeAnalysisTK.jar -T PrintReads \
		-R hg19.fa \
		-I realigned.bam \
		-BQSR recal.table \   #Original qualities can be retaind with OQ tag.
		-o recal.bam
		# Creates a new bam file using the input table generated previously which has exquisitely accurate base substitution, insertion, and deletion quality scoures

#c: To plot the before/after plot, we need to recalibrate the resultant recalibrated bam file.
	# so we will repeter the step #a once more with the recalibrated bam file.
> java -jar GenomeAnalysisTK.jar -T BaseRecalibrator \
		-R hg19.fa \
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
      -R hg19.fasta \
      -BQSR recal.table \
      -plots BQSR.pdf
	
	# Plot before (first pass) and after (second pass) recalibration tables to compare them 
> java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates \
		-R hg19.fa \
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
		-R hg19.fa \
		-I sample1.bam  # analysis ready samples
		-o output.raw.snps.indels.g.vcf \  #extension must be in .vcf, NOT in .gvcf
		-ERC GVCF \  # --emitRefConfidence GVCF\ , also available in BP_RESOLUTION
		--variant_index_type LINEAR \   #variant index arguments are related to file compression
		--variant_index_parameter 128000 \
		--dbsnp dbsnp137.vcf \ #optional
		-L targets.intervals_list \ #optional, URL:  http://gatkforums.broadinstitute.org/discussion/4133/when-should-i-use-l-to-pass-in-a-list-of-intervals#latest

#GenotypeGVCFs: performs joint genotyping on all samples together:
> java -jar GenomeAnalysisTK.jar -T GenotypeGVCFs \ 
		-R hg19.fa \
		-V sample1.g.vcf \
		-V sample2.g.vcf \
		-V sampleN.g.vcf \
		-o output.vcf
		# if >200 samples, combine in batches first using CombineGVCFs

......................................................................................

#Step10: variant Quality Score Recalibration: 
			#Assigning accurate confidance scores to each putative mutation call
			#Buiding a model of what true genetic variation looks like will allow us to rank-order variants based on their likelihood of being real.

	#a: VariantRecalibrator:
			#Note! SNPs and Indels must be recalibrated separately!
	> java -jar GenomeAnalysisTK.jar -T VariantRecalibrator \
			-R hg19.fa \
			-input raw_SNPs_Indels.vcf \  #its -input not -I
			-resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg19.sites.vcf  \
   			-resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.b37.sites.vcf \
   			-resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.vcf \ #issue with the file, contig not in same order.
			-resource:dbsnp,known=true,training=false,truth=false,prior=6.0 dbsnp_135.b37.vcf \
			-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR --maxGaussians 4  # -an InbreedingCoeff \ #The InbreedingCoeff statistic is a population-level calculation that is only available with 10 or more samples. If you have fewer samples you will need to omit that particular filter statement.
			-mode SNP \   #like wise for Indels
			-recalFile raw_SNPs.recal \  #likewise raw_INDELs.recal
			-tranchesFile raw_SNPs.tranches \ #likewise raw_INDELs.tranches
			-rscriptFile recal.plots.R

			#Indelspecific
			--maxGaussians 4 \   #set it to "2" if face error message "NaN LOD value assigned. Clustering with this few variants and these annotations is unsafe. 
   			-resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.b37.sites.vcf \
   			-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.b37.vcf \
   			-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \ #The InbreedingCoeff statistic is a population-level calculation that is only available with 10 or more samples. If you have fewer samples you will need to omit that particular filter statement.
   			-mode INDEL \
# http://gatkforums.broadinstitute.org/discussion/1259/what-vqsr-training-sets-arguments-should-i-use-for-my-specific-project
	#b: ApplyRecalibration:
	>java -jar  GenomeAnalysisTK.jar -T ApplyRecalibration \
			-R hg19.fa \
			-input raw_SNPs_Indels.vcf \   #its -input not -I
			-mode SNP \   # -mode INDEL 
			-recalFile raw_SNPs.recal \
			-tranchesFile raw_SNPs.tranches \
			-o recal_SNPs.vcf \
			-ts_filter_level 99.5    # for indel 99.0




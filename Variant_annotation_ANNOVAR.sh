# Annovar annotation:

#Step1: Preparing INPUT files:
> perl convert2annovar.pl -format vcf4 example/ex2.vcf -outfile ex2.avinput -includeinfo -withzyg -comment 
	#IF there are multiple samples in the VCF file (in cases of merging the different samples BAM files or Joint genotyping)
		> perl convert2annovar.pl -format vcf4 example/ex2.vcf -outfile ex2 -allsample
	#!Note: output files will be written to ex2.<samplename>.avinput
	-includeinfo
	-withzyg #If you need the zygosity, quality and read coverage information in the output line as well, add the -withzyg argument
	-withfreq #When -withfreq is set, it will print out the allele frequency of each SNP in the VCF file, based on all samples within the file
	#-withfreq and -withzyg are mutually exclusive
	-comment
	--coverage #-coverage requires an argument
---------------------------------------------------------------------------------------------
#Gene-Based annotation:
	#Before working on gene-based annotation, a gene definition file and associated FASTA file must be downloaded into a directory if they are not already downloaded. Let's call this directory as humandb/.
	# > annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene humandb/
The gene-based annotation can be issued by the following command (by default, --geneanno -dbtype refGene is assumed):

	> perl annotate_variation.pl -out ex1 -build hg19 example/ex1.avinput humandb/ -splicing_threshold 3
		# two output files: ex1.variant_function and ex1.exonic_variant_function
		# The first file ex1.variant_function contains annotation for all variants, by adding two columns to the beginning of each input line 
		# The second output file, ex1.exonic_variant_function, contains the amino acid changes as a result of the exonic variant. 
		'--separate' # If the users want to have all functional consequences printed out (rather than just the most important one defined by the precedence above), the --separate argument should be used. In this case, several output lines may be present for each variant, representing several possible functional consequences.
		'-precedence intronic,utr5,utr3' # The "-precedence" argument can be used to fine-tune the priority of variant function.
		' -splicing_threshold 3'  #by default it is 2: splicing" in ANNOVAR is defined as variant that is within 2-bp away from an exon/intron boundary by default
		# switching the database : '-dbtype ensGene'
		# Technical Notes: Technically, the RefSeq Gene and UCSC Gene are transcript-based gene definitions. They built gene model based on transcript data, and then map the gene model back to human genomes. In comparison, Ensemble Gene and Gencode Gene are assembly-based gene definitions that attempt to build gene model directly from reference human genome. They came from different angles, trying to do the same thing: define genes in human genome.
----------------------------------------------------------------------------------------------
#Region-based annotation:  
		'--regionanno'
#conserved genomic elements:
	#these are 'phastCons' conserved elements, which means specific genomic regions that are conserved. 
> perl annotate_variation.pl -regionanno -build hg19 -out ex1 -dbtype phastConsElements46way example/ex1.avinput humandb/
	# The output is saved in the ex1.hg19_phastConsElements46way file. The first column in the output is "phastConsElements46way" indicating the type of annotation. The second column contains two pieces of information: Score and Name
	# Score is the normalized score assigned by UCSC Genome Browser, and this score range from 0 to 1000 for the sole purpose of having a standard range of values to display in browser.
	# The "Name=lod=x" is used to tell the user a name for the region 

# Segmental Duplications:
	 	( # for downloading database # perl annotate_variation.pl -build hg19 -downdb genomicSuperDups humandb/)
 > perl annotate_variation.pl -regionanno -build hg19 -out ex1 -dbtype genomicSuperDups example/ex1.avinput humandb/

# Identifying variants in regions specified in BED files
|
|
----------------------------------------------------------------------------------------------
#Filter-based annotation:  
		'-filter'
	# 1000 Genomes Project (2015 Aug) annotations
		#getting Database 
			# > perl annotate_variation.pl -downdb 1000g2015aug humandb -buildver hg19
		> perl annotate_variation.pl -filter -dbtype 1000g2015aug_eur -buildver hg19 -out ex1 example/ex1.avinput humandb/
	# dbSNP annotations
		# GETTING DATABASE
			# > perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar snp138 humandb
		> perl annotate_variation.pl -filter -out ex1 -build hg19 -dbtype snp138 example/ex1.avinput 
		# Two output files are generated. The ex1.hg19_snp138_filtered file contains SNPs not in dbSNP. The ex1.hg19_snp130_dropped file contains variants that are annotated in dbSNP, and print out their rs identifiers (as the second column)
	
	# LJB* (dbNSFP) non-synonymous variants annotation.
		# GETTING DATABASE
			# > perl annotate_variation.pl -downdb -webfrom annovar -buildver hg19 dbnsfp30a humandb/
	-------------------------
		#not working
		Step1: 
		#building the scores:
	 	> table_annovar.pl ex1.avinput humandb/ -protocol dbnsfp30a -operation f -build hg19 -nastring . 
	 	# The command above takes an input file and generates 20 different scores and predictions for all the non-synonymous variants in the file.
	 	# The output include various scores in the following columns: SIFT_score SIFT_pred Polyphen2_HDIV_score Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred LRT_score LRT_pred MutationTaster_score MutationTaster_pred MutationAssessor_score MutationAssessor_pred FATHMM_score FATHMM_pred PROVEAN_score PROVEAN_pred VEST3_score CADD_raw CADD_phred DANN_score fathmm-MKL_coding_score fathmm-MKL_coding_pred MetaSVM_score MetaSVM_pred MetaLR_score MetaLR_pred integrated_fitCons_score integrated_confidence_value GERP++_RS phyloP7way_vertebrate phyloP20way_mammalian phastCons7way_vertebrate phastCons20way_mammalian SiPhy_29way_logOdds.
	 	# If a variant does not have a score (for example, for intronic variant), the corresponding position will be denoated by a period (which is specified by the '-nastring' argument above.
	-------------------------
	 	Step2:
	 	# SIFT: annotation based on score:
	 		> perl annotate_variation.pl -filter -dbtype ljb26_sift -buildver hg19 -out ex1 example/ex1.avin -otherinfo

		# - PolyPhen 2 annotation
			> perl annotate_variation.pl -filter -dbtype ljb26_pp2hvar -buildver hg19 -out ex1 example/ex1.avinput humandb/ -otherinfo
			# There are three possible predictions: "D" ("porobably damaging"), "P" ("possibly damaging") and "B" ("benign").

		# - MutationTaster annotation	
			> perl annotate_variation.pl -filter -dbtype ljb26_mt -buildver hg19 -out ex1 example/ex1.avinput humandb/ -otherinfo
			# There are four possible predictions: "A" ("disease_causing_automatic"), "D" ("disease_causing"), "N" ("polymorphism") or "P" ("polymorphism_automatic").
		# - PhyloP annotation
			>  perl annotate_variation.pl -filter -dbtype ljb26_phylop -buildver hg19 -out ex1 example/ex1.avinput humandb/ -otherinfo

	#ESP (Exome sequencing project) annotations:
		# getting Database
			# > perl annotate_variation.pl -downdb -webfrom annovar -build hg19 esp6500si_all humandb/
		> annotate_variation.pl -filter -dbtype esp6500si_all -build hg19 -out ex1 example/ex1.avinput humandb/
		# Note! To change to EA (European American) or AA (African American), just change the database name:esp6500si_ea or esp6500si_aa

	# ExAC (Exome Aggregation Consortium)annotations
		#DB download
			# > perl annotate_variation.pl -downdb -webfrom annovar -build hg19 exac03 humandb/
		> perl annotate_variation.pl -filter -build hg19 -dbtype exac03 example/ex1.avinput humandb/ -otherinfo

	# CLINVAR annotations:
		> annotate_variation.pl example/ex1.avinput humandb/ -filter -dbtype clinvar_20150629 -buildver hg19 -out ex1

	# COSMIC annotation:
		> annotate_variation.pl example/ex1.avinput humandb/ -filter -dbtype cosmic70 -buildver hg19 -out ex1



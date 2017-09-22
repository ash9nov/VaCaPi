#!/bin/bash

echo "Analysis starteded at:"; date +"%T" ; date +'%d/%m/%Y';
#writing the log file.........
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>/data/CNV_analysis/Norsk_sample/CNV_analysis.log 2>&1

# step1: Generating Sequencing-accessible regions
#cnvkit.py access /data/Data/hg19/ucsc.hg19.fasta -s 10000 -o /data/$dir1/Data/Intensities/BaseCalls/cnv_analysis/access-10kb.mn10.bed
#running all steps in single "Batch" command
cnvkit.py  batch /data/CNV_analysis/Norsk_sample/bam/*.bam --normal --targets /data/CNV_analysis/HereditaryCRC_genepanel_Final_selected.bed --fasta /data/Data/hg19/ucsc.hg19.fasta --split --access /data/Data/access-10kb.mn10.bed --output-reference /data/CNV_analysis/Norsk_sample/flat_reference.cnn -d /data/CNV_analysis/Norsk_sample
# Identifying the gender of pasients
cnvkit.py gender /data/CNV_analysis/Norsk_sample/*cns > /data/CNV_analysis/Norsk_sample/Norsk_sample_gender.txt
#changing the long names
sed -i "s/\/data\/CNV_analysis\/Norsk_sample\///g" Norsk_sample_gender.txt
sed -i "s/.cns//g" Norsk_sample_gender.txt
cut -f 1-2 Norsk_sample_gender.txt > Norsk_sample_gender_final.txt

# creating VCF files
for i in ./*.cns;
        do    
              i2=${i%.*}.CNV.vcf;
              cnvkit.py export vcf $i -i "SampleID" -o $i2
        done
# Estimate integer copy number of each segment
for i in ./*.cns;
        do
              i3=${i%.*}.call.cns;
              cnvkit.py call $i -o $i3
        done
#Show estimated integer copy number of all regions
for i in ./*.call.cns;
        do
              i4=${i%.*}.bed;
              cnvkit.py export bed $i --show all -o $i4
        done

echo "Analysis finished at:"; date +"%T" ; date +'%d/%m/%Y';

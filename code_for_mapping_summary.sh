#!/bin/bash

read dir1;
for i in /data/mapping_statistics/$dir1/*_sorted_dedup_RG_IndReAl_Baserecal.bam;
	do
		i2=${i%.*}_mapping_statistics.txt;
		samtools flagstat $i > $i2
	done

for i in /data/mapping_statistics/$dir1/*_mapping_statistics.txt;
	do
		echo "$i\n$(cat $i)" > $i # it adds the sample name on to of FLAGSTAT output
		sed -i "s/\/data\/mapping_statistics\/$dir1\///g" $i
		sed -i "s/_sorted_dedup_RG_IndReAl_Baserecal_mapping_statistics.txt//g" $i
	done
cp /data/Data/flag_property "$dir1.mapping_summary.txt"
for i in /data/mapping_statistics/$dir1/*_mapping_statistics.txt;
	do
		cut -d "+" -f 1 "$i" > "$dir1.temp" 
		paste -d "\t" "$dir1.mapping_summary.txt" "$dir1.temp" > "$dir1.map"
		mv "$dir1.map" "$dir1.mapping_summary.txt"
	done

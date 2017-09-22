#!bin/bash

read dir1;
# changing the long names
#sed -i "s/\/data\/$dir1\/Data\/Intensities\/BaseCalls\///g" *

echo "$dir1.sample_interval_summary" 
cut -f 1 "$dir1.sample_interval_summary" > "$dir1.sample_interval_mean_coverage"
cut -f 1 "$dir1.sample_interval_summary" > "$dir1.temp1"
read smpl

echo $smpl

i=1
while [ "$i" -le "$smpl" ] ; 
do
        jj=$(((6*$i)-1))
        cut -f $jj "$dir1.sample_interval_summary" > "$dir1.temp"      
	paste -d "\t" "$dir1.sample_interval_mean_coverage" "$dir1.temp" > "$dir1.op"
	mv "$dir1.op" "$dir1.sample_interval_mean_coverage"
	read -r FIRSTLINE < $dir1.temp
	paste -d "\t" "$dir1.temp1" "$dir1.temp" > "$FIRSTLINE.txt"
	awk 'FNR==NR { a[$1]=$2; next } $1 in a { print $1"\t"a[$1]"\t"$2"\t"$3 }' "$FIRSTLINE.txt" Target_ex_n_gn_mdfd_strtPlus1_merged_regions_for_CovRept_2120.bed > "$FIRSTLINE.annotated.txt" 
	awk -v x=30 -F "\t" '$2<=x {print}' "$FIRSTLINE.annotated.txt" > "$FIRSTLINE.low_coverage_region.txt"
	rm "$FIRSTLINE.txt"
        rm "$FIRSTLINE.annotated.txt"
	i=$(($i+1));
done
rm $dir1.temp*

sort "$dir1.sample_interval_mean_coverage" > bb1
sort Target_ex_n_gn_mdfd_strtPlus1_merged_regions_for_CovRept_2120.bed > bb2
join -j 1 bb1 bb2 > "$dir1.sample_interval_mean_coverage_annotated"
rm bb1
rm bb2
#join -j 1 <(sort "$dir1.sample_interval_mean_coverage") <(sort Target_ex_n_gn_mdfd_strtPlus1_merged_regions_for_CovRept_2120.bed) > "$dir1.sample_interval_mean_coverage_annotated"

#!/bin/bash

for i in ./*_genomeSummary.csv ;
	do 
		i2=${i%.*}_only_22_MMR_genes.csv;
		grep 'Gene\|MLH1\|MSH2\|MSH6\|PMS2\|MSH3\|PMS1\|MLH3\|EXO1\|RFC1\|RFC2\|RFC3\|RFC4\|RFC5\|PCNA\|LIG1\|RPA1\|RPA2\|RPA3\|POLD1\|POLD2\|POLD3\|POLD4' $i > $i2 ;
	done

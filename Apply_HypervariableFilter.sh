#!/bin/bash

for FILE in `ls *2HomSap_snps_phylop_200sp_autosomes.bed`; do
SP=`echo $FILE | sed 's/2HomSap_snps_phylop_200sp_autosomes.bed//g'`

awk -F'[\t@]' '{print $5"\t"$6-1"\t"$6"\t"$7"\t"$1"@"$2"@"$3"@"$4}' $FILE \
| sort-bed - > $SP'_mapfile'
sort-bed <(tail -n +2 '~/USS/aryn/Zoonomia/filter/'$SP'_mask_a0.1.txt') > \
'~/USS/aryn/Zoonomia/filter/'$SP'_mask_a0.1_sort.bed'

#filter with bedmaps:
bedmap --echo --echo-ref-name --delim "\t" --skip-unmapped \
$SP'_mapfile' <(grep FALSE '~/USS/aryn/Zoonomia/filter/'$SP'_mask_a0.1_sort.bed' | cut -f 1-3) | \
awk -F'[\t@]' '{print $5"\t"$6"\t"$7"\t"$8"\t"$1"@"$3"\t"$4}' \
> $SP'2HomSap_snps_phylop_200sp_autosomes_filtered.bed'

rm $SP'_mapfile'

done
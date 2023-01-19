#!/bin/bash

cd /USS/aryn/Zoonomia/snp_phylop_200


SP=$1
FILE=$SP'2HomSap_snps_filt_phylop_200sp_coding_hiqual.bed'

#make a bed of windows in human genome and sort
awk -F'[\t@]' '{print $7"\t"$8-1"\t"$8"\t"$6}' $FILE > $SP"2HomSap_snps_filt_phylop_200sp_coding_hiqual_humancoords.txt"
awk -F'[\t|]' '{print $6"\t"$7}' $FILE > $SP"2HomSap_snps_filt_phylop_200sp_coding_hiqual_humancoords.impact"
paste $SP"2HomSap_snps_filt_phylop_200sp_coding_hiqual_humancoords.txt" $SP"2HomSap_snps_filt_phylop_200sp_coding_hiqual_humancoords.impact" \
> $SP"2HomSap_snps_filt_phylop_200sp_coding_hiqual_humancoords.txt2"
mv $SP"2HomSap_snps_filt_phylop_200sp_coding_hiqual_humancoords.txt2" $SP"2HomSap_snps_filt_phylop_200sp_coding_hiqual_humancoords.txt"
rm $SP"2HomSap_snps_filt_phylop_200sp_coding_hiqual_humancoords.impact"

sort-bed $SP"2HomSap_snps_filt_phylop_200sp_coding_hiqual_humancoords.txt" | gzip > \
$SP"2HomSap_snps_filt_phylop_200sp_coding_hiqual_humancoords_sort.txt.gz"

####summarize coding substitutions in windows
#counts of missense substitutions
bedmap --echo --count --mean --delim "\t" /home/centos/USS/aryn/Zoonomia/roh/hg38_50kbwindows_sort.bed \
<(zcat $SP"2HomSap_snps_filt_phylop_200sp_coding_hiqual_humancoords_sort.txt.gz" \
| awk '$5=="missense_variant" {print $1"\t"$2"\t"$3"\t"$5"\t"$4}') > $SP'_200sp_autosomes_human50kb_miss_cnt_phylop.bed'

#counts of high-impact substitutions
bedmap --echo --count --mean --delim "\t" /home/centos/USS/aryn/Zoonomia/roh/hg38_50kbwindows_sort.bed \
<(zcat $SP"2HomSap_snps_filt_phylop_200sp_coding_hiqual_humancoords_sort.txt.gz" \
| awk '$6=="HIGH" {print $1"\t"$2"\t"$3"\t"$5"\t"$4}') > $SP'_200sp_autosomes_human50kb_high_cnt_phylop.bed'

#counts of synonymous substitutions
bedmap --echo --count --mean --delim "\t" /home/centos/USS/aryn/Zoonomia/roh/hg38_50kbwindows_sort.bed \
<(zcat $SP"2HomSap_snps_filt_phylop_200sp_coding_hiqual_humancoords_sort.txt.gz" \
| awk '$5=="synonymous_variant" {print $1"\t"$2"\t"$3"\t"$5"\t"$4}') > $SP'_200sp_autosomes_human50kb_syn_cnt_phylop.bed'

#counts of synonymous conserved substitutions
bedmap --echo --count --delim "\t" /home/centos/USS/aryn/Zoonomia/roh/hg38_50kbwindows_sort.bed \
<(zcat $SP"2HomSap_snps_filt_phylop_200sp_coding_hiqual_humancoords_sort.txt.gz" \
| awk '$5=="synonymous_variant" && $4>2.27 {print $1"\t"$2"\t"$3"\t"$5"\t"$4}') > $SP'_200sp_autosomes_human50kb_syn_cons_cnt.bed'

#counts of high-impact conserved substitutions
bedmap --echo --count --delim "\t" /home/centos/USS/aryn/Zoonomia/roh/hg38_50kbwindows_sort.bed \
<(zcat $SP"2HomSap_snps_filt_phylop_200sp_coding_hiqual_humancoords_sort.txt.gz" \
| awk '$6=="HIGH" && $4>2.27 {print $1"\t"$2"\t"$3"\t"$5"\t"$4}') > $SP'_200sp_autosomes_human50kb_high_cons_cnt.bed'

#counts of missense conserved substitutions
bedmap --echo --count --delim "\t" /home/centos/USS/aryn/Zoonomia/roh/hg38_50kbwindows_sort.bed \
<(zcat $SP"2HomSap_snps_filt_phylop_200sp_coding_hiqual_humancoords_sort.txt.gz" \
| awk '$5=="missense_variant" && $4>2.27 {print $1"\t"$2"\t"$3"\t"$5"\t"$4}') > $SP'_200sp_autosomes_human50kb_miss_cons_cnt.bed'

#combine into one bed file
paste $SP'_200sp_autosomes_human50kb_syn_cnt_phylop.bed' <(cut -f4 $SP'_200sp_autosomes_human50kb_syn_cons_cnt.bed') \
<(cut -f4-5 $SP'_200sp_autosomes_human50kb_miss_cnt_phylop.bed') <(cut -f4 $SP'_200sp_autosomes_human50kb_miss_cons_cnt.bed') \
<(cut -f4-5 $SP'_200sp_autosomes_human50kb_high_cnt_phylop.bed') <(cut -f4 $SP'_200sp_autosomes_human50kb_high_cons_cnt.bed') \
| gzip > $SP'_200sp_autosomes_human50kb_coding_phylop.bed.gz'

rm $SP'_200sp_autosomes_human50kb_syn_cnt_phylop.bed' $SP'_200sp_autosomes_human50kb_syn_cons_cnt.bed' \
$SP'_200sp_autosomes_human50kb_miss_cnt_phylop.bed' $SP'_200sp_autosomes_human50kb_miss_cons_cnt.bed' \
$SP'_200sp_autosomes_human50kb_high_cnt_phylop.bed' $SP'_200sp_autosomes_human50kb_high_cons_cnt.bed'

done


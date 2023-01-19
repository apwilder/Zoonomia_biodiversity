#!/bin/bash

SP=$1

#50kb windows with mean heterozygosity and ROH status (ROH=0, non-ROH=1) lifted to human genome (created by Ross Swofford @ Broad Inst.)
#eg:
#human_chr	human_start human_stop	sp_contig@sp_start@sp_stop@heterozygosity@ROH
#chr6	166709746	166709841	chr1@54650000@54699999@0.00157875852170793@0
#chr6	166709848	166709895	chr1@54650000@54699999@0.00157875852170793@0

zcat $SP".bed.gz" | sort-bed - > $SP"_sort.bed"

#assign phylop scores from human autosomes
for CHR in `seq 1 22`; do
bedmap --echo --max --mean --count --delim "\t" --skip-unmapped $SP"_sort.bed" \
"~/USS/aryn/Zoonomia/phylop/human_onlychr_v2_mdong_10Mb_241MAMMALS_MTDF_species_count/chr"$CHR"_positive_200sp.bed" \
> $SP"_200sp_chr"$CHR".bed"
done

cat $SP"_200sp_chr"*".bed" > $SP"_200sp_autosomes.bed"
rm $SP"_200sp_chr"*".bed"

sort-bed $SP"_200sp_autosomes.bed" | gzip > $SP'_200sp_autosomes_sort.bed.gz'
gzip $SP"_200sp_autosomes.bed" &

#summarize phylop and heterozygosity stats within 50kb windows of the human genome

#mean ROH assignment of bases in window
bedmap --echo --wmean --delim "\t" ../hg38_50kbwindows_sort.bed <(zcat $SP'_200sp_autosomes_sort.bed.gz' | \
awk -F'[\t@]' '{print $1"\t"$2"\t"$3"\t"$7"\t"$8}') | gzip > $SP'_200sp_autosomes_human50kb_wmeanROH.bed.gz'

#mean heterozygosity
bedmap --echo --wmean --delim "\t" ../hg38_50kbwindows_sort.bed <(zcat $SP'_200sp_autosomes_sort.bed.gz' | \
awk -F'[\t@]' '{print $1"\t"$2"\t"$3"\t"$8"\t"$7}') | gzip > $SP'_200sp_autosomes_human50kb_wmeanhet.bed.gz'

#number of sites that lifted over to human genome
bedmap --echo --sum --delim "\t" ../hg38_50kbwindows_sort.bed <(zcat $SP'_200sp_autosomes_sort.bed.gz' | \
awk -F'[\t@]' '{print $1"\t"$2"\t"$3"\t"$7"\t"$11}') | gzip > $SP'_200sp_autosomes_human50kb_numsites.bed.gz'

#mean phylop of window
bedmap --echo --wmean --delim "\t" ../hg38_50kbwindows_sort.bed <(zcat $SP'_200sp_autosomes_sort.bed.gz' | \
awk -F'[\t]' '{print $1"\t"$2"\t"$3"\t"$5"\t"$6}') | gzip > $SP'_200sp_autosomes_human50kb_wmeanphylop.bed.gz'

#max phylop of window
bedmap --echo --max --delim "\t" ../hg38_50kbwindows_sort.bed <(zcat $SP'_200sp_autosomes_sort.bed.gz' | \
awk -F'[\t]' '{print $1"\t"$2"\t"$3"\t"$6"\t"$5}') | gzip > $SP'_200sp_autosomes_human50kb_maxphylop.bed.gz'

#mean phylop of substitutions in window
sort-bed <(zcat "/home/centos/USS/aryn/Zoonomia/snp_phylop_200/"$SP"2HomSap_snps_phylop_200sp_autosomes_filtered.bed.gz") \
| gzip > "/home/centos/USS/aryn/Zoonomia/snp_phylop_200/"$SP"2HomSap_snps_phylop_200sp_autosomes_filtered_sort.bed.gz"

bedmap --echo --mean --count --delim "\t" ../hg38_50kbwindows_sort.bed \
<(zcat "/home/centos/USS/aryn/Zoonomia/snp_phylop_200/"$SP"2HomSap_snps_phylop_200sp_autosomes_filtered_sort.bed.gz" \
| awk -F'[\t]' '{print $1"\t"$2"\t"$3"\t"$5"\t"$6}') | gzip > $SP'_200sp_autosomes_human50kb_snpPhyloP.bed.gz'

paste <(zcat $SP'_200sp_autosomes_human50kb_wmeanROH.bed.gz') <(zcat $SP'_200sp_autosomes_human50kb_wmeanphylop.bed.gz' | cut -f 4) \
<(zcat $SP'_200sp_autosomes_human50kb_wmeanhet.bed.gz' | cut -f 4) \
<(zcat $SP'_200sp_autosomes_human50kb_maxphylop.bed.gz' | cut -f 4) \
<(zcat $SP'_200sp_autosomes_human50kb_snpPhyloP.bed.gz' | cut -f 4,5) \
<(zcat $SP'_200sp_autosomes_human50kb_numsites.bed.gz' | cut -f 4) \
| gzip > $SP'_200sp_autosomes_filtered_human50kb.bed.gz'

rm $SP'_200sp_autosomes_human50kb_wmeanROH.bed.gz' $SP'_200sp_autosomes_human50kb_wmeanphylop.bed.gz' $SP'_200sp_autosomes_human50kb_wmeanhet.bed.gz' \
$SP'_200sp_autosomes_human50kb_maxphylop.bed.gz' $SP'_200sp_autosomes_human50kb_numsites.bed.gz' $SP'_200sp_autosomes_human50kb_snpPhyloP.bed.gz'


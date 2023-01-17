#!/bin/bash
SP=$1 #focal species

#call substitutions, insertions and deletions in focal species' genome relative to most recent reconstructed ancestral node
halBranchMutations --refFile $SP'_ins.bed' --parentFile $SP'_del.bed' \
--snpFile $SP'_snps.txt' \
~/USS/aryn/Zoonomia/HAL/241-mammalian-2020v2.hal \
$SP \
>& $SP"_mut.nohup"

gzip $SP'_snps.txt' &

#liftover coordinates of substitutions in focal genome to human genome coordinates
halLiftover ~/USS/aryn/Zoonomia/HAL/241-mammalian-2020v2.hal --bedType 3 $SP --noDupes \
<(zcat '~/USS/aryn/Zoonomia/snps/'$SP'_snps.txt' | grep -v "#" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$1"@"$3}') Homo_sapiens \
$SP"2HomSap_snps.bed"


sort-bed $SP"2HomSap_snps.bed" > $SP"2HomSap_snps_sort.bed"

#assign phyloP scores (negative phylop scores=0) to substitutions at human positions that have at least 200 mammals represented in the alignment (autosomes only)
for CHR in `seq 1 22`; do

bedmap --echo --mean --count --delim "\t" --skip-unmapped $SP"2HomSap_snps_sort.bed" \
"~/USS/aryn/Zoonomia/phylop/human_onlychr_v2_mdong_10Mb_241MAMMALS_MTDF_species_count/chr"$CHR"_positive_200sp.bed" \
> $SP"2HomSap_snps_phylop_200sp_chr"$CHR".bed"
done

cat $SP"2HomSap_snps_phylop_200sp_chr"*".bed" > $SP"2HomSap_snps_phylop_200sp_autosomes.bed"
rm $SP"2HomSap_snps_phylop_200sp_chr"*".bed"

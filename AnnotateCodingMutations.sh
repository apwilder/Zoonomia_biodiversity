#!/bin/bash

SP=$1
GENOME=$2
GTF=$3
VCF=$4 #vcf of heterozygous sites in a genome, with phylop assigned by halLiftover
SNPEFFDIR=$5


##############build SnpEff database#############
mkdir $SNPEFFDIR'/data/'$SP
ln -sfn '~/USS/aryn/Zoonomia/gtfs/'$GTF $SNPEFFDIR'/data/'$SP'/genes.gtf.gz'
ln -sfn '~/USS/aryn/Zoonomia/genomes/'$GENOME $SNPEFFDIR'/data/'$SP'/sequences.fa.gz'

java -Xmx12g -jar $SNPEFFDIR'/snpEff.jar' \
build -gtf22 -v $SP

java -jar $SNPEFFDIR'/snpEff.jar' dump $SP | head -n 35 \
> $SP'_DB.stats'

##############Annotate vcfs of heterozygous sites#############
cd ~/USS/aryn/Zoonomia/human_phyloP_liftover #directory of vcfs 

java -Xmx12g -jar snpEff.jar $SP $VCF -o vcf \
-csvStats $SP"_variants_phylop_snpeff.stats" > $SP"_variants_phylop_snpeff.vcf"

gzip $SP"_variants_phylop_snpeff.vcf"

java -Xmx12g -jar SnpSift.jar \
filter "(ANN[*].IMPACT != 'MODIFIER') & ( GEN[*].GT = '0/1' )" $SP"_variants_snpeff.vcf.gz" |
~/USS/aryn/Zoonomia/snpeff/snpEff/scripts/vcfEffOnePerLine.pl | \
java -Xmx12g -jar ~/USS/aryn/Zoonomia/snpeff/snpEff/SnpSift.jar \
extractFields - CHROM POS "GEN[*].GT" "GEN[*].GQ" "GEN[*].DP" \
"ANN[*].EFFECT" "EFF[*].IMPACT" "EFF[*].GENE" "LOF[*].GENE" "LOF[*].NUMTR" "LOF[*].PERC" | gzip \
> $SP"_variants_snpeff_GQ_eff.txt.gz"

#############Annotate substitutions#############
cd /home/centos/USS/aryn/Zoonomia/snpeff

#create vcfs from lists of substitutions, adding basic vcf header
cp vcfheader $SP"_mutations.vcf"
grep -v "#" "~/USS/aryn/Zoonomia/snps/"$SP"_snps.txt" | awk \
'{print $1"\t"$3"\t.\t"substr($4,4,1)"\t"substr($4,3,1)"\t.\t.\t.\tGT\t0/1"}' \
| gzip >> "/home/centos/USS/aryn/Zoonomia/snpeff/"$SP"_mutations.vcf.gz"

#create list of genes to filter out (NEED TO FIGURE OUT HOW I GENERATED HIGHQUAL LIST)
zcat "~/USS/aryn/Zoonomia/gtfs2/"$SP"_highQual.gtf.gz" | \
awk '$3=="gene" {print $1"\t"($4 > 5000 ? $4 - 5000 : 0)"\t"$5+5000}' > \
"/home/centos/USS/aryn/Zoonomia/snpeff/"$SP"_coding.bed"

bedtools intersect -a "/home/centos/USS/aryn/Zoonomia/snpeff/"$SP"_coding.bed" -b \
<(zcat "~/USS/aryn/Zoonomia/filter/"$SP"_mask_a0.1_sort.bed.gz" | grep FALSE | cut -f1-3) > \
"/home/centos/USS/aryn/Zoonomia/snpeff/"$SP"_coding_filt.bed"
bedtools sort -i "/home/centos/USS/aryn/Zoonomia/snpeff/"$SP"_coding_filt.bed" > \
"/home/centos/USS/aryn/Zoonomia/snpeff/"$SP"_coding_filts.bed"
rm "/home/centos/USS/aryn/Zoonomia/snpeff/"$SP"_coding_filt.bed" "/home/centos/USS/aryn/Zoonomia/snpeff/"$SP"_coding.bed"


#annotate substitutions with predicted impact on protein coding genes
java -Djava.io.tmpdir=`pwd`/tmp -Xmx12g -jar ~/USS/aryn/Zoonomia/snpeff/snpEff/snpEff.jar \
$SP "/home/centos/USS/aryn/Zoonomia/snpeff/"$SP"_mutations.vcf.gz" -o vcf \
-fi "/home/centos/USS/aryn/Zoonomia/snpeff/"$SP"_coding_filts.bed" \
-csvStats "/home/centos/USS/aryn/Zoonomia/snpeff/"$SP"_mutations_filt_coding_snpeff.stats" \
| gzip > "/home/centos/USS/aryn/Zoonomia/snpeff/"$SP"_mutations_filt_coding_snpeff.vcf.gz"

###############Combine substitution annotations with phylop info###############
cd ~/USS/aryn/Zoonomia/snp_phylop_200

#create bed from vcf, reorder columns and sort
zcat '/home/centos/USS/aryn/Zoonomia/snpeff/'$SP'_mutations_filt_coding_snpeff.vcf.gz' \
| awk '{print $1"\t"$2-1"\t"$2"\t.\t"$8}' | sort-bed - > \
$SP'_mutations_filt_coding_snpeff_sorted.bed'
awk -F'[@\t]' '{print $5"\t"$6-1"\t"$6"\t"$1"@"$3"\t"$7}' $SP'2HomSap_snps_phylop_200sp_autosomes.bed' | \
sort-bed - > $SP'2HomSap_snps_phylop_200sp_autosomes_sorted.bed'

#combine protein-coding annotations with phylop scores
bedmap --echo --mean --echo-map-id --delim "\t" --skip-unmapped $SP'_mutations_snpeff_coding_sorted.bed' \
$SP'2HomSap_snps_phylop_200sp_autosomes_sorted.bed' \
> $SP'2HomSap_snps_phylop_200sp_coding.bed'

rm $SP'_mutations_snpeff_coding_sorted.bed'
rm $SP'2HomSap_snps_phylop_200sp_autosomes_sorted.bed'

#filter hypervariable windows
sort-bed $SP'2HomSap_snps_phylop_200sp_coding.bed' > $SP'coding_mapfile'

#filter with bedmaps:
bedmap --echo --echo-ref-name --delim "\t" --skip-unmapped \
$SP'coding_mapfile' <(grep FALSE '~/USS/aryn/Zoonomia/filter/'$SP'_mask_a0.1_sort.bed' | cut -f 1-3) \
> $SP'2HomSap_snps_phylop_200sp_coding_filtered.bed'

rm $SP'coding_mapfile'

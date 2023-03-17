#!/bin/bash

cd /home/centos/USS/aryn/Zoonomia/busco/busco_downloads/lineages/mammalia_odb10

for SITE in `awk -F'[\t]' '{print $1}' links_to_ODB10.txt`; do
wget --no-check-certificate -q 'https://www.orthodb.org/tab?query='$SITE'&species=9606_0&singlecopy=0.9' -O 'gene_desc/'$SITE'.out'
done

awk '
    FNR==1 && NR!=1 { while (/^<header>/) getline; }
    1 {print}
' gene_desc/*.out > busco_genes.txt

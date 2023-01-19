library(data.table)
library(ggplot2)
library(Polychrome)
library(RColorBrewer)
library(alphaOutlier)
library(ape)
library(phylolm)
library(MuMIn)
library(reshape)

setwd('/home/centos/USS/aryn/Zoonomia/snpeff_filt')
fusil=read.table('../genes_with_FUSIL_bins_20210922.txt',sep='\t',header=T)
buscogenes=read.table('/home/centos/USS/aryn/Zoonomia/datatables/buscogenes.txt')
buscosummary=read.table('/home/centos/USS/aryn/Zoonomia/datatables/buscoSummary.txt',header=T)
busco=as.character(buscosummary$gene)#[which(buscosummary$nSpecies>200)] (makes no difference)

files=list.files(path = ".", pattern = '_mutations_filt_coding_snpeff.stats.genes.txt')
sps=gsub('_mutations_filt_coding_snpeff.stats.genes.txt','',files)
coding=matrix(nrow=length(sps),ncol=33)
for (i in 1:length(sps)){
sp=sps[i]
genes=read.table(files[i])
headers=scan(files[i],nlines=2,what="character",sep='\t')
names(genes)=sub("#","",headers[2:length(headers)])
names(genes)=gsub("variants_impact_","",names(genes))
names(genes)=gsub("variants_effect_","",names(genes))
if(sum(c('GeneName','HIGH','LOW','MODERATE','missense_variant','synonymous_variant') %in% names(genes)==F)>0){
  genes$missing=0
  names(genes)[ncol(genes)]=c('GeneName','HIGH','LOW','MODERATE','missense_variant','synonymous_variant')[which(c('GeneName','HIGH','LOW','MODERATE','missense_variant','synonymous_variant') %in% names(genes)==F)]
}
genes=genes[,c('GeneName','HIGH','LOW','MODERATE','missense_variant','synonymous_variant')]
genes=genes[which(genes$GeneName %in% busco),]
genes$coding_vars=rowSums(genes[,which(names(genes) %in% c('missense_variant','start_lost','stop_gained','synonymous_variant'))])
genes=merge(genes,fusil[,c('human_symbol','Viability.Phenotype.HOMs.HEMIs','confidence','Achilles_gene_effect_mean','FUSIL_bin')],by.x='GeneName',by.y='human_symbol',all.x=T)

vars=c(sp,sum(genes$missense_variant),sum(genes$synonymous_variant),sum(genes$HIGH),sum(genes$coding_vars))
vars=c(vars,table(genes$FUSIL_bin))
vars=c(vars,aggregate(genes$synonymous_variant,by=list(genes$FUSIL_bin),sum,na.rm=T)$x)
vars=c(vars,aggregate(genes$missense_variant,by=list(genes$FUSIL_bin),sum,na.rm=T)$x)
vars=c(vars,aggregate(genes$HIGH,by=list(genes$FUSIL_bin),sum,na.rm=T)$x)
vars=c(vars,table(genes$Viability.Phenotype.HOMs.HEMIs))
vars=c(vars,aggregate(genes$synonymous_variant,by=list(genes$Viability.Phenotype.HOMs.HEMIs),sum,na.rm=T)$x)
vars=c(vars,aggregate(genes$missense_variant,by=list(genes$Viability.Phenotype.HOMs.HEMIs),sum,na.rm=T)$x)
vars=c(vars,aggregate(genes$HIGH,by=list(genes$Viability.Phenotype.HOMs.HEMIs),sum,na.rm=T)$x)
coding[i,1:33]=vars
}

coding=as.data.frame(coding)
coding[,2:ncol(coding)]=apply(coding[,2:ncol(coding)],2,function(x) as.numeric(as.character(x)))
names(coding)=c('Sp','missense','synonymous','high','coding_vars',paste('fusil_genes_',c('CL','DL','CV','V'),sep=""),paste('synonymous_fusil_',c('CL','DL','CV','V'),sep=""),paste('missense_fusil_',c('CL','DL','CV','V'),sep=""),paste('high_fusil_',c('CL','DL','CV','V'),sep=""),paste('impc_genes_',c('L','SV','V'),sep=""),paste('synonymous_impc_',c('L','SV','V'),sep=""),paste('missense_impc_',c('L','SV','V'),sep=""),paste('high_impc_',c('L','SV','V'),sep=""))
coding$Sp=as.character(coding$Sp)
write.table(coding,'~/USS/aryn/Zoonomia/datatables/Coding_mutations_IMPC_FUSIL_busco.txt',sep='\t',row.names=F,quote=F)


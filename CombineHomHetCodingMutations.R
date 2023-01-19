library(data.table)
library(ggplot2)
library(Polychrome)
library(RColorBrewer)
library(alphaOutlier)
library(ape)
library(phylolm)
library(MuMIn)

setwd('~/USS/aryn/Zoonomia/datatables')

######Read and filter data##########
my_tree <- read.tree('../Zoonomia_ChrX_lessGC40_241species_30Consensus.tree')

adjstats=read.table('Mutation-adjustedLoad_PSMC_het_241sp.txt',header=T,sep='\t')
adjstats$IUCN=factor(adjstats$IUCN,levels=c('LC','NT','VU','EN','CR','DD'))
stats=read.table('Load_PSMC_het_241sp.txt',header=T,sep='\t')
row.names(stats)=stats$Sp
row.names(adjstats)=adjstats$Sp
row.names(stats)=stats$Sp
stats$IUCN=factor(stats$IUCN,levels=c('LC','NT','VU','EN','CR','DD'))
lowsp=stats$Sp[which(stats$coding_vars<10000)]


######Summarize impacts of heterozygous variants on protein coding genes##########
coding_het=read.table('~/USS/aryn/Zoonomia/datatables/Coding_heterozygous_filtered_lof_noindels_IMPC_FUSIL_GQ80_busco.txt',sep='\t',header=T)#exclude loeuf coding_het=coding_het[,c(1:8,34:68)]
coding_het=merge(coding_het,adjstats[,c('Sp','Order','IUCN','Source','psmc_min','harm_mean_wt2','froh','NeNc','meanROH','froh_1MB','gw_het_mean','outbred_het_mode')],by='Sp',all.x=T)
coding_het$IUCN=factor(coding_het$IUCN,levels=c('LC','NT','VU','EN','CR','DD'))
row.names(coding_het)=as.character(coding_het$Sp)

#filter low quality species from heterozygous variatn analysis
coding_het=coding_het[which(coding_het$meanDP>=20),]#use only species with >=20x mean depth
coding_het=coding_het[which(coding_het$meanGQ>=80),]#use only species with meanGQ>=80
coding_het=coding_het[which(coding_het$coding_vars>=1000),]#exclude 1 sp with <1000 coding variants


##########counts of heterozygous var types############
lof_impc_L=which(names(coding_het)=='lof_impc_L')
synonymous_impc_L=which(names(coding_het)=='synonymous_impc_L')
missense_impc_L=which(names(coding_het)=='missense_impc_L')
lof_impc_SV=which(names(coding_het)=='lof_impc_SV')
missense_impc_SV=which(names(coding_het)=='missense_impc_SV')
synonymous_impc_SV=which(names(coding_het)=='synonymous_impc_SV')
lof_impc_V=which(names(coding_het)=='lof_impc_V')
missense_impc_V=which(names(coding_het)=='missense_impc_V')
synonymous_impc_V=which(names(coding_het)=='synonymous_impc_V')
lof_fusil_DL=which(names(coding_het)=='lof_fusil_DL')
lof_fusil_CL=which(names(coding_het)=='lof_fusil_CL')
missense_fusil_DL=which(names(coding_het)=='missense_fusil_DL')
missense_fusil_CL=which(names(coding_het)=='missense_fusil_CL')
synonymous_fusil_DL=which(names(coding_het)=='synonymous_fusil_DL')
synonymous_fusil_CL=which(names(coding_het)=='synonymous_fusil_CL')
lof=which(names(coding_het)=='lof')
synonymous=which(names(coding_het)=='synonymous')
missense=which(names(coding_het)=='missense')
high=which(names(coding_het)=='high')

#calculate proportions of variant impact types
val0=0
coding_het$highsyn_impc_L=apply(coding_het,1,function(x) ifelse(as.numeric(x[lof_impc_L])==0,val0,as.numeric(x[lof_impc_L])/(as.numeric(x[synonymous_impc_L])+as.numeric(x[missense_impc_L])+as.numeric(x[lof_impc_L]))))
coding_het$highsyn_impc_V=apply(coding_het,1,function(x) ifelse(as.numeric(x[lof_impc_V])==0,val0,as.numeric(x[lof_impc_V])/(as.numeric(x[synonymous_impc_V])+as.numeric(x[missense_impc_V])+as.numeric(x[lof_impc_V]))))
coding_het$highsyn_impc_SV=apply(coding_het,1,function(x) ifelse(as.numeric(x[lof_impc_SV])==0,val0,as.numeric(x[lof_impc_SV])/(as.numeric(x[synonymous_impc_SV])+as.numeric(x[missense_impc_SV])+as.numeric(x[lof_impc_SV]))))
coding_het$highsyn_fusil_DL=apply(coding_het,1,function(x) ifelse(as.numeric(x[lof_fusil_DL])==0,val0,as.numeric(x[lof_fusil_DL])/(as.numeric(x[synonymous_fusil_DL])+as.numeric(x[missense_fusil_DL])+as.numeric(x[lof_fusil_DL]))))
coding_het$highsyn_fusil_CL=apply(coding_het,1,function(x) ifelse(as.numeric(x[lof_fusil_CL])==0,val0,as.numeric(x[lof_fusil_CL])/(as.numeric(x[synonymous_fusil_CL])+as.numeric(x[missense_fusil_CL])+as.numeric(x[lof_fusil_CL]))))
coding_het$misssyn_impc_L=(coding_het$missense_impc_L)/(coding_het$synonymous_impc_L+coding_het$missense_impc_L+coding_het$high_impc_L)
coding_het$misssyn_impc_SV=(coding_het$missense_impc_SV)/(coding_het$synonymous_impc_SV+coding_het$missense_impc_SV+coding_het$high_impc_SV)
coding_het$misssyn_impc_V=(coding_het$missense_impc_V)/(coding_het$synonymous_impc_V+coding_het$missense_impc_V+coding_het$high_impc_V)
coding_het$misssyn_fusil_DL=(coding_het$missense_fusil_DL)/(coding_het$synonymous_fusil_DL+coding_het$missense_fusil_DL+coding_het$high_fusil_DL)
coding_het$misssyn_fusil_CL=(coding_het$missense_fusil_CL)/(coding_het$synonymous_fusil_CL+coding_het$missense_fusil_CL+coding_het$high_fusil_CL)
coding_het$misssyn_fusil_CV=(coding_het$missense_fusil_CV)/(coding_het$synonymous_fusil_CV+coding_het$missense_fusil_CV+coding_het$high_fusil_CV)
coding_het$misssyn_fusil_V=(coding_het$missense_fusil_V)/(coding_het$synonymous_fusil_V+coding_het$missense_fusil_V+coding_het$high_fusil_V)
coding_het$highsyn_fusil_CV=(coding_het$lof_fusil_CV)/(coding_het$synonymous_fusil_CV+coding_het$missense_fusil_CV+coding_het$high_fusil_CV)
coding_het$highsyn_fusil_V=(coding_het$lof_fusil_V)/(coding_het$synonymous_fusil_V+1)
coding_het$highsyn=apply(coding_het,1,function(x) ifelse(as.numeric(x[lof])==0,val0,as.numeric(x[lof])/(as.numeric(x[synonymous])+as.numeric(x[missense])+as.numeric(x[high]))))
coding_het$misssyn=(coding_het$missense)/(coding_het$synonymous+coding_het$missense+coding_het$high)

#remove estimates for species with <15 genes in a category
coding_het$highsyn_fusil_CV[which(coding_het$fusil_genes_CV<15)]=NA
coding_het$misssyn_fusil_CV[which(coding_het$fusil_genes_CV<15)]=NA
coding_het$highsyn_fusil_CL[which(coding_het$fusil_genes_CL<15)]=NA
coding_het$misssyn_fusil_CL[which(coding_het$fusil_genes_CL<15)]=NA
coding_het$highsyn_fusil_DL[which(coding_het$fusil_genes_DL<15)]=NA
coding_het$misssyn_fusil_DL[which(coding_het$fusil_genes_DL<15)]=NA
coding_het$highsyn_fusil_V[which(coding_het$fusil_genes_V<15)]=NA
coding_het$misssyn_fusil_V[which(coding_het$fusil_genes_V<15)]=NA
coding_het$highsyn_impc_V[which(coding_het$impc_genes_V<15)]=NA
coding_het$misssyn_impc_V[which(coding_het$impc_genes_V<15)]=NA
coding_het$highsyn_impc_SV[which(coding_het$impc_genes_SV<15)]=NA
coding_het$misssyn_impc_SV[which(coding_het$impc_genes_SV<15)]=NA
coding_het$highsyn_impc_L[which(coding_het$impc_genes_L<15)]=NA
coding_het$misssyn_impc_L[which(coding_het$impc_genes_L<15)]=NA

coding_het$threatened=ifelse(coding_het$IUCN=='LC','non-threatened','threatened')
coding_het$threatened=ifelse(coding_het$IUCN=='DD',NA,coding$threatened)

######Summarize impacts of homozygous substitutions on protein coding genes##########
coding=read.table('~/USS/aryn/Zoonomia/datatables/Coding_mutations_IMPC_FUSIL_busco.txt',sep='\t',header=T)


coding=coding[!is.na(coding$Sp),]
coding=merge(coding,adjstats[,c('Sp','IUCN','Order','psmc_min','harm_mean_wt2','froh','NeNc','meanROH','froh_1MB','mutations','distance')],by='Sp',all.x=T)
row.names(coding)=as.character(coding$Sp)

coding=coding[which(coding$coding_vars>10000),]#include only species with >10000 coding substitutions

coding$highsyn_impc_L=(coding$high_impc_L+1)/(coding$synonymous_impc_L+coding$missense_impc_L+coding$high_impc_L)
coding$highsyn_impc_V=(coding$high_impc_V+1)/(coding$synonymous_impc_V+coding$missense_impc_V+coding$high_impc_V)
coding$highsyn_impc_SV=(coding$high_impc_SV+1)/(coding$synonymous_impc_SV+coding$missense_impc_SV+coding$high_impc_SV)
coding$highsyn_fusil_DL=(coding$high_fusil_DL+1)/(coding$high_fusil_DL+coding$missense_fusil_DL+coding$synonymous_fusil_DL)
coding$highsyn_fusil_CL=(coding$high_fusil_CL+1)/(coding$high_fusil_CL+coding$missense_fusil_CL+coding$synonymous_fusil_CL)
coding$highsyn_fusil_CV=(coding$high_fusil_CV+1)/(coding$high_fusil_CV+coding$missense_fusil_CV+coding$synonymous_fusil_CV)
coding$highsyn_fusil_V=(coding$high_fusil_V+1)/(coding$high_fusil_V+coding$missense_fusil_V+coding$synonymous_fusil_V)
coding$misssyn_impc_L=(coding$missense_impc_L+1)/(coding$synonymous_impc_L+coding$missense_impc_L+coding$high_impc_L)
coding$misssyn_impc_SV=(coding$missense_impc_SV+1)/(coding$synonymous_impc_SV+coding$missense_impc_SV+coding$high_impc_SV)
coding$misssyn_impc_V=(coding$missense_impc_V+1)/(coding$synonymous_impc_V+coding$missense_impc_V+coding$high_impc_V)
coding$misssyn_fusil_DL=(coding$missense_fusil_DL+1)/(coding$high_fusil_DL+coding$missense_fusil_DL+coding$synonymous_fusil_DL)
coding$misssyn_fusil_CL=(coding$missense_fusil_CL+1)/(coding$high_fusil_CL+coding$missense_fusil_CL+coding$synonymous_fusil_CL)
coding$misssyn_fusil_CV=(coding$missense_fusil_CV+1)/(coding$high_fusil_CV+coding$missense_fusil_CV+coding$synonymous_fusil_CV)
coding$misssyn_fusil_V=(coding$missense_fusil_V+1)/(coding$synonymous_impc_V+coding$missense_impc_V+coding$high_impc_V)
coding$miss=(coding$missense+1)/(coding$synonymous+coding$missense+coding$high)

#remove categories with fewer than 15 genes
coding$highsyn_impc_V[which(coding$impc_genes_V<15)]=NA
coding$misssyn_impc_V[which(coding$impc_genes_V<15)]=NA
coding$highsyn_impc_SV[which(coding$impc_genes_SV<15)]=NA
coding$misssyn_impc_SV[which(coding$impc_genes_SV<15)]=NA
coding$highsyn_impc_L[which(coding$impc_genes_L<15)]=NA
coding$misssyn_impc_L[which(coding$impc_genes_L<15)]=NA

#adjust load estimates for number of mutations between ancestor and focal species
mutadj=function(x){
  llr=phylolm(log(x)~log(mutations),data=coding,phy=my_tree)
  return(residuals.phylolm(llr))
}
llr.res=apply(coding[,64:78],2,mutadj)
coding=merge(coding[,1:63],llr.res,by.x='Sp',by.y=0)
row.names(coding)=coding$Sp

coding$threatened=ifelse(coding$IUCN=='LC','non-threatened','threatened')
coding$threatened=ifelse(coding$IUCN=='DD',NA,coding$threatened)

write.table(coding,'~/USS/aryn/Zoonomia/datatables/Coding_mutations_IMPC_FUSIL_busco_adjusted_ppn.txt',sep='\t')

####Combine homozygous and heterozygous for plotting###
##Homozygous##
loadplot=reshape2::melt(coding[,c('Sp','threatened','Order','IUCN','harm_mean_wt2','psmc_min','froh','NeNc','highsyn_impc_L','misssyn_impc_L','highsyn_impc_SV','misssyn_impc_SV','highsyn_impc_V','misssyn_impc_V','miss')],id.vars=c('Sp','threatened','Order','IUCN','harm_mean_wt2','psmc_min','froh','NeNc'))
loadplot$threatenedbin=loadplot$threatened
loadplot$threatened=NA
loadplot$threatened[which(loadplot$threatenedbin==1)]='threatened'
loadplot$threatened[which(loadplot$threatenedbin==0)]='non-threatened'

loadplot$variant=NA
loadplot$variant[grep("high",loadplot$variable)]='LoF'
loadplot$variant[grep("miss",loadplot$variable)]='Missense'
loadplot$impact='All'
loadplot$impact[grep("_L",loadplot$variable)]='Lethal'
loadplot$impact[grep("_SV",loadplot$variable)]='Subviable'
loadplot$impact[grep("_V",loadplot$variable)]='Viable'

##Heterozygous##
loadplot_het= reshape2::melt(coding_het[,c('Sp','threatened','Order','IUCN','harm_mean_wt2','psmc_min','froh','NeNc','gw_het_mean','outbred_het_mode','highsyn_impc_L','misssyn_impc_L','highsyn_impc_SV','misssyn_impc_SV','highsyn_impc_V','misssyn_impc_V','highsyn','misssyn',"impc_genes_L",'impc_genes_SV','impc_genes_V')],id.vars=c('Sp','threatened','Order','IUCN','harm_mean_wt2','psmc_min','froh','NeNc','gw_het_mean','outbred_het_mode',"impc_genes_L",'impc_genes_SV','impc_genes_V'))
loadplot_het$threatened[which(loadplot_het$threatened==1)]='threatened'
loadplot_het$threatened[which(loadplot_het$threatened==0)]='non-threatened'

loadplot_het$variant=NA
loadplot_het$variant[grep("high",loadplot_het$variable)]='LoF'
loadplot_het$variant[grep("miss",loadplot_het$variable)]='Missense'
loadplot_het$impact='All'
loadplot_het$impact[grep("_L",loadplot_het$variable)]='Lethal'
loadplot_het$impact[grep("_SV",loadplot_het$variable)]='Subviable'
loadplot_het$impact[grep("_V",loadplot_het$variable)]='Viable'


##Combined
lp=rbind(loadplot[,c(1:10,12:13)],loadplot_het[,c(1:8,14:17)])
lp$zyg=NA
lp$zyg[1:nrow(loadplot)]='Homozygous'
lp$zyg[(nrow(loadplot)+1):nrow(lp)]='Heterozygous'
lp=lp[which((lp$variant=='LoF' & lp$zyg=='Homozygous')==F),]
#write.table(lp,'IMPC_CodingVars_Ne_4plot.txt',sep='\t',row.names=F)# without mutation-adjustment
write.table(lp,'IMPC_CodingVars_adjusted_Ne_4plot.txt',sep='\t',row.names=F)


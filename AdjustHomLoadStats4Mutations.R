library(ape)
library(phylolm)
library(MuMIn)
library(Polychrome)
library(RColorBrewer)
library(ggplot2)
library(reshape)

setwd('~/USS/aryn/Zoonomia/datatables')

##read in data#####
my_tree <- read.tree('~/USS/aryn/Zoonomia/Zoonomia_ChrX_lessGC40_241species_30Consensus.tree')
dists=my_tree$edge.length[sapply(1:length(my_tree$tip.label),function(x,y) which (y==x),y=my_tree$edge[,2])]
dists=as.data.frame(dists)
dists$Sp=my_tree$tip.label

load=read.table('LoadMatrix_241sp_pp200_filtered.txt',header=T,sep='\t')
meta=read.table('Metadata_14Sept2021.txt',header=T,na.string="",sep='\t')
load=merge(load,meta[,c(1:9,12)],by='Sp',all.x=T)
load$IUCN=factor(load$IUCN,levels=c('LC','NT','VU','EN','CR','DD'))

psmc=read.csv('zoonomia_psmc_metrics_12sept2021.csv',header=T)
psmc$genus_species=as.character(psmc$genus_species)
psmc=merge(load,psmc[,c(3,17:18)],by.y='genus_species',by.x='Sp',all=T)
psmc$Sp=as.character(psmc$Sp)

het=read.table('zoonomia_240spp_het_metrics_10sept2021.csv',header=T,sep=',')
names(het)[1]='Sp'
het$Sp=as.character(het$Sp)
froh=read.table('froh_human_209sp_manualfix.txt',header=T)
het=merge(het,froh,by='Sp',all=T)
het=het[,c(1:3,5:6)]
names(het)[4]='froh'
het$gw_het_mean[which(is.na(het$gw_het_mean) & !is.na(het$het))]=het$het[which(is.na(het$gw_het_mean) & !is.na(het$het))]
het=het[,1:4]
het=het[!duplicated(het$Sp),]
#add 11 species that were manually adjusted
outbredgaps=read.table('Outbred_het_mode_manualfix_11sp.txt',header=T)
outbredgaps=outbredgaps[!is.na(outbredgaps$V4),]
outbredgaps$Sp=as.character(outbredgaps$Sp)
het=merge(het,outbredgaps,by='Sp',all.x=T)
het$outbred_het_mode[!is.na(het$V4)]=het$V4[!is.na(het$V4)]
het=het[,1:4]

psmc=merge(psmc,het[!duplicated(het$Sp),],by='Sp',all.x=T)
psmc[which(psmc$Sp=="Oryctolagus_cuniculus"),15:19]=NA #remove european rabbit from analysis because short-reads were from pooled individuals

snpeff=read.table('Coding_Vars_and_phyloP_240sp.txt',header=T,sep='\t')#update
psmc=merge(psmc,snpeff[,c(1,3:ncol(snpeff))],by='Sp',all=T)
row.names(psmc)=as.character(psmc$Sp)
psmc$Source=factor(psmc$Source,levels=c("1. Zoonomia","2. Existing assembly"))
psmc$wild_status_ref=factor(psmc$wild_status_ref,levels=c('wild','captive','domestic'))
psmc$wild_status_reseq=factor(psmc$wild_status_reseq,levels=c('wild','captive','domestic'))
psmc$threatened=ifelse(psmc$IUCN %in% c('LC','DD'),0,1)
psmc$threatened=as.factor(psmc$threatened)

psmc[which(psmc$Sp=="Capra_hircus"),c(20:37)]=NA #remove because it had too few annotated genes
psmc=merge(psmc,dists,by='Sp')
row.names(psmc)=psmc$Sp

######log-log regression correction of homozygous load stats for number of mutations between focal species and ancestral sequence########
mutadj=function(x){
llr=phylolm(log(x)~log(mutations),data=psmc,phy=my_tree)
return(residuals.phylolm(llr))
}
llr.res=apply(psmc[,c(2:3,5)],2,mutadj)# apply correction to "mean_phylop","ppn_conserved","phylop_kurtosis"
llr.res=as.data.frame(llr.res,makes.names = T)
llr.res2=apply(psmc[,c(20:39)],2,mutadj) #apply correction to
llr.res2=as.data.frame(llr.res2,makes.names = T)
llr.res=merge(llr.res,llr.res2,by=0,all=T)
llr.res=merge(llr.res,psmc[,which(names(psmc) %in% names(llr.res)==F)],by.y=0,by.x='Row.names',all.x=T)
row.names(llr.res)=llr.res$Row.names
llr.res=llr.res[,2:ncol(llr.res)]

llr.res=llr.res[,names(psmc)]
write.table(llr.res,'Mutation-adjustedLoad_PSMC_het_241sp.txt',sep='\t',quote=F,row.names=F)
write.table(psmc,'Load_PSMC_het_241sp.txt',sep='\t',quote=F,row.names=F)

library(stringr)
library(ggplot2)
library(Polychrome)
library(data.table)
library(moments)
library(ape)
library(phylolm)
library(MuMIn)


setwd('~/USS/aryn/Zoonomia/snp_phylop_200')

filelist=list.files(path = ".", pattern ='2HomSap_snps_phylop_200sp_autosomes_filtered.bed')#list files of all species' substitutions with phylop scores
filelist=filelist[which(file.size(filelist)!=0)]

#create matrix for load stats for each species
load=matrix(nrow=length(filelist),ncol=6)
#matrix for creating hisotgram of phylop distribution across snps
histmat=matrix(nrow=51,ncol=length(filelist))

#summarize phylop stats for each species and add to matrix
for (i in 1:length(filelist)){
  file=filelist[i]
  sp=gsub('2HomSap_snps_phylop_200sp_autosomes_filtered.bed','',file)
  print(sp)
  bedsp=fread(file,tmpdir=getwd())
  names(bedsp)=c('chr','pos1','pos2','alleles','genome_coords','phylop')
  bedsp=bedsp[!(duplicated(bedsp[,1:3]) | duplicated(bedsp[,1:3], fromLast = TRUE)), ]#remove sites that liftover to multiple loci
  bedsp=bedsp[!duplicated(bedsp$genome_coords), ]
  p=ggplot(bedsp[which(bedsp$phylop>0)], aes(x = phylop)) + geom_histogram(breaks=0:50*0.18) + scale_y_log10()+ggtitle(sp) #create histogram
  pg <- ggplot_build(p)
  pg=pg$data[[1]]
  histmat[1:length(pg[,2]),i]=pg[,2] #add phylop distribution to species matrix
  
  load[i,1]=sp #species
  load[i,2]=mean(bedsp$phylop) #mean phylop across substitutions
  load[i,3]=sum(bedsp$phylop>2.27)/length(bedsp$phylop) #proportion of substitutions that are constrained (phylop>2.27)
  load[i,4]=length(bedsp$phylop) #number of substititions with a phylop scores
  load[i,5]=kurtosis(bedsp$phylop) #kurtosis of phylop distribution
  load[i,6]=mean(bedsp$phylop[which(bedsp$phylop>0)]) #mean phylop score of sites with positive scores
}

load=as.data.frame(load)
names(load)=c('Sp','mean_phylop','ppn_conserved','mutations','phylop_kurtosis')
load[,1:5]=apply(load[,1:5],2,function(x) as.character(x))
load[,2:5]=apply(load[,2:5],2,function(x) as.numeric(x))
load$Sp=gsub('2HomSap_snps_phylop_200sp_autosomes.bed','',load$Sp)

#add metadata to matrix
meta=read.table('SpeciesGenomeMetadata.txt',sep='\t',header=T)
dmat=read.table('species.closest.distances',header=T)
names(dmat)=c('Sp1','Sp2','distance')

load=merge(load,meta[,c('name','Order','IUCN')],all.x=T,by.x='Sp',by.y='name')
load$IUCN=factor(load$IUCN,levels=c('LC','NT','VU','EN','CR','DD'))
load=merge(load,dmat,by.x='Sp',by.y='Sp1',all.x=T)
names(load)[8:9]=c('closest_sp','distance_to_closest')

write.table(load,'LoadMatrix_241sp_pp200_filtered.txt',sep='\t',quote=F,row.names=F)

write.table(histmat,'~/USS/aryn/Zoonomia/Phylop_histmat_w0s.txt',row.names=F,col.names=F,sep='\t',quote=F)

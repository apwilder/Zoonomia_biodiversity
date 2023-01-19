library(stringr)
library(data.table)
library(reshape)
library(circlize)
library(tidyr)
library(ape)
library(phylolm)
library(modeest)

#throw out ones that failed AB or other tests in first paper (Table S3)
setwd('/home/centos/USS/aryn/Zoonomia/roh/')

###############read in metadata###############

meta=read.table('~/USS/aryn/Zoonomia/datatables/Metadata_14Sept2021.txt',header=T,sep='\t')
meta$IUCN=factor(meta$IUCN,levels=c('LC','NT','VU','EN','CR','DD'))

iucn=meta$IUCN
names(iucn)=meta$name
iucn=factor(iucn,levels=c('LC','NT','VU','EN','CR','DD'))

##read in phylogenetic data#####
my_tree <- read.tree('~/USS/aryn/Zoonomia/Zoonomia_ChrX_lessGC40_241species_30Consensus.tree')
dists=my_tree$edge.length[sapply(1:length(my_tree$tip.label),function(x,y) which (y==x),y=my_tree$edge[,2])]
dists=as.data.frame(dists)
dists$Sp=my_tree$tip.label

###############make matrices of human 50kb windows###############
rohmat=NULL #roh scores
phylopmat=NULL #phylop scores for sites that map to window (should be roughly equal across species)
hetmat=NULL #mean het
maxphylopmat=NULL #max phylop
sitesmat=NULL #number of sites that map
snpphylopmat=NULL #mean phylop of snps in window
snpmat=NULL #number of snps in window

#summarize ROH, het, phylop, etc. within 50kb windows lifted to human genome
#to manually call failed species, ran ManualReplaceSOHfails_50KB.sh
for (filei in list.files(path = "human_200sp", pattern = '_200sp_autosomes_filtered_human50kb.bed.gz',full.names=T)){
  sp=gsub('_200sp_autosomes_filtered_human50kb.bed.gz','',basename(filei))
  roh=fread(filei,sep='\t',nThread=6,showProgress=T,tmpdir='/home/centos/USS/aryn/Zoonomia/roh')
  names(roh)=c('chr','start','stop','roh','mean_phylop','het','max_phylop','snp_phylop','snps','sites')
  names(roh)[4:10]=paste(sp,names(roh)[4:10],sep=".")
  roh[which(roh[,10]<5000 | roh[,10]>=50000 | is.na(roh[,10])),4:10]=NA
  rohmat=cbind(rohmat,roh[,4])
  phylopmat=cbind(phylopmat,roh[,5])
  hetmat=cbind(hetmat,roh[,6])
  maxphylopmat=cbind(maxphylopmat,roh[,7])
  snpphylopmat=cbind(snpphylopmat,roh[,8])
  snpmat=cbind(snpmat,roh[,9])
  sitesmat=cbind(sitesmat,roh[,10])
}
rohmat=cbind(roh[,1:3],rohmat)
rohmat1=rohmat[,..idx]
phylopmat=cbind(roh[,1:3],phylopmat)
hetmat=cbind(roh[,1:3],hetmat)
maxphylopmat=cbind(roh[,1:3],maxphylopmat)
sitesmat=cbind(roh[,1:3],sitesmat)
snpphylopmat=cbind(roh[,1:3],snpphylopmat)
snpmat=cbind(roh[,1:3],snpmat)

#throw out Oryctolagus_cuniculus (multiple individuals in read data)
hetmat=hetmat[,c(1:202,204:210)]
rohmat=rohmat[,c(1:202,204:210)]
sitesmat=sitesmat[,c(1:202,204:210)]

write.table(snpmat,'~/USS/aryn/Zoonomia/datatables/snpmat_human50KB_240sp.txt',row.names=F,col.names=T,sep='\t',quote=F)
write.table(snpphylopmat,'~/USS/aryn/Zoonomia/datatables/snpphylopmat_human50KB_240sp.txt',row.names=F,col.names=T,sep='\t',quote=F)
write.table(rohmat,'~/USS/aryn/Zoonomia/roh/rohmat_human200_207sp.txt',row.names=F,col.names=T,sep='\t',quote=F)
write.table(hetmat,'~/USS/aryn/Zoonomia/roh/hetmat_human200_207sp.txt',row.names=F,col.names=T,sep='\t',quote=F)
write.table(sitesmat,'~/USS/aryn/Zoonomia/roh/sitesmat_human200_207sp.txt',row.names=F,col.names=T,sep='\t',quote=F)




###############make matrices of coding substitutions###############
setwd('/home/centos/USS/aryn/Zoonomia/snp_phylop_200')

missconmat=NULL 
missppmat=NULL 
missmat=NULL
highconmat=NULL 
highppmat=NULL 
highmat=NULL
synconmat=NULL 
synppmat=NULL 
synmat=NULL

for (filei in list.files(path = ".", pattern = '_200sp_autosomes_human50kb_coding_phylop.bed.gz',full.names=F)){
  sp=gsub('_200sp_autosomes_human50kb_coding_phylop.bed.gz','',basename(filei))
  roh=fread(filei,sep='\t',nThread=6,showProgress=T,tmpdir='.')
  names(roh)=c('chr','start','stop','syn','syn_phylop','syn_consvd','miss','miss_phylop','miss_consvd','high','high_phylop','high_consvd')
  names(roh)[4:12]=paste(sp,names(roh)[4:12],sep=".")
  synmat=cbind(synmat,roh[,4])
  synppmat=cbind(synppmat,roh[,5])
  synconmat=cbind(synconmat,roh[,6])
  missmat=cbind(missmat,roh[,7])
  missppmat=cbind(missppmat,roh[,8])
  missconmat=cbind(missconmat,roh[,9])
  highmat=cbind(highmat,roh[,10])
  highppmat=cbind(highppmat,roh[,11])
  highconmat=cbind(highconmat,roh[,12])
}

synmat=cbind(roh[,1:3],synmat)
synppmat=cbind(roh[,1:3],synppmat)
synconmat=cbind(roh[,1:3],synconmat)
missmat=cbind(roh[,1:3],missmat)
missppmat=cbind(roh[,1:3],missppmat)
missconmat=cbind(roh[,1:3],missconmat)
highmat=cbind(roh[,1:3],highmat)
highppmat=cbind(roh[,1:3],highppmat)
highconmat=cbind(roh[,1:3],highconmat)

write.table(synmat,'~/USS/aryn/Zoonomia/datatables/synonymous_counts_human50kb_240sp.txt',row.names=F,col.names=T,sep='\t',quote=F)
write.table(synppmat,'~/USS/aryn/Zoonomia/datatables/synonymous_phylop_human50kb_240sp.txt',row.names=F,col.names=T,sep='\t',quote=F)
write.table(synconmat,'~/USS/aryn/Zoonomia/datatables/synonymous_conserved_human50kb_240sp.txt',row.names=F,col.names=T,sep='\t',quote=F)
write.table(missmat,'~/USS/aryn/Zoonomia/datatables/missense_counts_human50kb_240sp.txt',row.names=F,col.names=T,sep='\t',quote=F)
write.table(missppmat,'~/USS/aryn/Zoonomia/datatables/missense_phylop_human50kb_240sp.txt',row.names=F,col.names=T,sep='\t',quote=F)
write.table(missconmat,'~/USS/aryn/Zoonomia/datatables/missense_conserved_human50kb_240sp.txt',row.names=F,col.names=T,sep='\t',quote=F)
write.table(highmat,'~/USS/aryn/Zoonomia/datatables/high_counts_human50kb_240sp.txt',row.names=F,col.names=T,sep='\t',quote=F)
write.table(highppmat,'~/USS/aryn/Zoonomia/datatables/high_phylop_human50kb_240sp.txt',row.names=F,col.names=T,sep='\t',quote=F)
write.table(highconmat,'~/USS/aryn/Zoonomia/datatables/high_conserved_human50kb_240sp.txt',row.names=F,col.names=T,sep='\t',quote=F)


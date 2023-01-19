setwd('~/USS/aryn/Zoonomia/vcfs/snpeff_files')

fusil=read.table('~/USS/aryn/Zoonomia/genes_with_FUSIL_bins_20210922.txt',sep='\t',header=T)#IMPC data
buscogenes=read.table('~/USS/aryn/Zoonomia/datatables/buscogenes.txt')

files=list.files(path = ".", pattern = '_variants_snpeff_GQ_eff.txt.gz')
sps=gsub('_variants_snpeff_GQ_eff.txt.gz','',files)
coding_het=matrix(nrow=length(sps),ncol=43)
for (i in 1:length(sps)){
  sp=sps[i]
  print(sp)
dat=read.table(files[i],sep='\t',header=T,fill=T)
names(dat)=c('chr','pos','GT','GQ','DP','effect','impact','gene',"lof_gene",'lof_ntr','lof_pct')
dat=dat[!duplicated(dat[,c("chr",'pos')]),]
mgq=mean(dat$GQ,na.rm=T)
mdp=mean(dat$DP,na.rm=T)
dat=dat[which(dat$GQ>=80),]
thresh=median(dat$DP,na.rm=T)+3*sd(dat$DP,na.rm=T)

svs=unique(c(grep("insertion",dat$effect),grep("deletion",dat$effect),grep("frameshift",dat$effect)))
dat=dat[(1:nrow(dat) %in% svs==F),]
dat=as.data.table(dat[which(dat$DP<thresh),])
dat=dat[which(dat$gene %in% buscogenes$V1),]

genes=dat[,list(coding_vars=.N,high=sum(impact=="HIGH"),lof=sum(lof_gene!=""),missense=length(grep('missense',effect)),synonymous=length(grep('synonymous',effect)),DP=mean(DP),GQ=mean(GQ)),by='gene']
genes=merge(genes,fusil[,c('human_symbol','Viability.Phenotype.HOMs.HEMIs','confidence','Achilles_gene_effect_mean','FUSIL_bin')],by.x='gene',by.y='human_symbol',all.x=T)

head(table(dat$effect[which(dat$impact=="HIGH")])[order(table(dat$effect[which(dat$impact=="HIGH")]),decreasing=T)])
sum(dat$lof_gene!="")
genes$FUSIL_bin=factor(genes$FUSIL_bin,levels=c('CL','DL','SV','V'))

vars=c(sp,sum(genes$missense),sum(genes$synonymous),sum(genes$high),sum(genes$lof),sum(genes$coding_vars),mdp,mgq)
vars=c(vars,table(genes$FUSIL_bin))
vars=c(vars,aggregate(genes$synonymous,by=list(genes$FUSIL_bin),sum,na.rm=T,drop=F)$x)
vars=c(vars,aggregate(genes$missense,by=list(genes$FUSIL_bin),sum,na.rm=T,drop=F)$x)
vars=c(vars,aggregate(genes$high,by=list(genes$FUSIL_bin),sum,na.rm=T,drop=F)$x)
vars=c(vars,aggregate(genes$lof,by=list(genes$FUSIL_bin),sum,na.rm=T,drop=F)$x)
vars=c(vars,table(genes$Viability.Phenotype.HOMs.HEMIs))
vars=c(vars,aggregate(genes$synonymous,by=list(genes$Viability.Phenotype.HOMs.HEMIs),sum,na.rm=T,drop=F)$x)
vars=c(vars,aggregate(genes$missense,by=list(genes$Viability.Phenotype.HOMs.HEMIs),sum,na.rm=T,drop=F)$x)
vars=c(vars,aggregate(genes$high,by=list(genes$Viability.Phenotype.HOMs.HEMIs),sum,na.rm=T,drop=F)$x)
vars=c(vars,aggregate(genes$lof,by=list(genes$Viability.Phenotype.HOMs.HEMIs),sum,na.rm=T,drop=F)$x)
coding_het[i,1:43]=vars
}


coding_het=as.data.frame(coding_het)
coding_het[,2:ncol(coding_het)]=apply(coding_het[,2:ncol(coding_het)],2,function(x) as.numeric(as.character(x)))
names(coding_het)=c('Sp','missense','synonymous','high','lof','coding_vars','meanDP','meanGQ',paste('fusil_genes_',c('CL','DL','CV','V'),sep=""),paste('synonymous_fusil_',c('CL','DL','CV','V'),sep=""),paste('missense_fusil_',c('CL','DL','CV','V'),sep=""),paste('high_fusil_',c('CL','DL','CV','V'),sep=""),paste('lof_fusil_',c('CL','DL','CV','V'),sep=""),paste('impc_genes_',c('L','SV','V'),sep=""),paste('synonymous_impc_',c('L','SV','V'),sep=""),paste('missense_impc_',c('L','SV','V'),sep=""),paste('high_impc_',c('L','SV','V'),sep=""),paste('lof_impc_',c('L','SV','V'),sep=""))
coding_het$Sp=as.character(coding_het$Sp)
coding_het[is.na(coding_het)]=0

write.table(coding_het,'~/USS/aryn/Zoonomia/datatables/Coding_heterozygous_busco_IMPC_FUSIL_GQ80.txt',sep='\t',row.names=F,quote=F)

library(data.table)
library(alphaOutlier)

setwd('~/USS/aryn/Zoonomia/snps')
files=list.files(path = ".", pattern = '_snps.txt.gz')#list of files of substitutions for each species

#for each species, get distribution of #snps/1KB window, then test which windows are outliers on the right tail of the distribution and write outlier windows to file
for (i in 1:length(files)){#
	filei=files[i]
	sp=gsub("_snps.txt.gz","",filei)
	print(sp)
	snps=fread(filei,skip=2,tmpdir=getwd())
	snps=snps[,1:3]
	names(snps)=c('scaf','pos1','pos')
	snps$win=ceiling(snps$pos/1000)*1000
	block=snps[, list(snps=.N), by=c('scaf','win')]
	outs=aout.pois(block$snps, median(block$snps), alpha = 0.1, hide.outliers = FALSE)
	sum(outs$is.outlier)/length(outs$is.outlier)
	thresh=max(outs$data[outs$is.outlier==F])
	block$throwout=ifelse(block$snps>thresh,T,F)
	block$win1=block$win-1000
	block=block[,c(1,5,2:4)]
	block[,2] = format(block[,2], scientific = FALSE,trim=T)
	block[,3] = format(block[,3], scientific = FALSE,trim=T)
	write.table(block,paste(c('~/USS/aryn/Zoonomia/filter/',sp,'_mask_a0.1.txt'),collapse=''),sep='\t',row.names=F,quote=F)
}


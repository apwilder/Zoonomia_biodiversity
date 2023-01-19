setwd('~/USS/aryn/Zoonomia/busco/busco_downloads/lineages/mammalia_odb10')

busco=read.table('busco_genes.txt',sep='\t',header=T,quote="\"")
busco$altID1=NA
busco$altID2=NA
busco$altID3=NA
busco$altID4=NA
busco$altID1=str_split_fixed(busco$pub_gene_id,";",Inf)[,1]
busco$altID2=str_split_fixed(busco$pub_gene_id,";",Inf)[,2]
busco$altID3=str_split_fixed(busco$pub_gene_id,";",Inf)[,3]
busco$altID4=str_split_fixed(busco$pub_gene_id,";",Inf)[,4]
busco[,10:12]=apply(busco[,10:12],2,function(x) ifelse(x=='',NA,x))
busco1=busco[which(rowSums(is.na(busco[10:12]))==3),]
busco2=busco[which(rowSums(is.na(busco[10:12]))==2),]
busco3=busco[which(rowSums(is.na(busco[10:12]))==1),]
busco4=busco[which(rowSums(is.na(busco[10:12]))==0),]
busco1$pub_gene_id=busco1$altID1
busco2$pub_gene_id=busco2$altID2
busco3$pub_gene_id=busco3$altID3
busco4$pub_gene_id=busco4$altID4
busco=rbind(busco1,busco2,busco3,busco4)
busco=busco[,1:8]
buscogenes=unique(busco$pub_gene_id)
write.table(buscogenes,'~/USS/aryn/Zoonomia/datatables/buscogenes.txt',quote=F,row.names=F,col.names=F)

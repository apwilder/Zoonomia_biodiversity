library(ape)
library(phylolm)
library(MuMIn)
library(Polychrome)
library(RColorBrewer)
library(ggplot2)
library(reshape)

setwd('/home/centos/USS/aryn/Zoonomia/datatables')

##read in compiled data#####
#phylogeny
my_tree <- read.tree('../Zoonomia_ChrX_lessGC40_241species_30Consensus.tree')

#mutation-adjusted load stats
adjstats=read.table('Mutation-adjustedLoad_PSMC_het_241sp.txt',header=T,sep='\t')
adjstats$IUCN=factor(adjstats$IUCN,levels=c('LC','NT','VU','EN','CR','DD'))
row.names(adjstats)=adjstats$Sp
#raw (non-adjusted) load stats
stats=read.table('Load_PSMC_het_241sp.txt',header=T,sep='\t')
row.names(stats)=stats$Sp
stats$ppn_missense=stats$missense_vars/stats$coding_vars
stats$IUCN=factor(stats$IUCN,levels=c('LC','NT','VU','EN','CR','DD'))

#exclude species with fewer than 10,000 coding variants
adjstats$ppn_missense[adjstats$Sp %in% adjstats$Sp[which(stats$coding_vars<10000)]]=NA
stats$ppn_missense[which(stats$coding_vars<10000)]=NA

lp=read.table('IMPC_CodingVars_Ne_4plot.txt',header=T,sep='\t')
lpa=read.table('IMPC_CodingVars_adjusted_Ne_4plot.txt',header=T,sep='\t')

#IUCN colors for plotting
iucncols=brewer.pal(9,'Set1')[c(1,5,2,3,7)]
iucncols=c(iucncols[1],'gray',iucncols[2:5])
names(iucncols)=c('CR','DD','EN','LC','NT','VU')
#Order colors for plotting
P36 = createPalette(nlevels(stats$Order),  iucncols[c(1,3:6)])
names(P36)=levels(stats$Order)
P36=c(brewer.pal(12,'Paired'),brewer.pal(8,'Dark2')[8],brewer.pal(8,'Set2')[1],'black',brewer.pal(8,'Pastel2')[8],brewer.pal(4,'Set3')[4],brewer.pal(9,'Greens')[9],brewer.pal(9,'Blues')[9])
P36=P36[c(1,3,5,14,7,9,12:13,16:17,6,2,4,8,10,11,13,18,19)]
names(P36)=names(table(stats$Order)[order(table(stats$Order),decreasing=T)])

####statistical tests in main text#####
####Ne is lower in threatened species (phylolm)####
hist(I(stats$harm_mean_wt2*1e4))#should be log-transformed to deal with skew across Orders
test=phylolm(log10(harm_mean_wt2*1e4)~threatened,data=stats,phy=my_tree)
summary(test)#signif
#means for each group
10^test$coefficient[1]
10^(test$coefficient[1]+test$coefficient[2])

#including within orders
#CARNIVORA
hist(I(stats$harm_mean_wt2[which(stats$Order=='CARNIVORA')]*1e4))#no need to transform
test=phylolm(harm_mean_wt2*1e4~threatened,data=stats[which(stats$Order=='CARNIVORA'),],phy=my_tree)
summary(test)#7.409e-06 #distrib looks better when not log-transformed
#means for each group
test$coefficient[1]
test$coefficient[1]+test$coefficient[2]

#CETARTIODACTYLA
test=phylolm(harm_mean_wt2*1e4~threatened,data=stats[which(stats$Order=='CETARTIODACTYLA'),],phy=my_tree)
summary(test)#0.02338
#means for each group
test$coefficient[1]
test$coefficient[1]+test$coefficient[2]

#but not primates
#PRIMATES
test=phylolm(harm_mean_wt2*1e4~threatened,data=stats[which(stats$Order=='PRIMATES'),],phy=my_tree)
summary(test)#0.30670
#means for each group
test$coefficient[1]
test$coefficient[1]+test$coefficient[2]

#Figure 1B
gg=ggplot(stats,aes(x=factor(threatened,labels=c('Non-threatened','Threatened')),y=harm_mean_wt2*1e4))+xlab('IUCN status')+ylab(expression(italic('N'[e])))
gg+geom_violin(draw_quantiles =  0.5,size=1)+ geom_jitter(shape=16,position=position_jitter(width=0.1),col=P36[as.character(stats$Order)],size=2.5)+theme(text = element_text(size=14),axis.text.x = element_text(vjust=-3),axis.title.x=element_text(vjust=-5),panel.background = element_blank(),rect = element_rect(fill = "transparent"),axis.line = element_line(colour = "black",size=0.75),plot.margin = margin(0.1, 0.1, 1, 0.2, "cm"))+ scale_y_log10(labels = scales::comma)

#Figure 1C
plain <- function(x,...) {
  format(x, ..., scientific = FALSE, drop0trailing = TRUE)
}

gg=ggplot(stats,aes(x=factor(threatened,labels=c('Non-threatened','Threatened')),y=NeNc*10))+xlab('IUCN status')+ylab(expression(italic('N'[e]*'/N'[c])))
gg+geom_violin(draw_quantiles =  0.5,size=1)+ geom_jitter(shape=16,position=position_jitter(width=0.1),col=P36[as.character(stats$Order)],size=2.5)+theme(text = element_text(size=14),axis.text.x = element_text(vjust=-3),axis.title.x=element_text(vjust=-5),panel.background = element_blank(),rect = element_rect(fill = "transparent"),axis.line = element_line(colour = "black",size=0.75),plot.margin = margin(0.1, 0.1, 1, 0.2, "cm"))+ scale_y_log10(labels = plain)

#legend
plot(1:nlevels(stats$Order),col=P36[levels(stats$Order)],pch=16,cex=0)
legend('topleft',legend=levels(stats$Order),col=P36[levels(stats$Order)],pch=16)


####Ne/Nc higher in threatened species (phylolm)####
#Nc is correlated
hist(log10(stats$NeNc))#log-transform
test=phylolm(log10(NeNc)~threatened,data=stats,phy=my_tree)
summary(test)#signif + if log-transformed
#means for each group
10^test$coefficient[1]
10^(test$coefficient[1]+test$coefficient[2])

#also works for glm
summary(phyloglm(threatened~log10(NeNc),data=stats,phy=my_tree))
summary(phylolm(log10(NeNc)~threatened,data=stats[which(stats$Order=="PRIMATES"),],phy=my_tree))#signif within primates

####Ne vs. genetic load ####
#Ne vs. constrained
summary(phylolm(ppn_conserved~log10(harm_mean_wt2*1e4),data=adjstats,phy=my_tree))#p=0.009652
test=phylolm(ppn_conserved~log10(harm_mean_wt2*1e4),data=stats,phy=my_tree)
summary(test)#p=0.1315
#slope
test$coefficient[2]#every 10-fold increase in Ne changes load by -0.001142767

#Ne vs. missense
summary(phylolm(ppn_missense~log10(harm_mean_wt2*1e4),data=adjstats,phy=my_tree))#p=7.763e-05
test=phylolm(ppn_missense~log10(harm_mean_wt2*1e4),data=stats,phy=my_tree)
summary(test)#p=5.917e-08
#slope
test$coefficient[2]#every 10-fold increase in Ne changes load by -0.0198274

#Ne vs. missense constrained
summary(phylolm(ppn_miss_conserved~log10(harm_mean_wt2*1e4),data=adjstats,phy=my_tree))#p=0.004283
test=phylolm(ppn_miss_conserved~log10(harm_mean_wt2*1e4),data=stats,phy=my_tree)
summary(test)#p=0.2727728
#slope
test$coefficient[2]#every 10-fold increase in Ne changes load by -0.004937023

####IMPC coding vars vs. Ne########
regres=matrix(nrow=6,ncol=7)
i=1
#test variable types*IMPC categories
for(impact in c('Lethal','Viable')){
  for (variant in unique(lp$variant)) {
    for (zyg in unique(lp$zyg)) {
      loadtest=lp[which(lp$impact==impact & lp$variant==variant & lp$zyg==zyg  & !is.na(lp$threatened)),]
      if(nrow(loadtest)>0){
        row.names(loadtest)=loadtest$Sp
        loadtest$threatened=ifelse(loadtest$threatened=='threatened',1,0)
        plm=phylolm(value~log10(harm_mean_wt2*10000),data=loadtest,phy=my_tree)
        regres[i,1:2]=plm$coefficients
        regres[i,3]=summary(plm)$coefficients[2,4]
        regres[i,4]=plm$n
        regres[i,5:7]=c(impact,variant,zyg)
        i=i+1
      }
    }
  }
}
regres=as.data.frame(regres)
regres[,1:4]=apply(regres[,1:4],2,function(x) as.numeric(as.character(x)))
names(regres)=c('intercept','slope','p','n','impact','variant','zyg')
regres$x1=min(log10(lp$harm_mean_wt2*10000),na.rm=T)
regres$x2=max(log10(lp$harm_mean_wt2*10000),na.rm=T)
regres$y1=regres$x1*regres$slope+regres$intercept
regres$y2=regres$x2*regres$slope+regres$intercept
regres$y2[which(regres$y2<0)]=0

plain <- function(x,...) {
  format(x, ..., scientific = FALSE, drop0trailing = TRUE)
}

#Fig 2
gg=ggplot(lp[which(lp$impact %in% c("Subviable",'All')==F),],aes(x=log10(harm_mean_wt2*10000),y=value,col=IUCN))+xlab(expression('log'[10]*'('*italic('N'[e])*')'))+ylab('Deleterious mutations/coding mutations')+scale_color_manual(values = iucncols)
gg+facet_grid(factor(zyg,levels=c('Homozygous','Heterozygous'))*factor(variant,levels=c('Missense','LoF'))~factor(impact,levels=c('Lethal','Viable')),scales='free_y')+geom_point()+theme(text = element_text(size=14),axis.text.x = element_text(),panel.background = element_blank(),axis.line = element_line(colour = "black",size=0.5),strip.text= element_text(size=12),legend.position='none')+annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1.5)+annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1.5)+geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),data = regres,col='black',size=0.75)#+geom_text(data=regres,aes(x = 5, y = max(y2)), label = p)

#results are qualitatively the same for lpa (adjusted stats)
regresa=matrix(nrow=6,ncol=7)
i=1
#test variable types*IMPC categories
for(impact in c('Lethal','Viable')){
  for (variant in unique(lpa$variant)) {
    for (zyg in unique(lpa$zyg)) {
      loadtest=lpa[which(lpa$impact==impact & lpa$variant==variant & lpa$zyg==zyg  & !is.na(lpa$threatened)),]
      if(nrow(loadtest)>0){
        row.names(loadtest)=loadtest$Sp
        loadtest$threatened=ifelse(loadtest$threatened=='threatened',1,0)
        plm=phylolm(value~log10(harm_mean_wt2*10000),data=loadtest,phy=my_tree)
        regresa[i,1:2]=plm$coefficients
        regresa[i,3]=summary(plm)$coefficients[2,4]
        regresa[i,4]=plm$n
        regresa[i,5:7]=c(impact,variant,zyg)
        i=i+1
      }
    }
  }
}
regresa=as.data.frame(regresa)
regresa[,1:4]=apply(regresa[,1:4],2,function(x) as.numeric(as.character(x)))
names(regresa)=c('intercept','slope','p','n','impact','variant','zyg')
regresa$x1=min(log10(lpa$harm_mean_wt2*10000),na.rm=T)
regresa$x2=max(log10(lpa$harm_mean_wt2*10000),na.rm=T)
regresa$y1=regresa$x1*regresa$slope+regresa$intercept
regresa$y2=regresa$x2*regresa$slope+regresa$intercept
regresa$y2[which(regresa$y2<0)]=0

#There were proportionally fewer missense mutations in IMPC lethal genes relative to IMPC viable genes
summary(lm(value~factor(impact,levels=c('Viable','Lethal')),data=lp[which(lp$variant=='Missense' & lp$zyg=='Homozygous'),]))#<2e-16
summary(lm(value~factor(impact,levels=c('Viable','Lethal')),data=lp[which(lp$variant=='Missense' & lp$zyg=='Heterozygous'),]))#4.42e-09

#Fig S5
gg=ggplot(lp[which(!is.na(lp$threatened) & lp$impact %in% c("Subviable",'All')==F),],aes(x=factor(impact,levels=c('Lethal','Viable')),y=value))+xlab('IMPC gene category')+ylab('Deleterious mutations/coding mutations')
gg+facet_grid(factor(zyg,levels=c('Homozygous','Heterozygous'))*factor(variant,levels=c('Missense','LoF'))~.,scales='free_y')+geom_jitter(shape=16,position=position_jitter(width=0.15),size=2,aes(col=Order))+geom_violin(position='dodge',draw_quantiles =  0.5,fill='NA',size=0.7)+scale_color_manual(values=P36)+theme(text = element_text(size=14),axis.text.x = element_text(),panel.background = element_blank(),axis.line = element_line(colour = "black",size=0.5),strip.text= element_text(size=12),legend.position='none')+annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

#Fig S4
regressup=matrix(nrow=12,ncol=7)
i=1
#test variable types*IMPC categories
for(impact in unique(lp$impact)){
  for (variant in unique(lp$variant)) {
    for (zyg in unique(lp$zyg)) {
      loadtest=lp[which(lp$impact==impact & lp$variant==variant & lp$zyg==zyg  & !is.na(lp$threatened)),]
      if(nrow(loadtest)>0){
        row.names(loadtest)=loadtest$Sp
        loadtest$threatened=ifelse(loadtest$threatened=='threatened',1,0)
        plm=phylolm(value~log10(harm_mean_wt2*10000),data=loadtest,phy=my_tree)
        regressup[i,1:2]=plm$coefficients
        regressup[i,3]=summary(plm)$coefficients[2,4]
        regressup[i,4]=plm$n
        regressup[i,5:7]=c(impact,variant,zyg)
        i=i+1
      }
    }
  }
}
regressup=as.data.frame(regressup)
regressup[,1:4]=apply(regressup[,1:4],2,function(x) as.numeric(as.character(x)))
names(regressup)=c('intercept','slope','p','n','impact','variant','zyg')
regressup$x1=min(log10(lp$harm_mean_wt2*10000),na.rm=T)
regressup$x2=max(log10(lp$harm_mean_wt2*10000),na.rm=T)
regressup$y1=regressup$x1*regressup$slope+regressup$intercept
regressup$y2=regressup$x2*regressup$slope+regressup$intercept
regressup$y2[which(regressup$y2<0)]=0

gg=ggplot(lp,aes(x=log10(harm_mean_wt2*10000),y=value,col=IUCN))+xlab(expression('log'[10]*'('*italic('N'[e])*')'))+ylab('Deleterious mutations/coding mutations')+scale_color_manual(values = iucncols)
gg+facet_grid(factor(zyg,levels=c('Homozygous','Heterozygous'))*factor(variant,levels=c('Missense','LoF'))~factor(impact,levels=c('All','Lethal','Subviable','Viable')),scales='free_y')+geom_point()+theme(text = element_text(size=14),axis.text.x = element_text(),panel.background = element_blank(),axis.line = element_line(colour = "black",size=0.5),strip.text= element_text(size=12),legend.position='none')+annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1.5)+annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1.5)+geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),data = regressup,col='black',size=0.75)

#for reporting p-values for mutation-adjusted load estimates
regressup=matrix(nrow=12,ncol=7)
i=1
#test variable types*IMPC categories
for(impact in unique(lp$impact)){
  for (variant in unique(lp$variant)) {
    for (zyg in unique(lp$zyg)) {
      loadtest=lpa[which(lpa$impact==impact & lpa$variant==variant & lpa$zyg==zyg  & !is.na(lpa$threatened)),]
      if(nrow(loadtest)>0){
        row.names(loadtest)=loadtest$Sp
        loadtest$threatened=ifelse(loadtest$threatened=='threatened',1,0)
        plm=phylolm(value~log10(harm_mean_wt2*10000),data=loadtest,phy=my_tree)
        regressup[i,1:2]=plm$coefficients
        regressup[i,3]=summary(plm)$coefficients[2,4]
        regressup[i,4]=plm$n
        regressup[i,5:7]=c(impact,variant,zyg)
        i=i+1
      }
    }
  }
}
regressup=as.data.frame(regressup)
regressup[,1:4]=apply(regressup[,1:4],2,function(x) as.numeric(as.character(x)))
names(regressup)=c('intercept','slope','p','n','impact','variant','zyg')
regressup$x1=min(log10(lp$harm_mean_wt2*10000),na.rm=T)
regressup$x2=max(log10(lp$harm_mean_wt2*10000),na.rm=T)
regressup$y1=regressup$x1*regressup$slope+regressup$intercept
regressup$y2=regressup$x2*regressup$slope+regressup$intercept
regressup$y2[which(regressup$y2<0)]=0
regressup[which(regressup$zyg=='Homozygous'),1:7]


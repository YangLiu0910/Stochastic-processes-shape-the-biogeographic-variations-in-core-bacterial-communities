###PCoA and three test methods for bacterial communities compositions
library(vegan)
library(ggplot2)
library(ggrepel)
library(devtools)
library(ape)

core_bray<-vegdist(core_otu)#or overall_otu
res<- pcoa(core_bray)#PCoA
res$values$Relative_eig
a=as.matrix(res$values$Relative_eig)[1,1]
b=as.matrix(res$values$Relative_eig)[2,1]
si<-res$vectors[,1:2]#PCoA
si<-as.data.frame(si)
tfa$Compartment<-factor(tfa$Compartment,levels=c("BS","RS","R","S","L"))
p<-ggplot(si,aes(x=Axis.1,y=Axis.2))
p=p+geom_point(aes(color=Compartment,shape=Site), size=4.5,alpha= 0.75)+theme_bw()
p+xlab(paste("PCoA1=",round(a,4)*100,"%",sep=""))+ylab(paste("PCoA2=",round(b,4)*100,"%",sep="")) + labs(title = "Core Bray-Curtis PCoA analysis/ 80% confidence ellipses are shown")
p+stat_ellipse(aes(colour =Type),type="norm",linetype = 2,level=.9,lwd=1,show.legend =NA)

##Three test methods
#adonis
adonis(core_bray~Compartment*Site)
#anosim
anosim(core_bray,Compartment)
#mrpp
mrpp(core_bray,Compartment)
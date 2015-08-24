setwd("~/proj/Battle_GC/etude")

library(ggplot2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(devtools)

coverage <- tbl_df(read.csv("coverage.csv"))

plotCov <- function(caseNum = 1, typeG = 1,  makeLog=FALSE){

p1<-ggplot()+theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(),
        legend.key = element_blank())

  # reads per locus
if(typeG == 1) {
readsPerLocus <- coverageM %>% 
  filter(case==caseNum) %>% 
  group_by(case,trial,locus) %>% 
  summarize(totReads=sum(numReads)+1) 
p1 <- p1+
  geom_line(data=readsPerLocus, stat="density", aes(x=totReads, color=factor(trial)), show_guide=FALSE)

if(makeLog) {p1 <- p1 + scale_x_log10()+xlab("Number of reads per locus (log scale)")}
else{p1 <- p1 + xlab("Number of reads per locus")}
}

else if (typeG == 2) {
readsPerIndiv <- coverageM %>%
  filter(case==caseNum) %>%
  group_by(case,trial,indiv) %>%
  summarize(totReads=sum(numReads)+1) 
p1 <- p1+
  geom_line(data=readsPerIndiv, stat="density", aes(x=totReads, color=factor(trial)), show_guide=FALSE)

if(makeLog) {p1 <- p1 + scale_x_log10()+xlab("Number of reads per indiv (log scale)")}
else{p1 <- p1 + xlab("Number of reads per indiv")}
}

else{
copyPerHaplo <- coverageM %>% 
  filter(case==caseNum) %>%
  group_by(case,trial,locus) %>% 
  summarize(frac=sum(numReads*haplotype)/sum(numReads)) %>%
  filter(frac>=0)
p1 <- p1+
  geom_line(data=copyPerHaplo, stat="density", aes(x=frac, color=factor(trial)), show_guide=FALSE)+
  xlab("Frac of haplotype 1 per locus")+
  xlim(c(0,1))
}
  return(p1)
}

grid.newpage()

p1 <- lapply(1:30, function(i) {j=i%%3; k=ceiling(i/3); if((k==9 && j==1)) {plotCov(k,j, TRUE)} else{plotCov(k,j)} })
#gridp <- arrangeGrob(p1, ncol=3, nrow=10)

#grid.arrange(plotCov(1,2), ncol=3, nrow=10)
#g<-arrangeGrob(p1, nrow=3)
pdf("distrib.pdf",width = 34, height = 22,paper="USr")# height = 1051, width = 818)# paper = "USr")
do.call(grid.arrange, c(p1, ncol=3))
dev.off()

#pushViewport(viewport(layout = grid.layout(2, 2)))
#print(plot5, vp = vplayout(1, 1:2))
#print(plot6, vp = vplayout(2, 1))


#read reference
refGeno <-tbl_df(read.table("ref_geno.tbl", fill=TRUE))
allGeno <-tbl_df(read.table("all_geno.tbl", fill=TRUE, as.is = TRUE))
colnames(allGeno) <- c("origin", "case", "trial", "contig", "pos", "anc", "der", paste("indiv", 1:48, sep=""))

MatchToRef <- function(d){
  matchCt <- c(0,0,0,0)
  flipCt <- c(0,0,0)
  if("ref" %in% d$origin) {
    refIndx = which("ref"==d$origin)
    i = 0;
    matchCt[4] = 1
    for (m in c("freeBayes", "angsdSam", "angsdGatk")) {
      i = i+1;
      if(m %in% d$origin) {
        targIndx = which(m==d$origin)
        flip = 0
        if(d$anc[refIndx]!=d$anc[targIndx] && d$der[refIndx]!=d$der[targIndx]){
          if (d$anc[refIndx]==d$der[targIndx] && d$der[refIndx] == d$anc[targIndx]) {
            flip = 1
            flipCt[i] = 1
          }
          else {
            next;
          }
        }
        for(p in 8:55){
          if (!is.na(d[refIndx,p]))
          {
            if((flip == 0 && d[refIndx,p] == d[targIndx,p]) ||
             (flip==1 && d[refIndx,p] == abs(2-d[targIndx,p]))){
              matchCt[i] = matchCt[i]+1;}
          }
        }
    }}}
  return(c(matchCt,flipCt))
}

#reload(inst("plyr"))
numMatchTbl <- ddply(allGeno, .(case, trial, contig, pos), function(i)MatchToRef(i))

numMatchDf <-tbl_df(numMatchTbl)
colnames(numMatchDf) <- c("case", "trial", "contig", "pos", "fb", "angS", "angG", "ref", "fbF", "angSF", "angGF")

#reload(inst("dplyr"))
predRateComponent <- numMatchDf %>% group_by(case, trial) %>% summarise(nPos = sum(ref),
                                                   nSnpfb = sum(fb>0),
                                                   nSnpangS = sum(angS>0),
                                                   nSnpangG = sum(angG>0),
                                                   nGenofb = sum(fb),
                                                   nGenoangS = sum(angS),
                                                   nGenoangG = sum(angG),
                                                   nFlipfb = sum(fbF),
                                                   nFlipangS =sum(angSF),
                                                   nFlipangG = sum(angGF))

getnPred <- allGeno %>% group_by(case, trial) %>% summarise(nfb = sum(origin=="freeBayes"),
                                                            nangS = sum(origin=="angsdSam"),
                                                            nangG = sum(origin=="angsdGatk")) 

predRateCompl <- left_join(predRateComponent, getnPred, by=c("case", "trial")) 


tem <- predRateCompl %>% group_by(case, trial) %>% summarise(method="FreeBayes",
                                                      metrics="Recall (%)",
                                                      val=nSnpfb*100/nPos)
tem <- rbind(tem,
             predRateCompl %>% group_by(case, trial) %>% summarise(method="ANGSD SAMtools",
                                                             metrics="Recall (%)",
                                                             val=nSnpangS*100/nPos))
tem <- rbind(tem,
             predRateCompl %>% group_by(case, trial) %>% summarise(method="ANGSD GATK",
                                                             metrics="Recall (%)",
                                                             val=nSnpangG*100/nPos))


tem <- rbind(tem,
             predRateCompl %>% group_by(case, trial) %>% summarise(method="FreeBayes",
                                                             metrics="Precision (%)",
                                                             val=nSnpfb*100/nfb))
tem <- rbind(tem,
             predRateCompl %>% group_by(case, trial) %>% summarise(method="ANGSD SAMtools",
                                                                   metrics="Precision (%)",
                                                                   val=nSnpangS*100/nangS))
tem <- rbind(tem,
             predRateCompl %>% group_by(case, trial) %>% summarise(method="ANGSD GATK",
                                                                   metrics="Precision (%)",
                                                                   val=nSnpangG*100/nangG))


tem <- rbind(tem,
             predRateCompl %>% group_by(case, trial) %>% summarise(method="FreeBayes",
                                                                   metrics="Genotype Accuracy (%)",
                                                                   val=nGenofb*100/((48*(case[1]!=2)+
                                                                                       10*(case[1]==2))*nSnpfb)))
tem <- rbind(tem,
             predRateCompl %>% group_by(case, trial) %>% summarise(method="ANGSD SAMtools",
                                                                   metrics="Genotype Accuracy (%)",
                                                                   val=nGenoangS*100/((48*(case[1]!=2)+
                                                                                         10*(case[1]==2))*nSnpangS)))
finalRate <- rbind(tem,
             predRateCompl %>% group_by(case, trial) %>% summarise(method="ANGSD GATK",
                                                                   metrics="Genotype Accuracy (%)",
                                                                   val=nGenoangG*100/((48*(case[1]!=2)+
                                                                                         10*(case[1]==2))*nSnpangG)))


finalRate$metrics = factor(finalRate$metrics, levels=c('Recall (%)','Precision (%)',"Genotype Accuracy (%)"))
rate.stat <- finalRate %>% group_by(case, method, metrics) %>% summarise(mean=round(mean(val),2),
                                                                         sd=round(sd(val),2))

save(finalRate,file="finalRate.Rda")
#reload(inst("ggplot2"))
ggplot(finalRate, aes(y=method, x=val))+
  xlab("")+
  geom_point(alpha=0.9, color="grey")+
  geom_text(data=rate.stat, aes(y = method, x = 65, label = format(paste0("", mean,"%\t (sd: ", sd, ")" ), nsmall = 1)), size = 4)+
  #geom_rect(aes(ymin = -2, ymax = 0, xmin=60, xmax = 100), fill="white")+
  facet_grid(case~metrics, scales="free_x")+
  theme_bw()+
  theme(panel.border = element_rect(colour = "light grey"),#element_blank(), 
        #axis.line = element_line(colour = "white"),
        axis.text.y = element_text(),
        legend.key = element_blank())

  

#  theme(plot.margin = unit(c(0.5,1, 0, 0.5), "lines")) + xlab(NULL) + ylab(NULL)
    




#mean = 80000
#variance = 400000000
#shape = mean*mean /variance
#scale = variance /mean
#hist(rgamma(96,shape=shape,scale=scale))

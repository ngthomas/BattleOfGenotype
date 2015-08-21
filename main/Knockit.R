setwd("~/proj/Battle_GC/etude")

library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

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
#grid.arrange(p1,p2,p3, ncol=3, nrow=1)
g<-arrangeGrob(p1, p2, p3, nrow=3)









#mean = 80000
#variance = 400000000
#shape = mean*mean /variance
#scale = variance /mean
#hist(rgamma(96,shape=shape,scale=scale))

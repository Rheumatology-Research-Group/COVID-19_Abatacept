library(survival)
library(edgeR)
library(limma)
###get Pheno

library(GEOquery)

GSE<-"GSE157103"
set1<- getGEO(GSE,GSEMatrix=T,destdir="./",AnnotGPL=F)
set1<-set1[[1]]

####get pheno

pheno<-pData(set1)
pheno$Group<-gsub("disease state: ","",pheno$characteristics_ch1)

###clean pheno


pheno<-pheno[,63:84]


##load scores

scores<-readRDS("Abatacept_Score/AbataceptScore.rds")
scores$X<-rownames(scores)
PHENO<-merge(pheno,scores,by="row.names")

#####split in two groups acording simple pos/neg score

high.scores<-as.character(scores[scores$score>0,"X"])
low.scores<-as.character(scores[scores$score<0,"X"])

##Suppl_Fig7A
###using quantiles extreme values
percent<-quantile(as.numeric(scores$score), probs = seq(0, 1, by= 0.1))##in 10%

###plot scores for all samples
h <- hist(scores$score, plot=F) # h$breaks and h$mids

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
cols <- rbPal(10)[as.numeric(cut(percent,breaks = 10))]

k <- cols[findInterval(h$mids,percent, rightmost.closed=T, all.inside=F) + 1]
# plot the histogram with the colours
plot(h, col=k,main = "Score Distribution",xlab = "Abatacept score" )
#legend(20, 28, h$mids, lwd=2, col=cols)
abline(v=c(as.numeric(percent[8]),as.numeric(percent[3])), col=c("red", "red"), lty=c(2,2), lwd=c(3, 3))


#####then try to deal with three

high.scores<-as.character(scores[scores$score>=percent[8],"X"])##8 works for mech intub as censor
low.scores<-as.character(scores[scores$score<=percent[3],"X"])##3 works for mech intub as censor

####subset in COVID19 pos and negative

my.threshold<-c(as.numeric(percent[8]),as.numeric(percent[3]))

scores.pos<-PHENO[PHENO$Group=="COVID-19",]

##Suppl_Fig7B

percent<-quantile(as.numeric(scores.pos$score), probs = seq(0, 1, by= 0.1))##in 10%

###plot scores for all samples
h <- hist(scores.pos$score, plot=F) # h$breaks and h$mids

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
cols <- rbPal(10)[as.numeric(cut(percent,breaks = 10))]

k <- cols[findInterval(h$mids,percent, rightmost.closed=T, all.inside=F) + 1]
# plot the histogram with the colours
plot(h, col=k,main = "Score Distribution COVID-19 patients",xlab = "Abatacept score" )
abline(v=my.threshold, col=c("red", "red"), lty=c(2,2), lwd=c(3, 3))

high.scores<-as.character(scores.pos[scores.pos$score>=percent[8],"X"])##8 works for mech intub as censor
low.scores<-as.character(scores.pos[scores.pos$score<=percent[3],"X"])##3 works for mech intub as censor


####same for non COVID-19 patients

my.threshold<-c(as.numeric(percent[8]),as.numeric(percent[3]))

scores.neg<-PHENO[PHENO$Group=="non-COVID-19",]


##Suppl_Fig7C

percent<-quantile(as.numeric(scores.neg$score), probs = seq(0, 1, by= 0.1))##in 10%

###plot scores for all samples
h <- hist(scores.neg$score, plot=F) # h$breaks and h$mids

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
cols <- rbPal(10)[as.numeric(cut(percent,breaks = 10))]

k <- cols[findInterval(h$mids,percent, rightmost.closed=T, all.inside=F) + 1]
# plot the histogram with the colours
plot(h, col=k,main = "Score Distribution non COVID-19 patients",xlab = "Abatacept score",
     xlim=c(-9,25))
abline(v=my.threshold, col=c("red", "red"), lty=c(2,2), lwd=c(3, 3))


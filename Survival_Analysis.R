##Survival figures script generator
library(survival)
library(edgeR)
library(limma)
library(survminer)

###get Pheno from published data

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

#####split in two groups acording simple pos/neg score and see how they behave

high.scores<-as.character(scores[scores$score>0,"X"])
low.scores<-as.character(scores[scores$score<0,"X"])

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

#####split ABT score according previous hist plots

high.scores<-as.character(scores[scores$score>=percent[8],"X"])##8 works for mech intub as censor
low.scores<-as.character(scores[scores$score<=percent[3],"X"])##3 works for mech intub as censor
###now, do the Survivalplots


###Fig4A
################Ventilation

toplot<-pheno
toplot$G<-ifelse(rownames(toplot) %in% low.scores,"LOW",ifelse(rownames(toplot) %in% high.scores,"HIGH",NA))###

###Ventilation
##days using ventilation
toplot$Vent_days<-as.numeric(toplot$`ventilator-free days:ch1`)##recuperation?
toplot<-toplot[!(toplot$Vent_days==0),]
toplot$vital_status<-ifelse(toplot$`mechanical ventilation:ch1`=="yes",1,0)

toplot$age<-as.numeric(toplot$`age (years):ch1`)
toplot$sex<-as.factor(toplot$`Sex:ch1`)
toplot$mech_vent<-as.factor(toplot$`mechanical ventilation:ch1`)

####remove samples without age or sex

toplot<-toplot[!(is.na(toplot$age) | is.na(toplot$sex)),]
toplot<-toplot[!(toplot$sex=="unknown"),]
toplot$G<-as.factor(toplot$G)
toplot<-toplot[!(is.na(toplot$G)),]
toplot$G<-droplevels(toplot$G)

toplot$sex<-droplevels(toplot$sex)

####subset only to patients

toplot<-toplot[toplot$Group=="COVID-19",]

####Plotting

####coxph adjusted with covariates
####sex as mean age

toplot$AGE<-ifelse(toplot$age>59,"1","0")###mean is 59 years old


fit2<- coxph(Surv(as.numeric(Vent_days), as.numeric(vital_status)) ~ G +AGE+sex, data = toplot)

sum.f<-summary(fit2)
pval<-round(as.numeric(sum.f$coefficients[,"Pr(>|z|)"][1]),digits = 4)
####plot coxph
surv.plot<-ggadjustedcurves(fit2, data = toplot,
                   variable  = "G",   # Variable of interest
                    palette = "npg",             # nature publishing group color palettes
                    curv.size = 2                # Change line size
)
p<-surv.plot

p<-surv.plot+ylab("% Ventilation free days")

p+theme(legend.position = "right")

##Fig 4B
####same with hospitalization days
####get pheno again, because we will merge with another label

pheno<-pData(set1)
pheno$Group<-gsub("disease state: ","",pheno$characteristics_ch1)

pheno<-pheno[,c(1,63:84)]
pheno$Sample_label<-ifelse(grepl("^COVID",pheno$title, perl=T),substring(pheno$title,1,9),substring(pheno$title,1,11))
pheno$Sample_label<-gsub("_$","",pheno$Sample_label)

toplot<-pheno
toplot$G<-ifelse(rownames(toplot) %in% low.scores,"LOW",ifelse(rownames(toplot) %in% high.scores,"HIGH",NA))###

###use here the survival data coming from pheno

surv.data<-read.delim("Source_data/Survival_data.csv")
##order based on Sample_Label
surv.data<-surv.data[order(surv.data$Sample_label),]

###merge

toplot<-merge(toplot,surv.data,by="Sample_label")

toplot$age<-as.numeric(toplot$`age (years):ch1`)
toplot$sex<-as.factor(toplot$`Sex:ch1`)
toplot$mech_vent<-as.factor(toplot$`mechanical ventilation:ch1`)

####remove samples without age or sex

toplot<-toplot[!(is.na(toplot$age) | is.na(toplot$sex)),]
toplot<-toplot[!(toplot$sex=="unknown"),]
toplot$G<-as.factor(toplot$G)
toplot<-toplot[!(is.na(toplot$G)),]
toplot$G<-droplevels(toplot$G)

toplot$sex<-droplevels(toplot$sex)
####subset only to patients

toplot<-toplot[toplot$Group=="COVID-19",]
####censors, herer HFD45
toplot$vital_status<-ifelse(toplot$`hospital-free days post 45 day followup (days):ch1`<25,1,0)##HFD45

####Plotting

####coxph adjusted with covariates
####sex as mean age

toplot$AGE<-ifelse(toplot$age>59,"1","0")###mean is 59 years old


fit2<- coxph(Surv(as.numeric(Stay), as.numeric(vital_status)) ~ G +AGE+sex, data = toplot)

sum.f<-summary(fit2)
pval<-round(as.numeric(sum.f$coefficients[,"Pr(>|z|)"][1]),digits = 4)
####plot coxph
surv.plot<-ggadjustedcurves(fit2, data = toplot,
                            variable  = "G",   # Variable of interest
                            palette = "npg",             # nature publishing group color palettes
                            curv.size = 2                # Change line size
)
p<-surv.plot+ylab("% Hospitalization days")

p+theme(legend.position = "right")

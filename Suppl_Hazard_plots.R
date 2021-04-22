####Suppl Fig 2A
###Load Hazards Mech info
mydf<-read.delim("Source_data/Hazards_Table_Mech.csv",na.strings = "")
mydf$SubgroupH<-as.character(mydf$SubgroupH)
mydf$Subgroup<-as.character(mydf$Subgroup)
mydf$Pvalue<-as.character(mydf$Pvalue)
####plot


rowseq <- seq(nrow(mydf),1)
par(mai=c(1,0,0,0))
plot(mydf$HazardRatio, rowseq, pch=15,
     xlim=c(-30,115), ylim=c(0,12),
     xlab='', ylab='', yaxt='n', xaxt='n',
     bty='n')
axis(1, seq(-30,115,by=10), cex.axis=.5)

segments(1,-1,1,12, lty=3)
segments(mydf$HazardLower, rowseq, mydf$HazardUpper, rowseq)

mtext('Benefit',1, line=2.5, at=-5, cex=1, font=2)
mtext('Risk',1.5, line=2.5, at=25, cex=1, font=2)

text(-25,10, "Subgroup", cex=1.25, font=2, pos=2)
t1h <- ifelse(!is.na(mydf$SubgroupH), mydf$SubgroupH, '')
text(-22,rowseq, t1h, cex=1.25, pos=4, font=3)
t1 <- ifelse(!is.na(mydf$Subgroup), mydf$Subgroup, '')
text(-21.5,rowseq, t1, cex=1.25, pos=4)

text(-4,12, "No. of\nPatients", cex=1.25, font=2, pos=2)
t2 <- ifelse(!is.na(mydf$NoOfPatients), format(mydf$NoOfPatients,big.mark=","), '')
text(-4, rowseq, t2, cex=1.25, pos=2)

#text(70,12, substitute(paste(italic("Coxph Mech Vent Days"))), cex=2, font=2, pos=2)

text(20,12, "Hazard Ratio (95%)", cex=1.25, font=2, pos=2)
t3 <- ifelse(!is.na(mydf$HazardRatio), with(mydf, paste(HazardRatio,' (',HazardLower,'-',HazardUpper,')',sep='')), '')
text(80,rowseq, t3, cex=1.25, pos=4)

text(110,12, "P Value", cex=1.25, font=4, pos=2)
t4 <- ifelse(!is.na(mydf$Pvalue), mydf$Pvalue, '')
text(108,rowseq, t4, cex=1.25, pos=4)

####Suppl Fig 2B
###Load Hazards HFD45 Info
mydf<-read.delim("Source_data/Hazards_Table.csv",na.strings = "")
mydf$SubgroupH<-as.character(mydf$SubgroupH)
mydf$Subgroup<-as.character(mydf$Subgroup)
mydf$Pvalue<-as.character(mydf$Pvalue)
####plot


rowseq <- seq(nrow(mydf),1)
par(mai=c(1,0,0,0))
plot(mydf$HazardRatio, rowseq, pch=15,
     xlim=c(-10,25), ylim=c(0,12),
     xlab='', ylab='', yaxt='n', xaxt='n',
     bty='n')
axis(1, seq(-10,25,by=10), cex.axis=.5)

segments(1,-1,1,12, lty=3)
segments(mydf$HazardLower, rowseq, mydf$HazardUpper, rowseq)

mtext('Benefit',1, line=2.5, at=-5, cex=1, font=2)
mtext('Risk',1.5, line=2.5, at=5, cex=1, font=2)

text(-7,10, "Subgroup", cex=1.25, font=2, pos=2)
t1h <- ifelse(!is.na(mydf$SubgroupH), mydf$SubgroupH, '')
text(-7,rowseq, t1h, cex=1.25, pos=4, font=3)
t1 <- ifelse(!is.na(mydf$Subgroup), mydf$Subgroup, '')
text(-5.75,rowseq, t1, cex=1.25, pos=4)

text(-2,12, "No. of\nPatients", cex=1.25, font=2, pos=2)
t2 <- ifelse(!is.na(mydf$NoOfPatients), format(mydf$NoOfPatients,big.mark=","), '')
text(-2, rowseq, t2, cex=1.25, pos=2)

#text(18,12, substitute(paste(italic("Coxph HFD-45"))), cex=2, font=2, pos=2)

text(8,12, "Hazard Ratio (95%)", cex=1.25, font=2, pos=2)
t3 <- ifelse(!is.na(mydf$HazardRatio), with(mydf, paste(HazardRatio,' (',HazardLower,'-',HazardUpper,')',sep='')), '')
text(16,rowseq, t3, cex=1.25, pos=4)

text(22,12, "P Value", cex=1.25, font=4, pos=2)
t4 <- ifelse(!is.na(mydf$Pvalue), mydf$Pvalue, '')
text(22,rowseq, t4, cex=1.25, pos=4)






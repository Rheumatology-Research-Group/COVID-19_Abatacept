library(corrplot)
library(ggplot2)
library(org.Hs.eg.db)
library(GO.db)
library(gridExtra)
library(fgsea)
library(prabclus)
library(ggrepel)


# Select curated BPs to analyze ==============================

bps1 <- c("GO:0002224","GO:0039528","GO:0060337","GO:0016197", # Virus sensing
          "GO:0035747","GO:0042267" # NKs
)

bps2 <- c("GO:0002456","GO:0002369","GO:0019882", # T cell and antigen processing and presentation
          "GO:0002573","GO:0042116","GO:0034341","GO:0071356", # Monocytes and macrophages activation and differentiation (& response to IFNG i TNF)
          "GO:0032612","GO:0032635","GO:0032640","GO:0032637",  # Monocyte and macrophage - produced cytokines 
          "GO:0098760","GO:0032613","GO:0006956", # Other associated cytokines + complement activation
          "GO:0019724", # B cell mediated immunity
          "GO:0030193" # Blood coagulation
          
)

# Determine genes in BP
dictGO <- as.list(org.Hs.egGO2ALLEGS)
dictTERMS <- as.list(GOTERM)
dictGO <- lapply(dictGO, function(x) {unique(x)})
dictGO <- dictGO[c(bps1,bps2)]
bps1 <- sapply(bps1, function(x) Term(dictTERMS[[x]]))
bps2 <- sapply(bps2, function(x) Term(dictTERMS[[x]]))
names(dictGO) <- c(bps1,bps2)
rm(dictTERMS)

genes_dict <- readRDS("Source_data/symbol2entrez.RDS")
genes_dict <- genes_dict[!is.na(genes_dict$ENTREZ_ID),]
rownames(genes_dict) <- genes_dict$ENTREZ_ID
dictGO <- sapply(dictGO, function(x) {genes_dict$SYMBOL[genes_dict$ENTREZ_ID %in% x]})

# Plot overlap among selected pathways (Jaccard similarity)
genes <- unique(unlist(dictGO))
genes_in_bp_mat <- matrix(0,nrow=length(genes),ncol=length(dictGO),dimnames=list(genes,names(dictGO)))
for (path in names(dictGO)){
  path_genes <- dictGO[[path]]
  genes_in_bp_mat[path_genes,path] <- 1
}
Sim <- 1 - jaccard(genes_in_bp_mat)
row.names(Sim)[2] <- colnames(Sim)[2] <- "cytoplasmic pattern recognition receptor \nsignaling pathway in response to virus"

pdf("Figures/SupplementaryFigure6.pdf")
corrplot::corrplot(Sim, is.corr=T, tl.cex = 0.7, pch.cex=0.7, type = "upper",
         tl.col = "black", cl.lim = c(0,1), method = "pie", diag = F)
dev.off()
rm(genes_in_bp_mat,Sim)

# Datasets 
dt <- gsea <- gsea_tb <- list(earlyCOVID=NULL,lateCOVID=NULL,Abatacept=NULL)
dt[["earlyCOVID"]] <- readRDS("DEG/Xiong_COVID19_CtrlDisease.RDS")
dt[["lateCOVID"]] <- readRDS("DEG/COVIDSevere_Mild_HFD45_25cut.RDS")
dt[["lateCOVID"]]$symbols <- row.names(dt[["lateCOVID"]])
dt[["Abatacept"]] <- readRDS("DEG/RA_Abatacept.RDS")
dt[["Abatacept"]]$symbols <- row.names(dt[["Abatacept"]])

# GSEA of selected BPs
for (dt_name in names(dt)){
  vals2test <- -log10(dt[[dt_name]]$PValue) * sign(dt[[dt_name]]$logFC) 
  names(vals2test) <- dt[[dt_name]]$symbols
  vals2test <- vals2test[order(vals2test,decreasing = T)]
  
  gsea[[dt_name]] <- fgsea(pathways=dictGO,stats=vals2test,nperm=1000000,minSize=1,maxSize=300,nproc=3)
  gsea_tb[[dt_name]] <- as.data.frame(gsea[[dt_name]])[,c('pathway','NES','pval','padj','leadingEdge')]
  gsea[[dt_name]] <- gsea[[dt_name]][,c('pathway','NES','pval','padj')]
  names(gsea[[dt_name]]) <- c("Description","NES","pvalue","p.adjust")
  gsea[[dt_name]] <- as.data.frame(gsea[[dt_name]])
  names(gsea_tb[[dt_name]]) <- paste(dt_name,names(gsea_tb[[dt_name]]),sep="_")
  row.names(gsea[[dt_name]]) <- gsea[[dt_name]]$Description
}

# Supplementary table
gsea_tb <- do.call(cbind, list(gsea_tb$Abatacept, gsea_tb$earlyCOVID[,2:4], gsea_tb$lateCOVID[,2:4]))
gsea_tb$Abatacept_leadingEdge <- sapply(gsea_tb$Abatacept_leadingEdge, function(x) paste(x, collapse=", "))
names(gsea_tb)[1] <- "Pathway"
row.names(gsea_tb) <- c(names(bps1),names(bps2))
write.table(gsea_tb, "Tables/SupplementaryTable2.csv", sep=";")

# FIGURE: Curated BPs modulation by Abatacept summary ====================================
earlyCOVID <- gsea[["earlyCOVID"]]
lateCOVID <- gsea[["lateCOVID"]]
Abatacept <- gsea[["Abatacept"]]
nm1 <- "COVID-19 vs Control"
nm2 <- "COVID-19 Severe vs Moderate"
nm3 <- "Abatacept in RA"
nm_ref <- "COVID-19 literature"

Descriptions1 <- earlyCOVID[bps1,"Description"]
Descriptions2 <- earlyCOVID[bps2,"Description"]
nes1_ref <- c(3.5,3.5,3.5,3.5,3.5,-3.5)
nes2_ref <- c(3.5,3.5,3.5,3.5,3.5,
              3.5,3.5,3.5,3.5,3.5,
              3.5,3.5,3.5,3.5,3.5,3.5)
pval1_ref <- rep(0.5,length(bps1))
pval2_ref <- rep(0.5,length(bps2))

summary_nes1 <- cbind(nes1_ref, earlyCOVID[bps1,"NES"], lateCOVID[bps1,"NES"], Abatacept[bps1,"NES"])
summary_pvals1 <- cbind(pval1_ref, earlyCOVID[bps1,"pvalue"], lateCOVID[bps1,"pvalue"], Abatacept[bps1,"pvalue"])
row.names(summary_nes1) <- row.names(summary_pvals1) <- Descriptions1

summary_nes2 <- cbind(nes2_ref, earlyCOVID[bps2,"NES"], lateCOVID[bps2,"NES"], Abatacept[bps2,"NES"])
summary_pvals2 <- cbind(pval2_ref, earlyCOVID[bps2,"pvalue"], lateCOVID[bps2,"pvalue"], Abatacept[bps2,"pvalue"])
row.names(summary_nes2) <- row.names(summary_pvals2) <- Descriptions2

colnames(summary_nes1) <- colnames(summary_pvals1) <- colnames(summary_nes2) <- colnames(summary_pvals2) <- c(nm_ref,nm1,nm2,nm3)
row.names(summary_nes1)[2] <- row.names(summary_pvals1)[2] <- "cytoplasmic pattern recognition receptor \nsignaling pathway in response to virus"

pdf("Figures/figure_1A.pdf", width = 8, height = 7.2)
corrplot::corrplot(as.matrix(summary_nes1), is.corr=F, p.mat=as.matrix(summary_pvals1), 
         sig.level=c(0.00001,0.0001,0.001,0.01,0.05), tl.cex = 1, pch.cex=0.7, 
         insig="label_sig", tl.col = "black", cl.lim = c(-4,4), cl.align.text='l')
dev.off()

pdf("Figures/figure_1B.pdf")
corrplot::corrplot(as.matrix(summary_nes2), is.corr=F, p.mat=as.matrix(summary_pvals2), 
         sig.level=c(0.00001,0.0001,0.001,0.01,0.05), tl.cex = 1, pch.cex=0.7, cl.ratio = 0.5, 
         insig="label_sig", tl.col = "black", cl.lim = c(-4,4),  cl.align.text='l')
dev.off()

# Curated BPs modulation by Abatacept volcano plots =================================
# (base plot with a fraction of the points)
data <- dt[["Abatacept"]]
data$SIG <- ifelse(data$logFC < 0, 'NEG', 'POS')
data$PValue <- -log10(data$PValue)

n_parts <- round(nrow(data)/5)
del_smp <- c(sample(50:n_parts,(n_parts/2)),
             sample((n_parts+1):(n_parts*2),(n_parts/1.5)),
             sample(((n_parts*2)+1):(n_parts*3),(n_parts/1.3)),
             sample(((n_parts*3)+1):(n_parts*4),(n_parts/1.1)),
             sample(((n_parts*4)+1):(n_parts*5),(n_parts/1.05)))
base_data <- data[-del_smp, ]

# FIGURE: BP1 volcano plots ---------------
p <- vector(mode="list",length=3)
bps1[2] <- names(dictGO)[2] <- "cytoplasmic pattern recognition receptor \nsignaling pathway in response to virus"

for (line in 0:2){
  bp_line <- bps1[(1:2)+line*2]
  bp_line <- bp_line[!is.na(bp_line)]
  
  base_data_line <- NULL
  bp_data_line <- NULL
  for (bp in bp_line){
    base_data_line <- rbind(base_data_line,base_data)
    genes <- dictGO[[bp]]
    bp_data_line[[bp]] <- data[row.names(data) %in% genes, ]
    bp_data_line[[bp]]$BP <- bp
  }
  bp_data_line <- do.call(rbind,bp_data_line)
  bp_data_line$BP <- factor(bp_data_line$BP, levels=bp_line)
  base_data_line$BP <- factor(rep(bp_line,each=nrow(base_data)),levels=bp_line)
  
  
  p[[line+1]] <- ggplot() +
    geom_point(aes(x = logFC, y = PValue), data = base_data_line, color =  "darkgray", size = 3) +
    geom_point(aes(x = logFC, y = PValue, color = SIG), data = bp_data_line, size = 3, show.legend = F) +
    facet_grid(. ~ BP) +
    scale_color_manual(values=c("POS"="#0000cc","NEG"="#cc0000")) +
    scale_shape_identity() +
    geom_hline(yintercept=-log10(0.05),linetype="twodash",col="seagreen4") +
    xlab("Effect size") + 
    ylab("-log10 P Value") +
    theme(legend.position="none", 
          panel.border=element_rect(colour = "darkgray", fill=NA)) +
    theme_classic(base_size = 29, base_line_size = 1.25, base_rect_size = 1)
}

jpeg(file="Figures/SupplementaryFigure1A.jpeg",width = 1080, height = 1350)
grid.arrange(p[[1]], p[[2]], p[[3]], ncol = 1, nrow = 3)
dev.off()

# FIGURE: BP2 volcano plots ---------------------
bps2[7] <- names(dictGO)[13] <- "cellular response to TNF"
p <- vector(mode="list",length=4)
for (line in 0:3){
  bp_line <- bps2[(1:4)+line*4]
  bp_line <- bp_line[!is.na(bp_line)]
  
  base_data_line <- NULL
  bp_data_line <- NULL
  for (bp in bp_line){
    base_data_line <- rbind(base_data_line,base_data)
    genes <- dictGO[[bp]]
    bp_data_line[[bp]] <- data[row.names(data) %in% genes, ]
    bp_data_line[[bp]]$BP <- bp
  }
  bp_data_line <- do.call(rbind,bp_data_line)
  bp_data_line$BP <- factor(bp_data_line$BP, levels=bp_line)
  base_data_line$BP <- factor(rep(bp_line,each=nrow(base_data)),levels=bp_line)
  
  
  p[[line+1]] <- ggplot() +
    geom_point(aes(x = logFC, y = PValue), data = base_data_line, color =  "darkgray", size = 3) +
    geom_point(aes(x = logFC, y = PValue, color = SIG), data = bp_data_line, size = 3, show.legend = F) +
    facet_grid(. ~ BP) +
    scale_color_manual(values=c("POS"="#0000cc","NEG"="#cc0000")) +
    scale_shape_identity() +
    geom_hline(yintercept=-log10(0.05),linetype="twodash",col="seagreen4") +
    xlab("Effect size") + 
    ylab("-log10 P Value") +
    theme(legend.position="none", 
          panel.border=element_rect(colour = "darkgray", fill=NA)) +
    theme_classic(base_size = 29, base_line_size = 1.25, base_rect_size = 1)
  
}

jpeg(file="Figures/SupplementaryFigure1B.jpeg",width = 1600, height = 1440)
grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], ncol = 1, nrow = 4)
dev.off()

library(org.Hs.eg.db)
library(GO.db)
library(metap)
library(corrplot)
library(eulerr)
library(prabclus)

# Choose between early or late COVID cohort
analysis <- "early_cohort" # late_cohort / early_cohort

# Find BPs clusters to reduce redundancy =========================
# Load BP GO database
if (! file.exists("Source_data/GO3.10_h06_BP_clusters.rds")){
  dictGO <- as.list(org.Hs.egGO2ALLEGS)
  dictGO <- lapply(dictGO, function(x) {unique(x)})
  dictGO <- dictGO[sapply(names(dictGO),function(x) length(dictGO[[x]])>9&length(dictGO[[x]])<301)]
  dictTERMS <- as.list(GOTERM)
  dictGO <- dictGO[sapply(names(dictGO), function(x){Ontology(dictTERMS[[x]])=='BP'})]
  
  # Jaccard similarity matrix
  genes <- unique(unlist(dictGO))
  genes_in_bp_mat <- matrix(0,nrow=length(genes),ncol=length(dictGO),dimnames=list(genes,names(dictGO)))
  for (path in names(dictGO)){
    path_genes <- dictGO[[path]]
    genes_in_bp_mat[path_genes,path] <- 1
  }
  Sim <- 1 - jaccard(genes_in_bp_mat)
  rm(genes_in_bp_mat)
  
  # Cluster BPs
  tree <- hclust(d=as.dist(1-Sim))
  cut <- cutree(tree, h=0.60)
  table(table(cut))
  
  bps_cluster <- data.frame(GOID = names(dictGO),
                            Description = sapply(names(dictGO), function(x) Term(dictTERMS[[x]])),
                            Cluster = cut[names(dictGO)],
                            Size = sapply(dictGO, function(x) length(x), USE.NAMES=F),
                            stringsAsFactors = FALSE)
} else {
  bps_cluster <- readRDS("Source_data/GO3.10_h06_BP_clusters.rds")
}


# GSEA results antagonism =========================
# Load datasets and write figure / table paths
if (analysis == "late_cohort"){
  # Figure 3 datasets
  COVID <- readRDS("GSEA/Overmeyer_COVID19_SevereModerate.RDS")
  Abatacept <- readRDS("GSEA/RA_Abatacept.RDS")  
  
  tables_path <- c("Tables/SupplementaryTable6.csv","Tables/SupplementaryTable7.csv","Tables/SupplementaryTable8.csv")
  figures_path <- c("Figures/figure_3A.jpeg","Figures/figure_3B.pdf")
  corplot_names <- c("COVID-19 Severe\nvs Moderate","\nAbatacept in RA")
  
} else if (analysis == "early_cohort"){
  # Figure 2 datasets
  COVID <- readRDS("GSEA/Xiong_COVID19_CtrlDisease.RDS")
  Abatacept <- readRDS("GSEA/RA_Abatacept_GranulocytesAdj.RDS")  
  
  tables_path <- c("Tables/SupplementaryTable3.csv","Tables/SupplementaryTable4.csv","Tables/SupplementaryTable5.csv")
  figures_path <- c("Figures/figure_2A.jpeg","Figures/figure_2B.pdf")
  corplot_names <- c("COVID-19 vs CTRL","Abatacept in RA")
}

# Remove redundant BPs (keep the cluster-BP most strongly associated to COVID19 as representative)
COVID <- merge(bps_cluster[,2:4],COVID[,4:11],by=0)
Abatacept <- merge(bps_cluster[,2:4],Abatacept[,4:11],by=0)
names(COVID)[1] <- names(Abatacept)[1] <- "ID"
row.names(COVID) <- COVID$ID
row.names(Abatacept) <- Abatacept$ID
 
COVID <- COVID[order(COVID$pvalue), ]
Abatacept <- Abatacept[order(Abatacept$pvalue), ]
COVID <- COVID[COVID$ID %in% Abatacept$ID, ]
COVID <- COVID[!duplicated(COVID$Cluster), ]
Abatacept <- Abatacept[Abatacept$ID %in% COVID$ID, ]
n <- nrow(COVID)

# Significant and overlapping BPs
COVID <- COVID[COVID$p.adjust < 0.05, ]
Abatacept <- Abatacept[Abatacept$p.adjust < 0.05, ]
overlap <- merge(COVID[,c(1:4,6:8,10,11)],Abatacept[,c(6:8,10,11)],by=0,suffixes=c("_COVID","_Abatacept"))

# Write SUPPLEMENTARY TABLES 3 to 8
write.csv(file = tables_path[1], row.names = F,
          COVID[,c("ID","Description","Size","NES","pvalue","p.adjust","ratio")])

write.csv(file = tables_path[2], row.names = F,
          Abatacept[,c("ID","Description","Size","NES","pvalue","p.adjust","ratio")])

write.csv(file = tables_path[3], row.names = F,
          overlap[,c("ID","Description","Size","NES_COVID","pvalue_COVID","p.adjust_COVID","ratio_COVID",
                     "NES_Abatacept","pvalue_Abatacept","p.adjust_Abatacept","ratio_Abatacept")])

# Calculate antagonism
cov_pos_coherent <- sum(overlap$NES_COVID > 0 & overlap$NES_Abatacept < 0)
cov_neg_coherent <- sum(overlap$NES_COVID < 0 & overlap$NES_Abatacept > 0)
cov_pos_incoherent <- sum(overlap$NES_COVID > 0 & overlap$NES_Abatacept > 0)
cov_neg_incoherent <- sum(overlap$NES_COVID < 0 & overlap$NES_Abatacept < 0)

# Calculate antagonism probability
p0 <- (nrow(COVID[COVID$NES > 0, ])/n)*(nrow(Abatacept[Abatacept$NES < 0, ])/n) + (nrow(COVID[COVID$NES < 0, ])/n)*(nrow(Abatacept[Abatacept$NES > 0, ])/n)

p <- binom.test(x = cov_pos_coherent + cov_neg_coherent, n = n , p = p0, alternative = "greater")$p.value

print(paste0("n coherent = ",cov_pos_coherent + cov_neg_coherent, "; n incoherent = ", 
             cov_pos_incoherent + cov_neg_incoherent,"; analyzed BPs = ", n, "; binomial test probability = ", p))

# FIGURE 1A and 2A: Euler diagram
cov_pos <- sum(COVID$NES > 0) - cov_pos_coherent - cov_pos_incoherent
cov_neg <- sum(COVID$NES < 0) - cov_neg_coherent - cov_neg_incoherent
abatacept_pos <- sum(Abatacept$NES > 0) - cov_neg_coherent - cov_pos_incoherent
abatacept_neg <- sum(Abatacept$NES < 0) - cov_pos_coherent - cov_neg_incoherent

if (analysis == "late_cohort"){
  fit1 <- euler(c("COVID-19\nSEVERE UP\n" = cov_pos, "COVID-19\nSEVERE DOWN\n" = cov_neg,
                  "\n\n\n\n\n\n\n\nABATACEPT\nUP" = abatacept_pos, "ABATACEPT\nDOWN\n" = abatacept_neg, 
                  "COVID-19\nSEVERE UP\n&ABATACEPT\nDOWN\n" = cov_pos_coherent,
                  "COVID-19\nSEVERE UP\n&\n\n\n\n\n\n\n\nABATACEPT\nUP" = cov_pos_incoherent,
                  "COVID-19\nSEVERE DOWN\n&\n\n\n\n\n\n\n\nABATACEPT\nUP" = cov_neg_coherent,
                  "COVID-19\nSEVERE DOWN\n&ABATACEPT\nDOWN\n" = cov_neg_incoherent)
  )
  
} else if (analysis == "early_cohort"){
  fit1 <- euler(c("COVID-19 UP" = cov_pos, "COVID-19 DOWN" = cov_neg,
                  "ABATACEPT UP" = abatacept_pos, "ABATACEPT DOWN" = abatacept_neg, 
                  "COVID-19 UP&ABATACEPT DOWN" = cov_pos_coherent,
                  "COVID-19 UP&ABATACEPT UP" = cov_pos_incoherent,
                  "COVID-19 DOWN&ABATACEPT UP" = cov_neg_coherent,
                  "COVID-19 DOWN&ABATACEPT DOWN" = cov_neg_incoherent)
  )  
}

jpeg(figures_path[1], width=1350, height=900)
  plot(fit1, shape="ellipse",
       labels = list(fontsize = 25),
       edges ="azure4",
       quantities = list(fontsize = 30))
dev.off()


# FIGURE 1B and 2B: Most antagonistic BPs
# Determine most different BPs by combining the P-value of both exposures (sum of the logs Fisher's method)
overlap$sum_pval <- apply(overlap[,c("pvalue_COVID","pvalue_Abatacept")], MARGIN = 1, function(x){y<-sumlog(x);y$p})
overlap <- overlap[order(overlap$sum_pval), ]
overlap <- overlap[1:20, ]

pvals <- overlap[,c("pvalue_COVID","pvalue_Abatacept")]
nes <- overlap[,c("NES_COVID","NES_Abatacept")]

row.names(pvals) <- row.names(nes) <- overlap$Description
colnames(nes) <- colnames(pvals) <- corplot_names

pdf(figures_path[2])
  corrplot::corrplot(as.matrix(nes), is.corr=F, p.mat=as.matrix(pvals), 
                     sig.level=c(0.00001,0.0001,0.001,0.01,0.05), tl.cex = 0.8, pch.cex=0.7, 
                     insig="label_sig", tl.col = "black", cl.lim = c(-4,4), cl.pos = "r",
                     cl.align.text="l", cl.ratio = 1)
dev.off()


# DEG overlap: Supplementary Table 9 =========================
COVIDlate <- readRDS("/media/IMID/Projects/COVID19/Github_repo/COVID-19_Abatacept_good/DEG/Overmeyer_COVID19_SevereModerate.RDS")
Abatacept <- readRDS("/media/IMID/Projects/COVID19/Github_repo/COVID-19_Abatacept_good/DEG/RA_Abatacept.RDS")
genes <- row.names(Abatacept[Abatacept$FDR < 0.05, ])
dt <- COVIDlate[row.names(COVIDlate) %in% genes, ]
dt <- cbind(Abatacept=Abatacept[row.names(dt), ], Covid=dt)
dt <- dt[order(dt$Covid.PValue), ]
write.csv(dt, file = "Tables/SupplementaryTable9.csv", row.names = F)








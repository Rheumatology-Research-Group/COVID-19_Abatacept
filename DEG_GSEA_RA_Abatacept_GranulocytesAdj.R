library(edgeR)
library(limma)
library(org.Hs.eg.db)
library(GO.db)
library(fgsea)

# DIFFERENTIAL GENE EXPRESSION ANALYSIS ==================================================
load("Source_data/RA_Abatacept_cts.RData")
sampleTable <- sampleTable[!is.na(sampleTable$Granulocytes), ]
cts <- cts[ , row.names(sampleTable)]
cts <- DGEList(counts=cts,group=sampleTable$week)

# Filter genes with low counts
cts <- cts[rowSums(cpm(cts) > 0.6) >= 20,]
dim(cts)
cts$samples$lib.size <- colSums(cts$counts)

# TMM normalization
cts <- calcNormFactors(cts, method="TMM")

#Model design
attach(sampleTable)
design <- model.matrix(~ week + seq.plate + Sex + age + Granulocytes)
detach(sampleTable)

# Transformation to logCPM with voom
# Calculate within sample correlation for blocking, running corfit twice
v <- voom(cts,design,plot=F)
corfit <- duplicateCorrelation(v,design,block=sampleTable$Code)
v <- voom(cts, design, block=sampleTable$Code, correlation=corfit$consensus.correlation,plot=T)
corfit <- duplicateCorrelation(v, design, block=sampleTable$Code)

# Fit linear model
fit <- lmFit(v, design, block=sampleTable$Code, correlation=corfit$consensus.correlation)
fit <- eBayes(fit, trend=TRUE)

# Differentially expressed genes between w12 and baseline
res <- topTable(fit, coef="weekw12", number=nrow(cts),adjust.method="BH",sort.by="P")
colnames(res) <- c('logFC','AveExpr','t','PValue','FDR','B')
saveRDS(res, file=paste0("DEG/RA_Abatacept_GranulocytesAdj_DEG.RDS"))

# GENE SET ENRICHMENT ANALYSIS ========================================================================
# GO BPs to test ----------
dictGO <- as.list(org.Hs.egGO2ALLEGS)
dictGO <- lapply(dictGO, function(x) {unique(x)})
dictGO <- dictGO[sapply(names(dictGO),function(x) length(dictGO[[x]])>9&length(dictGO[[x]])<301)]

dictTERMS <- as.list(GOTERM)
dictGO <- dictGO[sapply(names(dictGO), function(x){Ontology(dictTERMS[[x]])=='BP'})]
GO_description <- sapply(names(dictGO), function(x) Term(dictTERMS[[x]]))
rm(dictTERMS)

# DEG results setup --------
#Translate symbol to entrez and exclude genes without translation
genes_dict <- readRDS("Source_data/symbol2entrez.RDS")
genes_dict <- genes_dict[!is.na(genes_dict$ENTREZ_ID), ]
rownames(genes_dict) <- genes_dict$SYMBOL

res$symbols <- rownames(res)
res <- res[res$symbols %in% genes_dict$SYMBOL, ]
res$entrez <- genes_dict[res$symbols,"ENTREZ_ID"]

# Rank genes based on sign(logFC) * -log10(PValue)
res$sLogP <- -log10(res$PValue) * sign(res$logFC) 
vals2test <- res$sLogP
names(vals2test) <- res$entrez
vals2test<-vals2test[order(vals2test,decreasing = T)]

# GSEA ---------------------
edo<-fgsea(pathways=dictGO, stats=vals2test, nperm=1000000, minSize=10, maxSize=300, nproc=3)

# Calculate enrichment ratio: Number of genes in core enrichment/set size
edo$counts <- sapply(edo$leadingEdge, function(x) length(x), USE.NAMES=F) 
edo$ratio <- edo$counts/edo$size 
edo <- edo[order(edo$pval),]

# Format
edo$Description <- GO_description[edo$pathway]
edo$leadingEdge <- sapply(edo$leadingEdge, function(x) paste(x, collapse='/'))
names(edo) <- c('ID','pvalue','p.adjust','enrichmentScore','NES','nMoreExtreme','setSize','core_enrichment','counts','ratio','Description')
edo <- edo[,c('ID','Description','setSize','enrichmentScore','NES','pvalue','p.adjust','core_enrichment','counts','ratio','nMoreExtreme')]
edo <- as.data.frame(edo)
row.names(edo) <- edo$ID

saveRDS(edo, file=paste0("GSEA/RA_Abatacept_GranulocytesAdj_GSEA.RDS"))
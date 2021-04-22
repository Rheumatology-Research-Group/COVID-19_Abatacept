library(edgeR)
library(limma)
library(org.Hs.eg.db)
library(GO.db)
library(fgsea)
library(ggplot2)
library(shiny)

# DIFFERENTIAL GENE EXPRESSION ANALYSIS ==================================================
# Read data ---------
gene_info <- read.table("Source_data/hg38_gencode.v32.info.tsv", sep="\t", header=T, stringsAsFactors=F)[, c("Gid", "GeneName")]
gene_info <- gene_info[!duplicated(gene_info), ]
names(gene_info) <- c("Gid", "Name")
row.names(gene_info) <- gene_info$Gid

raw_cts <- read.table("Source_data/Xiong_COVID19early_cts.tsv", sep="\t", header=T, stringsAsFactors=F)
rownames(raw_cts) <- raw_cts$Gid
raw_cts <- as.matrix(raw_cts[gene_info$Gid, -1])
sampleTable <- data.frame(ID=colnames(raw_cts), condition=factor(c(rep("N", 3), rep("P", 3)), levels = c("N", "P")), stringsAsFactors = F)

# Remove contaminating globin transcripts
raw_cts <- raw_cts[!gene_info$Name %in% c("HBA1","HBA2","HBB"), ]
gene_info <- gene_info[!gene_info$Name %in% c("HBA1","HBA2","HBB"), ]

# Determine cell composition (ABIS) --------
# Transform to TPM
gene_length <- readRDS("Source_data/gencode.v32.genelength_as_mean_transcript_length.RDS")
row.names(gene_length) <- gene_length$Gid
table(gene_info$Gid %in% gene_length$Gid)
gene_length <- gene_length[gene_info$Gid, ]
gene_info <- gene_length

tpm <- raw_cts/gene_info$Gene.Length # Gene length correction
tpm <- t(t(tpm)*1e6/colSums(tpm)) # Library size correction

# Summarize at symbol level (mean)
tpm <- aggregate(tpm, list(gene_info$Name), mean, na.rm=T)
row.names(tpm) <- tpm$Group.1
tpm <- tpm[,-1]
write.table(tpm,"Source_data/temp_COVID19_TPM_cts.txt",quote=T,sep="\t")

# Calculate cell composition using ABIS
runGitHub("ABIS", user="giannimonaco")
cells <- read.csv("Source_data/temp_COVID19_cellprop.txt",sep="\t",stringsAsFactors = F)
colSums(cells)
sampleTable$Neutrophils <- unlist(cells["Neutrophils LD", ])
system("rm Source_data/temp_COVID19*.txt")

# FIGURE: Neutrophil contamination
cells <- data.frame(Sample = colnames(cells), Neutrophils = t(cells["Neutrophils LD",]),stringsAsFactors = F)
names(cells)[2] <- "Neutrophils"
cells$Sample <- factor(cells$Sample, levels=cells$Sample)
cells$Group <- factor(rep(c("Healthy","COVID-19"),each=3),levels=c("Healthy","COVID-19"))

pdf("Figures/SupplementaryFigure5.pdf")
  ggplot(data=cells, aes(x = Sample, y = Neutrophils, fill = Group)) +
    geom_bar(stat="identity", color="black") +
    scale_fill_manual(values=c("#999999", "#56B4E9")) +
    theme_minimal(base_size = 20) +
    ylab("Estimated neutrophils") +
    theme(axis.title.x = element_blank())
dev.off()

# Differential gene expression analysis ------------
cts <- DGEList(counts=raw_cts,group=sampleTable$condition)

# Filter genes with low counts
cts <- cts[rowSums(cpm(cts)>2) >= 2,]
dim(cts)
gene_info <- gene_info[gene_info$Gid %in% row.names(cts), ]
cts$samples$lib.size <- colSums(cts$counts) 

# TMM normalization
cts <- calcNormFactors(cts, method="TMM")

# Model design
design <- model.matrix(formula(~Neutrophils+condition), sampleTable)

# Fit glmQL model
y <- estimateDisp(cts, design)
fit <- glmQLFit(y, design, robust=T)

# Differentially expressed genes between patients and controls
fit <- glmQLFTest(fit, coef = "conditionP") 
res <- topTags(fit, n = dim(fit$table)[1], adjust.method = "BH", sort.by = "PValue")@.Data[[1]]
res$symbols <- gene_info[row.names(res),"Name"]
res <- res[!duplicated(res$symbols), ]
saveRDS(res, file=paste0("DEG/Xiong_COVID19_CtrlDisease.RDS"))


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

res <- res[res$symbols %in% genes_dict$SYMBOL, ]
res$entrez <- genes_dict[res$symbols,"ENTREZ_ID"]

# Rank genes based on sign(logFC) * -log10(PValue)
res$sLogP <- -log10(res$PValue) * sign(res$logFC) 
vals2test <- res$sLogP
names(vals2test) <- res$entrez
vals2test <- vals2test[order(vals2test,decreasing = T)]

# GSEA ---------------------
edo <- fgsea(pathways=dictGO, stats=vals2test, nperm=1000000, minSize=10, maxSize=300, nproc=3)

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

saveRDS(edo, file=paste0("GSEA/Xiong_COVID19_CtrlDisease.RDS"))


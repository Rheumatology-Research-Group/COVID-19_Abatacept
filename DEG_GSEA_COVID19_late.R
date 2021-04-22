library(edgeR)
library(limma)

# DIFFERENTIAL GENE EXPRESSION ANALYSIS ==================================================
# Prepare phenotypes table -----------------
###get Pheno
library(GEOquery)
GSE<-"GSE157103"
set1<- getGEO(GSE,GSEMatrix=T,destdir="./",AnnotGPL=F)
set1<-set1[[1]]

####clean dataframe
pheno<-pData(set1)
pheno$Group<-gsub("disease state: ","",pheno$characteristics_ch1)
pheno<-pheno[,63:84]
pheno[pheno=="unknown"] <- NA
pheno[pheno==":"] <- NA
pheno[pheno==">89"] <- 90
pheno<-pheno[,c(1,12,13,15,18,22)]

######set group based on HFD-45
pheno$days <- as.numeric(pheno$`hospital-free days post 45 day followup (days):ch1`)
hist(pheno$days)

######More than 25 days is less severe
pheno$Severity<-ifelse(pheno$days<25,"HIGH","LOW")

# Prepare expression matrix -----------------
###load RSEM estimated Gene cnts
cnts<-read.delim("GSE157103_genes.ec.tsv",head=T)
rownames(cnts)<-cnts$X.symbol
cnts<-cnts[,-1]

###filter out non-expressed genes
y<-cnts
keep <- rowSums(cpm(y)>1) >= min(table(pheno$Group))#### change group depending replicates
y<- y[keep, ]
y <- DGEList(y,group=pheno$Group)
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y, method="TMM")

##exclude incomplete cases
pheno$title<-colnames(y)
pheno<-pheno[complete.cases(pheno),]

####subset to only COVID samples
pheno<-pheno[pheno$Group=="COVID-19",]
y<-y[,colnames(y) %in% pheno$title]

# Differential expression analysis -----------------
group<-as.factor(pheno$Severity)
group<-relevel(group,"LOW")##Health always below!

design <- model.matrix(~group+as.numeric(pheno$`age (years):ch1`)+
                         as.factor(pheno$`Sex:ch1`)
)
colnames(design)[2]<-"Comp"

####use limma voom for speed and 
###low number of samples
v <- voom(y, design=design, plot=TRUE)
fit<-lmFit(v,design)
fit<-eBayes(fit,robust=TRUE)
results <- topTable(fit,coef="Comp",num=dim(fit)[1],sort.by="P",adjust.method = "bonferroni") 
colnames(results)[4:5] <- c("PValue","FDR")
saveRDS(results,file="DEG/Overmeyer_COVID19_SevereModerate.RDS")

# GENE SET ENRICHMENT ANALYSIS ==================================================
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

results$symbols <- row.names(results)
results <- results[results$symbols %in% genes_dict$SYMBOL, ]
results$entrez <- genes_dict[results$symbols,"ENTREZ_ID"]

# Rank genes based on sign(logFC) * -log10(PValue)
results$sLogP <- -log10(results$PValue) * sign(results$logFC) 
vals2test <- results$sLogP
names(vals2test) <- results$entrez
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

saveRDS(edo, file=paste0("GSEA/Overmeyer_COVID19_SevereModerate.RDS"))



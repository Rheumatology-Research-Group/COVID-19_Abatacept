library(edgeR)

# Note: set path to working directory

#### Abatacept differential expression ####
DiffExp_all <- readRDS("DEG/RA_Abatacept.RDS")

# Get signature genes (FDR significant)
Resp_genes_all <- row.names(DiffExp_all[DiffExp_all$FDR < 0.05, ])

#### COVID dataset ####
# Phenodata
pheno <- readRDS("Source_data/GSE157103_COVID19late_pheno.rds")
pheno$Group <- as.factor(gsub("disease state: ", "", pheno$characteristics_ch1))
pheno$Group <- relevel(pheno$Group,"non-COVID-19")
pheno <- pheno[, 63:84]

# Load COVID dataset count matrix
cnts <- read.delim("Source_data/GSE157103_COVID19late_RSEMcts.tsv",head=T)

rownames(cnts) <- cnts$X.symbol
cnts <- cnts[,-1] 
colnames(cnts) <- rownames(pheno)
y <- cnts

# Filter expressed genes, based on minimum group frequency
keep <- rowSums(cpm(y)>1) >= min(table(pheno$Group)) 
y <- y[keep, ]
y <- DGEList(y, group = pheno$Group)
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y, method="TMM")
cpm_log <- cpm(y,log = T, prior.count = 2)


#### Compute Score ####
# Abatacept signature genes
AGenes <- Resp_genes_all[Resp_genes_all %in% rownames(cpm_log)]
saveRDS(AGenes, file = "Abatacept_Score/AbataceptScore_genes.rds")

# Standardize the covid dataset expression matrix
cpm_log <- scale(t(cpm_log[AGenes,]))

# Compute score based on the linear combination given by the logFC of the Abatacept results
scores <- data.frame(score = cpm_log %*% DiffExp_all[AGenes, "logFC"])
saveRDS(scores, file = "Abatacept_Score/AbataceptScore.rds")

# load the required libraries
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
library(curatedTCGAData)
library(pheatmap)
library(SummarizedExperiment)
library(TCGAbiolinks)

# get a list of projects
gdcprojects <- getGDCprojects()
gdcprojects$disease_type
summarypro <- getProjectSummary('TCGA-SKCM')

#tcga_ids <- read.csv("/home/admin/ids.csv",stringsAsFactors = FALSE)

######## 1 make an TCGA object #######
# build a query to retrieve gene expression data ------------
query_TCGA <- GDCquery(project = 'TCGA-SKCM',data.category = 'Transcriptome Profiling',experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts', access = 'open')
GDCquery()

# download data - GDCdownload
####### 2 download the object ########
GDCdownload(query_TCGA, directory = "/home/balqees/Melanoma-Analysis/")

####### 3 prepare data #########
tcga_skin_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE,directory = "/home/balqees/Melanoma-Analysis/")
# Assuming your data frame is named 'phdata'
tcga_skin_data <- tcga_skin_data[,tcga_skin_data@colData@listData[["definition"]] != "Solid Tissue Normal"]
table(tcga_skin_data@colData@listData[["definition"]])
######## 4 get the assay ########

table(tcga_skin_data@colData@listData[["primary_diagnosis"]]) # subset Malignant melanoma, NOS 422
table(tcga_skin_data@colData@listData[["tumor_descriptor"]]) # Exclude Not Applicable
table(tcga_skin_data@colData@listData[["ajcc_pathologic_stage"]]) # exclude Not Reported 14
table(tcga_skin_data@colData@listData[["gender"]])
table(tcga_skin_data@colData@listData[["race"]]) # subset the white 450
table(tcga_skin_data@colData@listData[["site_of_resection_or_biopsy"]]) # for batch correction
#================================================================================
# Subset Malignant melanoma, NOS 422
malignant_melanoma <- tcga_skin_data[tcga_skin_data@colData@listData[["primary_diagnosis"]] == "Malignant melanoma, NOS"]

# Exclude Not Applicable from tumor_descriptor
applicable_tumors <- malignant_melanoma[malignant_melanoma@colData@listData[["tumor_descriptor"]] != "Not Applicable"]

# Exclude Not Reported and NA from ajcc_pathologic_stage
reported_stage <- applicable_tumors[!is.na(applicable_tumors@colData@listData[["ajcc_pathologic_stage"]]) & 
                                      applicable_tumors@colData@listData[["ajcc_pathologic_stage"]] != "Not Reported"]


skin_matrix <- assay(tcga_skin_data, 'stranded_first')

skin_matrix <- as.data.frame(skin_matrix)

######## 5 convert Ensemble Ids to gene symbols ########

# Install and load org.Hs.eg.db package
library(org.Hs.eg.db)

# Convert Ensembl IDs to gene symbols
ensembl_ids <- rownames(skin_matrix)

# Remove version numbers for compatibility
ensembl_ids <- sub("\\..*", "", ensembl_ids)

# Map IDs
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Remove NA values and keep only mapped genes
mapped_genes <- !is.na(gene_symbols)
skin_matrix <- skin_matrix[mapped_genes, ]
gene_symbols <- gene_symbols[mapped_genes]

unique_gene_symbols <- make.unique(gene_symbols, sep = "_")
length(unique_gene_symbols)

# Assign gene symbols to prostate_matrix
rownames(skin_matrix) <- unique_gene_symbols



############### Deseq2 Normalizations ####################

# my phenodata 
phdata <- colnames(skin_matrix)
phdata <- as.data.frame(phdata)

# Assuming your data frame is called 'df'
phdata$disease_state <- tcga_skin_data@colData@listData[["definition"]]
# Assuming your data frame is named 'df' and the column is named 'disease_state'
phdata$disease_state <- gsub("Primary solid Tumor", "Primary", phdata$disease_state)
phdata$disease_state <- gsub("Additional Metastatic", "Metastatic", phdata$disease_state)
table(phdata$disease_state)
names(phdata)[1] <- "ids"
rownames(phdata) <- phdata$ids

# making the rownames and column names identical
all(rownames(phdata) %in% colnames(skin_matrix))
all(rownames(phdata) == colnames(skin_matrix))
library(DESeq2)
# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = round(skin_matrix),
                              colData = phdata,
                              design = ~ disease_state)

# Pre-filtering
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Set factor levels
dds$disease_state <- relevel(dds$disease_state, ref = "Primary")
levels(dds$disease_state)
# Run DESeq
dds <- DESeq(dds)
# Assuming 'res' is the result object from DESeq2
res <- results(dds)

# Subset for p-value < 0.05 and |log2FoldChange| > 2
filtered_res <- res[which(res$pvalue < 0.05 & abs(res$log2FoldChange) > 2), ]
summary(filtered_res)



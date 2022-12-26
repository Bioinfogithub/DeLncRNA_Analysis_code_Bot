#install DESeq2 using Bioconductor
if (!require('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

#BiocManager::install('DESeq2')

#Install EnhancedVolcano using Bioconductor
if (!require('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

#BiocManager::install('EnhancedVolcano')

#converting the gene id into gene symbol and gene name
# library(org.Hs.eg.db)# loading homosapiens org db library
# setwd("D:\\Mtech project data\\GSE63979_Read_Counts") 
# input <- read.csv("GSE63979_ReadCount.csv",TRUE,",") #read input file
# 
# annots <- select(org.Hs.eg.db, keys=rownames(input),
#                  columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")# to getting symbol and genename
# resultTable <- merge(input, annots, by.x=0, by.y="ENTREZID")#entrezid
# head(resultTable)# results
# 
# write.csv(resultTable, file = "GSE63979_Read_Counts1.csv") #write the list of gene_symbols and genename to a CSV file


library(DESeq2)
library(EnhancedVolcano)
#Load the file containing the experiment design
#'path'=location of the file
 
#Load raw counts data
#'path'=location of the file
setwd("D:\\Mtech project data\\GSE63979_Read_Count")
data <- read.csv("GSE63979_ReadCount.csv", sep = ",", header = TRUE)

#View(design)
#Now lets view the raw counts dataset and make sure the counts data is there.
View(data)
summary(data)
dim(data)
data<- na.omit(data)

# load the package
#library(dplyr)
 
 

#Now that all the data is confirmed, validated and examined, the data wrangling part can be started.
#This file is important because it contatin columns relating samples with phenotypes
#design=design[1:37,]

#converting first column as a row in original dataset so that ncol(deg_data) == nrow(meta)
data<- data.frame(data[,-1], row.names = data[,1])
 

View(data)
dim(data)

#Removing rows having all zeros ???
data<- data[rowSums(data[])>0,]
head(data)



#step 1: preparing the count data ......
#read in count data as a deg_data
# deg_data<- read.csv('FPKM data.csv', header = TRUE, sep = ",")
# head(deg_data)
# dim(deg_data)

condition <- factor(c(rep("healthy",8),rep("psoriatic",7)))

#converting first column as a row in original dataset so that ncol(deg_data) == nrow(meta)
#deg_data<- data.frame(deg_data[,-1], row.names = deg_data[,1])
library(DESeq2)
#read in samples info....
meta <- data.frame(condition)

ncol(data) == nrow(meta)

library(DESeq2)
#Load the genetable
#Select the variables of interest data and add gene names
#dataset=raw.counts[,c('Gene.Name',design$Run)]
#genetable=raw.counts[,design$Run]

#step2: Construct DESEQDataSet Object......
dds <- DESeqDataSetFromMatrix(
  countData = data,
  design = ~ condition,
  colData = meta)

data = ceiling(data) #converting some values in assay that are not integers
normCounts<-rlog(dds,blind = FALSE)#regularized log' transformation

#Now we're ready to run DESEQ function
dds <- DESeq(dds)

#Create DESeq results object using Benjamini-Hochberg correction
# res=results(object = dds, contrast = c('type','cancer','normal'),
#             pAdjustMethod = 'BH',alpha = 0.000001)
# row.names(res)=dataset$Gene.Name
# summary(res)

#Create DESeq results object using Holmr correction
# res=results(object = dds, contrast = c('type','cancer','normal'),
#             pAdjustMethod = 'holm', alpha = 0.000001)
# row.names(res)=dataset$Gene.Name
res <- results(dds)
head(results(dds, tidy=TRUE)) #let's look at the results table

#Summary of differential gene expression
summary(res) #summary of results

#Sort summary list by p-value
res <- res[order(res$padj),]
head(res)
#Create the final dataframe consisting of ordered deseq results based on log2fc
resord=as.data.frame(res)
finaltable1=resord[order(resord$padj),]
write.table(finaltable1, file = 'deseq2_result.csv', sep = ',',
            col.names = NA)



#Plot the most basic volcano plot of deseq2 result
data1<- read.csv("deseq2_result.csv", sep = ",", header = TRUE)
EnhancedVolcano(data1,
                lab = data1$GENE_SYMBOL,
                x = 'log2FoldChange',
                y = 'pvalue')


#Create publication grade volcano plot with marked genes of interest
EnhancedVolcano(data1,
                lab = data1$GENE_SYMBOL,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 10e-5,
                FCcutoff = 1.333,
                xlim = c(-5.7, 5.7),
                ylim = c(0, -log10(10.2e-12)),
                pointSize = 1.3,
                labSize = 2.6,
                title = 'The results',
                subtitle = 'Differential expression analysis',
                caption = 'log2fc cutoff=1.333; p value cutof=10e-5',
                legendPosition = "right",
                legendLabSize = 14,
                col = c('lightblue', 'orange', 'blue', 'red2'),
                colAlpha = 0.6,
                drawConnectors = TRUE,
                hline = c(10e-8),
                widthConnectors = 0.5)

#Create publication grade volcanoplot with marked genes of interest
EnhancedVolcano(data1,
                lab = data1$GENE_SYMBOL,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-7,
                FCcutoff = 2.5,
                xlim = c(-5.7, 5.7),
                ylim = c(0, -log10(10.2e-12)),
                pointSize = 1.3,
                labSize = 2.6,
                title = 'The results',
                subtitle = 'Differential expression analysis',
                caption = 'log2fc cutoff=1.333; p value cutof=10e-6',
                legendPosition = "right",
                legendLabSize = 14,
                col = c('lightblue', 'orange', 'blue', 'red2'),
                colAlpha = 0.6,
                drawConnectors = TRUE,
                hline = c(10e-8),
                widthConnectors = 0.5)

## DESeq2 ##
library(DESeq2)
#loading libraries,edge R,limma,Biobase,deseq2
library(BiocManager)
library(edgeR)       
library(limma)
library(Biobase)
library(DESeq2)
#DE with edge R
dge <- DGEList(counts = data, group=meta$condition)
head(dge)

# Normalize by total count
dge <- calcNormFactors(dge)

# Create the contrast matrix
design.mat <- model.matrix(~ 0 + dge$samples$group)
colnames(design.mat) <- levels(dge$samples$group)

# Estimate dispersion parameter for GLM
#Generalized Linear Model(GLM) for differential expression analysis in RNA-Seq data, and in the model two indicator 
#variables Group and Module are adopted to fit the GLM
dge <- estimateGLMCommonDisp(dge, design.mat)
dge <- estimateGLMTrendedDisp(dge, design.mat, method="power")
dge<- estimateGLMTagwiseDisp(dge,design.mat)

# Plot mean-variance
plotBCV(dge)#Value of variation of biological variation (BCV)

# limma-voom
# Create design matrix
design <- model.matrix(~meta$condition)

# Apply voom transformation
nf <- calcNormFactors(data)
v <- voom((data), design, lib.size=colSums(data)*nf, 
          normalize.method="quantile", plot=FALSE)
##DEG analysis using different tools 
# set the threshold
p.threshold <- 0.05

## edgeR ##
# Design matrix
design.mat <- model.matrix(~ 0 + dge$samples$group)
colnames(design.mat) <- c("healthy", "psoriatic")

# Model fitting
fit.edgeR <- glmFit(dge, design.mat)
# Differential expression
contrasts.edgeR <- makeContrasts(healthy - psoriatic, levels=design.mat)
lrt.edgeR <- glmLRT(fit.edgeR, contrast=contrasts.edgeR)
names(fit.edgeR)

# Access results tables
edgeR_results <- lrt.edgeR$table
View(edgeR_results)
sig.edgeR <- decideTestsDGE(lrt.edgeR, adjust.method="BH", p.value =0.05)
genes.edgeR <- row.names(edgeR_results)[which(sig.edgeR != 0)]#inequality operator(!=)
View(sig.edgeR)
head(genes.edgeR)
summary(edgeR_results)

############
## DESeq2 ##
library(DESeq2)
contrast.deseq2 <- as.list(resultsNames(dds)[2:3])
contrast.deseq2 <- list(c("condition_psoriatic_vs_healthy")) 
# contrast.deseq2 <- list("strainC57BL.6J", "strainDBA.2J")
deseq2_results <- results(dds, contrast=contrast.deseq2)
deseq2_results
summary(deseq2_results)

deseq2_results$threshold <- as.logical(deseq2_results$padj < 0.05)
genes.deseq <- row.names(deseq2_results)[which(deseq2_results$threshold)]
length(genes.deseq)
summary(deseq2_results, alpha = 0.05)
View(deseq2_results)
View(genes.deseq)

# Usual limma pipeline
fit.voom <- lmFit(v, design)
fit.voom <- eBayes(fit.voom)
voom_results <- topTable(fit.voom, coef=2,  adjust="BH", number = nrow(data))
View(voom_results)
voom_results$threshold <- as.logical(voom_results$adj.P.Val < 0.05)
genes.voom <- row.names(voom_results)[which(voom_results$threshold)]

#how many genes are significant? 
length(genes.voom)
length(genes.deseq)
length(genes.edgeR)

### Overlapping genes (We will use a venn diagram to depict the overlaps). 
#install.packages('gplots')
library(gplots)
venn(list(edgeR = genes.edgeR, DESeq2 = genes.deseq, voom = genes.voom))

#write.table(edgeR_results, file = "edgeR_results.txt", sep = "\t" ,col.names = NA, quote = F)
write.csv(edgeR_results,file = "edgeR_results.csv",
          quote = F)
write.table(deseq2_results, file = "deseq2_results.csv", sep = "\t" ,
            col.names = NA, quote = F)
write.table(voom_results, file = "voom_results.csv", sep = "\t" ,
            col.names = NA, quote = F)
new_results <- results(dds, contrast = contrast.deseq2,alpha = 0.05)



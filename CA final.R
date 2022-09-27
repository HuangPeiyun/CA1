## CA1 from Peiyun Huang and Xuerui Zhong
library(GEOquery)
library(limma)
library(tidyverse)
library(oligo)
library(oligoClasses)
library(hgu133plus2.db)
library(org.Hs.eg.db)

## Retrieving the whole GSE as a list of ExpressionSets
gse50737 <- getGEO('GSE50737')
gse50737 <- gse50737[[1]]
length(gse50737)
class(gse50737)
varLabels(gse50737)

## Importing raw (unprocessed) Affymetrix microarray data
gse50737_celdata <- read.celfiles(list.celfiles('/Users/apple/Desktop/BL5631/GSE50737_RAW',full.names=TRUE,listGzipped=TRUE))
image(gse50737_celdata[,1])

## Appending the supplementary_file and re-reading out CEL files
library(dplyr)
library(purrr)
library(stringr)
pd <- pData(gse50737)
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)
gse50737_celdata <- read.celfiles(paste0('/Users/apple/Desktop/BL5631/GSE50737_RAW/',pd$cel_file),phenoData=phenoData(gse50737))
pData(gse50737_celdata)[,c("data_row_count","tissue:ch1","treatment:ch1")]

## Normalisation
gse50737_eset <- rma(gse50737_celdata)

## To identify differentially expressed genes using linear models
~ treat
gse50737<- getGEO('gse50737',GSEMatrix=TRUE,AnnotGPL=TRUE,destdir="./")
gse50737 <- gse50737[[1]]
pData(gse50737)[['treatment:ch1']]

# Generate a model matrix specifying the design
design <- model.matrix( ~ gse50737_eset[['treatment:ch1']])
design

design <- model.matrix( ~ 0 + gse50737_eset [['treatment:ch1']])
design

colnames(design) <- c("benzene_exposed","benzene_poisoning","control")
design

# Fit the model to the design
fit <- lmFit(gse50737_eset,design)

# Empirical Bayes correction in limma
fitted.ebayes <- eBayes(fit)

# Extracting differentially expressed genes
topTable(fitted.ebayes)
summary(decideTests(fitted.ebayes[,"benzene_exposed"],lfc=1))
summary(decideTests(fitted.ebayes[,"benzene_poisoning"],lfc=1))
summary(decideTests(fitted.ebayes[,"control"],lfc=1))

## Remove the intercept and use contrasts
design <- model.matrix( ~ 0 + gse50737_eset[['treatment:ch1']])
colnames(design) <- c("benzene_exposed","benzene_poisoning","control")
contrast_matrix <- makeContrasts(benzene_exposed-control, benzene_exposed-benzene_poisoning, benzene_poisoning-control,levels=design)
contrast_matrix
fit2 <- contrasts.fit(fit,contrasts=contrast_matrix)
fit2 <- eBayes(fit2)
topTable(fit2)
summary(decideTests(fit2,lfc=1))

## Annotation of genomic data in AnnotationDb
ps <- rownames(topTable(fit2))
library(hugene20sttranscriptcluster.db)
ls('package:hugene20sttranscriptcluster.db')
columns(hugene20sttranscriptcluster.db)
keytypes(hugene20sttranscriptcluster.db)
head(keys(hugene20sttranscriptcluster.db,keytype="PROBEID"))
AnnotationDbi::select(hugene20sttranscriptcluster.db,ps,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
ps2 <- topTable(fit2,number=Inf,p.value = 0.05,lfc=1)
## Since there is no logFC in our dataset, so the following lines actually would not work properly for our data.
## We just wrote them down for your information.
# ps2_up <- rownames(ps2[ps2$logFC > 0,]) 
# df <- AnnotationDbi::select(hugene20sttranscriptcluster.db,ps2_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
# dplyr::mutate(df,GENENAME=stringr::str_trunc(GENENAME,30))

## Volcano plots for visualization
volcanoplot(fit2)
interesting_genes <- topTable(fit2,number=Inf, p.value = 0.05,lfc = 2)
interesting_genes
volcanoplot(fit2, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes)))
# points(interesting_genes[['logFC']],-log10(interesting_genes[['P.Value']]),col='red')
# This data do not have 'logFC' in the topTable, so we can not run line 90 for circling interesting genes in red.

## Heatmaps
eset_of_interest <- gse50737_eset[rownames(interesting_genes),]
heatmap(exprs(eset_of_interest))
library(RColorBrewer)
heatmap(exprs(eset_of_interest),
        labCol=gse50737_eset[['treatment:ch1']] ,labRow=NA,
        col       = rev(brewer.pal(10, "RdBu")),
        distfun   = function(x) as.dist(1-cor(t(x))))

legend(x = "topright", legend = c("high", "medium", "low"), cex = 0.6, fill = colorRampPalette(brewer.pal(10, "RdBu"))(3))

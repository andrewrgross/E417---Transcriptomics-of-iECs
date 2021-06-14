### DESeq2 Pipeline - Andrew R Gross - 2018-08-06
### A standard exectution of the DESeq2 pipeline.  
### INPUT: This script requires a counts table.  Normalized TPM data and sample data is recommended.
### OUTPUT: This script generates a table normalized expression with fold change and p-values between sample groups.

####################################################################################################################################################
### Header
library("pasilla")
library("DESeq2")
library("biomaRt")
library("VennDiagram")

ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
listMarts(host="www.ensembl.org")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

round.DESeq.results <- function(dataframe) {
  dataframe$baseMean <- round(dataframe$baseMean, 2)
  dataframe$log2FoldChange <- round(dataframe$log2FoldChange, 2)
  dataframe$lfcSE <- round(dataframe$lfcSE, 3)
  dataframe$stat <- round(dataframe$stat, 2)
  #dataframe$pvalue <- formatC(dataframe$pvalue, format = "e", digits = 2)
  #dataframe$padj <- formatC(dataframe$padj, format = "e", digits = 2)
  return(dataframe)
}
convertIDs <- function(dataframe) {
  ensemblIDs <- c()
  for (rowName in row.names(dataframe)) {
    ensemblID <- strsplit(rowName,"\\.")[[1]][1]
    ensemblIDs <- c(ensemblIDs,ensemblID)
  }
  row.names(dataframe)<-ensemblIDs
  return(dataframe)
}
addGene <- function(dataframe) {
  genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
  genes <- genes[match(row.names(dataframe),genes[,1]),]
  Gene <- c()
  for (rowNumber in 1:length(genes[,1])) {
    newGene <- genes[rowNumber,][,2]
    Gene <- c(Gene, newGene)
  }
  dataframe[length(dataframe)+1] <- Gene
  names(dataframe)[ncol(dataframe)] <- "Gene"
  return(dataframe)
}
add.description <- function(dataframe) {
  descr <- getBM(attributes=c('ensembl_gene_id','description'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
  descr <- descr[match(row.names(dataframe),descr[,1]),]
  descriptions <- c()
  for (rowNumber in 1:length(descr[,1])) {
    newDescr <- descr[rowNumber,][,2]
    newDescr <- strsplit(newDescr, " \\[")[[1]][1]
    descriptions <- c(descriptions, newDescr)
  }
  dataframe[length(dataframe)+1] <- descriptions
  names(dataframe)[ncol(dataframe)] <- "Description"
  return(dataframe)
}
####################################################################################################################################################
### Input

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/RNAseq Data/Motor Neurons/E099 - PM - ALS v CTRL/')
counts.als <- read.csv('PM-5119--07--02--2018_COUNTS.csv', row.names = 1)
tpm.als <- read.csv('PM-5119--07--02--2018_TPM.csv', row.names = 1)
sample.names.als <- c('02iCTR','03iCTR','159iALS','172iCTR','372iALS_n1','372iALS_n2','372iALS_n3','395iCTR')

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/RNAseq Data/Motor Neurons/E0x - PP - D10 KO v CTRL/')
counts.ko <- read.csv('COUNTS.csv', row.names = 1)
tpm.ko <- read.csv('TPM.csv', row.names = 1)
sample.names.ko <- read.csv('Sample names.txt', header = FALSE)

####################################################################################################################################################
### Format
##########################################################################
### Reassign names
names(counts.als) <- sample.names.als
names(tpm.als) <- sample.names.als

names(counts.ko) <- sample.names.ko[,1]
names(tpm.ko) <- sample.names.ko[,1]

##########################################################################
### Reorder columns
counts.als <- counts.als[c(1,2,4,8,5,6,7,3)]
tpm.als <- tpm.als[c(1,2,4,8,5,6,7,3)]

### Make a column data data frame
columnData.als <- data.frame(rownames = names(counts.als), disease = c("ctrl","ctrl","ctrl","ctrl","ALS","ALS","ALS","ALS"))
columnData.ko <- data.frame(rownames= names(counts.ko), disease = c('KO', 'KO', 'KO', 'Iso', 'Iso', 'Iso'))

# Convert counts into a matrix
counts.als <- round(as.matrix(counts.als), 0)
counts.ko <- round(as.matrix(counts.ko), 0)

head(counts.als)
head(counts.ko)

####################################################################################################################################################
### Differential Expression
### Make our DESeq data sets
dds.als <- DESeqDataSetFromMatrix(countData = counts.als, colData = columnData.als, design = ~ disease)
dds.ko <- DESeqDataSetFromMatrix(countData = counts.ko, colData = columnData.ko, design = ~ disease)
### Run DESeq
dds.als <- DESeq(dds.als)
dds.ko <- DESeq(dds.ko)
#results(dds.als)

####################################################################################################################################################
### Format Results
results.als <- as.data.frame(results(dds.als))
results.ko <- as.data.frame(results(dds.ko))

### Drop NAs
na.check <- is.na(results.als$padj)
results.als <- results.als[!is.na(results.als$padj),]
na.check <- is.na(results.ko$padj)
results.ko <- results.ko[!is.na(results.ko$padj),]



##########################################################################
### Filter results
### Round figures
pvalue.cutoff <- 0.05

### Keep p-adjusted
results.als <- round.DESeq.results(results.als)[-c(3,4,5)]
results.ko <- round.DESeq.results(results.ko)[-c(3,4,5)]

### Keep raw p-value
results.als <- round.DESeq.results(results.als)[-c(3,4,6)]
#results.ko <- round.DESeq.results(results.ko)[-c(3,4,6)]

### Apply cutoff
results.als <- results.als[results.als$padj <= pvalue.cutoff,]
results.ko <- results.ko[results.ko$padj <= pvalue.cutoff,]

results.als <- results.als[results.als$pvalue <= pvalue.cutoff,]
#results.ko <- results.ko[results.ko$pvalue <= pvalue.cutoff,]

### Reorder by p-value
results.als <- results.als[order(results.als$padj),]
results.ko <- results.ko[order(results.ko$padj),]

results.als <- results.als[order(results.als$pvalue),]
results.ko <- results.ko[order(results.ko$pvalue),]

### Reorder by Fold-change
#results.als <- results.als[order(results.als$log2FoldChange),]
#results.ko <- results.ko[order(results.ko$log2FoldChange),]

##########################################################################
### Apply cutoffs
summary(counts.als)
counts.cutoff = 1
results.als <- results.als[results.als$basemean >= counts.cutoff,]

#summary(counts.ko)
#counts.cutoff = 1
#results.ko <- results.ko[results.ko$basemean >= counts.cutoff]

##########################################################################
### Add in normalized expression values
tpm.als <- tpm.als[row.names(results.als),]
results.als <- cbind(results.als, tpm.als)
tpm.ko <- tpm.ko[row.names(results.ko),]
results.ko <- cbind(results.ko, tpm.ko)

### Annotate genes
results.als <- convertIDs(results.als)
results.als <- addGene(results.als)
results.als <- add.description(results.als)

results.ko <- convertIDs(results.ko)
results.ko <- addGene(results.ko)
results.ko <- add.description(results.ko)

##########################################################################
### Output the data

getwd()
setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/Input files/")

write.csv(results.als, 'DEGs - ALS 0.05 pv.csv')
#write.csv(results.ko, 'DEGs - KO 0.05 pv.csv')
write.csv(results.als, 'DEGs - ALS 0.05 padj.csv')
write.csv(results.ko, 'DEGs - KO 0.05 padj.csv')


##########################################################################
##########################################################################
### Compare KO data to ALS line data

setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/Input files/")

results.als <- read.csv('DEGs - ALS 0.05 pv.csv', row.names= 1)
results.ko <- read.csv('DEGs - KO 0.05 pv.csv', row.names = 1)

### Find the overlap

shared.genes <- intersect(results.als$Gene, results.ko$Gene)

length(shared.genes)

### Draw Venn Diagram
draw.pairwise.venn(area1 = nrow(results.als), area2 = nrow(results.ko), cross.area = length(shared.genes), 
                   category = c("Mutant vs. Wildtype", "Knockout vs. Wildtype"),
                   lty = rep("blank", 2), fill = c("light blue", "pink"), 
                   alpha = rep(0.5, 2), cat.pos = c(0, 0), 
                   cat.dist = rep(0.025, 2))
                                                                                                                                                        0), cat.dist = rep(0.025, 2))

### Find expression values
shared.pos <- match(shared.genes, results.als$Gene)
results.als.shrd <- results.als[shared.pos,]

shared.pos <- match(shared.genes, results.ko$Gene)
results.ko.shrd <- results.ko[shared.pos,]

### Create data frame of shared genes
names(results.als.shrd)[1:3] <- c('Mut: Mean', 'Mut: lFC', 'Mut: pv')
names(results.ko.shrd)[1:3] <- c('KO: Mean', 'KO: lFC', 'KO: pv')

results.als.shrd <- results.als.shrd[order(row.names(results.als.shrd)),]
results.ko.shrd <- results.ko.shrd[order(row.names(results.ko.shrd)),]

results.shared.genes <- cbind(results.als.shrd, results.ko.shrd)
results.shared.genes <- results.shared.genes[c(12,13,1,2,3,14,15,16)]
shared.genes.order <- order(results.shared.genes$`Mut: pv`,decreasing = TRUE)
results.shared.genes <- results.shared.genes[shared.genes.order,]
dir.check <- sign(results.shared.genes$`Mut: lFC`*results.shared.genes$`KO: lFC`)==1
results.shared.genes$same.dir <- dir.check

setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/DEG analyses/")

write.csv(results.shared.genes, 'Shared genes - 0.05 padj.csv')

### Sort by direction

up.reg <- results.ALS.shared$log2FoldChange > 0
als.up <- results.ALS.shared[up.reg,]
als.down <- results.ALS.shared[!up.reg,]

up.reg <- results.ko.shared$log2FoldChange > 0
ko.up <- results.ko.shared[up.reg,]
ko.down <- results.ko.shared[!up.reg,]

### Find Shared by direction

shared.up <- intersect(als.up$Gene, ko.up$Gene)
shared.down <- intersect(als.down$Gene, ko.down$Gene)

shared.up.pos <- match(shared.up, als.up$Gene)
shared.up.df <- als.up[shared.up.pos,]

shared.down.pos <- match(shared.down, als.down$Gene)
shared.down.df <- als.down[shared.down.pos,]

### Save results

setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Prasanthi/")

write.csv(shared.up.df, "Shared UP reg - ALS-KO.csv")
write.csv(shared.down.df, "Shared DOWN reg - ALS-KO.csv")
write.csv(results.ALS.shared, 'Shared DEGs - ALS-KO.csv')


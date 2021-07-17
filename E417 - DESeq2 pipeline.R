### DESeq2 Pipeline - Andrew R Gross - 21JUN21
### Identification of genes perturbed in iEC differentiation with various levels of ANG1  
### INPUT: This script requires a counts table.  Normalized TPM data and sample data is recommended.
### OUTPUT: This script generates a table normalized expression with fold change and p-values between sample groups.

####################################################################################################################################################
### 1- Header
####################################################################################################################################################
library("pasilla")
library("DESeq2")
library("biomaRt")
library(gplots)
#library("VennDiagram")

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
convert.ids <- function(dataframe, add.gene.name.column = TRUE) {
  ### This function will convert a row name consisting of a contactenated ensembl ID and gene to one or the other,
  ### based on the users instruction (2018-10-04)
  ensemblIDs <- c()                                           # Empty lists are initialized to receive IDs as they're created
  gene.names <- c()
  for (rowName in row.names(dataframe)) {                     # Loops through all rows in the data frame
    ensemblID <- strsplit(rowName,"\\.")[[1]][1]                 # Splits the row name and declares the ensembl ID
    gene.name <- strsplit(rowName,"\\_")[[1]][2]                 # Splits the row name, declares the gene name
    ensemblIDs <- c(ensemblIDs, ensemblID)                       # Adds ensembl ID and gene name to appropriate lists
    gene.names <- c(gene.names, gene.name)
  }
  row.names(dataframe) <- make.unique(ensemblIDs)                          # assigns the new row names
  if(add.gene.name.column == TRUE) {
    dataframe$Gene <- gene.names
  }
  return(dataframe)                                           # Returns the data frame with new rows
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
filter.low.rows <- function(dataframe, rm.rows.w.max.below.quant = 0.25) {
  rm.rows.w.max.below.quant <- as.numeric(rm.rows.w.max.below.quant)
  print(paste('Initial number of genes:', nrow(dataframe)))
  row.max <- apply(dataframe,1, max)
  empty.rows <- which(row.max == 0)
  print(paste('Genes with 0 expression:', length(empty.rows)))
  dataframe <- dataframe[-empty.rows,]
  row.max <- apply(dataframe,1, max)
  quantile.25 = quantile(row.max, rm.rows.w.max.below.quant)[[1]]
  print(paste(rm.rows.w.max.below.quant * 100, 'th percentile:', quantile.25))
  below.quant <- which(row.max <= quantile.25 )
  print(paste('Genes below the 25th percentile:', length(below.quant)))
  dataframe <- dataframe[-below.quant]
  print(paste('Final number of genes:', nrow(dataframe)))
  return(dataframe)
}

segregate.unique.genes <- function(dataframe) {
  gene.uniqueness <- isUnique(dataframe$Gene)
  dataframe.u <- dataframe[which(gene.uniqueness==TRUE),]
  dataframe.nu <- dataframe[which(gene.uniqueness==FALSE),]
  return(list(dataframe.u, dataframe.nu))
}

make.df.w.unique.gene.names <- function(dataframe.pair) {
  dataframe <- dataframe.pair[[2]]
  non.unique.genes <- names(table(dataframe$Gene)>1)
  
  non.unique.gene.rows <- dataframe[0,]
  for(gene in non.unique.genes) {
    rows.with.gene <- dataframe[which(dataframe$Gene == gene),][-ncol(dataframe)]
    rows.with.gene[1,] <- apply(rows.with.gene,2,max)
    row.names(rows.with.gene)[1] <- gene
    non.unique.gene.rows <- rbind(non.unique.gene.rows, rows.with.gene[1,])
  }
  dataframes.joined <- dataframe.pair[[1]][-ncol(dataframe.pair[[1]])]
  row.names(dataframes.joined) <- dataframe.pair[[1]]$Gene
  dataframes.joined <- rbind(dataframes.joined, non.unique.gene.rows)
  return(dataframes.joined)
}
remove.na.and.filter.adj.pval <- function(results.df, adj.pval.cutoff = 0.005) {
  results.df <- results.df[!is.na(results.df$padj),]
  results.df <- results.df[results.df$padj <= adj.pval.cutoff,]
  results.df <- results.df[order(results.df$padj, decreasing = FALSE),]
  return(results.df)
}
####################################################################################################################################################
### 2- Input
############################################################################################
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/RNAseq Data/iECs/')
counts.data <- read.csv('RDSS-12511--04--28--2021_COUNTS.csv', row.names = 1)
tpm.data <- read.csv('RDSS-12511--04--28--2021_TPM.csv', row.names = 1 )
metadata.data <- read.csv('metadata.csv', fileEncoding="UTF-8-BOM")

####################################################################################################################################################
### 3 - Format
##########################################################################
### Reassign sample names (columns)
names(counts.data) <- metadata.data$Shortname
names(tpm.data)    <- metadata.data$Shortname

### Filter out the lowest 25th percentile after removing empty rows
counts.data <- filter.low.rows(counts.data, rm.rows.w.max.below.quant = 0.25)
### Convert transcript names(rows)
counts.data <- convert.ids(counts.data, add.gene.name.column = FALSE)

### Convert TPM ID labels
tpm.data <- convert.ids(tpm.data, add.gene.name.column = FALSE)
### Review the levels of non-unique genes:
#gene.uniqueness <- isUnique(counts.data2$Gene)
#summary(gene.uniqueness)
#counts.data2.u <- counts.data2[which(gene.uniqueness==TRUE),]
#counts.data2.nu <- counts.data2[which(gene.uniqueness==FALSE),]



#write.csv(counts.data3, 'RDSS-12511--04--28--2021_COUNTS-gene-names.csv')

##########################################################################
### 3.1 - Format for DESeq

### Make a column data data frame
(columnData.v.ipsc   <- data.frame(row.names = metadata.data$Shortname, condition = metadata.data$Celltype))
(columnData.50.ctrl  <- data.frame(row.names = metadata.data$Shortname, condition = metadata.data$Condition)[4:9, , drop = FALSE])
(columnData.150.ctrl <- data.frame(row.names = metadata.data$Shortname, condition = metadata.data$Condition)[c(1,2,3,7,8,9), , drop = FALSE])
(columnData.50.150   <- data.frame(row.names = metadata.data$Shortname, condition = metadata.data$Condition)[1:6, , drop = FALSE])
(columnData.ang.ctrl <- data.frame(row.names = metadata.data$Shortname, condition = c(rep('Ang',6), rep('Ctrl',6)))[1:9, , drop = FALSE])

####################################################################################################################################################
### 4 - Differential Expression
############################################################################################
### Make our DESeq data sets

dds.iec.v.ipsc <- DESeqDataSetFromMatrix(countData = as.matrix(counts.data), colData = columnData.v.ipsc, design = ~ condition) ; title = 'iECs vs iPSCs'
dds.50.v.ctrl <- DESeqDataSetFromMatrix(countData = as.matrix(counts.data[row.names(columnData.50.ctrl)]), colData = columnData.50.ctrl, design = ~ condition) ; title = 'ANG1-50 vs Ctrl'
dds.150.v.ctrl <- DESeqDataSetFromMatrix(countData = as.matrix(counts.data[row.names(columnData.150.ctrl)]), colData = columnData.150.ctrl, design = ~ condition) ; title = 'ANG1-150 vs Ctrl'
dds.50.v.150 <- DESeqDataSetFromMatrix(countData = as.matrix(counts.data[row.names(columnData.50.150)]), colData = columnData.50.150, design = ~ condition) ; 
dds.ang.v.ctrl<-DESeqDataSetFromMatrix(countData = as.matrix(counts.data[row.names(columnData.ang.ctrl)]), colData = columnData.ang.ctrl, design = ~ condition) ; 


### Run DESeq
dds.iec.v.ipsc <- DESeq(dds.iec.v.ipsc)
dds.50.v.ctrl <- DESeq(dds.50.v.ctrl)
dds.150.v.ctrl <- DESeq(dds.150.v.ctrl)
dds.50.v.150 <- DESeq(dds.50.v.150)
dds.ang.v.ctrl <- DESeq(dds.ang.v.ctrl)

### Define results table
results.iec.v.ipsc <- as.data.frame(results(dds.iec.v.ipsc))
results.50.v.ctrl  <- as.data.frame(results(dds.50.v.ctrl))
results.150.v.ctrl <- as.data.frame(results(dds.150.v.ctrl))
results.50.v.150   <- as.data.frame(results(dds.50.v.150))
results.ang.v.ctrl <- as.data.frame(results(dds.ang.v.ctrl))


####################################################################################################################################################
### 4.1 - Format Results

### Remove NAs, filter by adj. pval., and sort by adj. pval
results.iec.v.ipsc <- remove.na.and.filter.adj.pval(results.iec.v.ipsc, adj.pval.cutoff = 0.005)
results.50.v.ctrl  <- remove.na.and.filter.adj.pval(results.50.v.ctrl, adj.pval.cutoff = 0.005)
results.150.v.ctrl <- remove.na.and.filter.adj.pval(results.150.v.ctrl, adj.pval.cutoff = 0.005)
results.50.v.150   <- remove.na.and.filter.adj.pval(results.50.v.150, adj.pval.cutoff = 0.005)
results.ang.v.ctrl <- remove.na.and.filter.adj.pval(results.ang.v.ctrl, adj.pval.cutoff = 0.005)


### Add back in the expression counts
results.iec.v.ipsc <- cbind(results.iec.v.ipsc, round(tpm.data[row.names(results.iec.v.ipsc),],0))
results.50.v.ctrl  <- cbind(results.50.v.ctrl,  round(tpm.data[row.names(results.50.v.ctrl),],0))
results.150.v.ctrl <- cbind(results.150.v.ctrl, round(tpm.data[row.names(results.150.v.ctrl),],0))
results.50.v.150   <- cbind(results.50.v.150,   round(tpm.data[row.names(results.50.v.150),],0))
results.ang.v.ctrl <- cbind(results.ang.v.ctrl, round(tpm.data[row.names(results.ang.v.ctrl),],0))

### Annotate genes
results.all  <- convert.ids(results.iec.v.ipsc)
results.all  <- add.description(results.all)

results.50   <- convert.ids(results.50.v.ctrl)
results.50   <- add.description(results.50)

results.150  <- convert.ids(results.150.v.ctrl)
results.150  <- add.description(results.150)

results.ava  <- convert.ids(results.50.v.150)
results.ava  <- add.description(results.ava)

results.ang  <- convert.ids(results.ang.v.ctrl)
results.ang  <- add.description(results.ang)
############################################################################################
### 5 - Output the data
############################################################################################
getwd()
setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E417 - RNAseq analysis of iECs/DE")

write.csv(results.all, 'DEGs - iECs v. iPSCs.csv')
write.csv(results.50 , 'DEGs - ANG1-50 v. Ctrl.csv')
write.csv(results.150, 'DEGs - ANG1-150 v. Ctrl.csv')
write.csv(results.ava, 'DEGs - ANG1-50 v. ANG1-150.csv')
write.csv(results.ang, 'DEGs - ANG1-50&150 v. Ctrl.csv')


############################################################################################
### 6 - Identify all the genes which went up relative to iPSCs
############################################################################################
setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E417 - RNAseq analysis of iECs/DE")
results.all  <- read.csv('DEGs - iECs v. iPSCs.csv')
results.50   <- read.csv('DEGs - ANG1-50 v. Ctrl.csv')
results.150  <- read.csv('DEGs - ANG1-150 v. Ctrl.csv')
results.ang  <- read.csv()

up.in.iec   <- results.all[which(results.all$log2FoldChange <0),]
down.in.iec <- results.all[which(results.all$log2FoldChange >0),]

up.in.iec   <- up.in.iec[order(up.in.iec$log2FoldChange, decreasing = FALSE),]
down.in.iec <- down.in.iec[order(down.in.iec$log2FoldChange, decreasing = TRUE),]

### Here are the genes that went from absent in iPSC to present:
median.in.ipscs <- apply(up.in.iec[c(16,17,18)], 1, median)
up.in.iec.no.exp.in.ipsc <- up.in.iec[which(median.in.ipscs <= 1),]

### Retain only genes expressed in the top 80th percentile in iECs
(percentile.8 = quantile(as.matrix(results.all[7:15]), 0.8))
min.in.iec <- apply(up.in.iec.no.exp.in.ipsc[7:15], 1, min)
up.in.iec.no.exp.in.ipsc <- up.in.iec.no.exp.in.ipsc[min.in.iec >= percentile.8[[1]],]

### Here are the genes with the highest expression among iECs
(percentile.95 = quantile(as.matrix(results.all[7:15]), 0.95))
median.in.iecs <- apply(up.in.iec[7:15], 1, median)
up.in.iec.w.high.exp <- up.in.iec[median.in.iecs >= percentile.95[[1]],]
up.in.iec.w.high.exp <- up.in.iec.w.high.exp[order(up.in.iec.w.high.exp$baseMean, decreasing = TRUE),]

### Separate out genes with no expression in iEC
# Later


############################################################################################
### 7 - Output the data
############################################################################################
getwd()
setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E417 - RNAseq analysis of iECs/DE")

write.csv(up.in.iec.no.exp.in.ipsc, 'DEGs - iEC genes activated from zero.csv')
write.csv(up.in.iec.w.high.exp , 'DEGs - iEC genes activated to high expr.csv')


############################################################################################
###  8 - Identify all the genes went up with expression of ANG1
############################################################################################
med.in.ctrl      <- apply(up.in.iec[13:15], 1, median)
med.in.ang.50    <- apply(up.in.iec[10:12], 1, median)
med.in.ang.150   <- apply(up.in.iec[7:9], 1, median)
med.across.conditions <- cbind(up.in.iec, med.in.ctrl, med.in.ang.50, med.in.ang.150)
med.across.conditions$del.1 <- med.across.conditions$med.in.ang.50  - med.across.conditions$med.in.ctrl
med.across.conditions$del.2 <- med.across.conditions$med.in.ang.150 - med.across.conditions$med.in.ang.50

### Filter to keep only those that went up each time
med.across.conditions <- med.across.conditions[med.across.conditions$del.1 >0,]
med.across.conditions <- med.across.conditions[med.across.conditions$del.2 >0,]

med.across.conditions$del.mult <- med.across.conditions$del.1 * med.across.conditions$del.2
med.across.conditions <- med.across.conditions[c(6,7,8,9,10,11,12,13,14,15,16,18,19,21,22,23,24,25,26,19,20)]
med.across.conditions <- med.across.conditions[order(med.across.conditions$del.mult, decreasing = TRUE),]


############################################################################################
### 8.1 - Add in the pvalue for the t-test between treatments

results.50.reordered <- results.50[row.names(med.across.conditions),]
results.150.reordered <- results.150[row.names(med.across.conditions),]
med.across.conditions <- cbind(med.across.conditions, results.50.reordered$padj, results.150.reordered$padj)
med.across.conditions <- med.across.conditions[!is.na(med.across.conditions$`results.50.reordered$padj`),]
med.across.conditions <- med.across.conditions[!is.na(med.across.conditions$`results.150.reordered$padj`),]
med.across.conditions$sig.padj <- apply(med.across.conditions[c(1,21,22)],1, sum)
med.across.conditions <- med.across.conditions[order(med.across.conditions$sig.padj, decreasing = FALSE),]


getwd()
setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E417 - RNAseq analysis of iECs/DE")

write.csv(med.across.conditions, 'Genes that increased with ANG1.csv')

####################################################################################################################################################
### 9 - Reimport the data
############################################################################################

setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E417 - RNAseq analysis of iECs/DE")

up.in.iec.no.exp.in.ipSC <- read.csv('DEGs - iEC genes activated from zero.csv', row.names = 1)
up.in.iec.w.high.exp <- read.csv('DEGs - iEC genes activated to high expr.csv', row.names = 1)
med.across.conditions <- read.csv('Genes that increased with ANG1.csv', row.names = 1)
results.all <- read.csv('DEGs - iECs v. iPSCs.csv', row.names = 1)
results.50 <-  read.csv('DEGs - ANG1-50 v. Ctrl.csv', row.names = 1)
results.150 <-read.csv('DEGs - ANG1-150 v. Ctrl.csv', row.names = 1)
results.ava <- read.csv('DEGs - ANG1-50 v. ANG1-150.csv', row.names = 1)
results.ang <- read.csv('DEGs - ANG1-50&150 v. Ctrl.csv', row.names = 1)

genes.that.went.up <- read.csv('Diff Exp - Genes that went up with ANG1.csv', row.names = 1)

############################################################################################
### Select the genes to plot

genes.to.plot <- genes.that.went.up[1:10,]
expr.to.plot <- tpm.data[row.names(genes.to.plot),]


geneNumber = 10
(df.to.plot <- data.frame(Condition = factor(ang.metadata$Condition, levels = c('iPSC', 'Ctrl', 'Ang1_50', 'Ang1_150')), Expression = unlist(expr.to.plot[geneNumber,])))
(currentGene = genes.to.plot$Description[geneNumber])

##########################################################################
### Plot bar graphs

currentPlot <- ggplot(data = df.to.plot, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_boxplot(varwidth = FALSE, size = 0.5) +
  scale_y_continuous(limits = c(0,NA), breaks = seq(0,10000,1000)) +
  scale_x_discrete(labels = c('iPSC', 'Control iECs', '+50 ng/mL ANG1', '+150 ng/mL ANG1')) +
  labs(title = currentGene, 
       x = 'Condition', 
       y = 'Expression [TPM]') +
  theme(plot.title = element_text(color="black", face="bold", size=18, hjust = 0.5),
        axis.text = element_text(size = 14),
        axis.title.x = element_text(face="bold", size=12,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
        legend.position = 'none',
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(color = "black", fill = NA)) +
  scale_fill_brewer(palette = 'BuGn') 
currentPlot
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E417 - RNAseq analysis of iECs/Boxplots/')
png(filename=paste0(currentGene,'-',strftime(Sys.time(),"%a%b%d%H%M"),".png"), 
    type="cairo",
    units="in", 
    width=14, 
    height=14, 
    pointsize=12, 
    res=300)
currentPlot

dev.off()

##########################################################################
### Heatmaps
pal <- colorRampPalette(c('yellow','black','blue'))(100)
expression.df <- results.ang
expression.df <- expression.df[expression.df$padj < 0.0001,]; (nrow(expression.df))
expression.df <- expression.df[order(expression.df$log2FoldChange),]
expression.m <- as.matrix(expression.df[7:18])
expression.m <- as.matrix(expression.df[7:15])

heatmap(expression.m, col = pal, Rowv = NA, Colv = NA, scale = "row", margins = c(8,6))

heatmap.2(expression.m, scale = 'row', col = pal, trace = 'none', dendrogram="none", 
          Rowv=FALSE, symm=TRUE, density.info='none', labRow=NA)


heatmap.2(expression.m, scale = 'row', col = pal, trace = 'none', dendrogram="none", 
          Rowv=FALSE, symm=TRUE, density.info='none', labRow=NA,
          lmat=rbind(c(4, 2), c(1, 3)), lhei=c(2, 8), lwid=c(4, 1), margins = c(4,0))





















generate.bargraphs <- function(dataframe) {
  plot.list <- list()
  for(gene.column in 1:(ncol(dataframe)-2)){
    current.plot <- ggplot(data = dataframe, 
                           aes_string(x = 'Sample', y = names(dataframe)[gene.column], fill = "Type")) +
      geom_bar(stat = 'identity') +
      labs(title = names(dataframe[gene.column])) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x = element_blank(),
            axis.title.y = element_blank(), legend.position = 'none',
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            panel.border = element_rect(color = "black", fill = NA)) +
      scale_fill_brewer(palette = 'Set1')
    plot.list[[length(plot.list)+1]] <- current.plot
  }
  print(paste("List contains",length(plot.list),"plots"))
  return(plot.list)
}
generate.boxplots <- function(dataframe) {
  plot.list <- list()
  for(gene.column in 1:(ncol(dataframe)-2)){
    current.plot <- ggplot(data = dataframe, 
                           aes_string(x = 'Type', y = names(dataframe)[gene.column], fill = "Type")) +
      geom_boxplot(varwidth = FALSE) +
      scale_y_continuous(limits = c(0,NA)) + geom_jitter(width = 0, size = 2) +
      labs(title = names(dataframe[gene.column])) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x = element_blank(),
            axis.title.y = element_blank(), legend.position = 'none',
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            panel.border = element_rect(color = "black", fill = NA)) +
      scale_fill_brewer(palette = 'Set1')
    plot.list[[length(plot.list)+1]] <- current.plot
  }
  print(paste("List contains",length(plot.list),"plots"))
  return(plot.list)
}
########################################################################
### Prepare plot data
########################################################################

#genes.df$Reference <- ht.reference$Brain...Hypothalamus[1:nrow(genes.df)]
plot.df <- data.frame(t(genes.df))
plot.df$Type <- factor(metaData$Type[match(row.names(plot.df),row.names(metaData))], levels = c('aHT','CTR','OBS','iMN'))
plot.df$Sample <- factor(row.names(plot.df), levels = names(genes.df))

########################################################################
### Generate plots
########################################################################
###### Bars

bar.list <- generate.bargraphs(plot.df)

###### Boxes
box.list <- generate.boxplots(plot.df)

########################################################################
### Select data to plot
########################################################################

plot.list <- bar.list

plot.list <- box.list#[13:24]

########################################################################
### Display tiled plots
########################################################################

(tiled.figure <- grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]],
                              ncol = 3, top = title))

(tiled.figure <- grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]],
                              plot.list[[7]],plot.list[[8]],plot.list[[9]],plot.list[[10]],plot.list[[11]],plot.list[[12]],
                              ncol = 4, top = title))

(tiled.figure <- grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],ncol = 4, top = title))

(tiled.figure <- grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]],
                              plot.list[[7]],plot.list[[8]],plot.list[[9]],
                              ncol = 3, top = title))


########################################################################
### Export plot
########################################################################

setwd("z:/Uthra/HT paper/Bioinformatics figures/Plots of gene expression/")

tiff(filename=paste0(substr(title,1,6),'_2-',strftime(Sys.time(),"%a%b%d%H%M"),".tiff"), 
     type="cairo",
     units="in", 
     width=14, 
     height=14, 
     pointsize=12, 
     res=300)
tiled.figure <- grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]],
                             plot.list[[7]],plot.list[[8]],plot.list[[9]],plot.list[[10]],plot.list[[11]],plot.list[[12]],
                             ncol = 4, top = title)
dev.off()

png(filename=paste0(substr(title,1,6),'_02-',strftime(Sys.time(),"%a%b%d%H%M"),".png"), 
    type="cairo",
    units="in", 
    width=14, 
    height=14, 
    pointsize=12, 
    res=300)
tiled.figure <- grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]],
                             plot.list[[7]],plot.list[[8]],plot.list[[9]],plot.list[[10]],plot.list[[11]],plot.list[[12]],
                             ncol = 4, top = title)
dev.off()

########################################################################
########################################################################
### Generate Single Plots
########################################################################
########################################################################




########################################################################
### Select geneset
########################################################################

genes.of.interest <- uthras.genes           ;title <- "Selected Hypothalamus genes"

genes.of.interest <- ht.genes_lit           ;title <- "Selected Hypothalamus genes"

genes.of.interest <- ht.genes_0.0001.df     ;title <- "Hypothalamus genes, pSI = 1e-4"

########################################################################
### Subsample
########################################################################
### Subsample rows

gene.positions <- match(genes.of.interest$Ensembl.ID,row.names(TPMdata))  # Declare the row numbers in the rnaseq data which correspond to genes in or list of genes of interest
genes.df <- TPMdata[gene.positions,]  # Make a dataframe containing just the rows of RNAseq data corresponding to genes of interest
### Rename the rows with genes
row.names(genes.df) <- genes.of.interest$Gene  # Replace the row names with the gene names
genes.df <- genes.df[complete.cases(genes.df),]

genes.with.med <- addMedSD(genes.df)
genes.df <- genes.df[genes.with.med$median>5,]


########################################################################
### Format data
########################################################################

#genes.df$Reference <- ht.reference$Brain...Hypothalamus[1:nrow(genes.df)]
plot.df <- data.frame(t(genes.df))
plot.df$Type <- as.character(metaData$Type[match(row.names(plot.df),row.names(metaData))], levels = c('aHT','CTR','OBS','iMN'))
plot.df$Sample <- factor(row.names(plot.df), levels = names(genes.df))

### Subsample columns
plot.df <- plot.df[plot.df$Type != 'iMN',]

### Reclassify type
plot.df$Type[plot.df$Type == 'CTR'] <- 'iHT'
plot.df$Type[plot.df$Type == 'OBS'] <- 'iHT'

########################################################################
### Loop through columns and output plots
########################################################################

setwd("z:/Uthra/HT paper/Bioinformatics figures/Plots of gene expression/Individual plots/")

for(col.numb in 1:(ncol(plot.df)-2)){
  gene <- names(plot.df)[col.numb]
  png(filename=paste0(gene,'-HT-',strftime(Sys.time(),"%a%b%d%H%M"),".png"), 
      type="cairo",
      units="in", 
      width=10, 
      height=10, 
      pointsize=20, 
      res=300)
  
  plot <-ggplot(data = plot.df, 
                aes_string(x = 'Type', y = gene, fill = "Type")) +
    geom_boxplot(varwidth = FALSE) +
    scale_y_continuous(limits = c(0,NA)) + geom_jitter(width = 0, size = 2) +
    labs(title = gene) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 20),axis.title.x = element_blank(),
          axis.text.y = element_text(size = 20), axis.title.y = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          panel.border = element_rect(color = "black", fill = NA)) +
    scale_fill_brewer(palette = 'Set1')
  print(plot)
  dev.off()
  print(paste(gene,'exported'))
}


########################################################################
### Export
########################################################################





####################################################################################################################################################
### Format Results


##########################################################################
### Heatmaps
pal <- colorRampPalette(c('yellow','black','blue'))(100)

heatmap.2(expression.m, scale = 'row', col = pal, trace = 'none', dendrogram="none", 
          Rowv=FALSE, symm=TRUE, density.info='none', labRow=NA,
          lmat=rbind(c(4, 2), c(1, 3)), lhei=c(2, 8), lwid=c(4, 1), margins = c(9,0))

heatmap(expression.m, col = pal, Rowv = NA, Colv = NA, scale = "row", margins = c(8,6))

##########################################################################
### Volcano Plots
volcano.data.cov <- iec.de
p.less.than = volcano.data.cov$padj <0.5
p.less.than[is.na(p.less.than)] = FALSE
volcano.data.cov <- volcano.data.cov[p.less.than,]
volcano.data.cov$Gene = volcano.data.cov[,1]
volcano.data.cov$Gene = rownames(volcano.data.cov)
subtitle = 'iEC vs HUVECs'
genes.ds <- volcano.data.cov$Gene[1:20]    # Genes of interest

### Converting FC from log2 to log10
volcano.data.cov$log10FC <- log10(2^volcano.data.cov$log2FoldChange)

volcano.data.cov.sig <- volcano.data.cov[volcano.data.cov$padj<0.005,]
volcano.data.cov.sig.up <- volcano.data.cov.sig[volcano.data.cov.sig$log2FoldChange>0.5,]
volcano.data.cov.sig.down <- volcano.data.cov.sig[-volcano.data.cov.sig$log2FoldChange>0.5,]
volcano.data.cov.sig <- volcano.data.cov.sig[order(volcano.data.cov.sig$baseMean),]

### Selecting genes to list
volcano.data.cov.sig$Edge.up <- -log(volcano.data.cov.sig$padj) * volcano.data.cov.sig$log2FoldChange
volcano.data.cov.sig$Edge.down <- -log(volcano.data.cov.sig$padj) * -volcano.data.cov.sig$log2FoldChange

volcano.text.up <- volcano.data.cov.sig[volcano.data.cov.sig$Edge.up > 9,]; nrow(volcano.text.up)
volcano.text.down <- volcano.data.cov.sig[volcano.data.cov.sig$Edge.down > 9,]; nrow(volcano.text.down)

### Genes of interest
volcano.genes.of.interest <- volcano.data.cov.sig[volcano.data.cov.sig$Gene %in% genes.ds,]
volcano.genes.of.interest.null <- volcano.data.cov[volcano.data.cov$Gene %in% genes.ds,]

### Manually filtering
text.size = 2


ggplot( ) +
  geom_point(data = volcano.data.cov, aes(x=log2FoldChange, y = log10(padj) ),color = 'grey', size = 0.9) +
  geom_point(data = volcano.data.cov.sig.up, aes(x=log2FoldChange, y = log10(padj), size = log10(baseMean)), color = 'red') +
  geom_point(data = volcano.data.cov.sig.down, aes(x=log2FoldChange, y = log10(padj), size = log10(baseMean)), color = 'blue') +
  #  geom_text(data = volcano.text.up, aes(x=log2FoldChange + 0.2, y = log10(padj) - 0, label = Gene), 
  #            hjust = 0, vjust = 0.1, size = 3, check_overlap = TRUE) +
  #  geom_text(data = volcano.text.down, aes(x=log2FoldChange - 0.2, y = log10(padj) - 0, label = Gene), 
  #            hjust = 1, vjust = 0.1, size = 3, check_overlap = TRUE) +
  geom_point(data = volcano.genes.of.interest, aes(x=log2FoldChange, y = log10(padj), size = log10(baseMean)), color = 'yellow') +
  geom_text(data = volcano.genes.of.interest, aes(x=log2FoldChange + 0.2, y = log10(padj) - 0, label = Gene), hjust = 0, vjust = 0.1, size = 3, check_overlap = TRUE) +  
  scale_color_gradient(low="pink", high="red") +
  scale_size('Log10 Expression', range = c(0.5,4)) +
  ylim(c(0, min(log10(volcano.data.cov$padj)))) +
  xlim(c(min(volcano.data.cov$log2FoldChange-0.5), max(volcano.data.cov$log2FoldChange)+3)) +
  labs(title="P-Value vs. Log Fold Change Volcano Plot",
       subtitle = subtitle,
       x = 'Log2 Fold Change', 
       y = 'Log10 P-value (Adjusted)') +
  theme(plot.title = element_text(hjust = 0.5, color="black", face="bold", size=18, margin=margin(0,0,5,0)),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.x = element_text(face="bold", size=13,margin =margin(5,0,5,0)),
        axis.title.y = element_text(face="bold", size=13,margin =margin(0,5,0,5)),
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12),
        legend.position = 'none') 

volcano.text.down 
volcano.text.up





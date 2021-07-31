### E417 - Volcano plots -- Andrew R Gross -- 29JUL21
### A Volcano plot generation script  
### INPUT: This script requires a differential expression table 
### OUTPUT: This script outputs volcano plots

####################################################################################################################################################
### 1 - Header
############################################################################################
library (ggplot2)



####################################################################################################################################################
### 2 - Import differential expression tables
############################################################################################
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/RNAseq Data/iECs/')
tpm.data <- read.csv('RDSS-12511--04--28--2021_TPM.csv', row.names = 1 )
metadata.data <- read.csv('metadata.csv', fileEncoding="UTF-8-BOM")


results.50  <-  read.csv('DEGs - ANG1-50 v. Ctrl.csv', row.names = 1)
results.150 <-  read.csv('DEGs - ANG1-150 v. Ctrl.csv', row.names = 1)
results.ava <-  read.csv('DEGs - ANG1-50 v. ANG1-150.csv', row.names = 1)


genes.that.went.up <- read.csv('Diff Exp - Genes that went up with ANG1.csv', row.names = 1)


####################################################################################################################################################
### Format
##########################################################################

quantile(results.test$padj, c(0, 0.001, 0.002,  0.01, 0.05, 0.075, 0.1, 0.2, 0.3, 0.5))

test <- head(results.50, 20)
t.test( test[1,][c('Ang-1_50-1', 'Ang-1_50-2', 'Ang-1_50-3')], test[1,][c('Ctrl-1', 'Ctrl-2', 'Ctrl-3')])
tTestTest <- t.test(test)

volcano.data <- results.150[c('log2FoldChange')]

### Find the spread of each gene

counts.data
findRange <- function(set) {
  setRange = max(set)-min(set)
  return(setRange)
}
test <- counts.data[order(apply(counts.data, 1, max), decreasing = TRUE),]
test$std = apply(counts.data, 1, sd)

test$A150_range = apply(counts.data[1:3], 1, sd)
test$A50_range = apply(counts.data[4:6], 1, sd)
test$ctrl_range = apply(counts.data[7:9], 1, sd)
test$ipsc_range = apply(counts.data[10:12], 1, sd)
test$iec_range = apply(counts.data[1:9], 1, sd)

counts

##########################################################################
### Volcano Plots

### Format #############################################################################################
volcano.data <- results.150                                               # Assign the table to plot
p.less.than = volcano.data$padj <0.5                                      #
p.less.than[is.na(p.less.than)] = FALSE
volcano.data <- volcano.data[p.less.than,]
volcano.data$Gene = volcano.data[,1]
volcano.data$Gene = rownames(volcano.data)
genes.ds <- volcano.data$Gene[1:20]    # Genes of interest

### Converting FC from log2 to log10
volcano.data$log10FC <- log10(2^volcano.data$log2FoldChange)

volcano.data.sig <- volcano.data[volcano.data$padj<0.005,]
volcano.data.sig.up <- volcano.data.sig[volcano.data.sig$log2FoldChange>0.5,]
volcano.data.sig.down <- volcano.data.sig[-volcano.data.sig$log2FoldChange>0.5,]

volcano.data.sig <- volcano.data.sig[order(volcano.data.sig$baseMean),]


### Selecting genes to list
volcano.data.sig$Edge.up <- -log(volcano.data.sig$padj) * volcano.data.sig$log2FoldChange
volcano.data.sig$Edge.down <- -log(volcano.data.sig$padj) * -volcano.data.sig$log2FoldChange

volcano.text.up <- volcano.data.sig[volcano.data.sig$Edge.up > 9,]; nrow(volcano.text.up)
volcano.text.down <- volcano.data.sig[volcano.data.sig$Edge.down > 9,]; nrow(volcano.text.down)

### Genes of interest
volcano.genes.of.interest <- volcano.data.sig[volcano.data.sig$Gene %in% genes.ds,]
volcano.genes.of.interest.null <- volcano.data[volcano.data$Gene %in% genes.ds,]

### Manually filtering
#genes.to.omit <- c('ANKRD36', 'AC022400.7', 'LINC00342' )
#volcano.text.up <- volcano.text.up[!volcano.text.up$Gene %in% genes.to.omit,]
#volcano.text.down <- volcano.text.down[!volcano.text.down$Gene %in% genes.to.omit,]
text.size = 2


volcano <- ggplot( ) +
  geom_point(data = volcano.data, aes(x=log2FoldChange, y = log10(padj) ),color = 'grey', size = 0.9) +
  geom_point(data = volcano.data.sig.up, aes(x=log2FoldChange, y = log10(padj), size = log10(baseMean)), color = 'red') +
  geom_point(data = volcano.data.sig.down, aes(x=log2FoldChange, y = log10(padj), size = log10(baseMean)), color = 'blue') +
  #  geom_text(data = volcano.text.up, aes(x=log2FoldChange + 0.2, y = log10(padj) - 0, label = Gene), 
  #            hjust = 0, vjust = 0.1, size = 3, check_overlap = TRUE) +
  #  geom_text(data = volcano.text.down, aes(x=log2FoldChange - 0.2, y = log10(padj) - 0, label = Gene), 
  #            hjust = 1, vjust = 0.1, size = 3, check_overlap = TRUE) +
  geom_point(data = volcano.genes.of.interest, aes(x=log2FoldChange, y = log10(padj), size = log10(baseMean)), color = 'yellow') +
  geom_text(data = volcano.genes.of.interest, aes(x=log2FoldChange + 0.2, y = log10(padj) - 0, label = Gene), hjust = 0, vjust = 0.1, size = 5, check_overlap = TRUE) +  
  scale_color_gradient(low="pink", high="red") +
  scale_size('Log10 Expression', range = c(0.5,4)) +
  ylim(c(0, min(log10(volcano.data$padj))-2)) +
  xlim(c(min(volcano.data$log2FoldChange-0.1), max(volcano.data$log2FoldChange)+2)) +
  labs(title= title,
       x = 'Log2 Fold Change', 
       y = 'Log10 P-value (Adjusted)') +
  theme(plot.title = element_text(color="black", face="bold", hjust = 0.5, size=22, margin=margin(10,0,20,0)),
        axis.title.x = element_text(face="bold", size=18,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=18,margin =margin(0,20,0,10)),
        panel.background = element_rect(fill = 'white', color = 'black', size = 3),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 16),
        legend.position = 'none') 


volcano

############################################################################################
### Write to folder
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E401 - Analysis of proteomics data/Volcano Plots/')


### Save plot
tiff(filename= paste0('Volcano - ',title, '.tiff'), width = 2000, height = 1600, units = "px", pointsize = 20, res = 250)
volcano
dev.off()



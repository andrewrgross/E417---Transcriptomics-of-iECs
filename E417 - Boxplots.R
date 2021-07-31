### E417 - Andrew R Gross - 29JUL21
### Import lists of genes that went up with administration of ANG1 and their expression levels and plot
### INPUT: This script requires a table listing the expression of genes that went up with increase in the adminstration of ANG1
### OUTPUT: This script generates boxplots.

####################################################################################################################################################
### 1- Header
####################################################################################################################################################
library("pasilla")
library("DESeq2")
library("biomaRt")
library(gplots)

####################################################################################################################################################
### 2 - Reimport filtered expression data
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
### 3 - Select the genes to plot
############################################################################################
genes.to.plot <- genes.that.went.up[1:10,]
expr.to.plot <- tpm.data[row.names(genes.to.plot),]


geneNumber = 10
(df.to.plot <- data.frame(Condition = factor(ang.metadata$Condition, levels = c('iPSC', 'Ctrl', 'Ang1_50', 'Ang1_150')), Expression = unlist(expr.to.plot[geneNumber,])))
(currentGene = genes.to.plot$Description[geneNumber])

##########################################################################
### 4 - Plot box plots
############################################################################################
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

########################################################################################################################################################################################
### 5 - Save plot
############################################################################################

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E417 - RNAseq analysis of iECs/Boxplots/')

tiff(filename= paste0('PCA 1v2.tiff'), width = 2000, height = 1600, units = "px", pointsize = 12, res = 250)
pca.1.v.2
dev.off()


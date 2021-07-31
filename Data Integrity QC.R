### Data Integrity Check - Andrew R Gross - 30JUL21
### This script performs statistical checks to look for any irregularities in data
### INPUT: This script requires a counts table.  
### OUTPUT: This script returns data comparisons

####################################################################################################################################################
### 1- Header
####################################################################################################################################################
library(reshape2)
library(ggplot2)

filter.low.rows <- function(dataframe, rm.rows.w.max.below.quant = 0.25) {
  rm.rows.w.max.below.quant <- as.numeric(rm.rows.w.max.below.quant)
  print(paste('Initial number of genes:', nrow(dataframe)))
  row.max <- apply(dataframe,1, max)
  empty.rows <- which(row.max == 0)
  print(paste('Genes with 0 expression:', length(empty.rows)))
  ### Now filter out based on quantile:
  dataframe <- dataframe[-empty.rows,]
  row.max <- apply(dataframe,1, max)
  quantile.25 = quantile(row.max, rm.rows.w.max.below.quant)[[1]]
  print(paste(rm.rows.w.max.below.quant * 100, 'th percentile:', quantile.25))
  below.quant <- which(row.max <= quantile.25 )
  print(paste('Genes below the', rm.rows.w.max.below.quant*100, 'th percentile:', length(below.quant)))
  dataframe <- dataframe[-below.quant,]
  print(paste('Final number of genes:', nrow(dataframe)))
  return(dataframe)
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

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/RNAseq Data/COVID tissues and cell cultures/For submission')
tpm.2 <- read.csv('VW-10228--08--04--2020_TPM.csv', row.names = 1)
sample.names <- c('mock-d1-r1', 'mock-d1-r2', 'mock-d1-r3', 
                  'cov-d1-r1', 'cov-d1-r2', 'cov-d1-r3', 
                  'mock-d3-r1', 'mock-d3-r2', 'mock-d3-r3', 
                  'cov-d3-r1', 'cov-d3-r2', 'cov-d3-r3')


####################################################################################################################################################
### 3 - Format
##########################################################################
### Reassign sample names (columns)
names(tpm.data)    <- metadata.data$Shortname
names(tpm.2)       <- sample.names

### Filter out the lowest 25th percentile after removing empty rows
tpm.data <- filter.low.rows(tpm.data, rm.rows.w.max.below.quant = 0.1)
quantile(as.matrix(tpm.data), seq(0,1,0.1))

tpm.2    <- filter.low.rows(tpm.2, rm.rows.w.max.below.quant = 0.1)
quantile(as.matrix(tpm.2), seq(0,1,0.1))

####################################################################################################################################################
### 4 - Compose a data frame of Std. Dev.
##########################################################################
### Compose a dataframe of SD
sd.df <- data.frame(Ang1_150 = apply(tpm.data[1:3], 1, sd),
                    Ang1_50  = apply(tpm.data[4:6], 1, sd),
                    ctrl     = apply(tpm.data[7:9], 1, sd),
                    ipsc     = apply(tpm.data[10:12], 1, sd) )

sd.df2 <- data.frame(mockD1 = apply(tpm.2[1:3], 1, sd),
                     covD1  = apply(tpm.2[4:6], 1, sd),
                     mockD3     = apply(tpm.2[7:9], 1, sd),
                     covD3     = apply(tpm.2[10:12], 1, sd) )

sd.df3 <- data.frame(SA = apply(tpm.2[2:4], 1, sd),
                     SB  = apply(tpm.2[5:7], 1, sd),
                     SC     = apply(tpm.2[8:10], 1, sd),
                     SD     = apply(tpm.2[c(1,11,12)], 1, sd) )

### Examine the distribution of std. dev
sd.m = as.matrix(sd.df)
sd.m2 = as.matrix(sd.df2)
sd.m3 = as.matrix(sd.df3)

quantile((sd.m), seq(0,1,0.1))
quantile((sd.m2), seq(0,1,0.1))
quantile((sd.m3), seq(0,1,0.1))


test <- melt(c(sd.df))
test <- melt(c(sd.df2))
test <- melt(c(sd.df3))

test <- melt(c(sd.df,sd.df2))
test <- melt(c(sd.df2, sd.df3))
test <- melt(c(sd.df, sd.df2, sd.df3))
### This method doesn't discern it well. This is a pretty tricky challenge. I guess I'd need several good data sets, and I'd need to 
## ... determine the range of values that are common. All I know right now is that the data looks fishy.

test$value = log(sqrt(test$value))


ggplot(data = test, aes(x = value, color = L1)) + geom_density()


quantile(sqrt(sd.m), seq(0,1,0.1))
quantile(sqrt(sd.m2), seq(0,1,0.1))

hist(log(sd.m))

hist(log(as.matrix(sd.df[1])))
hist(log(as.matrix(sd.df[2])))
hist(log(as.matrix(sd.df[3])))

hist(log(sd.m2))

hist(log(as.matrix(sd.df2[1])))
hist(log(as.matrix(sd.df2[2])))
hist(log(as.matrix(sd.df2[3])))



median(sd.m)
median(as.matrix(sd.df[1]))
median(as.matrix(sd.df[2]))
median(as.matrix(sd.df[3]))
median(as.matrix(sd.df[4]))

median(sd.m2)
median(as.matrix(sd.df2[1]))
median(as.matrix(sd.df2[2]))
median(as.matrix(sd.df2[3]))
median(as.matrix(sd.df2[4]))

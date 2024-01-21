library(edgeR)
library(limma)
library(org.Hs.eg.db)
library(singleCellNet)

# Primary Author: Yuqi Tan 
#filtering lowly expressed genes
#ref: https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html
expDatGood = singleCellNet::utils_loadObject("../input/TC71_alignment/expDatGood_091319.rda")
sampTabGood = singleCellNet::utils_loadObject("../input/TC71_alignment/sampTabGood_091319.rda")
myCPM <- cpm(expDatGood)
head(myCPM)
thresh <- myCPM > 0.5
table(rowSums(thresh))
keep <- rowSums(thresh) >= 1
counts.keep <- expDatGood[keep,]
summary(keep)
dim(counts.keep)

#convert counts into DGElist object
y <- DGEList(counts.keep)

#look at library size and distribution
barplot(y$samples$lib.size,names=colnames(y),las=2)
title("Barplot of library sizes")

# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)

# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

#multidimensional scaling plot
col.cell <- c("#1B9E77", "#7570B3", "orangered3")[as.factor(sampTabGood$description1)]
plotMDS(y,col=col.cell)
legend("bottomright",fill=c("#1B9E77", "#7570B3","orangered3"),legend=levels(as.factor(sampTabGood$description1)))
title("multidimensional scaling plot")

#try a different dimensions
#still the same
plotMDS(y,dim.plot = c(1,3),col=col.cell)
legend("bottomright",fill=c("#1B9E77", "#7570B3","orangered3"),legend=levels(as.factor(sampTabGood$description1)))
title("multidimensional scaling plot")

#estimate variable genes
var_genes <- apply(logcounts, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
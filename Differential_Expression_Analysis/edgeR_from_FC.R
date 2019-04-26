rm(list=ls())

library(edgeR)
library(gplots)
library(RColorBrewer)
library(dplyr)
library(GO.db)
library(org.Mm.eg.db)

#Note: Classically you would expect to normalize for gene length, but we are only concerned with the differnces between samples. Gene length is expected to have the same effect on read counts across samples. 

#Note: We should go back and make sure of these results, weirdly the library size for sample6-3 is 1/3 the other library sizes

#Read in featureCounts input and modify it for use with edgeR
fc <- read.table('/Users/MacProMatt/Desktop/mouse/sorted/features_count_all/total_file_seq.count', header = T)
row.names(fc) <- fc$Geneid
fc_min <- subset(fc, select = -c(Geneid,Chr,Start,End,Strand,Length))
group <- c("naive", "naive", "naive", "mCherry Neg", "mCherry Low", "mCherry High", "mCherry Neg", "mCherry Low", "mCherry High", "mCherry Neg", "mCherry Low", "mCherry High")
dge <- DGEList(counts = fc_min, group = group)

#Filter out genes that are lowly expressed
keep <- rowSums(cpm(dge)>1) >= 2
dge <- dge[keep, , keep.lib.sizes=FALSE]

#adjust for RNA composition levels
dge <- calcNormFactors(dge)

#Estimate the dispersion within the samples
dge <- estimateDisp(dge)

#Output to a csv file (this one is for low versus high)
et <- exactTest(dge, pair=c("mCherry Low","mCherry High"))
csv_out <- topTags(et, n=Inf)
write.csv(csv_out, file = '/Users/MacProMatt/Desktop/edgeR_out/low_versus_high.csv')

#high vs naive
et <- exactTest(dge, pair=c('naive','mCherry High'))
csv_out <- topTags(et, n=Inf)
write.csv(csv_out, file = '/Users/MacProMatt/Desktop/edgeR_out/high_versus_naive.csv')

#low vs naive
et <- exactTest(dge, pair=c('naive','mCherry Low'))
csv_out <- topTags(et, n=Inf)
write.csv(csv_out, file = '/Users/MacProMatt/Desktop/edgeR_out/low_versus_naive.csv')

#high vs negative
et <- exactTest(dge, pair=c('mCherry Neg','mCherry High'))
csv_out <- topTags(et, n=Inf)
write.csv(csv_out, file = '/Users/MacProMatt/Desktop/edgeR_out/high_versus_negative.csv')

#low vs negative
et <- exactTest(dge, pair=c('mCherry Neg','mCherry Low'))
csv_out <- topTags(et, n=Inf)
write.csv(csv_out, file = '/Users/MacProMatt/Desktop/edgeR_out/low_versus_negative.csv')


#Merge in gene identifiers from ensembl gene id
#Table generated from http://www.informatics.jax.org/faq/GM_batch.shtml
#Need to save the above table as a csv rather than excel file
#This works by changing the name of table1 column one to 'Input'
files <- list.files(path="/Users/MacProMatt/Desktop/edgeR_out", pattern="*.csv", full.names=TRUE, recursive=FALSE)
gene_table = read.csv('/Users/MacProMatt/Desktop/Langlois_Lab/MGIBatchReport_20190226_114333.csv')
for (i in files){
  table1 = read.csv(i)
  colnames(table1) = c('Input', 'logFC', 'logCPM', 'PValue', 'FDR')
  new_out = merge(table1, gene_table)
  new_out$Symbol = toupper(new_out$Symbol) #Switch gene column to all caps
  write.csv(new_out, file = i)
}

pch = c(1,1,1,15,16,17,15,16,17,15,16,17)
col = c(1,1,1,1,12,26,1,12,26,1,12,26)

all = plotMDS(dge, top = Inf, pch = pch, col = col, ylab = "logFC dim 2", xlab = 'logFC dim 1')
legend("topleft", legend= c("naive", "mCherry Neg", "mCherry Low", "mCherry High"), pch= c(1,15,16,17), col= c(1,1,12,26), ncol = 2, pt.cex = 1, cex = 0.8)

jpeg(filename = '/Users/MacProMatt/Desktop/all_grey.jpeg', height = 4, width = 5, units = "in", res = 800, pointsize = 12)
top_500 = plotMDS(dge, top = 500, pch = pch, col = col, xlab = 'logFC Dim1', ylab = 'logFC Dim2', cex = 1.2)
legend("topleft", legend= c("naive", "mCherry Neg", "mCherry Low", "mCherry High"), pch= c(1,15,16,17), col= c(1,1,12,26), ncol = 2, pt.cex = 1, cex = 0.8)
dev.off()

#Get output for making a heatmap
logcpm <- cpm(dge, log = T)
write.csv(logcpm, file = "/Users/MacProMatt/Desktop/logcpm_out.csv")

#This is probably not the right way to filter out some genes
ISGList <- read.csv("/Users/MacProMatt/Desktop/Langlois_Lab/ISG_list_R.csv")
in_file = read.csv("/Users/MacProMatt/Desktop/logcpm_out.csv")
names(in_file) = c('Gene', 'Sample1.1', 'Sample2.1', 'Sample3.1', 'Sample4.1', 'Sample4.2', 'Sample4.3', 'Sample5.1', 'Sample5.2', 'Sample5.3', 'Sample6.1', 'Sample6.2', 'Sample6.3')
ISG_dat = filter(in_file, Gene %in% ISGList$genes)
rnames = ISG_dat[,1]
ISG_dat = as.matrix(ISG_dat)
mat_ISG_dat = data.matrix(ISG_dat[,2:ncol(ISG_dat)])
rownames(mat_ISG_dat) <- rnames

mode(mat_ISG_dat) = "numeric"

##Re-order rows/columns by mean, use 1-Pearson's correlation distance, and complete linkage
jpeg(filename = '/Users/MacProMatt/Desktop/all_grey.jpeg', height = 4, width = 5, units = "in", res = 800, pointsize = 12)
my_plot = heatmap.2(mat_ISG_dat, scale='row', density.info="none", trace='none', col = colorRampPalette(c( 'blue', 'orange'))(4), reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean), distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x) hclust(x, method="complete"), labCol = F, labRow = FALSE, ColSideColors = c('black','black', 'black', 'grey','blue','red','grey','blue','red','grey','blue','red'), main = "ISG Expression 12 hpi")
dev.off()



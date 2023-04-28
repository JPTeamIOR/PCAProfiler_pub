#!/usr/bin/Rscript

neededlibraries <- c("jsonlite","plotly", "edgeR", "tximport","DESeq2", "ggplot2", "slingshot", "sva", "tidymodels", "biomaRt",  "BiocParallel", "dplyr", "slingshot", "SingleCellExperiment",  "RColorBrewer", "GSVA" )
loaded_libraries = lapply(neededlibraries, require, character.only = TRUE)

source("./00_functions.R")
pcaPlot <- function(indata=NULL,  colorRef=NULL,  ngenes=500, pca=NULL){
  if ( is.null(pca) ){
    rv_all= rowVars(indata)
    pcagenes_all = rownames(indata)[order(rv_all, decreasing = TRUE)]
    selected_all = pcagenes_all[1:ngenes]
    pca = prcomp( t(indata[selected_all,]))  
  }
  percentVar <- setNames(object = round(100* pca$sdev^2/sum(pca$sdev^2)), nm = colnames(pca$x))
  print(ggplot(data.frame(PC1=pca$x[,1], PC2=pca$x[,2], color = colorRef ))+geom_point(aes(x=PC1, y=PC2, color=color))+ylab(paste("PC2",percentVar[2],"%"))+xlab(paste("PC1", percentVar[1],"%")))
  return(pca)
}
mdata = read.csv("./metadata/metadata.tsv", sep="\t")
nrow(mdata)
mdata= mdata[mdata$LibraryType != "Hybrid Capture" , ]


tx2gene = read.table("./metadata/salmon_tx2gene.tsv")
colnames(tx2gene) = c("transcript_id","gene_id", "gene_name")

str(mdata)
base_dir="./export/"
files <- file.path(base_dir, mdata$Source.UID, "quant.sf.gz")
names(files)= mdata$Name
files = files[file.exists(files)]
length(files)
rownames(mdata)=mdata$Name

if ( ! all(file.exists(files)) ) {
  stop("Error! some files don't exist") 
}
dir.create("./objects", showWarnings = F)
if ( file.exists("./objects/txi.gene.RDS") ){
  txi.tx = readRDS("./objects/txi.tx.RDS")
  txi.gene = readRDS("./objects/txi.gene.RDS")
} else {
  txi.tx <- tximport(files, type = "salmon", txOut = TRUE)
  txi.gene <- summarizeToGene(txi.tx, tx2gene)
  saveRDS(txi.gene, "./objects/txi.gene.RDS")
  saveRDS(txi.tx, "./objects/txi.tx.RDS")
}

saveRDS(txi.gene$abundance, "./objects/txi.gene.abundance.RDS")
saveRDS(txi.gene$counts, "./objects/txi.gene.counts.RDS")
saveRDS(txi.gene$length, "./objects/txi.gene.length.RDS")
saveRDS(txi.tx$abundance, "./objects/txi.tx.abundance.RDS")
saveRDS(txi.tx$counts,"./objects/txi.tx.counts.RDS")
saveRDS(txi.tx$length,"./objects/txi.tx.length.RDS")


mdata = mdata[colnames(txi.tx$counts), ]


dir.create("./plots/", showWarnings = F)
if ( file.exists("./objects/ncount.RDS")){
  ncount = readRDS("./objects/ncount.RDS")
} else {
  ncount = log2(sweep(x=txi.gene$counts, MARGIN= 2, STATS = colSums(txi.gene$counts)/1000000 , FUN="/" ) +1 )
  pca=pcaPlot(indata=ncount,  colorRef = mdata$PC.Type, ngenes=2000)
  p<- plot_ly(x=pca$x[,1], y=pca$x[,2], z=pca$x[,3], type="scatter3d", mode="markers", color=mdata$PC.Type)
  p
  htmlwidgets::saveWidget(p, "./plots/GenesPCA.html")
  saveRDS(ncount, "./objects/ncount.RDS")
}
## Removed single end for batch effect ( tried combat, ComBat-seq and limma before )

paired_mask= mdata$LibraryLayout == "PAIRED"
if ( file.exists("./objects/ncount.txi.RDS")){
  ncount.txi = readRDS("./objects/ncount.txi.RDS")
}else {
  ncount.txi = log2(sweep(x=txi.tx$counts, MARGIN= 2, STATS = colSums(txi.tx$counts)/1000000 , FUN="/" ) +1 )
  pca=pcaPlot(ncount.txi, mdata$PC.Type, 2000)
  p<- plot_ly(x=pca$x[,1], y=pca$x[,2], z=pca$x[,3], type="scatter3d", mode="markers", color=mdata$PC.Type)
  p
  htmlwidgets::saveWidget(p, "./plots/Transcripts.html")
  pca=pcaPlot(ncount.txi[,paired_mask], mdata$PC.Type[paired_mask], 2000)
  p<- plot_ly(x=pca$x[,1], y=pca$x[,2], z=pca$x[,3], type="scatter3d", mode="markers", color=mdata$PC.Type[paired_mask])
  p
  htmlwidgets::saveWidget(p, "./plots/Transcripts_nosingle.html")
  
  saveRDS(ncount.txi, "./objects/ncount.txi.RDS")
  gc()
}

pdf("./plots/PCAs.pdf")
pca=pcaPlot(indata=ncount,  colorRef =  mdata$PC.Type, ngenes=2000)
pca=pcaPlot(pca=pca, colorRef = mdata$LibraryType)
pca=pcaPlot(pca=pca, colorRef =  mdata$LibraryLayout)
pca=pcaPlot(pca=pca, colorRef = mdata$Dataset)
pca=pcaPlot(indata= ncount.txi, colorRef =mdata$PC.Type, ngenes=2000)
pca=pcaPlot(pca=pca,  colorRef =mdata$Dataset)
dev.off()


pca = pcaPlot(ncount, mdata$PC.Type,2000)
pca$sdev = NULL
pca$x = NULL
tmp = predict(pca, t(ncount))
saveRDS(pca, "./objects/pca.RDS")

pca = pcaPlot(ncount, mdata$PC.Type,2000)


set.seed(123)
keep <- rowSums(txi.gene$counts > 3) >= (ncol(txi.gene$counts)*0.05)
## start pseudoptime calculation

sim2d <- SingleCellExperiment(assays = List(counts = as.matrix(txi.gene$counts[keep,]) ))

assays(sim2d)$norm <- as.matrix(ncount[keep,])
reducedDims(sim2d) <- SimpleList(PCA = pca$x[,1:2])

cl1= kmeans(pca$x[,1:2] , centers = 3)$cluster
rd1 = pca$x[,1:2]

plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

mdata$clus= cl1
pdf("./plots/clusters.pdf")
ggplot(mdata) + geom_bar(aes(x=clus, fill=PC.Type))
dev.off()
colData(sim2d)$kmeans <- cl1
sim2d_2 <- slingshot(sim2d, reducedDim = 'PCA',clusterLabels = "kmeans",allow.breaks=F, start.clus=3 , end.clus = 2)

pdf("./plots/slingshot.pdf")
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sim2d_2$slingPseudotime_1, breaks=100)]
plot(pca$x[,1:2], col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sim2d_2), lwd=2,  col='black')
ggplot(data.frame(x=sim2d_2$slingPseudotime_1, y=mdata$PC.Type)) + geom_point(aes(x=x, y=y))
ggplot(data.frame(x=sim2d_2$slingPseudotime_1, y=mdata$PC.Type))+geom_histogram(aes(x=x, fill=y))
dev.off()
curva = slingCurves(sim2d_2)[[1]]$s
ordine = slingCurves(sim2d_2)[[1]]$ord
curve_i <- curva[ordine, ]
colnames(curve_i) <- c("PCX", "PCY")
curve_i = as.data.frame(curve_i)

write(toJSON(list(curve_i$PCX, curve_i$PCY)), "./objetcs/slingshot_trajectory_curve.json" )

saveRDS(curve_i,"./objects/curve_i.RDS")

### To predict pseudotime for new data: predict(sds, pca$x)
sds = SlingshotDataSet(sim2d_2)
saveRDS(sim2d_2, "./objects/sim2d.RDS")
saveRDS(sds,"./objects/slingshot.RDS")
tmp = slingshot::predict(sds, pca$x[,1:2])
pseudotimes = slingPseudotime(tmp)[,1]

mdata$PC1 = pca$x[,1]
mdata$PC2 = pca$x[,2]
mdata$pseudotime = sim2d_2$slingPseudotime_1

ggplot(mdata)+geom_point(aes(x=pseudotime, y=PC.Type))

saveRDS(mdata , "./objects/mdata.RDS")


txi_ids = data.frame(ensembl = rownames(ncount.txi), pvalue = 1 , corr = 0 , stat=0) 

for ( i in 1:nrow(txi_ids)){
  ct = cor.test(as.numeric(ncount.txi[i,paired_mask]), mdata$pseudotime[paired_mask])
  txi_ids$corr[i]=ct$estimate
  txi_ids$pvalue[i]=ct$p.value
  txi_ids$stat[i]=ct$statistic
}
txi_ids$pvalue[is.na(txi_ids$pvalue)]=1
txi_ids$stat[is.na(txi_ids$stat)]=0
txi_ids$corr[is.na(txi_ids$corr)]=0
txi_ids$padj = p.adjust(txi_ids$pvalue, method = "fdr")
saveRDS(txi_ids, "./objetcs/tx.corr.RDS")

g_ids = data.frame(ensembl = rownames(ncount), pvalue = 1 , corr = 0 , stat = 0) 
for ( i in 1:nrow(g_ids)){
  ct = cor.test(as.numeric(ncount[i,paired_mask]), mdata$pseudotime[paired_mask])
  g_ids$corr[i]=ct$estimate
  g_ids$pvalue[i]=ct$p.value
  g_ids$stat[i]=ct$statistic
}
g_ids$pvalue[is.na(g_ids$pvalue)]=1
g_ids$corr[is.na(g_ids$corr)]=0
g_ids$padj = p.adjust(g_ids$pvalue, method = "fdr")
saveRDS(g_ids, "./objetcs/gene.corr.RDS")

dat = g_ids
dat$mlogpadj = round(-log10(dat$padj), digits=4)
dat$mlogpadj[is.na(dat$mlogpadj)]=1
dat$mlogpadj[is.infinite(dat$mlogpadj)]=max(dat$mlogpadj[is.finite(dat$mlogpadj)])+1
dat$stat[is.na(dat$stat)]=0
dat$stat = round(dat$stat, digits=4)
dat$corr = round(dat$corr, digits=4)
write.csv(dat[,c("ensembl","corr","stat","mlogpadj","padj", "pvalue")], "./objects/gene.corr.csv", row.names=F, quote=F)

dat = txi_ids
dat$mlogpadj = round(-log10(dat$padj), digits=4)
dat$mlogpadj[is.na(dat$mlogpadj)]=1
dat$mlogpadj[is.infinite(dat$mlogpadj)]=max(dat$mlogpadj[is.finite(dat$mlogpadj)])
dat$stat[is.na(dat$stat)]=0
dat$stat = round(dat$stat, digits=4)
dat$corr = round(dat$corr, digits=4)
write.csv(dat[,c("ensembl","corr","stat","mlogpadj","padj", "pvalue")], "./objects/tx.corr.csv", row.names=F, quote=F)


#!/usr/bin/Rscript

neededlibraries <- c("jsonlite","plotly", "edgeR", "tximport","DESeq2", "ggplot2", "slingshot", "sva",  "tximport",  "BiocParallel", "dplyr", "slingshot", "SingleCellExperiment",  "RColorBrewer", "GSVA" )
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
mdata = read.csv("./metadata/metadata_final.tsv", sep="\t")
tx2gene = read.table("./metadata/salmon_tx2gene.tsv")
colnames(tx2gene) = c("transcript_id","gene_id", "gene_name")
outdir="./results_03_PCA/"
dir.create(outdir, showWarnings = F)

dir.create("./objects", showWarnings = F)
if ( file.exists("./objects/txi.gene.RDS") ){
  txi.tx = readRDS("./objects/txi.tx.RDS")
  txi.gene = readRDS("./objects/txi.gene.RDS")
  ncount = readRDS("./objects/ncount.RDS")
  ncount.tx = readRDS("./objects/ncount.tx.RDS")
  mdata = readRDS("./objects/mdata.RDS")
} else {
  base_dir="./export/"
  mdata$file = file.path(base_dir, mdata$Source.UID, "quant.sf.gz")
  mdata = mdata[file.exists(mdata$file),]
  summary(factor(mdata$LibraryType))
  files <- file.path(base_dir, mdata$Source.UID, "quant.sf.gz")
  names(files)= mdata$Name
  sum(file.exists(files)) == nrow(mdata)
  rownames(mdata)=mdata$Name
  mdata$file = NULL
  txi.tx <- tximport(files, type = "salmon", txOut = TRUE)
  txi.gene <- summarizeToGene(txi.tx, tx2gene)
  ncount = log2(sweep(x=txi.gene$counts, MARGIN= 2, STATS = colSums(txi.gene$counts)/1000000 , FUN="/" ) +1 )
  ncount.tx = log2(sweep(x=txi.tx$counts, MARGIN= 2, STATS = colSums(txi.tx$counts)/1000000 , FUN="/" ) +1 )
  ### Batch effect for Hybrid-seq
  
  pca=pcaPlot(indata=ncount,  colorRef = mdata$LibraryType, ngenes=2000)
  tmp = data.frame(PC1= pca$x[,1], PC2=pca$x[,2], PC3 = pca$x[,3] ,LibraryType = mdata$LibraryType)
  tmp$LibraryType = factor(tmp$LibraryType, levels=c("TotalRNA", "PolyA", "Hybrid Capture"))
  p<- plot_ly(data =  tmp, x=~PC1, y=~PC2, z=~PC3, type="scatter3d", mode="markers", color=~LibraryType, colors= c("#008C8F", "#99EDCC", "#CB958E"))
  p
  htmlwidgets::saveWidget(p, paste0(outdir, "/HybridCapturePCA.html"))
  pdf(paste0(outdir, "/HybridCapturePCA.pdf"), width = 8, height = 6)
    percentVar <- setNames(object = round(100* pca$sdev^2/sum(pca$sdev^2)), nm = colnames(pca$x))
    ggplot(tmp)+geom_point(aes(x=-PC1, y=PC2, fill=LibraryType),shape= 21, size=3)+ylab(paste("PC2",percentVar[2],"%"))+xlab(paste("PC1", percentVar[1],"%"))+
      theme_classic()+scale_fill_manual(values = c("#008C8F", "#99EDCC", "#CB958E"))
    ggplot(tmp)+geom_point(aes(x=-PC1, y=PC3, fill=LibraryType),shape= 21, size=3)+ylab(paste("PC3",percentVar[3],"%"))+xlab(paste("PC1", percentVar[1],"%"))+
      theme_classic()+scale_fill_manual(values = c("#008C8F", "#99EDCC", "#CB958E"))
    ggplot(tmp)+geom_point(aes(x=PC2, y=PC3, fill=LibraryType),shape= 21, size=3)+ylab(paste("PC3",percentVar[3],"%"))+xlab(paste("PC2", percentVar[2],"%"))+
      theme_classic()+scale_fill_manual(values = c("#008C8F", "#99EDCC", "#CB958E"))
  dev.off()
  
  
  ### remove the hybrid capture
  mdata= mdata[mdata$LibraryType != "Hybrid Capture",]
  colnames(ncount)
  ncount = ncount[,mdata$Name]
  ncount.tx = ncount.tx[,mdata$Name]
  for ( name in c("abundance", "counts", "length")){
    txi.gene[[name]] = txi.gene[[name]][,mdata$Name]
    txi.tx[[name]] = txi.tx[[name]][,mdata$Name]
  }
  
  pca=pcaPlot(indata=ncount,  colorRef = mdata$LibraryType, ngenes=2000)
  tmp = data.frame(PC1= pca$x[,1], PC2=pca$x[,2], PC3 = pca$x[,3] ,LibraryType = mdata$LibraryType, LibraryLayout = mdata$LibraryLayout)
  tmp$LibraryType = factor(tmp$LibraryType, levels=c("TotalRNA", "PolyA"))
  tmp$LibraryLayout = factor(tmp$LibraryLayout, levels=c("SINGLE", "PAIRED"))
  p<- plot_ly(data =  tmp, x=~PC1, y=~PC2, z=~PC3, type="scatter3d", mode="markers", color=~LibraryType, colors= c("#008C8F", "#99EDCC"))
  p
  htmlwidgets::saveWidget(p, paste0(outdir, "/NoHybridCapturePCA.html"))
  p<- plot_ly(data =  tmp, x=~PC1, y=~PC2, z=~PC3, type="scatter3d", mode="markers", color=~LibraryLayout, colors= c("#EE4266", "#2E5077"))
  p
  htmlwidgets::saveWidget(p, paste0(outdir, "/NoHybridCapturePCA_Layout.html"))
  pdf(paste0(outdir, "/NoHybridCapturePCA.pdf"), width = 8, height = 6)
  percentVar <- setNames(object = round(100* pca$sdev^2/sum(pca$sdev^2)), nm = colnames(pca$x))
  ggplot(tmp)+geom_point(aes(x=-PC1, y=PC2, fill=LibraryType),shape= 21, size=3)+ylab(paste("PC2",percentVar[2],"%"))+xlab(paste("PC1", percentVar[1],"%"))+
    theme_classic()+scale_fill_manual(values = c("#008C8F", "#99EDCC", "#CB958E"))
  ggplot(tmp)+geom_point(aes(x=-PC1, y=PC3, fill=LibraryType),shape= 21, size=3)+ylab(paste("PC3",percentVar[3],"%"))+xlab(paste("PC1", percentVar[1],"%"))+
    theme_classic()+scale_fill_manual(values = c("#008C8F", "#99EDCC", "#CB958E"))
  ggplot(tmp)+geom_point(aes(x=PC2, y=PC3, fill=LibraryType),shape= 21, size=3)+ylab(paste("PC3",percentVar[3],"%"))+xlab(paste("PC2", percentVar[2],"%"))+
    theme_classic()+scale_fill_manual(values = c("#008C8F", "#99EDCC", "#CB958E"))
  ggplot(tmp)+geom_point(aes(x=-PC1, y=PC2, fill=LibraryLayout),shape= 21, size=3)+ylab(paste("PC2",percentVar[2],"%"))+xlab(paste("PC1", percentVar[1],"%"))+
    theme_classic()+scale_fill_manual(values = c("#EE4266", "#2E5077"))
  dev.off()
  
  pca=pcaPlot(ncount.tx, mdata$PC.Type, 2000)
  tmp = data.frame(PC1= pca$x[,1], PC2=pca$x[,2], PC3 = pca$x[,3] ,LibraryLayout = mdata$LibraryLayout, LibraryType = mdata$LibraryType)
  tmp$LibraryLayout = factor(tmp$LibraryLayout, levels=c("SINGLE", "PAIRED"))
  p<- plot_ly(data =  tmp, x=~PC1, y=~PC2, z=~PC3, type="scatter3d", mode="markers", color=~LibraryLayout, colors= c("#EE4266", "#2E5077"))
  p
  htmlwidgets::saveWidget(p, paste0(outdir, "/TX_LayoutPCA.html"))
  p<- plot_ly(data =  tmp, x=~PC1, y=~PC2, z=~PC3, type="scatter3d", mode="markers", color=~LibraryType, colors=  c("#008C8F", "#99EDCC"))
  p
  htmlwidgets::saveWidget(p, paste0(outdir, "/TX_TypePCA.html"))
  pdf(paste0(outdir, "/TX_LayoutPCA.pdf"), width = 8, height = 6)
  percentVar <- setNames(object = round(100* pca$sdev^2/sum(pca$sdev^2)), nm = colnames(pca$x))
  ggplot(tmp)+geom_point(aes(x=-PC1, y=PC2, fill=LibraryLayout),shape= 21, size=3)+ylab(paste("PC2",percentVar[2],"%"))+xlab(paste("PC1", percentVar[1],"%"))+
    theme_classic()+scale_fill_manual(values =  c("#EE4266", "#2E5077"))
  ggplot(tmp)+geom_point(aes(x=-PC1, y=PC3, fill=LibraryLayout),shape= 21, size=3)+ylab(paste("PC3",percentVar[3],"%"))+xlab(paste("PC1", percentVar[1],"%"))+
    theme_classic()+scale_fill_manual(values =  c("#EE4266", "#2E5077"))
  ggplot(tmp)+geom_point(aes(x=PC2, y=PC3, fill=LibraryLayout),shape= 21, size=3)+ylab(paste("PC3",percentVar[3],"%"))+xlab(paste("PC2", percentVar[2],"%"))+
    theme_classic()+scale_fill_manual(values =  c("#EE4266", "#2E5077"))
  ggplot(tmp)+geom_point(aes(x=-PC1, y=PC2, fill=LibraryType),shape= 21, size=3)+ylab(paste("PC2",percentVar[2],"%"))+xlab(paste("PC1", percentVar[1],"%"))+
    theme_classic()+scale_fill_manual(values = c("#008C8F", "#99EDCC"))
  ggplot(tmp)+geom_point(aes(x=-PC1, y=PC3, fill=LibraryType),shape= 21, size=3)+ylab(paste("PC3",percentVar[3],"%"))+xlab(paste("PC1", percentVar[1],"%"))+
    theme_classic()+scale_fill_manual(values = c("#008C8F", "#99EDCC"))
  ggplot(tmp)+geom_point(aes(x=PC2, y=PC3, fill=LibraryType),shape= 21, size=3)+ylab(paste("PC3",percentVar[3],"%"))+xlab(paste("PC2", percentVar[2],"%"))+
    theme_classic()+scale_fill_manual(values = c("#008C8F", "#99EDCC"))
  dev.off()
  
  saveRDS(ncount, "./objects/ncount.RDS")
  saveRDS(ncount.tx, "./objects/ncount.tx.RDS")
  saveRDS(txi.gene, "./objects/txi.gene.RDS")
  saveRDS(txi.tx, "./objects/txi.tx.RDS")
  saveRDS(mdata, "./objects/mdata.RDS")
  
  saveRDS(txi.gene$abundance, "./objects/txi.gene.abundance.RDS")
  saveRDS(txi.gene$counts, "./objects/txi.gene.counts.RDS")
  saveRDS(txi.gene$length, "./objects/txi.gene.length.RDS")
  saveRDS(txi.tx$abundance, "./objects/txi.tx.abundance.RDS")
  saveRDS(txi.tx$counts,"./objects/txi.tx.counts.RDS")
  saveRDS(txi.tx$length,"./objects/txi.tx.length.RDS")
}





pca = pcaPlot(ncount, mdata$PC.Type,2000)

pca$sdev = NULL
pca$x = NULL
tmp = predict(pca, t(ncount))
saveRDS(pca, "./objects/pca.RDS")

pca = pcaPlot(ncount, mdata$PC.Type,2000)

pca$x[,1] = -pca$x[,1]

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
pdf(paste0(outdir,"/clusters.pdf"))
ggplot(mdata) + geom_bar(aes(x=clus, fill=PC.Type))
dev.off()
colData(sim2d)$kmeans <- cl1
sim2d_2 <- slingshot(sim2d, reducedDim = 'PCA',clusterLabels = "kmeans",allow.breaks=F, start.clus=3 , end.clus = 2)

pdf(paste0(outdir,"/slingshot.pdf"))
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

write(toJSON(list(curve_i$PCX, curve_i$PCY)), "./objects/slingshot_trajectory_curve.json" )

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

paired_mask = mdata$LibraryLayout == "PAIRED"
txi_ids = data.frame(ensembl = rownames(ncount.tx), pvalue = 1 , corr = 0 , stat=0) 

for ( i in 1:nrow(txi_ids)){
  ct = cor.test(as.numeric(ncount.tx[i,paired_mask]), mdata$pseudotime[paired_mask])
  txi_ids$corr[i]=ct$estimate
  txi_ids$pvalue[i]=ct$p.value
  txi_ids$stat[i]=ct$statistic
}
txi_ids$pvalue[is.na(txi_ids$pvalue)]=1
txi_ids$stat[is.na(txi_ids$stat)]=0
txi_ids$corr[is.na(txi_ids$corr)]=0
txi_ids$padj = p.adjust(txi_ids$pvalue, method = "fdr")
saveRDS(txi_ids, "./objects/tx.corr.RDS")

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
saveRDS(g_ids, "./objects/gene.corr.RDS")

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


### Compare the current pseudotimes with the ones of the previous publication

mdata = readRDS("./objects/mdata.RDS")
prev_ps = read.csv("./metadata/prev_pseudotimes.csv")
prev_md = read.table("./metadata/prev_cohort_annot.tsv", sep="\t", header = T)
nrow(prev_md)
nrow(prev_ps)
prev_ps = merge(prev_ps, prev_md[,c("BASENAME", "GROUP","SAMPLE_ID")], by.x="name", by.y="BASENAME")

mdata$common_name = mdata$Source.UID
mdata$common_name[mdata$Dataset == "TCGA-PRAD"] = gsub("[A|B]$", "",mdata$Sample[mdata$Dataset == "TCGA-PRAD"])
prev_ps$common_name = prev_ps$name
prev_ps$common_name[prev_ps$SAMPLE_ID %in% mdata$common_name] = prev_ps$SAMPLE_ID[prev_ps$SAMPLE_ID %in% mdata$common_name]
colnames(prev_ps)[1:2] = c("name", "prev_pseudotime")
comp_ps = merge(prev_ps, mdata[,c("common_name", "pseudotime", "PC1", "PC2")], by="common_name", all.x = T, all.y=F)
comp_ps = comp_ps[!is.na(comp_ps$pseudotime),]
comp_ps$GROUP = factor(comp_ps$GROUP, levels=c("NORMAL", "PRIMARY", "CRPC", "NEPC"))

pdf(paste0(outdir, "/previous_pseudotime.pdf"), width = 8, height = 6)
ggplot(comp_ps) + geom_point(aes(x=prev_pseudotime, y=pseudotime, fill=GROUP), size=3, shape=21)+ 
  theme_classic() + scale_fill_manual(values=c("#25681E", "#F63A21","#4889C7",  "#93355A")) + 
  ggtitle(paste0("PCC= ", cor(comp_ps$prev_pseudotime, comp_ps$pseudotime)))

ggplot(comp_ps, aes(x=PC1, y=PC2, fill=GROUP)) + geom_point(alpha=0.2, size=3, shape=21)+ 
  geom_point(data = comp_ps[abs(comp_ps$prev_pseudotime - comp_ps$pseudotime) > 100,], size=3, shape=21, alpha=1)+
  theme_classic() + scale_fill_manual(values=c("#25681E", "#F63A21","#4889C7",  "#93355A"))
dev.off()




#!/usr/bin/Rscript
library("ggplot2")
library(cluster)
library(clusterProfiler)
library("enrichplot")
library(org.Hs.eg.db)
library(msigdbr)

source("./00_functions.R")
ncount = readRDS("./objects/ncount.RDS")
ncount.txi = readRDS("./objects/ncount.txi.RDS")
mdata = readRDS("./objects/mdata.RDS")
curve_i = readRDS("./objects/curve_i.RDS")
tx2gene = read.table("./metadata/salmon_tx2gene.tsv")
txi.gene  = readRDS("./objects/txi.gene.RDS")
colnames(tx2gene) = c("transcript_id", "gene_id", "gene_name")
tx_info = read.csv("./metadata//transcripts_info.csv")
tx_info = merge(tx_info, tx2gene, by=c("transcript_id","gene_id"))
nrow(mdata)
### Signatures
g_info = data.frame(gene_name = tx_info$gene_name, gene_id = tx_info$gene_id)
nrow(g_info)
g_info = g_info[!duplicated(g_info$gene_id),]
nrow(g_info)

keep <- rowSums(txi.gene$counts > 3) >= (ncol(txi.gene$counts)*0.05)
saveRDS(keep, "./objects/keep.RDS")
CS <- list()
# AR-SCORE NELSON (Bluemn et al., 2017) ---->
CS$AR_GSET <- c("KLK3","KLK2","TMPRSS2","FKBP5","NKX3-1","PLPP1","PMEPA1","PART1","ALDH1A3","STEAP4")
CS$AR_GSET %in% tx_info$gene_name
# NE-SCORE NELSON (Bluemn et al., 2017) ---->
CS$NE_GSET <- c("SYP","CHGA","CHGB","ENO2","CHRNB2","SCG3","SCN3A","PCSK1","ELAVL4","NKX2-1")
CS$NE_GSET %in% tx_info$gene_name


for(cat in names(CS)){
  CS[[cat]] = unique(tx_info$gene_id[tx_info$gene_name %in% CS[[cat]] ])
}

nelson.ssgsea <- gsva(txi.gene$abundance[keep,], CS, method = "ssgsea", parallel.sz=4,  mx.diff=TRUE, ssgsea.norm	=F)

tmp = as.data.frame(t(nelson.ssgsea))
tmp$Name = rownames(tmp)
mdata = merge(mdata, tmp, by="Name")
rownames(mdata) = mdata$Name
mdata = mdata[colnames(ncount),]
pca = pcaPlot(ncount, mdata$PC.Type, 2000, 1, 2)
mdata$PC1 = pca$x[,1]
mdata$PC2 = pca$x[,2]
mdata$PC3 = pca$x[,3]
mdata$AR = ncount[g_info$gene_id[g_info$gene_name == "AR"],]
dir.create("./SignaturesGraphs", showWarnings = F)

pdf("./SignaturesGraphs/AR_NE_scores.pdf", width = 10)
ggplot(mdata)+ geom_point(aes(x=PC1, y=PC2, color=NE_GSET))
ggplot(mdata)+ geom_point(aes(x=PC1, y=PC2, color=AR_GSET))
ggplot(mdata)+ geom_point(aes(x=PC1, y=PC2, color=AR))

NE_thr=max(mdata$NE_GSET[mdata$PC.Type != "CRPC"])
ggplot(mdata)+ geom_point(aes(x=PC.Type, y=NE_GSET), position = position_jitter(0.1))+geom_hline(yintercept = NE_thr)
ggplot(mdata)+ geom_point(aes(x=PC.Type, y=AR_GSET), position = position_jitter(0.1))
ggplot(mdata)+ geom_point(aes(x=PC.Type, y=AR), position = position_jitter(0.1))
ggplot(mdata)+ geom_point(aes(x=AR_GSET, y=AR, color=PC.Type))

myScale = function(x){
  ((2*(x-min(x)) )/(max(x)-min(x)))-1
}

tmp = mdata[mdata$PC.Type == "CRPC", ]
for ( g in c("KLK3","KLK2","TMPRSS2","FKBP5","NKX3-1","PLPP1","PMEPA1","PART1","ALDH1A3","STEAP4", "SYP","CHGA","CHGB","ENO2","CHRNB2","SCG3","SCN3A","PCSK1","ELAVL4","NKX2-1") ){
  if ( sum(g_info$gene_name == g ) == 1 ){
    tmp[,g]=ncount[g_info$gene_id[g_info$gene_name == g], tmp$Name]
  } else {
    tmp[,g]=colMeans(ncount[g_info$gene_id[g_info$gene_name == g], tmp$Name])
  }
  tmp[,g]=myScale(tmp[,g])
}
tmp$AR = myScale(tmp$AR)
tmp$NE_GSET = myScale(tmp$NE_GSET)
tmp$AR_GSET = myScale(tmp$AR_GSET)

library(ComplexHeatmap)
tmp = tmp[order(tmp$AR_GSET, decreasing = T),]
tmp$Subtype = "ARPC"
tmp$Subtype[tmp$NE_GSET > 0.5]="NEPC"
tmp$Subtype[tmp$NE_GSET < 0.5 & tmp$AR_GSET < 0] = "DNPC"
sum(tmp$NE_GSET < 0.5 & tmp$AR_GSET < 0)
tmp$Subtype = factor(tmp$Subtype, levels=c("ARPC",  "DNPC","NEPC"))
tmp = tmp[order(tmp$Subtype),]
rownames(tmp) = NULL

Heatmap(top_annotation = HeatmapAnnotation(Subtype = tmp$Subtype, col=list(Subtype=c("ARPC"= "darkgreen", "NEPC" = "orange", "DNPC"= "blue")) ),
  t(tmp[,c( "AR", "AR_GSET", "KLK3","KLK2","TMPRSS2","FKBP5","NKX3-1","PLPP1","PMEPA1","PART1","ALDH1A3","STEAP4", "NE_GSET","SYP","CHGA","CHGB","ENO2","CHRNB2","SCG3","SCN3A","PCSK1","ELAVL4","NKX2-1"  )]),
        row_split = c( "AR", "AR score", rep("AR score genes", 10), "NE score", rep("NE score genes", 10)),
        column_split = tmp$Subtype,
        cluster_rows = F, cluster_columns = F, use_raster = F, 
)


mdata$Subtype=as.character(mdata$PC.Type)
mdata$Subtype[mdata$Name %in% tmp$Name[tmp$Subtype == "ARPC"]]="ARPC"
mdata$Subtype[mdata$Name %in% tmp$Name[tmp$Subtype == "DNPC"]]="DNPC"
mdata$Subtype[mdata$Name %in% tmp$Name[tmp$Subtype == "NEPC"]]="NEPC"
ggplot(mdata)+ geom_point(aes(x=NE_GSET, y=AR_GSET, color=Subtype))


ggplot(mdata)+ geom_point(aes(y=AR_GSET, x=NE_GSET, color=Subtype))
ggplot(mdata)+ geom_point(aes(y=AR, x=NE_GSET, color=Subtype))

ggplot(mdata)+ geom_point(aes(x=PC1, y=PC2, color=Subtype))
dev.off()


### Estimate score
library(estimate)

out_f = "./estimate.gct"
library(org.Hs.eg.db)
hs <- org.Hs.eg.db

mapper = biomaRt::select(hs, 
                keys = g_info$gene_name,
                columns = c("ENTREZID", "SYMBOL"),
                keytype = "SYMBOL")
colnames(mapper) = c("gene_name", "entrez")
mapper = merge(mapper, g_info, by="gene_name")
mapper = mapper[!is.na(mapper$entrez),]
nrow(mapper)
tmp = matrix(0, nrow=length(unique(mapper$entrez)), ncol = ncol(ncount))
rownames(tmp)=unique(mapper$entrez)
colnames(tmp) = colnames(ncount)

for (eid in unique(mapper$entrez) ){
  if ( sum(mapper$entrez == eid) == 1 ){
    tmp[eid,]=as.numeric(txi.gene$abundance[mapper$gene_id[mapper$entrez == eid],])
  } else {
    tmp[eid,]=colMeans(txi.gene$abundance[mapper$gene_id[mapper$entrez == eid],])  
  }
}


write.table(as.data.frame(tmp), file = "temp.txt", quote = FALSE, col.names = NA, sep = "\t")
filterCommonGenes(input.f="temp.txt", output.f=out_f, id="EntrezID")

estimateScore(out_f, "estimate_score.gct", platform="illumina")

scores <- read.table("estimate_score.gct", sep = "\t", header = T, skip = 2, row.names = 1) 
rownames(scores)
scores = as.data.frame(t(scores[,-1]))

all(mdata$CountColName == rownames(scores))

mdata$ESTIMATE.Stromal.Score = scores$StromalScore
mdata$ESTIMATE.Immune.Score = scores$ImmuneScore
mdata$ESTIMATE.Score = scores$ESTIMATEScore




pdf("./SignaturesGraphs/Estimate.pdf", width = 10)
ggplot(mdata) + geom_point(aes(x=PC1, y=PC2, color=ESTIMATE.Immune.Score))
ggplot(mdata) + geom_point(aes(x=PC1, y=PC2, color=ESTIMATE.Stromal.Score))
ggplot(mdata) + geom_point(aes(x=PC1, y=PC2, color=ESTIMATE.Score))

center = data.frame( x =(min(mdata$PC1) + max(mdata$PC1) )/2, y= (min(mdata$PC2) + max(mdata$PC2) )/2)
center$y=center$y+15
mdata$distances = apply(mdata[,c("PC1","PC2")],1,function(x,cnt) {(sqrt((x["PC1"] - center$x[1])^2+(x["PC2"]-center$y[1])^2))},center)
ggplot(mdata)+ geom_point(aes(x=PC1, y=PC2, fill=distances), colour="black",pch=21, stroke = 0.3, size=3.5)+geom_point(data=center, aes(x=x, y=y), fill="yellow", colour="black",pch=21, stroke = 0.3, size=3.5)+scale_fill_gradient()
ggplot(mdata) + geom_point(aes(x=distances, y=ESTIMATE.Stromal.Score, fill=PC.Type), colour="black",pch=21, stroke = 0.3, size=3.5)+ ggtitle(paste("All, PCC=", cor(mdata$ESTIMATE.Stromal.Score, mdata$distances,use = "pairwise.complete.obs")))
ggplot(mdata) + geom_point(aes(x=distances, y=ESTIMATE.Immune.Score, fill=PC.Type), colour="black",pch=21, stroke = 0.3, size=3.5)+ ggtitle(paste("All, PCC=", cor(mdata$ESTIMATE.Immune.Score, mdata$distances,use = "pairwise.complete.obs")))
ggplot(mdata) + geom_point(aes(x=distances, y=ESTIMATE.Score, fill=PC.Type), colour="black",pch=21, stroke = 0.3, size=3.5)+ ggtitle(paste("All, PCC=", cor(mdata$ESTIMATE.Score, mdata$distances,use = "pairwise.complete.obs")))
mask = mdata$PC.Type == "CRPC" 
ggplot(mdata[mask,]) + geom_point(aes(x=distances, y=ESTIMATE.Stromal.Score, fill=PC.Type), colour="black",pch=21, stroke = 0.3, size=3.5) + ggtitle(paste("CRPC, PCC=", cor(mdata[mask,]$ESTIMATE.Stromal.Score, mdata[mask,]$distances,use = "pairwise.complete.obs")))
ggplot(mdata[mask,]) + geom_point(aes(x=distances, y=ESTIMATE.Immune.Score, fill=PC.Type), colour="black",pch=21, stroke = 0.3, size=3.5)+ ggtitle(paste("CRPC, PCC=", cor(mdata[mask,]$ESTIMATE.Immune.Score, mdata[mask,]$distances,use = "pairwise.complete.obs")))
ggplot(mdata[mask,]) + geom_point(aes(x=distances, y=ESTIMATE.Score, fill=PC.Type), colour="black",pch=21, stroke = 0.3, size=3.5)+ ggtitle(paste("CRPC, PCC=", cor(mdata[mask,]$ESTIMATE.Score, mdata[mask,]$distances,use = "pairwise.complete.obs")))
dev.off()

saveRDS(mdata, "./objects/mdata_signatures.RDS")


###  Plot signatures

m_df = msigdbr(species = "Homo sapiens", category="H")

for ( gsn in unique(m_df$gs_name)){
  target_gene_list = rownames(txi.gene$abundance)[gsub("\\..*" , "", rownames(txi.gene$abundance)) %in% m_df$human_ensembl_gene[m_df$gs_name == gsn] ];
  h.ssgsea <- gsva(txi.gene$abundance[keep,], list(gsn = target_gene_list), method = "ssgsea", parallel.sz=4,  mx.diff=TRUE, ssgsea.norm	=T)
  mdata[,gsn]=h.ssgsea[1,]
}

pdf("./SignaturesGraphs/PCA_individual_Hallmarks.pdf")
for ( gsn in unique(m_df$gs_name)){
  if ( gsn %in% colnames(mdata)){
    tmp = data.frame(PC1=mdata$PC1, PC2=mdata$PC2, h=mdata[,gsn])
    tmp$h = (2*((tmp$h - min(tmp$h)) / (max(tmp$h)- min(tmp$h)))) -1
    p<-ggplot(tmp) + geom_point(aes(x=PC1, y=PC2, color=h))+ggtitle(gsn)+
      scale_color_gradient2(low = "blue", high = "red", mid = "grey")
    print(p)  
  }
}
dev.off()



pdf("./SignaturesGraphs/signatures_PCA.pdf")
for ( gsn in c("AR", "AR_GSET", "NE_GSET", "ESTIMATE.Immune.Score", "ESTIMATE.Stromal.Score", "ESTIMATE.Score")){
  if ( gsn %in% colnames(mdata)){
    tmp = data.frame(PC1=mdata$PC1, PC2=mdata$PC2, h=mdata[,gsn])
    tmp$h = (2*((tmp$h - min(tmp$h)) / (max(tmp$h)- min(tmp$h)))) -1
    p<-ggplot(tmp) + geom_point(aes(x=PC1, y=PC2, color=h))+ggtitle(gsn)+
      scale_color_gradient2(low = "blue", high = "red", mid = "grey")
    print(p)  
  }
}
dev.off()

library(factoextra)

pca = pcaPlot(ncount, mdata$Subtype, ngenes=nrow(ncount[keep,]))
percentVar <- setNames(object = round(100* pca$sdev^2/sum(pca$sdev^2)), nm = colnames(pca$x))
pca.var = get_pca_var(pca)

pca.contrib = pca.var$contrib[,1:3]
pca.contrib = round((pca.contrib / colSums(pca.contrib))*100, 1)
ggplot(as.data.frame(pca.var$coord))+geom_point(aes(x=Dim.1, y=Dim.2))
pca.coord = pca.var$coord[,1:3]
colnames(pca.coord) = c("PC1","PC2","PC3" )

write.csv(pca.contrib, "./objects/GenesContributions.csv")
gene_list = list(PC1 = gsub("\\..*", "" , rownames(pca.contrib)[pca.contrib[,1]>0.01]),
                 PC2 = gsub("\\..*", "" ,rownames(pca.contrib)[pca.contrib[,2]>0.01]),
                 PC3 = gsub("\\..*", "" ,rownames(pca.contrib)[pca.contrib[,3]>0.01])
)
universe=gsub("\\..*", "",rownames(ncount[keep,]))
summary(gene_list)
base_name="PCA_Importances"
base_dir="./GO/"
dir.create(base_dir, recursive = T,showWarnings =F)
categories = c("H")
for ( cat in categories ){
  m_df = msigdbr(species = "Homo sapiens", category=cat)
  em_list = list()
  if ( length(unique(m_df$gs_subcat)) == 1 ){
    for ( k in names(gene_list)){
      em_list[[k]]= enricher(gene_list[[k]],universe =universe, TERM2GENE = m_df[,c("gs_name", "human_ensembl_gene")])
    }
    dotPlotList(em_list, out_file = paste0(base_name, cat), plot_dir = base_dir)  
  } else {
    for ( subcat in unique(m_df$gs_subcat) ){
      for ( k in names(gene_list)){
        em_list[[k]]= enricher(gene_list[[k]],universe = universe, TERM2GENE = m_df[m_df$gs_subcat==subcat,c("gs_name", "human_ensembl_gene")])
      }
      dotPlotList(em_list, out_file = paste0(base_name, cat, "_", subcat), plot_dir = base_dir)  
    }
  }
}
categories=c("H")
base_name="PCA_Importances"
base_dir="./GSEA/"
dir.create(base_dir, recursive = T,showWarnings =F)
for ( cat in categories){
  m_df = msigdbr(species = "Homo sapiens", category=cat)
  for ( subcat in unique(m_df$gs_subcat) ){
    DE_gene_list = list()
    for (contrast in c("PC1", "PC2", "PC3") ){
      stat = pca.coord[!grepl("_PAR_", rownames(pca.coord)),contrast]
      names(stat) = gsub("\\..*", "", names(stat))
      DE_gene_list[[contrast]]=GSEA(sort(stat, decreasing = T),TERM2GENE = m_df[m_df$gs_subcat == subcat,c("gs_name", "human_ensembl_gene")] )
    }
    if ( subcat == "" ){
      dotPlotListGSEA(DE_gene_list, out_file = paste0("ALL_DE_", cat), plot_dir = base_dir)
    } else {
      dotPlotListGSEA(DE_gene_list, out_file = paste0("ALL_DE_", cat, "_" , subcat), plot_dir = base_dir)
    }
  }
}


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

mdata$AR = ncount[g_info$gene_id[g_info$gene_name == "AR"],mdata$Name]
outdir= "./results_04_SignaturesGraphs/"
dir.create(outdir, showWarnings = F)

pdf(paste0(outdir,"/AR_NE_scores.pdf"), width = 10)
ggplot(mdata)+ geom_point(aes(x=PC1, y=PC2, fill=NE_GSET), size=3, shape=21) + theme_classic() + scale_fill_gradient2(low="blue", high="red", midpoint = (max(mdata$NE_GSET) + min(mdata$NE_GSET) )/2 )
ggplot(mdata)+ geom_point(aes(x=PC1, y=PC2, fill=AR_GSET), size=3, shape=21)+ theme_classic() + scale_fill_gradient2(low="blue", high="red", midpoint = (max(mdata$AR_GSET) + min(mdata$AR_GSET) )/2)
ggplot(mdata)+ geom_point(aes(x=PC1, y=PC2, fill=AR), size=3, shape=21)+ theme_classic()+ scale_fill_gradient2(low="blue", high="red", midpoint = (max(mdata$AR) + min(mdata$AR) )/2)

NE_thr=max(mdata$NE_GSET[mdata$PC.Type != "CRPC"])
ggplot(mdata)+ geom_point(aes(x=PC.Type, y=NE_GSET), position = position_jitter(0.1))+geom_hline(yintercept = NE_thr)+ theme_classic()
ggplot(mdata)+ geom_point(aes(x=PC.Type, y=AR_GSET), position = position_jitter(0.1))+ theme_classic()
ggplot(mdata)+ geom_point(aes(x=PC.Type, y=AR), position = position_jitter(0.1))+ theme_classic()
ggplot(mdata)+ geom_point(aes(x=AR_GSET, y=AR, color=PC.Type))+ theme_classic()

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

Heatmap(top_annotation = HeatmapAnnotation(Subtype = tmp$Subtype, col=list(Subtype=c("ARPC"= "#4889C7", "NEPC" = "#93355A", "DNPC"= "#1010C0")) ),
  t(tmp[,c( "AR", "AR_GSET", "KLK3","KLK2","TMPRSS2","FKBP5","NKX3-1","PLPP1","PMEPA1","PART1","ALDH1A3","STEAP4", "NE_GSET","SYP","CHGA","CHGB","ENO2","CHRNB2","SCG3","SCN3A","PCSK1","ELAVL4","NKX2-1"  )]),
        row_split = c( "AR", "AR score", rep("AR score genes", 10), "NE score", rep("NE score genes", 10)),
        column_split = tmp$Subtype,
        cluster_rows = F, cluster_columns = F, use_raster = F, 
)

dev.off()
palette_subtypes = c("#25681E", "#F63A21","#4889C7", "#1010C0", "#93355A")
mdata$Subtype=as.character(mdata$PC.Type)
mdata$Subtype[mdata$Name %in% tmp$Name[tmp$Subtype == "ARPC"]]="ARPC"
mdata$Subtype[mdata$Name %in% tmp$Name[tmp$Subtype == "DNPC"]]="DNPC"
mdata$Subtype[mdata$Name %in% tmp$Name[tmp$Subtype == "NEPC"]]="NEPC"
mdata$Subtype = factor(mdata$Subtype, levels=c("NORMAL", "PRIMARY","ARPC", "DNPC", "NEPC"))


pdf(paste0(outdir, "/Subtypes.pdf"), width = 8, height = 6)
ggplot(mdata)+ geom_point(aes(x=PC1, y=PC2, fill=Subtype), shape=21, size=3) + theme_classic() +
  scale_fill_manual(values = palette_subtypes)
ggplot(mdata)+ geom_point(aes(x=NE_GSET, y=AR_GSET, fill=Subtype), shape=21, size=3) + theme_classic() +
  scale_fill_manual(values = palette_subtypes)
ggplot(mdata)+ geom_point(aes(y=AR, x=NE_GSET, fill=Subtype), shape=21, size=3) + theme_classic() +
  scale_fill_manual(values = palette_subtypes)

dev.off()

### Heatmap of the top 2000 genes by variability

rv_all= rowVars(ncount)
names(rv_all) = rownames(ncount)
genes_all = names(rv_all)[order(rv_all, decreasing = TRUE)]
selected_all = genes_all[1:2000]
library(ComplexHeatmap)
mat = as.matrix(ncount[selected_all,])
mat = t(scale(t(mat)))
tmp = mdata
tmp$DatasetLabel = tmp$Dataset
tmp$DatasetLabel[tmp$Dataset %in% c("ANTE", "EIMB RAS", "Lim Y et al.")]= "Other"

summary(factor(tmp$DatasetLabel))
palette_databases = palette.colors(n=length(unique(tmp$DatasetLabel)), palette = palette.pals()[7] )
all_databases = unique(tmp$DatasetLabel)
test_pal = c()
for ( i in 1:length(all_databases)){
  test_pal[all_databases[i]]=palette_databases[i]
}
sub_pal = c()
for ( i in 1:length(levels(tmp$Subtype))){
  sub_pal[levels(tmp$Subtype)[i]] = palette_subtypes[i]
}
column_ha = HeatmapAnnotation(Subtype = tmp$Subtype, Dataset = tmp$DatasetLabel, LibraryType = tmp$LibraryType, libraryLayout = tmp$LibraryLayout,
                              col = list(
                                Subtype =sub_pal,
                                Dataset = test_pal
                              ))
pdf(paste0(outdir, "/Heatmap.pdf"))
Heatmap(mat, show_row_names = F, show_column_names = F, show_row_dend = F, show_column_dend = F, top_annotation = column_ha)

rv_all= rowVars(ncount)
names(rv_all) = rownames(ncount)
genes_all = names(rv_all)[order(rv_all, decreasing = TRUE)]
selected_all = genes_all[1:500]
mat = as.matrix(ncount[selected_all,])
mat = t(scale(t(mat)))

Heatmap(mat, show_row_names = F, show_column_names = F, show_row_dend = F, show_column_dend = F, top_annotation = column_ha)
dev.off()
### Estimate score
library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
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




pdf(paste0(outdir, "/Estimate.pdf"), width = 8, height = 6)
ggplot(mdata) + geom_point(aes(x=PC1, y=PC2, fill=ESTIMATE.Immune.Score), size=3, shape=21) +
  theme_classic() + scale_fill_gradient2(low="blue", high="red", midpoint = (max(mdata$ESTIMATE.Immune.Score) + min(mdata$ESTIMATE.Immune.Score) )/2 )
ggplot(mdata) + geom_point(aes(x=PC1, y=PC2, fill=ESTIMATE.Stromal.Score), size=3, shape=21) +
  theme_classic() + scale_fill_gradient2(low="blue", high="red", midpoint = (max(mdata$ESTIMATE.Stromal.Score) + min(mdata$ESTIMATE.Stromal.Score) )/2 )
ggplot(mdata) + geom_point(aes(x=PC1, y=PC2, fill=ESTIMATE.Score), size=3, shape=21) +
  theme_classic() + scale_fill_gradient2(low="blue", high="red", midpoint = (max(mdata$ESTIMATE.Score) + min(mdata$ESTIMATE.Score) )/2 )

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
GS = list()

for ( gsn in unique(m_df$gs_name)){
  GS[[gsn]] =rownames(txi.gene$abundance)[gsub("\\..*" , "", rownames(txi.gene$abundance)) %in% m_df$human_ensembl_gene[m_df$gs_name == gsn] ];
}

h.ssgsea <- gsva(txi.gene$abundance[keep,mdata$Name], GS, method = "ssgsea", parallel.sz=4,  mx.diff=TRUE, ssgsea.norm	=T)
h.ssgsea = t(h.ssgsea)
all(rownames(h.ssgsea) == mdata$Name)
mdata[,colnames(h.ssgsea)]=h.ssgsea


pdf(paste0(outdir,"/PCA_individual_Hallmarks.pdf"), width = 8, height = 6)
for ( gsn in names(GS)){
    tmp = data.frame(PC1=mdata$PC1, PC2=mdata$PC2, h=mdata[,gsn])
    tmp$h = (2*((tmp$h - min(tmp$h)) / (max(tmp$h)- min(tmp$h)))) -1
    p<-ggplot(tmp) + geom_point(aes(x=PC1, y=PC2, fill=h), shape=21, size=3)+ggtitle(gsn)+ theme_classic()+
      scale_fill_gradient2(low = "blue", high = "red", mid = "grey")
    print(p)  
}
dev.off()



pdf(paste0(outdir,"/signatures_PCA.pdf"), width = 8, height = 6)
for ( gsn in c("AR", "AR_GSET", "NE_GSET", "ESTIMATE.Immune.Score", "ESTIMATE.Stromal.Score", "ESTIMATE.Score")){
  if ( gsn %in% colnames(mdata)){
    tmp = data.frame(PC1=mdata$PC1, PC2=mdata$PC2, h=mdata[,gsn])
    tmp$h = (2*((tmp$h - min(tmp$h)) / (max(tmp$h)- min(tmp$h)))) -1
    p<-ggplot(tmp) + geom_point(aes(x=PC1, y=PC2, fill=h), shape=21, size=3)+ggtitle(gsn)+ theme_classic()+
      scale_fill_gradient2(low = "blue", high = "red", mid = "grey")
    print(p)  
  }
}
dev.off()
install.packages("factoextra")
library(factoextra)
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
pca = pcaPlot(ncount, mdata$Subtype, ngenes=2000)
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
base_dir=paste0(outdir,"/GO/")
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
base_dir=paste0(outdir,"/GSEA/")
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


### Databases division

library(treemapify)

mdata.datasets = mdata %>% group_by(Dataset, Subtype) %>% summarise(count = n())
mdata.datasets$Subtype = factor(mdata.datasets$Subtype, levels=c("NORMAL", "PRIMARY", "ARPC", "DNPC", "NEPC"))
mdata.datasets$label = paste0(mdata.datasets$Dataset, "\nN=",mdata.datasets$count)
mdata.datasets$size = mdata.datasets$count
mdata.datasets$size[mdata.datasets$size < 10]=10
pdf(paste0(outdir,"/datasets.pdf"))
ggplot(mdata.datasets, aes(area=size, group=Dataset, fill=Subtype, subgroup=Subtype)) + geom_treemap() +  geom_treemap_subgroup_border(colour="black", size=3)+
  geom_treemap_text(aes(label=label), place="centre", color="white", min.size = 0)+
  scale_fill_manual(values=palette_subtypes)
dev.off()

  



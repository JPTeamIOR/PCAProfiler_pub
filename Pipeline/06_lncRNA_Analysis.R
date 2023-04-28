#!/usr/bin/Rscript
library("ggplot2")
library(ggrepel)
library(tidyverse)
source("./00_functions.R")


mdata=readRDS("./objects/mdata_signatures.RDS")
txi.tx = readRDS("./objects/txi.tx.RDS")
txi.gene = readRDS("./objects/txi.gene.RDS")
ncount = readRDS("./objects/ncount.RDS")
ncount.txi = readRDS("./objects/ncount.txi.RDS")

tx2gene = read.table("./metadata/salmon_tx2gene.tsv")
colnames(tx2gene) = c("transcript_id","gene_id", "gene_name")
tx_info = read.csv("./metadata/transcripts_info.csv")
tx_names = read.table("./metadata/transcript_protein.tsv", col.names = c("transcript_id", "protein_id", "transcript_name"), sep="\t")
g_info= read.csv("./metadata/genes_info.csv")
rownames(mdata)=mdata$Name
mdata = mdata[colnames(txi.tx$counts), ]

### Long Non coding RNAs gene level analysis
onlyPE = mdata$LibraryLayout == "PAIRED"
ptime = mdata$pseudotime[onlyPE]

keep.gene = rowSums(txi.gene$counts > 3) >= (ncol(txi.gene$counts)*0.05) & 
  rownames(ncount) %in% g_info$gene_id[g_info$gene_type == "lncRNA"]


dat = t(sapply(rownames(ncount)[keep.gene], 
               function(x){ 
                 ct = cor.test(ptime, ncount[x, mdata$Name[onlyPE]]) ; 
                 c(id= x, corr=as.numeric(ct$estimate), pval = as.numeric(ct$p.value), stat = as.numeric(ct$statistic))
               }))
dat = as.data.frame(dat)
dat$corr = as.numeric(dat$corr)
dat$pval = as.numeric(dat$pval)
dat$stat = as.numeric(dat$stat)
if ( any(is.na(dat$corr))){
  dat$corr[is.na(dat$corr)]=0
  dat$pval[is.na(dat$pval)]=1
  dat$stat[is.na(dat$stat)]=0
}
dat$padj = p.adjust(dat$pval,method = "fdr")
dat = merge(dat, data.frame(id=g_info$gene_id, gene_type = g_info$gene_type, name = g_info$gene_name), by="id", all.x=T, all.y=F)

dir.create("./lncRNA/",showWarnings = F)
columns_to_write= c("id", "name", "gene_type", "corr", "pval", "padj", "stat")
write.csv(dat[order(dat$padj),columns_to_write], file = "./lncRNA/Gene_Level_Correlations.csv", row.names = F)

to_plot = dat

to_plot$color="grey"
to_plot$color[to_plot$corr < -0.5]="blue"
to_plot$color[to_plot$corr > 0.5]="red"
to_show = c(to_plot$name[order(to_plot$corr)][1:5], to_plot$name[order(to_plot$corr, decreasing = T)][1:5])
to_plot$to_show=""
to_plot$to_show[to_plot$name %in% to_show]=to_plot$name[to_plot$name %in% to_show]
pdf("./lncRNA/Genes_correlations_to_pseudotime.pdf", width = 10, height = 7)
ggplot(to_plot, aes(x=corr, y=stat, label=to_show))+
  geom_point(data=to_plot[to_plot$color == "grey", ], color="grey", size=3 ) + 
  geom_point(data=to_plot[to_plot$color == "red", ],  color="red", size=4 )+
  geom_point(data=to_plot[to_plot$color == "blue", ], color="blue", size=4 )+
  geom_label_repel(max.overlaps = Inf, seed=123, size=4,  nudge_y = 1, nudge_x = .5, force = 5) + theme_bw()+ 
  xlab("Pearson Correlation Coefficient") + ylab("Pearson's product moment statistic")
dev.off()

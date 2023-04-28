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
g_info= unique(tx2gene[,c("gene_id", "gene_name")])
rownames(mdata)=mdata$Name
mdata = mdata[colnames(txi.tx$counts), ]


### For protein coding transcripts, compute the correlation coefficients between transcripts normalized counts
### and pseudotime and the percentage of contribution of each transcript to the gene counts. 
### Consider only Paired End samples

onlyPE = mdata$LibraryLayout == "PAIRED"
ptime = mdata$pseudotime[onlyPE]
keep.tx = rowSums(txi.tx$counts > 3) >= (ncol(txi.tx$counts)*0.05) & 
  rownames(ncount.txi) %in% tx_info$transcript_id[tx_info$transcript_type == "protein_coding"]

perc_tx = t(sapply(rownames(txi.tx$counts)[keep.tx], function(x){txi.tx$counts[x,mdata$Name[onlyPE]] / txi.gene$counts[tx_info$gene_id[tx_info$transcript_id == x], mdata$Name[onlyPE]] } ))
perc_tx[is.na(perc_tx)] = 0


### consider only the Transcript that at least in one sample contribute to more than 10% to the gene count
keep.tx = keep.tx & rownames(ncount.txi) %in% rownames(perc_tx)[apply(perc_tx, 1, max) > 0.10]
dat = t(sapply(rownames(ncount.txi)[keep.tx], 
            function(x){ 
              ct = cor.test(ptime, ncount.txi[x, mdata$Name[onlyPE]]) ; 
              c(tx= x, corr=as.numeric(ct$estimate), pval = as.numeric(ct$p.value), stat = as.numeric(ct$statistic))
              }))
dat = as.data.frame(dat)
dat$corr = as.numeric(dat$corr)
dat$pval = as.numeric(dat$pval)
dat$stat = as.numeric(dat$stat)
if ( sum(is.na(dat$pval)) > 0 ){
  dat[is.na(dat$pval),]$pval = 1
  dat[is.na(dat$stat),]$stat = 0
  dat[is.na(dat$corr),]$corr = 0
}
dat$padj = p.adjust(dat$pval)
dat = merge(dat, data.frame(tx=tx_info$transcript_id, gid= tx_info$gene_id, tx_type = tx_info$transcript_type ), by="tx", all.x=T, all.y=F)
dat = dat %>% group_by(gid) %>% mutate(gCorr = cor(ptime, ncount[unique(gid), mdata$Name[onlyPE]]))

dat = dat[dat$tx_type == "protein_coding", ]
dat = merge(dat, data.frame(gid = g_info$gene_id, name = g_info$gene_name), by="gid", all.x=T, all.y=F)
dat = merge(dat, data.frame(tx = tx_names$transcript_id, protein = tx_names$protein_id, tx_name = tx_names$transcript_name), by="tx")

dat_bkup = dat
dat = dat[abs(dat$gCorr) < 0.4 ,]
dir.create("./isoforms_analysis/", showWarnings = F)
columns_to_write= c("tx", "tx_name", "protein", "gid", "name", "corr", "stat", "padj", "pval", "gCorr", "tx_type")
write.csv(dat_bkup[order(dat_bkup$corr),columns_to_write], "./isoforms_analysis/Correlations.csv", row.names = F)
write.csv(dat_bkup[order(dat_bkup$corr),columns_to_write][1:100,], "./isoforms_analysis/Correlations_bottom_100_all.csv", row.names = F)
write.csv(dat[order(dat$corr),columns_to_write][1:100,], "./isoforms_analysis/Correlations_bottom_100_LowGeneCorr.csv", row.names = F)
write.csv(dat_bkup[order(dat_bkup$corr, decreasing = T),columns_to_write][1:100,], "./isoforms_analysis/Correlations_top_100_all.csv", row.names = F)
write.csv(dat[order(dat$corr, decreasing = T),columns_to_write][1:100,], "./isoforms_analysis/Correlations_top_100_LowGeneCorr.csv", row.names = F)

to_plot = dat
to_plot$color = "grey"
to_plot$color[abs(to_plot$corr) > 0.50 & abs(to_plot$gCorr) < 0.4]="red"
to_plot$color[to_plot$name == "AR"]="blue"
to_plot$label = paste(to_plot$tx, to_plot$tx_name, to_plot$name,sep="\n")
to_show= c( to_plot$tx[to_plot$name == "AR"], to_plot$tx[to_plot$color == "red" & 
                                                           to_plot$name %in% c("ERF","MXI1", "CYHR1", "SPTAN1")])
to_plot$to_show=to_plot$label
to_plot$to_show[! to_plot$tx %in% to_show]=""
pdf("./isoforms_analysis/Isoforms_correlations_Filtered.pdf", width = 10, height = 7)
ggplot(to_plot, aes(x=corr, y=stat, label=to_show))+
  geom_point(data=to_plot[to_plot$color == "grey", ], color="grey", size=3 ) + 
  geom_point(data=to_plot[to_plot$color == "red", ], aes(x=corr, y=stat), color="red", size=4 )+
  geom_point(data=to_plot[to_plot$color == "blue", ], aes(x=corr, y=stat), color="blue", size=4 )+
  geom_label_repel(max.overlaps = Inf, seed=123, size=4,  nudge_y = 1, nudge_x = .5, force = 5) + theme_bw()+ 
  xlab("Pearson Correlation Coefficient") + ylab("Pearson's product moment statistic")
dev.off()


to_show = to_plot$tx[to_plot$name %in% c("ERF")]
to_plot$color = "grey"
to_plot$color[abs(to_plot$corr) > 0.50 & abs(to_plot$gCorr) < 0.4]="red"
to_plot$color[to_plot$tx %in% to_show] = "blue"
to_plot$to_show=to_plot$label
to_plot$to_show[! to_plot$tx %in% to_show]=""

pdf("./isoforms_analysis/Isoforms_correlations_ERF.pdf", width = 10, height = 7)
ggplot(to_plot, aes(x=corr, y=stat, label=to_show))+
  geom_point(data=to_plot[to_plot$color == "grey", ], color="grey", size=3 ) + 
  geom_point(data=to_plot[to_plot$color == "red", ], aes(x=corr, y=stat), color="red", size=4 )+
  geom_point(data=to_plot[to_plot$color == "blue", ], aes(x=corr, y=stat), color="blue", size=4 )+
  geom_label_repel(max.overlaps = Inf, seed=123, size=4,  nudge_y = 1, nudge_x = 0.5, force=5) + theme_bw()+ 
  xlab("Pearson Correlation Coefficient") + ylab("Pearson's product moment statistic")
dev.off()


to_show = dat$tx[dat$name %in% c("MXI1")]
to_show = rowSums(ncount.txi[to_show,])
to_show = to_show[order(to_show, decreasing = T)]
to_show = names(to_show)[1:4]
to_plot$color = "grey"
to_plot$color[abs(to_plot$corr) > 0.50 & abs(to_plot$gCorr) < 0.4]="red"
to_plot$color[to_plot$tx %in% to_show] = "blue"
to_plot$to_show=to_plot$label
to_plot$to_show[! to_plot$tx %in% to_show]=""

pdf("./isoforms_analysis/Isoforms_correlations_MXI1.pdf", width = 10, height = 7)
ggplot(to_plot, aes(x=corr, y=stat, label=to_show))+
  geom_point(data=to_plot[to_plot$color == "grey", ], color="grey", size=3 ) + 
  geom_point(data=to_plot[to_plot$color == "red", ], aes(x=corr, y=stat), color="red", size=4 )+
  geom_point(data=to_plot[to_plot$color == "blue", ], aes(x=corr, y=stat), color="blue", size=4 )+
  geom_label_repel(max.overlaps = Inf, seed=123, size=4,  nudge_y = 1, nudge_x = 0.5, force=5) + theme_bw()+ 
  xlab("Pearson Correlation Coefficient") + ylab("Pearson's product moment statistic")
dev.off()



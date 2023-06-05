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
tx_names = read.table("./metadata/transcript_protein.tsv", col.names = c("transcript_id", "protein_id", "transcript_name", "CDS"), sep="\t")
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
outdir="./results_05_Isoform_analysis/"

dir.create(outdir, showWarnings = F)
columns_to_write= c("tx", "tx_name", "protein", "gid", "name", "corr", "stat", "padj", "pval", "gCorr", "tx_type")

write.csv(dat_bkup[order(dat_bkup$corr),columns_to_write], paste0(outdir,"Correlations.csv"), row.names = F)
write.csv(dat_bkup[abs(dat_bkup$corr) > 0.5 & abs(dat_bkup$gCorr) < 0.4 ,columns_to_write], paste0(outdir,"Correlations_Filtered.csv"), row.names = F)
write.csv(dat_bkup[order(dat_bkup$corr),columns_to_write][1:100,], paste0(outdir,"Correlations_bottom_100_all.csv"), row.names = F)
write.csv(dat[order(dat$corr),columns_to_write][1:100,], paste0(outdir,"Correlations_bottom_100_LowGeneCorr.csv"), row.names = F)
write.csv(dat_bkup[order(dat_bkup$corr, decreasing = T),columns_to_write][1:100,], paste0(outdir,"Correlations_top_100_all.csv"), row.names = F)
write.csv(dat[order(dat$corr, decreasing = T),columns_to_write][1:100,], paste0(outdir,"Correlations_top_100_LowGeneCorr.csv"), row.names = F)
dat.filtered = dat_bkup[abs(dat_bkup$corr) > 0.5 & abs(dat_bkup$gCorr) < 0.4 ,columns_to_write]
nrow(dat.filtered)
sum(dat.filtered$corr < 0)
sum(dat.filtered$corr > 0)
to_plot = dat
to_plot$color = "grey"
to_plot$color[abs(to_plot$corr) > 0.50 & abs(to_plot$gCorr) < 0.4]="red"
to_plot$color[to_plot$name == "AR"]="blue"
to_plot$label = paste(to_plot$tx, to_plot$tx_name, to_plot$name,sep="\n")
to_show= c( to_plot$tx[to_plot$name == "AR"], to_plot$tx[to_plot$color == "red" & 
                                                           to_plot$name %in% c("ERF","MXI1", "CYHR1", "SPTAN1")])
to_plot$to_show=to_plot$label
to_plot$to_show[! to_plot$tx %in% to_show]=""
pdf(paste0(outdir,"Isoforms_correlations_Filtered.pdf"), width = 10, height = 7)
ggplot(to_plot, aes(x=corr, y=stat, label=to_show))+
  geom_point(data=to_plot[to_plot$color == "grey", ], color="grey", size=3 ) + 
  geom_point(data=to_plot[to_plot$color == "red", ], aes(x=corr, y=stat), color="red", size=4 )+
  geom_point(data=to_plot[to_plot$color == "blue", ], aes(x=corr, y=stat), color="blue", size=4 )+
  geom_label_repel(max.overlaps = Inf, seed=123, size=4,  nudge_y = 1, nudge_x = .5, force = 5) + theme_classic()+ 
  xlab("Pearson Correlation Coefficient") + ylab("Pearson's product moment statistic")
dev.off()


to_show = to_plot$tx[to_plot$name %in% c("ERF")]
to_plot$color = "grey"
to_plot$color[abs(to_plot$corr) > 0.50 & abs(to_plot$gCorr) < 0.4]="red"
to_plot$color[to_plot$tx %in% to_show] = "blue"
to_plot$to_show=to_plot$label
to_plot$to_show[! to_plot$tx %in% to_show]=""

pdf(paste0(outdir,"Isoforms_correlations_ERF.pdf"), width = 10, height = 7)
ggplot(to_plot, aes(x=corr, y=stat, label=to_show))+
  geom_point(data=to_plot[to_plot$color == "grey", ], color="grey", size=3 ) + 
  geom_point(data=to_plot[to_plot$color == "red", ], aes(x=corr, y=stat), color="red", size=4 )+
  geom_point(data=to_plot[to_plot$color == "blue", ], aes(x=corr, y=stat), color="blue", size=4 )+
  geom_label_repel(max.overlaps = Inf, seed=123, size=4,  nudge_y = 1, nudge_x = 0.5, force=5) + theme_classic()+ 
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

pdf(paste0(outdir,"Isoforms_correlations_MXI1.pdf"), width = 10, height = 7)
ggplot(to_plot, aes(x=corr, y=stat, label=to_show))+
  geom_point(data=to_plot[to_plot$color == "grey", ], color="grey", size=3 ) + 
  geom_point(data=to_plot[to_plot$color == "red", ], aes(x=corr, y=stat), color="red", size=4 )+
  geom_point(data=to_plot[to_plot$color == "blue", ], aes(x=corr, y=stat), color="blue", size=4 )+
  geom_label_repel(max.overlaps = Inf, seed=123, size=4,  nudge_y = 1, nudge_x = 0.5, force=5) + theme_bw()+ 
  xlab("Pearson Correlation Coefficient") + ylab("Pearson's product moment statistic")
dev.off()


### AR V-7 analysis
palette_subtypes = c("#25681E", "#F63A21","#4889C7", "#1010C0", "#93355A")

to_plot = NULL
variants = list("AR-V7"="ENST00000504326.5", "AR-FL"="ENST00000374690.9")
for ( variant in names(variants)){
  tmp = mdata[onlyPE,c("Name", "PC.Type", "Subtype", "pseudotime")]
  tmp$isoform = variant
  tmp$ncount = ncount.txi[variants[[variant]], tmp$Name]
  tmp$perc = perc_tx[variants[[variant]], tmp$Name]
  to_plot = rbind(to_plot, tmp)
}
library(ggsignif)
library(gridExtra)

pdf(paste0(outdir, "AR_v7.pdf"))
to_plot$Subtype = factor(to_plot$Subtype, levels= c("NORMAL", "PRIMARY", "ARPC", "DNPC", "NEPC"))
ggplot(to_plot) + geom_boxplot(aes(x=Subtype, fill=isoform, y=ncount)) + theme_classic() + scale_fill_manual(values = c("blue", "red"))+ylab("transcript CPM")+
  facet_grid(rows=vars(isoform))
p1 = ggplot(to_plot[to_plot$isoform == "AR-FL",], aes(x=Subtype, fill=isoform, y=ncount)) + geom_boxplot()  + geom_signif(comparisons = list(c("PRIMARY", "ARPC"), c("NORMAL", "PRIMARY"), c("DNPC", "NEPC")), test = "t.test")+
  theme_classic() + scale_fill_manual(values = c("red"))+ylab("transcript CPM") + 
  facet_grid(rows=vars(isoform))
p2 = ggplot(to_plot[to_plot$isoform == "AR-V7",], aes(x=Subtype, fill=isoform, y=ncount)) + geom_boxplot()  + geom_signif(comparisons = list(c("PRIMARY", "ARPC"), c("NORMAL", "PRIMARY"), c("DNPC", "NEPC")), test = "t.test")+
  theme_classic() + scale_fill_manual(values = c("blue"))+ylab("transcript CPM") + 
  facet_grid(rows=vars(isoform))
grid.arrange(p1, p2, nrow=2)

ggplot(to_plot) + geom_point(aes(x=pseudotime, y=ncount, fill=Subtype), shape=21, size=3) + theme_classic() + 
  scale_fill_manual(values = palette_subtypes) + facet_grid(rows = vars(isoform)) + xlab("Disease Progression") + ylab("transcript CPM")
ggplot(to_plot) + geom_point(aes(x=pseudotime, y=perc, fill=Subtype), shape=21, size=3) + theme_classic() + 
  scale_fill_manual(values = palette_subtypes) + facet_grid(rows = vars(isoform)) + xlab("Disease Progression") + ylab("transcript/gene count")
dev.off()


### PKM1 -> adult tissues, PKM2 -> embryonic/tumor
pkm = tx_info[tx_info$gene_id == g_info$gene_id[g_info$gene_name == "PKM"] & tx_info$transcript_type == "protein_coding",]
pkm$maxExp = rowMax(ncount.txi[pkm$transcript_id,mdata$Name[onlyPE]])
pkm = pkm[order(pkm$maxExp, decreasing = T),]
pkm = pkm[pkm$length > 2000,]

to_plot = NULL
variants = list("PKM1"="ENST00000319622.10", "PKM2"="ENST00000335181.10")
for ( variant in names(variants)){
  tmp = mdata[onlyPE,c("Name", "PC.Type", "Subtype", "pseudotime")]
  tmp$isoform = variant
  tmp$ncount = ncount.txi[variants[[variant]], tmp$Name]
  tmp$perc = perc_tx[variants[[variant]], tmp$Name]
  to_plot = rbind(to_plot, tmp)
}

pdf(paste0(outdir, "PKM.pdf"))
to_plot$Subtype = factor(to_plot$Subtype, levels= c("NORMAL", "PRIMARY", "ARPC", "DNPC", "NEPC"))
ggplot(to_plot) + geom_boxplot(aes(x=Subtype, fill=isoform, y=ncount)) + theme_classic() + scale_fill_manual(values = c("blue", "red"))+ylab("transcript CPM")
p1 = ggplot(to_plot[to_plot$isoform == "PKM1",], aes(x=Subtype, fill=isoform, y=ncount)) + geom_boxplot()  + geom_signif(comparisons = list(c("PRIMARY", "ARPC"), c("NORMAL", "PRIMARY"), c("DNPC", "NEPC")), test = "t.test")+
  theme_classic() + scale_fill_manual(values = c("red"))+ylab("transcript CPM") + 
  facet_grid(rows=vars(isoform))
p2 = ggplot(to_plot[to_plot$isoform == "PKM2",], aes(x=Subtype, fill=isoform, y=ncount)) + geom_boxplot()  + geom_signif(comparisons = list(c("PRIMARY", "ARPC"), c("NORMAL", "PRIMARY"), c("DNPC", "NEPC")), test = "t.test")+
  theme_classic() + scale_fill_manual(values = c("blue"))+ylab("transcript CPM") + 
  facet_grid(rows=vars(isoform))
grid.arrange(p1, p2, nrow=2)

ggplot(to_plot) + geom_point(aes(x=pseudotime, y=ncount, fill=Subtype), shape=21, size=3) + theme_classic() + 
  scale_fill_manual(values = palette_subtypes) + facet_grid(rows = vars(isoform)) + xlab("Disease Progression") + ylab("transcript CPM")
ggplot(to_plot) + geom_point(aes(x=pseudotime, y=perc, fill=Subtype), shape=21, size=3) + theme_classic() + 
  scale_fill_manual(values = palette_subtypes) + facet_grid(rows = vars(isoform)) + xlab("Disease Progression") + ylab("transcript/gene count")
dev.off()



### AR V-7 analysis
palette_subtypes = c("#25681E", "#F63A21","#4889C7", "#1010C0", "#93355A")
gene_id = g_info$gene_id[g_info$gene_name == "AR"]
pkm = tx_info[tx_info$gene_id == gene_id ,]
pkm$maxExp = rowMax(ncount.txi[pkm$transcript_id,mdata$Name[onlyPE]])
pkm$mean = rowMeans(ncount.txi[pkm$transcript_id,mdata$Name[onlyPE]])
pkm = pkm[order(pkm$maxExp, decreasing = T),]
pkm

to_plot = NULL
variants = list("AR-dLBD"=c("ENST00000504326.5", "ENST00000613054.2","ENST00000514029.5"), "AR-LBD"=c("ENST00000374690.9", "ENST00000396043.3", "ENST00000612452.5"))
for ( variant in names(variants)){
  tmp = mdata[onlyPE,c("Name", "PC.Type", "Subtype", "pseudotime")]
  tmp$isoform = variant
  if ( length(variants[[variant]]) == 1 ){
    tmp$ncount = ncount.txi[variants[[variant]], tmp$Name]
    tmp$perc = perc_tx[variants[[variant]], tmp$Name]
  } else {
    tmp$ncount = colSums(txi.tx$counts[variants[[variant]], tmp$Name])
    tmp$perc = as.numeric(tmp$ncount/txi.gene$counts[gene_id, tmp$Name])
    tmp$ncount = log2( (tmp$ncount / (colSums(txi.tx$counts[,tmp$Name]) / 1000000))+1 )
  }
  to_plot = rbind(to_plot, tmp)
}

pdf(paste0(outdir, "AR_v7_grouped.pdf"))
to_plot$Subtype = factor(to_plot$Subtype, levels= c("NORMAL", "PRIMARY", "ARPC", "DNPC", "NEPC"))
ggplot(to_plot) + geom_boxplot(aes(x=Subtype, fill=isoform, y=ncount)) + theme_classic() + scale_fill_manual(values = c("blue", "red"))+ylab("transcript CPM")
ggplot(to_plot) + geom_point(aes(x=pseudotime, y=ncount, fill=Subtype), shape=21, size=3) + theme_classic() + 
  scale_fill_manual(values = palette_subtypes) + facet_grid(rows = vars(isoform)) + xlab("Disease Progression") + ylab("transcript CPM")
ggplot(to_plot) + geom_point(aes(x=pseudotime, y=perc, fill=Subtype), shape=21, size=3) + theme_classic() + 
  scale_fill_manual(values = palette_subtypes) + facet_grid(rows = vars(isoform)) + xlab("Disease Progression") + ylab("transcript/gene count")
dev.off()


### group similar isoforms PKM1 -> adult tissues, PKM2 -> embryonic/tumor
gene_id = g_info$gene_id[g_info$gene_name == "PKM"]
pkm = tx_info[tx_info$gene_id == gene_id & tx_info$transcript_type == "protein_coding",]
pkm$maxExp = rowMax(ncount.txi[pkm$transcript_id,mdata$Name[onlyPE]])
pkm = pkm[order(pkm$maxExp, decreasing = T),]
pkm = pkm[pkm$length > 2000,]

to_plot = NULL
variants = list("PKM-EX9"=c("ENST00000319622.10", "ENST00000565184.6", "ENST00000389093.7"), "PKM-EX10"=c("ENST00000335181.10"))
variant="PKM1"
for ( variant in names(variants)){
  tmp = mdata[onlyPE,c("Name", "PC.Type", "Subtype", "pseudotime")]
  tmp$isoform = variant
  if ( length(variants[[variant]]) == 1 ){
    tmp$ncount = ncount.txi[variants[[variant]], tmp$Name]
    tmp$perc = perc_tx[variants[[variant]], tmp$Name]
  } else {
    tmp$ncount = colSums(txi.tx$counts[variants[[variant]], tmp$Name])
    tmp$perc = as.numeric(tmp$ncount/txi.gene$counts[gene_id, tmp$Name])
    tmp$ncount = log2( (tmp$ncount / (colSums(txi.tx$counts[,tmp$Name]) / 1000000))+1 )
  }
  to_plot = rbind(to_plot, tmp)
}

pdf(paste0(outdir, "PKM_grouped.pdf"))
to_plot$Subtype = factor(to_plot$Subtype, levels= c("NORMAL", "PRIMARY", "ARPC", "DNPC", "NEPC"))
ggplot(to_plot) + geom_boxplot(aes(x=Subtype, fill=isoform, y=ncount)) + theme_classic() + scale_fill_manual(values = c("blue", "red"))+ylab("transcript CPM")+ 
  facet_grid(rows=vars(isoform))
ggplot(to_plot) + geom_boxplot(aes(x=Subtype, fill=isoform, y=perc)) + theme_classic() + scale_fill_manual(values = c("blue", "red"))+ ylab("transcript/gene count")+ 
  facet_grid(rows=vars(isoform))
ggplot(to_plot) + geom_point(aes(x=pseudotime, y=ncount, fill=Subtype), shape=21, size=3) + theme_classic() + 
  scale_fill_manual(values = palette_subtypes) + facet_grid(rows = vars(isoform)) + xlab("Disease Progression") + ylab("transcript CPM")
ggplot(to_plot) + geom_point(aes(x=pseudotime, y=perc, fill=Subtype), shape=21, size=3) + theme_classic() + 
  scale_fill_manual(values = palette_subtypes) + facet_grid(rows = vars(isoform)) + xlab("Disease Progression") + ylab("transcript/gene count")
dev.off()


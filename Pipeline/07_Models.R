#!/usr/bin/Rscript
library("ggplot2")
library(msigdbr, quietly = T)
library(tidyverse)
library("clusterProfiler")
library(ggrepel)
library(DESeq2)
library("GSVA", quietly = T)
library(sva, quietly = T)
library(tximport)
library(ggrepel)
library(slingshot)
source("./00_functions.R")

normalize.logCPM<- function(dat, lib.size){
  return(log2(sweep(x=dat, MARGIN= 2, STATS = lib.size/1000000 , FUN="/" ) +1))
}
outdir = "./results_07_Models/"
dir.create(outdir, showWarnings = F)
tx2gene = read.table("./metadata/salmon_tx2gene.tsv")
colnames(tx2gene) = c("transcript_id","gene_id", "gene_name")
tx_info = read.csv("./metadata/transcripts_info.csv")
tx_to_protein = read.table("./metadata/transcript_protein.tsv",  col.names = c("transcript_id", "protein_id", "transcript_name", "CDS"))
tx_info = merge(tx_info, tx_to_protein, by="transcript_id", all.x=T)
g_info= read.csv("./metadata/genes_info.csv")
atlas.txi.counts = readRDS("./objects/txi.tx.counts.RDS")
atlas.mdata=readRDS("./objects/mdata_signatures.RDS")

mdata = read.csv("./metadata/models_metadata.csv")
mdata$file = paste0("./models/", mdata$name, "/quant.sf")

files = mdata$file
names(files) = mdata$name
all(file.exists(files))
txi.tx <- tximport(files, type = "salmon", txOut = TRUE)
for ( nam in c("abundance", "length", "counts")){
  txi.tx[[nam]] = txi.tx[[nam]][rownames(atlas.txi.counts), mdata$name]
}

txi.gene <- summarizeToGene(txi.tx, tx2gene)

pca.model = readRDS("./objects/pca.RDS")
pseudotime.model = readRDS("./objects/slingshot.RDS")

ncount = normalize.logCPM(txi.gene$counts, colSums(txi.gene$counts))
colnames(ncount) = mdata$name
ncount.txi = normalize.logCPM(txi.tx$counts, colSums(txi.tx$counts))
colnames(ncount.txi) = mdata$name


pca = predict(pca.model, t(ncount))
pca[,1]= -pca[,1]
mdata$PC1 = pca[,1]
mdata$PC2 = pca[,2]

pseudotimes = predict(pseudotime.model,rbind(atlas.mdata[,c("PC1","PC2")] , pca[,1:2, drop=F] ))
mdata$pseudotime = slingPseudotime(pseudotimes)[mdata$name,1]


atlas.mdata$Subtype = factor(atlas.mdata$Subtype, levels= c("NORMAL", "PRIMARY", "ARPC", "DNPC", "NEPC"))
mdata$Subtype = factor(mdata$Subtype, levels= c("NORMAL", "PRIMARY", "ARPC", "DNPC", "NEPC"))

pdf(paste0(outdir, "./models_PCA.pdf"),width = 10, height = 7)
p<- ggplot()+geom_point(data=atlas.mdata, aes(x=PC1, y=PC2, fill=Subtype), alpha=0.5, size=2, shape=21) + 
  geom_point(data=mdata, aes(x=PC1, y=PC2, fill=Subtype), size=4, shape=21) + scale_fill_manual(values=palette_subtypes) + theme_classic()
print(p)
tmp =mdata
tmp$label = ""
tmp$label[!duplicated(tmp$type)]=tmp$type[!duplicated(tmp$type)]
p<-ggplot()+geom_hline(yintercept = 0)+geom_hline(yintercept = 1) + geom_point(data=atlas.mdata, aes(x=pseudotime, fill=Subtype, y=1), size=3, position = position_jitter(height = 0.05), pch=21)+ 
  geom_point(data=tmp, aes(x=pseudotime, y=0, fill=Subtype), size=3, pch=21, position=position_jitter(height = 0.05))+
  ylim(-1,2) + xlim(150, 250)+theme_classic() + scale_fill_manual(values=palette_subtypes) + geom_text_repel(data=tmp, aes(x=pseudotime, y=0, label=label), force = 20, force_pull = 1, max.overlaps = Inf)
print(p)
dev.off()




### Save them
saveRDS(list(
  ncount = ncount,
  ncount.txi = ncount.txi,
  txi.gene = txi.gene,
  txi.tx = txi.tx,
  mdata = mdata
), "./objects/models.RDS")
###

tx_info = tx_info[tx_info$transcript_type == "protein_coding",]
if ( file.exists("./objects/tx_contributions_models.RDS")){
  tx.contributions = readRDS("./objects/tx_contributions_models.RDS")  
} else {
  tx.contributions = t(sapply(tx_info$transcript_id, function(x){ txi.tx$counts[x,mdata$name ] / txi.gene$counts[tx_info$gene_id[tx_info$transcript_id == x], mdata$name] } ))
  tx.contributions[is.na(tx.contributions)] = 0
  saveRDS(tx.contributions, "./objects/tx_contributions_models.RDS")
}

### Include the ssGSEA
if ( file.exists("./objects/models.mdata.ssgsea.RDS")){
  mdata.ssgsea = readRDS("./objects/models.mdata.ssgsea.RDS")
} else {
  hallmarks = as.data.frame(msigdbr(species = "Homo sapiens", category = "H")) 
  hallmarks = rbind(hallmarks, as.data.frame(msigdbr(species = "Homo sapiens", category = "C2", subcategory = "PID")))
  hallmarks = rbind(hallmarks, as.data.frame(msigdbr(species = "Homo sapiens", category = "C6")))
  CS = list()
  # AR-SCORE NELSON (Bluemn et al., 2017) ---->
  CS$AR_GSET <- g_info$gene_id[g_info$gene_name %in% c("KLK3","KLK2","TMPRSS2","FKBP5","NKX3-1","PLPP1","PMEPA1","PART1","ALDH1A3","STEAP4")]
  # NE-SCORE NELSON (Bluemn et al., 2017) ---->
  CS$NE_GSET <- g_info$gene_id[g_info$gene_name %in% c("SYP","CHGA","CHGB","ENO2","CHRNB2","SCG3","SCN3A","PCSK1","ELAVL4","NKX2-1")]
  for ( h in unique(hallmarks$gs_name)){
    CS[[h]]=g_info$gene_id[g_info$gene_name %in% hallmarks$gene_symbol[hallmarks$gs_name == h]]
  }
  length(CS)
  keep = rowSums(ncount[,mdata$name] > 0) > nrow(mdata)*0.05
  sum(keep) / length(keep)
  THREADS=8
  res = gsva(as.matrix(ncount[keep,mdata$name,drop=F]), CS, mx.diff=TRUE, method = "ssgsea", parallel.sz=THREADS, ssgsea.norm=F)
  res=t(res)
  res = as.data.frame(res)
  res$name = rownames(res)
  nrow(mdata)
  nrow(res)
  mdata.ssgsea = merge(mdata, res, by="name")
  saveRDS(mdata.ssgsea, "./objects/models.mdata.ssgsea.RDS")
}

summary(factor(mdata$Subtype))

all_res = read.csv("./results_06_IsoformSwitch/IsoformSwitchAll/all_results.csv")  
all_res.bygene = read.csv("./results_06_IsoformSwitch/IsoformSwitchAll/all_results_by_gene.csv")
gene_list = unique(all_res.bygene$gene_name)

nrow(all_res)


odir = paste0(outdir, "/validation/")
cmdata = mdata.ssgsea
dir.create(odir, showWarnings = F)
all_res = all_res[all_res$target %in% colnames(cmdata),]
nrow(all_res)

candidates = read.table("./metadata/candidates.txt")$V1
all_res.bygene = all_res.bygene[order(all_res.bygene$delta_corr, decreasing = T),]
genes_to_print = unique(c(candidates,  all_res.bygene$gene_name[1:10], all_res$gene_name[all_res$target == "pseudotime"], "ERF"))
length(genes_to_print)

thr = 0.5

models.validation = NULL
i=1
for ( i in 1:nrow(all_res.bygene) ){
  plots = list()
  curr_gene_id = all_res.bygene$gene_id[i]
  curr_gene_name = all_res.bygene$gene_name[i]
  targets = unique(all_res$target[all_res$gene_id == curr_gene_id])
  if ( !is.null(models.validation)){
    print(paste(i, length(unique(models.validation$gene_name))))  
  }
  
  for ( target in targets ){
    tmp = all_res[all_res$gene_id ==curr_gene_id & all_res$target == target ,]
    iso_up = tmp$transcript_id[which.max(tmp$corr)]
    iso_down =  tmp$transcript_id[which.min(tmp$corr)]
    if ( sd(tx.contributions[iso_up, cmdata$name]) == 0  ){
      corr_a = 0
    } else {
      corr_a = cor(tx.contributions[iso_up, cmdata$name], cmdata[,target])  
    }
    if ( sd(tx.contributions[iso_down, cmdata$name]) == 0  ){
      corr_b = 0 
    } else {
      corr_b = cor(tx.contributions[iso_down, cmdata$name], cmdata[,target])  
    }
    models.validation = rbind(models.validation,
                              data.frame(
                                gene_id = curr_gene_id,
                                gene_name= curr_gene_name,
                                target = target,
                                iso_up = iso_up,
                                iso_down = iso_down,
                                corr_iso_up = corr_a,
                                corr_iso_down = corr_b,
                                valid = corr_a > thr & corr_b < -thr
                              ))
    if ( corr_a > thr & corr_b < -thr ){
      tmp = rbind(data.frame(name = iso_up,
                             counts = ncount.txi[iso_up, cmdata$name],
                             perc = tx.contributions[iso_up, cmdata$name],
                             sample = cmdata$name,
                             sampleType = cmdata$Subtype,
                             target = cmdata[,target]
                             ), 
                  data.frame(name = iso_down,
                             counts = ncount.txi[iso_down, cmdata$name], 
                             perc = tx.contributions[iso_down, cmdata$name],
                             sample = cmdata$name,
                             sampleType = cmdata$Subtype,
                             target = cmdata[,target]
                  ))
      if ( curr_gene_name %in% genes_to_print){
        tmp$sampleType = factor(tmp$sampleType, levels=c("NORMAL","PRIMARY", "ARPC", "DNPC", "NEPC"))
        plots[[target]] =
          ggplot(tmp, aes(x=target, y=perc)) + geom_point(aes(fill=sampleType), size=3, shape=21)+geom_smooth()+theme_bw()+xlab(target) +
          facet_grid(rows=vars(name)) + scale_fill_manual(values = palette_subtypes[2:length(palette_subtypes)])+
          ggtitle(paste0("Gene = " ,curr_gene_name, "\nTarget = " , target, "\n", iso_up, " (up) = ", corr_a, "\n", iso_down, " (down) = ", corr_b))
      }
    }
  }
  tmp = models.validation %>% filter(valid & gene_id == curr_gene_id ) %>% arrange(desc(corr_iso_up - corr_iso_down))
  if ( nrow(tmp) > 0 && curr_gene_name %in% genes_to_print ){
    pdf(paste0(odir, "/", curr_gene_name, ".pdf"))
    for ( tname in tmp$target[1:min(10,nrow(tmp))]){
      print(plots[[tname]])
    }
    dev.off()
  }
}

view(models.validation.by_gene[models.validation.by_gene$gene_name %in% candidates,])

nrow(all_res.bygene)
length(unique(models.validation$gene_name))

write.csv(models.validation, paste0(outdir, "/all_results.csv"), row.names = F)
models.validation.by_gene = models.validation %>% group_by(gene_id, gene_name) %>% 
  summarise(perc_valid = sum(valid) / n(), valid_targets = sum(valid), tot_target = n() , valid_targets_names = paste(target[valid], sep=",", collapse = ",") ) 
nrow(models.validation.by_gene)
models.validation.by_gene.small = models.validation.by_gene %>% select(gene_id, gene_name, valid_targets ,tot_target) %>% group_by(gene_id, gene_name,tot_target) %>%
  mutate(avg_targets = mean(valid_targets / tot_target) ) %>% arrange(desc(avg_targets))
write.csv(models.validation.by_gene, paste0(outdir, "/all_results_by_gene.csv"), row.names = F)
write.csv(models.validation.by_gene.small, paste0(outdir, "/all_results_by_gene_small.csv"), row.names = F)
write.csv(models.validation.by_gene.small[models.validation.by_gene.small$gene_name %in% genes_to_print,], paste0(outdir, "/all_results_by_gene_small_candidates.csv"), row.names = F)

all_res.bygene$valid = all_res.bygene$gene_name %in% models.validation$gene_name[models.validation$valid]
colnames(all_res.bygene)

pdf(paste0(outdir, "/validated_in_PDX.pdf"), width = 8, height = 6)
tmp = all_res.bygene
tmp$label = ""
mask = tmp$valid & ( tmp$gene_name %in% candidates | tmp$delta_corr > 1.5 | tmp$n_targets > 20) 
tmp$label[mask ] = tmp$gene_name[mask] 
ggplot(tmp, aes(x=n_targets, y=delta_corr)) + geom_point(data = tmp[tmp$valid,], size=3, shape=21, fill="red")+
  geom_point(data = tmp[!tmp$valid,], size=3, shape=21, fill="grey", alpha=0.2) + 
  geom_text_repel(aes(label=label), max.overlaps = 10000, size=2)+
  theme_classic() 
tmp = all_res.bygene
tmp$label = ""
mask = tmp$valid & ( tmp$gene_name %in% candidates | tmp$delta_corr > 1.5 | tmp$n_targets > 50) 
tmp$label[mask ] = tmp$gene_name[mask] 
ggplot(tmp, aes(x=n_targets, y=delta_corr)) + geom_point(data = tmp[tmp$valid,], size=3, shape=21, fill="red")+
  geom_point(data = tmp[!tmp$valid,], size=3, shape=21, fill="grey", alpha=0.2) + 
  geom_text_repel(aes(label=label), max.overlaps = 10000, size=3, min.segment.length = 0.001, force_pull = 0.3)+
  theme_classic() 

dev.off()

sum(all_res.bygene$valid)
sum(all_res.bygene$valid) / nrow(all_res.bygene)

###
pdf(paste0(outdir, "/pseudotime_isoswitch.pdf"))
tmp = unique(all_res[all_res$target == "pseudotime", c("gene_name", "delta_corr", "gcorr")])
tmp = tmp[order(tmp$delta_corr),]
tmp$gene_name = factor(tmp$gene_name, levels=tmp$gene_name)
ggplot(tmp) + geom_bar(aes(x=delta_corr, y=gene_name, fill=gcorr), stat="identity") + scale_fill_gradient2(low="blue", mid="grey", high = "red") + theme_classic()
dev.off()
###

curr_gene_name="ERF"
gene_id=g_info$gene_id[g_info$gene_name=="ERF"]
tx_info[tx_info$gene_id == gene_id,]
iso_up = "ENST00000222329.9"
iso_down = "ENST00000440177.6"
target="pseudotime"
tmp = rbind(data.frame(name = iso_up,
                       counts = ncount.txi[iso_up, cmdata$name],
                       perc = tx.contributions[iso_up, cmdata$name],
                       sample = cmdata$name,
                       sampleType = cmdata$Subtype,
                       target = cmdata[,target]
), 
data.frame(name = iso_down,
           counts = ncount.txi[iso_down, cmdata$name], 
           perc = tx.contributions[iso_down, cmdata$name],
           sample = cmdata$name,
           sampleType = cmdata$Subtype,
           target = cmdata[,target]
))

tmp$sampleType = factor(tmp$sampleType, levels=c("NORMAL","PRIMARY", "ARPC", "DNPC", "NEPC"))
ggplot(tmp, aes(x=target, y=perc)) + geom_point(aes(fill=sampleType), size=3, shape=21)+geom_smooth()+theme_bw()+xlab(target) +
  facet_grid(rows=vars(name)) + scale_fill_manual(values = palette_subtypes)
###



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
mask = tmp$valid & ( tmp$gene_name %in% candidates | tmp$delta_corr > 1.5 | tmp$n_targets > 50) 
tmp$label[mask ] = tmp$gene_name[mask] 
ggplot(tmp, aes(x=n_targets, y=delta_corr)) + geom_point(data = tmp[tmp$valid,], size=3, shape=21, fill="red")+
  geom_point(data = tmp[!tmp$valid,], size=3, shape=21, fill="grey", alpha=0.2) + 
  geom_text_repel(aes(label=label), max.overlaps = 10000, size=3, min.segment.length = 0, force_pull = 0.03, force = 10)+
  theme_classic() 
dev.off()

pdf(paste0(outdir, "/validated_in_PDX_bar.pdf"), width = 8, height = 6)
tmp = models.validation
all_res.bygene = all_res.bygene[order(all_res.bygene$n_targets, decreasing = T),]
tmp$gene_name = factor(tmp$gene_name, levels = all_res.bygene$gene_name)
ggplot(tmp[tmp$gene_name %in% all_res.bygene$gene_name[1:50], ]) + geom_bar(aes(x=gene_name, fill=valid), position = "stack") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ scale_fill_manual(values = c("grey", "red"))
ggplot(tmp[tmp$gene_name %in% models.validation.by_gene$gene_name[models.validation.by_gene$valid_targets > 0]  & tmp$gene_name %in% all_res.bygene$gene_name[1:100], ]) + 
  geom_bar(aes(x=gene_name, fill=valid), position = "stack") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values = c("grey", "red"))
models.validation.by_gene=models.validation.by_gene[order(models.validation.by_gene$valid_targets, decreasing = T),]
tmp$gene_name = factor(tmp$gene_name, levels=models.validation.by_gene$gene_name)
ggplot(tmp[ tmp$gene_name %in% models.validation.by_gene$gene_name[1:50], ]) + 
  geom_bar(aes(x=gene_name, fill=valid), position = "stack") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values = c("grey", "red"))
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

### Check ERF and MXI1
tx_names = read.table("./metadata/transcript_protein.tsv", col.names = c("transcript_id", "protein_id", "transcript_name", "CDS"), sep="\t")

gene_name = "ERF"
gene_id = g_info$gene_id[g_info$gene_name == gene_name]
isoforms = tx_names[tx_names$transcript_name %in% c("ERF-201", "ERF-202", "ERF-203"),]
target_tx = isoforms$transcript_id[isoforms$transcript_name == "ERF-202"]


to_plot= data.frame(
  id = gene_id,
  name = gene_name,
  ncount = ncount[gene_id, mdata$name],
  sample = mdata$name,
  subtype = mdata$Subtype,
  pseudotime = mdata$pseudotime
)

for ( i in 1:nrow(isoforms)){
  to_plot = rbind(to_plot, data.frame(
    id = isoforms$transcript_id[i],
    name = isoforms$transcript_name[i],
    ncount = ncount.txi[isoforms$transcript_id[i], mdata$name],
    sample = mdata$name,
    subtype = mdata$Subtype,
    pseudotime = mdata$pseudotime
  ))
}

to_plot = merge(to_plot, mdata[,c("name", "type")], by.x="sample" , by.y="name")

pdf(paste0(outdir, "ERF.pdf"))
ggplot(to_plot, aes(x=pseudotime, y=ncount)) + geom_point(aes(fill=subtype), size=3, shape=21) + 
  facet_grid(rows=vars(name))+ theme_bw() +scale_fill_manual(values = palette_subtypes)
ggplot(to_plot, aes(x=pseudotime, y=ncount)) + geom_point(aes(fill=subtype), size=3, shape=21) + 
  facet_grid(rows=vars(name))+ theme_bw() +scale_fill_manual(values = palette_subtypes) + geom_text_repel(aes(label=type), size=3)
dev.off()


gene_name = "MXI1"
gene_id = g_info$gene_id[g_info$gene_name == gene_name]
isoforms = tx_names[tx_names$transcript_name %in% c("MXI1-201", "MXI1-202", "MXI1-214"),]
target_tx = isoforms$transcript_id[isoforms$transcript_name == "ERF-214"]


to_plot= data.frame(
  id = gene_id,
  name = gene_name,
  ncount = ncount[gene_id, mdata$name],
  sample = mdata$name,
  subtype = mdata$Subtype,
  pseudotime = mdata$pseudotime
)

for ( i in 1:nrow(isoforms)){
  to_plot = rbind(to_plot, data.frame(
    id = isoforms$transcript_id[i],
    name = isoforms$transcript_name[i],
    ncount = ncount.txi[isoforms$transcript_id[i], mdata$name],
    sample = mdata$name,
    subtype = mdata$Subtype,
    pseudotime = mdata$pseudotime
  ))
}

to_plot = merge(to_plot, mdata[,c("name", "type")], by.x="sample" , by.y="name")

pdf(paste0(outdir, "MXI1.pdf"))
ggplot(to_plot, aes(x=pseudotime, y=ncount)) + geom_point(aes(fill=subtype), size=3, shape=21) + 
  facet_grid(rows=vars(name))+ theme_bw() +scale_fill_manual(values = palette_subtypes)
ggplot(to_plot, aes(x=pseudotime, y=ncount)) + geom_point(aes(fill=subtype), size=3, shape=21) + 
  facet_grid(rows=vars(name))+ theme_bw() +scale_fill_manual(values = palette_subtypes) + geom_text_repel(aes(label=type), size=3)
dev.off()




### Other models

omdata = read.csv("./metadata/models_metadata_others.csv")
omdata$file = paste0("./models/", omdata$name, "/quant.sf")

### remove PRECLH since it should be normal but the pseudotime is really high
omdata = omdata[omdata$name != "PRECLH_PROSTATE",]


files = omdata$file[omdata$AlignmentType == "Human"]
names(files) = omdata$name[omdata$AlignmentType == "Human"]
all(file.exists(files))
txi.tx <- tximport(files, type = "salmon", txOut = TRUE)
files = omdata$file[omdata$AlignmentType == "Human-Murine"]
names(files) = omdata$name[omdata$AlignmentType == "Human-Murine"]
all(file.exists(files))
txi.tx_hm <- tximport(files, type = "salmon", txOut = TRUE)

for ( nam in c("abundance", "length", "counts")){
  txi.tx_hm[[nam]] = txi.tx_hm[[nam]][rownames(txi.tx[[nam]]),]
  txi.tx[[nam]] = cbind(txi.tx[[nam]], txi.tx_hm[[nam]])
  txi.tx[[nam]]= txi.tx[[nam]][,omdata$name]
}

txi.gene <- summarizeToGene(txi.tx, tx2gene)

oncount = normalize.logCPM(txi.gene$counts, colSums(txi.gene$counts))
colnames(oncount) = omdata$name
oncount.txi = normalize.logCPM(txi.tx$counts, colSums(txi.tx$counts))
colnames(oncount.txi) = omdata$name


pca = predict(pca.model, t(oncount))
pca[,1]= -pca[,1]
omdata$PC1 = pca[,1]
omdata$PC2 = pca[,2]

pseudotimes = predict(pseudotime.model,rbind(atlas.mdata[,c("PC1","PC2")] , pca[,1:2, drop=F] ))
omdata$pseudotime = slingPseudotime(pseudotimes)[omdata$name,1]

omdata$type
pdf(paste0(outdir, "./other_models_PCA.pdf"),width = 10, height = 7)
for ( cohort in unique(omdata$cohort )){
  p<- ggplot()+geom_point(data=atlas.mdata, aes(x=PC1, y=PC2, fill=Subtype), alpha=0.2, size=2, shape=21) + 
    geom_point(data=omdata[omdata$cohort == cohort,], aes(x=PC1, y=PC2, color=type), size=4) + scale_fill_manual(values=palette_subtypes) + theme_classic() + ggtitle(cohort)
  print(p)
}
tmp =omdata
tmp$label = ""
tmp$label[!duplicated(tmp$type)]=tmp$type[!duplicated(tmp$type)]
tmp$y = as.numeric(factor(tmp$cohort))
p<-ggplot()+ geom_point(data=atlas.mdata, aes(x=pseudotime, fill=Subtype, y=0), size=3, position = position_jitter(height = 0.05), pch=21)+ 
  geom_point(data=tmp, aes(x=pseudotime, y=y, color=type), size=3,  position=position_jitter(height = 0.05))+
   xlim(min(tmp$pseudotime)-10, 260)+theme_classic() + scale_fill_manual(values=palette_subtypes) + geom_text_repel(data=tmp, aes(x=pseudotime, y=y, label=label), force = 20, force_pull = 1, max.overlaps = Inf)
print(p)
dev.off()

if ( file.exists("./objects/tx_contributions_other_models.RDS")){
  tx.contributions = readRDS("./objects/tx_contributions_other_models.RDS")
} else {
  
  tx.contributions = t(sapply(tx_info$transcript_id, function(x){ txi.tx$counts[x,omdata$name ] / txi.gene$counts[tx_info$gene_id[tx_info$transcript_id == x], omdata$name] } ))
  tx.contributions[is.na(tx.contributions)] = 0
  saveRDS(tx.contributions, "./objects/tx_contributions_other_models.RDS")
}




### Include the ssGSEA
if ( file.exists("./objects/other.models.mdata.ssgsea.RDS")){
  omdata.ssgsea = readRDS("./objects/other.models.mdata.ssgsea.RDS")
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
  keep_gene = rowSums(oncount[,omdata$name] > 0) > nrow(omdata)*0.05
  sum(keep_gene) / length(keep_gene)
  THREADS=8
  res = gsva(as.matrix(oncount[keep_gene,omdata$name,drop=F]), CS, mx.diff=TRUE, method = "ssgsea", parallel.sz=THREADS, ssgsea.norm=F)
  res=t(res)
  res = as.data.frame(res)
  res$name = rownames(res)
  nrow(omdata)
  nrow(res)
  omdata.ssgsea = merge(omdata, res, by="name")
  saveRDS(omdata.ssgsea, "./objects/other.models.mdata.ssgsea.RDS")
}

other.models.validation = NULL
odir = paste0(outdir, "/other_validations/")
dir.create(odir, showWarnings = F)
for ( curr_gene_name in  unique(models.validation$gene_name[models.validation$valid])){
  plots = list()
  curr_gene_id = unique(models.validation$gene_id[models.validation$gene_name == curr_gene_name])
  tmp = models.validation[models.validation$gene_id == curr_gene_id & models.validation$valid,]
  for ( i in 1:nrow(tmp) ){
    iso_up = tmp$iso_up[i]
    iso_down =  tmp$iso_down[i]
    target = tmp$target[i]
    other.models.validation.partial =data.frame(
      gene_id = curr_gene_id,
      gene_name= curr_gene_name,
      target = target,
      iso_up = iso_up,
      iso_down = iso_down
    )
    for ( cohort in unique(omdata$cohort)){
      cmdata = omdata.ssgsea[omdata.ssgsea$cohort == cohort,]
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
      other.models.validation.partial[1,paste0(cohort, c("_corr_iso_up", "_corr_iso_down"))] = c(corr_a, corr_b)
      other.models.validation.partial[1,paste0(cohort, "_valid")]=corr_a > thr & corr_b < -thr
    }
    
    other.models.validation = rbind(other.models.validation,other.models.validation.partial)
                              
    if ( curr_gene_name %in% genes_to_print ){
      toplot = rbind(data.frame(name = iso_up,
                             counts = oncount.txi[iso_up, omdata.ssgsea$name],
                             perc = tx.contributions[iso_up, omdata.ssgsea$name],
                             sample = omdata.ssgsea$name,
                             sampleType = omdata.ssgsea$type,
                             cohort = omdata.ssgsea$cohort,
                             target = omdata.ssgsea[,target]
      ), 
      data.frame(name = iso_down,
                 counts = oncount.txi[iso_down, omdata.ssgsea$name],
                 perc = tx.contributions[iso_down, omdata.ssgsea$name],
                 sample = omdata.ssgsea$name,
                 sampleType = omdata.ssgsea$type,
                 cohort = omdata.ssgsea$cohort,
                 target = omdata.ssgsea[,target]
      ))
      
      for ( cohort in unique(toplot$cohort)) {
        plots[[paste0(target, "_", cohort)]] =
          ggplot(toplot[toplot$cohort == cohort,], aes(x=target, y=perc)) + geom_point(aes(fill=sampleType), size=3, shape=21)+geom_smooth()+theme_bw()+xlab(target) +
          facet_grid(rows=vars(name)) + 
          ggtitle(paste0(cohort,"\nGene = " ,curr_gene_name, "\nTarget = " , target, "\n", iso_up, " (up) = ", other.models.validation.partial[1,paste0(cohort, "_corr_iso_up")], "\n", iso_down, " (down) = ", other.models.validation.partial[1,paste0(cohort, "_corr_iso_down")]))
      }
      
    }
  }
  
  if (curr_gene_name %in% genes_to_print ){
    pdf(paste0(odir, "/", curr_gene_name, ".pdf"))
    for ( tname in names(plots)){
      print(plots[[tname]])
    }
    dev.off()
  }
}



write.csv(other.models.validation, paste0(outdir, "/other_models_results.csv"), row.names = F)
colnames(other.models.validation)
other.models.validation.by_gene = other.models.validation %>% group_by(gene_id, gene_name) %>% 
  summarise(tot_target = n() ,
            perc_valid_CellLine = sum(CellLine_valid) / n(), 
            valid_targets_CellLine = sum(CellLine_valid), 
            valid_targets_CellLine_names = paste(target[CellLine_valid], sep=",", collapse = ","),
            perc_valid_LTL_longitudinal = sum(LTL_longitudinal_valid) / n(), 
            valid_targets_LTL_longitudinal = sum(LTL_longitudinal_valid), 
            valid_targets_LTL_longitudinal_names = paste(target[LTL_longitudinal_valid], sep=",", collapse = ","),
            perc_valid_PDX_longitudinal = sum(PDX_longitudinal_valid) / n(), 
            valid_targets_PDX_longitudinal = sum(PDX_longitudinal_valid), 
            valid_targets_PDX_longitudinal_names = paste(target[PDX_longitudinal_valid], sep=",", collapse = ",")
            ) 

write.csv(other.models.validation.by_gene, paste0(outdir, "/other_models_by_gene.csv"), row.names = F)
write.csv(other.models.validation.by_gene[other.models.validation.by_gene$gene_name %in% genes_to_print,], paste0(outdir, "/other_models_by_gene_candidates.csv"), row.names = F)



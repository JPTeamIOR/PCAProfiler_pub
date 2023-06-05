#!/usr/bin/Rscript


neededlibraries <- c( "tximport",  "reshape2", "readr", "DESeq2", "ggplot2", "dplyr", "VennDiagram", "plotly", "msigdbr", "limma", "GSVA")
loaded_libraries = lapply(neededlibraries, require, character.only = TRUE)
palette_subtypes = c("#25681E", "#F63A21","#4889C7", "#1010C0", "#93355A")
heatmapSamples = function(indata, ngenes=500){
  rv_all= rowVars(indata)
  pcagenes_all = rownames(indata)[order(rv_all, decreasing = TRUE)]
  selected_all = pcagenes_all[1:ngenes]
  corr_mat = round(cor(indata[selected_all,]),2)
  heatmap(corr_mat)
}

plotCount= function(mat, rowN, fnames){
  ggplot(data.frame(mRNA=mat[rowN,], grp=fnames))+geom_boxplot(aes(x=grp, y=mRNA))
}


volcanoPlot = function(res, levelNames){
  tmp = as.data.frame(res)
  tmp$padj[is.na(tmp$padj)] = 1
  tmp$mLog10Padj=-log10(tmp$padj)
  tmp$DE="-"
  tmp$DE[tmp$padj < 0.05 & tmp$log2FoldChange < -1  ]=levelNames[2]
  tmp$DE[tmp$padj < 0.05 & tmp$log2FoldChange > 1 ]=levelNames[1]
  tmp$DE = factor(tmp$DE, levels=c(levelNames, "-") )
  ggplot(tmp) + geom_point(aes(y=mLog10Padj, x=log2FoldChange, col=DE))+ theme_minimal() +ggtitle(paste(levelNames, collapse = " VS "))    
  
}


qOverlap = function(a,b){
  qA = quantile(a)
  qB = quantile(b)
  
  numA = rep(F, length(a))
  numB = rep(F, length(b))
  
  if (qA[2] >=  qB[2]){
    numB = b >= qA[2]
  } else {
    numA = a >= qB[2]
  }
  if (qB[4] <= qA[4]){
    numA = a <= qB[4] | numA
  } else {
    numB = b <= qA[4] | numB
  }
  return( (sum(numA) + sum(numB) ) / (length(a) + length(b)))
}


wilcoxTestDGE = function(count_norm, conditions, factorNames){
  ## factor 1 : target, factor 2: control
  conditions = factor(as.character(conditions), levels = factorNames)
  conditionsLevel<-levels(conditions)
  condA= conditions == conditionsLevel[1]
  condB= conditions == conditionsLevel[2]
  pvalues <- sapply(1:nrow(count_norm),function(i){
    return(wilcox.test(count_norm[i,condA],count_norm[i,condB])$p.value)
  })
  
  fdr=p.adjust(pvalues,method = "fdr")
  tot_l = length(conditions)
  log2FoldChanges=log2((2^rowMeans(count_norm[,condA]))/(2^rowMeans(count_norm[,condB])))
  outRst<-data.frame(log2FoldChange=log2FoldChanges, pValues=pvalues, padj = fdr)
  rownames(outRst)=rownames(count_norm)
  outRst=na.omit(outRst)
  return(outRst)  
}



myFormat = function(v){
  if ( abs(v) < 0.005 ){
    return(formatC(v, digits=2, format="e"))
  } else {
    return(round(v, digits=3))
  }
}

compareDW = function (a,b, FDRThr=0.05, FC_thr=1 ){
  if ( FC_thr > 0 ){
    x = list(
      rownames(a)[a$log2FoldChange > FC_thr & a$padj < FDRThr], 
      rownames(b)[b$log2FoldChange > FC_thr & b$padj < FDRThr] )  
  } else {
    x = list(
      rownames(a)[a$log2FoldChange < FC_thr & a$padj < FDRThr], 
      rownames(b)[b$log2FoldChange < FC_thr & b$padj < FDRThr] )  
  }
  grid.newpage()
  draw.pairwise.venn(area1 = length(x[[1]]), area2= length(x[[2]]), cross.area = sum(x[[2]] %in% x[[1]]))
  missing=x[[1]][! x[[1]] %in% x[[2]]]
  tmp=a[missing,]
  tmp=tmp[order(tmp$padj),]
  missing= rownames(tmp[1:min(10, nrow(tmp)),])
  
  g=missing[2]
  for (g in missing){
    dat=data.frame(mRNA= ncount[g,], type=mdata$PC.Type)
    dat2 <- dat %>% group_by(type) %>% summarise(mRNA=mean(mRNA))
    p<-ggplot(dat)+geom_boxplot(aes(x=type, y=mRNA))+ geom_point(data =dat2 , aes(x=type, y=mRNA), color="red")+
      ggtitle(paste0(g,  "\nFC Deseq= ",myFormat(a[g,]$log2FoldChange), "\nFC W = ",  myFormat(b[g,]$log2FoldChange) ,    "\npajd DEseq= ",myFormat(a[g,]$padj), " \npadj Wilcox= ",myFormat(b[g,]$padj) ))
    print(p)
  }
  missing=x[[2]][! x[[2]] %in% x[[1]]]
  tmp=b[missing,]
  tmp=tmp[order(tmp$padj),]
  missing= rownames(tmp[1:min(10, nrow(tmp)),])
  grid.newpage()
  for (g in missing){
    dat=data.frame(mRNA= ncount[g,], type=mdata$PC.Type)
    dat2 <- dat %>% group_by(type) %>% summarise(mRNA=mean(mRNA))
    p<-ggplot(dat)+geom_boxplot(aes(x=type, y=mRNA))+ geom_point(data =dat2 , aes(x=type, y=mRNA), color="red")+
      ggtitle(paste0(g,  "\nFC Deseq=",myFormat(a[g,]$log2FoldChange), "\nFC W = ",  myFormat(b[g,]$log2FoldChange) , "\npajd DEseq= ",myFormat(a[g,]$padj), " \npadj Wilcox= ",myFormat(b[g,]$padj) ))
    print(p)
  }
}



find_Interesting_DTE <- function (tx_data, g_data, tx_a_data,g_a_data, classes, tx_info ,ov_thr= 0.25,FC_thr=2, FDR_thr = 0.01, out_dir = "./out_interesting"){
  tmp = tx_data
  dir.create(out_dir, showWarnings = F, recursive = T)
  colnames(tx2gene)=c("transcript_id", "gene_id", "gene_name")
  tx_info = merge(tx_info, tx2gene, by=c("transcript_id", "gene_id"))
  tmp$transcript_id = rownames(tmp)
  tmp = merge(tmp, tx_info, by=c("transcript_id"),all.x=T, all.y=F)
  tmp = merge(tmp, data.frame(gene_id = rownames(g_data), padj_gene = g_data$padj, log2FC_gene = g_data$log2FoldChange), by=c("gene_id"), all.x=T, all.y=F)
  tmp$TX_DE = abs(tmp$log2FoldChange) > FC_thr & tmp$padj < FDR_thr  & ! is.na(tmp$padj)
  tmp$G_DE = abs(tmp$log2FC_gene) > FC_thr & tmp$padj_gene < FDR_thr & ! is.na(tmp$padj_gene)
  tmp$TX_Dir = "UP"
  tmp$TX_Dir[tmp$log2FoldChange < 0 ] = "DOWN"
  pdf(paste0(out_dir, "/Global.pdf"))
  p<- ggplot(tmp[tmp$TX_DE ,])+geom_bar(aes(x=transcript_type, fill=TX_Dir), position = position_dodge())+coord_flip()
  print(p)
  p<- ggplot(tmp[tmp$TX_DE & ! tmp$G_DE ,])+geom_bar(aes(x=transcript_type, fill=TX_Dir), position = position_dodge())+coord_flip()
  print(p)
  dev.off()
  tmp = tmp[order(tmp$padj),]
  out = list(TX_NO_G = c(), G_NO_TX = c(), DISCORDANT_TX = c() )
  
  for ( gid in unique(tmp$gene_id[tmp$TX_DE | tmp$G_DE])){
    tt = tmp[tmp$gene_id == gid, ]
    if ( nrow(tt) > 1 & "UP" %in% tt$TX_Dir[tt$TX_DE] & "DOWN" %in% tt$TX_Dir[tt$TX_DE] ){
      out[["DISCORDANT_TX"]]=c(out[["DISCORDANT_TX"]], gid)
    }
    if ( any(tt$TX_DE) & ( ! any(tt$G_DE) ) ){
        out[["TX_NO_G"]]= c(out[["TX_NO_G"]] ,gid)
    }
    if ( any(tt$g_DE) & ( ! any(tt$TX_DE) ) ){
      out[["G_NO_TX"]]= c(out[["G_NO_TX"]] ,gid)
    }
  }
  for ( k in names(out)){
    if (length(out[[k]])> 0){
      pdf(paste0(out_dir,"/", k, ".pdf"))
        for ( gid in out[[k]]){
          tt = tmp[tmp$gene_id == gid, ]
          if ( k == "DISCORDANT_TX" ||  k == "TX_NO_G"){
            tt = tt[tt$TX_DE,]
          } 
          dat = data.frame(grp = classes, count=g_a_data[gid,], name=gid)
          for ( tx in tt$transcript_id){
            dat=rbind(dat, data.frame(grp=classes, count = tx_a_data[tx,], name=tx))
          }
          p<-ggplot(dat)+geom_boxplot(aes(x=name, y=count, fill=grp))+ ggtitle(unique(tt$gene_name))
          print(p)
        }
      dev.off()
    }
  }
  
  return(out)
}

getSign <- function(df, f_name){
  out=NULL
  if ( f_name %in% rownames(df)){
    if ( df[f_name,]$padj < 0.05 & abs(df[f_name,]$log2FoldChange) > 1) {
      out="*"
      if (df[f_name,]$padj < 0.005 ){
        out="**"  
      } 
    }  
  }
  
  return(out)
}


dotPlot = function(res, out_file, plot_dir="./plots/"){
  tmp = res@result[res@result$p.adjust < 0.05,]
  if ( nrow(tmp) > 0 ){
    dir.create(plot_dir, showWarnings = F)
    tmp= tmp[order(tmp$p.adjust, decreasing = T),]
    tmp$Description = gsub(pattern = ".*\r: " ,replacement = "" , tmp$Description)
    tmp$Description = gsub(pattern = "_" ,replacement = " " , tmp$Description)
    tmp$Description = gsub(pattern = "^(.{30,40}) " ,replacement = "\\1\n" , tmp$Description)
    tmp$Description = gsub(pattern = "\n(.{30,40}) " ,replacement = "\n\\1\n" , tmp$Description)
    tmp$GeneRatio = unlist(lapply(tmp$GeneRatio, function(x){ as.numeric(strsplit(x, "/")[[1]][[1]]) /as.numeric(strsplit(x, "/")[[1]][[2]])  }))
    tmp = tmp[order(tmp$Count, decreasing = T),]
    tmp$Description = factor(tmp$Description, levels = tmp$Description)
    pdf(paste0(plot_dir, out_file, ".pdf"))
    for ( i in seq(1 , nrow(tmp), by=10)){
      p <- ggplot(tmp[i:min((i+10), nrow(tmp)),])+ geom_point(
        aes(x=GeneRatio, y=Description , color=p.adjust, size=Count))+ scale_y_discrete(limits=rev)+
        scale_color_gradient(low = "blue", high = "red", limits = c(min(tmp$p.adjust), max(tmp$p.adjust))) + xlim(c(0, max(tmp$GeneRatio)))
      print(p)
    }
    dev.off()
    write.csv(tmp, paste0(plot_dir, out_file, ".csv"), row.names = F)
  } else {
    print("No enrichment terms found")
  }
}



dotPlotList = function(res_list, out_file, plot_dir="./plots/", width_factor=2){
  all_res = NULL
  for ( k in names(res_list)){
    tmp = res_list[[k]]@result[res_list[[k]]@result$p.adjust < 0.05,]
    if ( nrow(tmp) > 0 ){
      dir.create(plot_dir, showWarnings = F)
      tmp= tmp[order(tmp$p.adjust, decreasing = T),]
      tmp$Description = gsub(pattern = ".*\r: " ,replacement = "" , tmp$Description)
      tmp$Description = gsub(pattern = "_" ,replacement = " " , tmp$Description)
      tmp$Description = gsub(pattern = "^(.{10,30}) " ,replacement = "\\1\n" , tmp$Description)
      tmp$Description = gsub(pattern = "\n(.{10,30}) " ,replacement = "\n\\1\n" , tmp$Description)
      tmp$Description = gsub(pattern = "\n(.{10,30}) " ,replacement = "\n\\1\n" , tmp$Description)
      tmp$GeneRatio = unlist(lapply(tmp$GeneRatio, function(x){ as.numeric(strsplit(x, "/")[[1]][[1]]) /as.numeric(strsplit(x, "/")[[1]][[2]])  }))
      tmp$Name = k
      all_res = rbind(all_res, tmp)
    } 
  }
  if ( !is.null(all_res) ){
    all_res = all_res[order(all_res$p.adjust),]
    all_res$Name = factor(all_res$Name, levels=names(res_list))
    IDS= unique(all_res$ID)
    pdf(paste0(plot_dir, out_file, ".pdf"), width =  round(length(names(res_list))*width_factor) )
    ydiv = 20
    for ( i in seq(1 , length(IDS), by=ydiv)){
      to_use = IDS[i:min((i+ydiv-1), length(IDS))]
      p <- ggplot(all_res[all_res$ID %in% to_use,])+ geom_point(
        aes(x=Name, y=Description , color=-log10(p.adjust), size=GeneRatio))+ scale_y_discrete(limits=rev)+scale_x_discrete(drop=F)+
        scale_color_gradient(low = "blue", high = "red", limits = c(min(-log10(all_res$p.adjust)), max(-log10(all_res$p.adjust))))
      print(p)
    }
    dev.off()
    write.csv(all_res, paste0(plot_dir, out_file, ".csv"), row.names = F)  
  }
  
}



dotPlotList2 = function(res_list){
  all_res = NULL
  out=list()
  for ( k in names(res_list)){
    tmp = res_list[[k]]@result
    if ( nrow(tmp) > 0 ){
      tmp= tmp[order(tmp$p.adjust, decreasing = T),]
      tmp$Description = gsub(pattern = ".*\r: " ,replacement = "" , tmp$Description)
      tmp$Description = gsub(pattern = "_" ,replacement = " " , tmp$Description)
      tmp$Description = gsub(pattern = "^(.{10,30}) " ,replacement = "\\1\n" , tmp$Description)
      tmp$Description = gsub(pattern = "\n(.{10,30}) " ,replacement = "\n\\1\n" , tmp$Description)
      tmp$Description = gsub(pattern = "\n(.{10,30}) " ,replacement = "\n\\1\n" , tmp$Description)
      tmp$GeneRatio = unlist(lapply(tmp$GeneRatio, function(x){ as.numeric(strsplit(x, "/")[[1]][[1]]) /as.numeric(strsplit(x, "/")[[1]][[2]])  }))
      tmp$Name = k
      all_res = rbind(all_res, tmp)
    } 
  }
  if ( !is.null(all_res) ){
    all_res = all_res[order(all_res$p.adjust),]
    all_res$Name = factor(all_res$Name, levels=names(res_list))
    IDS= unique(all_res$ID)
    for(grp in unique(h_groups$group) ){
      all_res$signif = ifelse(all_res$p.adjust < 0.05, 1, 0.2)
      p <- ggplot(all_res[all_res$ID %in% h_groups$hallmark[h_groups$group==grp],])+ geom_point(
        aes(x=Name, y=Description , color=-log10(p.adjust), size=GeneRatio,alpha=signif))+ scale_y_discrete(limits=rev)+scale_x_discrete(drop=F)+
        scale_color_gradient(low = "blue", high = "red", limits = c(min(-log10(all_res$p.adjust)), max(-log10(all_res$p.adjust))))+ ggtitle(grp)
      out[[grp]]=p
    }
    out[["all_res"]]=all_res
    
  }
  return(out)
  
}


dotPlotListGSEA = function(res_list, out_file, plot_dir="./plots/"){
  all_res = NULL
  for ( k in names(res_list)){
    tmp = res_list[[k]]@result[res_list[[k]]@result$p.adjust < 0.05,]
    if ( nrow(tmp) > 0 ){
      dir.create(plot_dir, showWarnings = F)
      tmp= tmp[order(tmp$p.adjust, decreasing = T),]
      tmp$Description = gsub(pattern = ".*\r: " ,replacement = "" , tmp$Description)
      tmp$Description = gsub(pattern = "_" ,replacement = " " , tmp$Description)
      tmp$Description = gsub(pattern = "^(.{10,30}) " ,replacement = "\\1\n" , tmp$Description)
      tmp$Description = gsub(pattern = "\n(.{10,30}) " ,replacement = "\n\\1\n" , tmp$Description)
      tmp$Description = gsub(pattern = "\n(.{10,30}) " ,replacement = "\n\\1\n" , tmp$Description)
      tmp$Name = k
      all_res = rbind(all_res, tmp)
    } 
  }
  if ( !is.null(all_res) ){
    all_res = all_res[order(all_res$p.adjust),]
    all_res$Name = factor(all_res$Name, levels=names(res_list))
    IDS= unique(all_res$ID)
    pdf(paste0(plot_dir, out_file, ".pdf"), width = max(10, length(names(res_list))))
    ydiv = 20
    for ( i in seq(1 , length(IDS), by=ydiv)){
      to_use = IDS[i:min((i+ydiv-1), length(IDS))]
      
      p <- ggplot(all_res[all_res$ID %in% to_use,])+ geom_point(
        aes(x=Name, y=Description , color=enrichmentScore, size=-log10(p.adjust) ))+ scale_y_discrete(limits=rev)+scale_x_discrete(drop=F)+
        scale_color_gradient2(low = "blue", high = "red", limits = c(min(all_res$enrichmentScore), max(all_res$enrichmentScore)))
      print(p)
    }
    dev.off()
    write.csv(all_res, paste0(plot_dir, out_file, ".csv"), row.names = F)  
  }
}


h_groups = rbind(
  data.frame(group = "Inflammatory" , hallmark=c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_COAGULATION"
)  ), data.frame(group = "Damage Response" , hallmark=c(
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_UV_RESPONSE_UP",
  "HALLMARK_UV_RESPONSE_DN"
)))

h_groups = rbind(h_groups, data.frame(group = "Cell Proliferation" , hallmark=c(
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2"
)))

h_groups = rbind(h_groups, data.frame(group = "Cell Metabolism" , hallmark=c(
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_HYPOXIA",
  "HALLMARK_PEROXISOME",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_XENOBIOTIC_METABOLISM",
  "HALLMARK_BILE_ACID_METABOLISM",
  "HALLMARK_HEME_METABOLISM"
)))

h_groups = rbind(h_groups, data.frame(group = "Cell Junction" , hallmark=c(
  "HALLMARK_APICAL_JUNCTION",
  "HALLMARK_APICAL_SURFACE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)))

h_groups = rbind(h_groups, data.frame(group = "Differentiation" , hallmark=c(
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_ADIPOGENESIS",
  "HALLMARK_MYOGENESIS",
  "HALLMARK_SPERMATOGENESIS",
  "HALLMARK_PANCREAS_BETA_CELLS"
)))

h_groups = rbind(h_groups, data.frame(group = "Signaling" , hallmark=c(
  "HALLMARK_ANDROGEN_RESPONSE",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_ESTROGEN_RESPONSE_LATE",
  "HALLMARK_NOTCH_SIGNALING",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_HEDGEHOG_SIGNALING",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_KRAS_SIGNALING_DN",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_PROTEIN_SECRETION"
)))


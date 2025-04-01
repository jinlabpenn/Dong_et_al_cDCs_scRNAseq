library(clusterProfiler)
library(org.Mm.eg.db)
library(scProgram)
library(viridis)
library(tidyverse)
#load rds file
#The rds can download from https://www.dropbox.com/scl/fi/8yufebal0uae30t5ftqwn/scRICA_qc_integration_allCond_auto.rds?rlkey=hdlpurnvbvjulg83dy0y86603&st=vq49hww5&dl=0
rds <- readRDS("scRICA_qc_integration_allCond_auto.rds")
#filter DC
DC.proj <- subset(rds,subset = seurat_clusters %in% c(0:9,11:14))
#re-order cluster ids
dat <- data.frame(
  seurat_cluster = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14),
  new_cluster    = c(11, 1, 6, 7, 9, 13, 8, 2, 4, 12, 3, 10, 5, 14)
)
DC.proj@meta.data$seurat_cluster <- as.numeric(as.character(DC.proj@meta.data$seurat_clusters))
table(DC.proj@meta.data$seurat_cluster)
tmp <- DC.proj@meta.data
tmp$cell <- rownames(tmp)
tmp <- merge(tmp, dat, by = 'seurat_cluster')
rownames(tmp) <- tmp$cell
DC.proj@meta.data$new_cluster <- as.factor(tmp[rownames(DC.proj@meta.data), "new_cluster"])

DefaultAssay(DC.proj) <- 'RNA'
DC.proj@meta.data$celltype <- 'N'
DC.proj@meta.data$celltype[DC.proj@meta.data$new_cluster %in% c(1:5)] <- 'cDC1'
DC.proj@meta.data$celltype[DC.proj@meta.data$new_cluster %in% c(6:11)] <- 'cDC2'
DC.proj@meta.data$celltype[DC.proj@meta.data$new_cluster==12] <- 'migDC1'
DC.proj@meta.data$celltype[DC.proj@meta.data$new_cluster==13] <- 'migDC2'
DimPlot(DC.proj,group.by = 'celltype',label = T)
table(DC.proj@meta.data$expCond1)

marker.dc.cond <- data.frame()
for (i in c('cDC1','cDC2','migDC1','migDC2')){
  tmp.seu <- subset(DC.proj,subset=celltype==i)
  DefaultAssay(tmp.seu) <- 'RNA'
  tmp.seu@meta.data$expCond1 <- factor(tmp.seu@meta.data$expCond1,levels = c('SPF','GF'))
  Idents(tmp.seu) <- 'expCond1'
  marker.all.dc <- FindAllMarkers(tmp.seu)
  marker.all.dc$condtion <- marker.all.dc$cluster
  marker.all.dc$celltype <- i
  marker.dc.cond <- rbind(marker.dc.cond,marker.all.dc)
}
table(marker.dc.cond$celltype)
topG <- marker.dc.cond[marker.dc.cond$p_val_adj < 0.05&
                         marker.dc.cond$avg_log2FC>0.75&marker.dc.cond$pct.1>0.1,] %>%
  as.data.frame()
table(topG$cluster,topG$celltype)
#load gmt file
kegmt <- read.gmt("msigdb.v2023.2.Mm.symbols.gmt")
#GSEA
ls <- list()
pqvalueCutoff <- 0.99#show all terms
for (i in c('cDC1','migDC1')){
  for (j in c('GF','SPF')){
    topG <- marker.dc.cond[marker.dc.cond$celltype==i,]
    gene <- topG[topG$condtion==j,]
    gene <- gene[order(gene$avg_log2FC,decreasing = T),]
    tmp <- gene$avg_log2FC
    names(tmp) <- gene$gene
    gsea <- GSEA(tmp,TERM2GENE = kegmt,pvalueCutoff = 0.99)
    ls[[paste0(i,'-',j)]] <- gsea
  }
}
datGSEA <- data.frame()
names(ls)
for (i in 1:length(ls)){
  tmp <- ls[[i]]@result
  tmp$cluster <- names(ls)[i]
  datGSEA <- rbind(datGSEA,tmp)
}
colnames(datGSEA)
datGSEA <- datGSEA[datGSEA$p.adjust<0.1,]
datGSEA <- datGSEA[datGSEA$enrichmentScore>0.5,]
datGSEA$Description <- tolower(datGSEA$Description)
datGSEA$Description <- gsub('_',' ',datGSEA$Description)
datGSEA$celltype <- gsub('-.*','',datGSEA$cluster)
datGSEA$condition <- gsub('.*-','',datGSEA$cluster)
#GO
ls <- list()
pqvalueCutoff <- 0.99
for (i in c('cDC1','migDC1')){
  for (j in c('GF','SPF')){
    marker.all.dc <- marker.dc.cond[marker.dc.cond$celltype==i,]
    topG <- marker.all.dc[marker.all.dc$p_val_adj < 0.05&
                            marker.all.dc$avg_log2FC>0.75&marker.all.dc$pct.1>0.1,] %>% as.data.frame()
    gene <- topG[topG$cluster==j,]$gene
    trans <- bitr(gene, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
    ego <- enrichGO(gene  = trans$ENTREZID,keyType = "ENTREZID", 
      OrgDb = org.Mm.eg.db,ont = "BP",pAdjustMethod = "BH",minGSSize = 3,
      pvalueCutoff = pqvalueCutoff, qvalueCutoff = pqvalueCutoff,readable = TRUE)
    ls[[paste0(i,'-',j)]] <- ego 
  }
}
datGO <- data.frame()
for (i in 1:length(ls)){
  tmp <- ls[[i]]@result
  tmp$cluster <- names(ls)[i]
  datGO <- rbind(datGO,tmp)
}
datGO <- datGO[datGO$p.adjust<0.1,]
datGO$Description <- tolower(datGO$Description)
datGO$Description <- gsub('_',' ',datGO$Description)
datGO$celltype <- gsub('-.*','',datGO$cluster)
datGO$condition <- gsub('.*-','',datGO$cluster)
#Figure3D selected terms to show####
dat <- as.data.frame(fread('gene_functional_terms.txt'))
dat$p <- -log10(dat$p.adjust)
dat[dat$condition=='SPF',]$p <- -(dat[dat$condition=='SPF',]$p)
dat <- dat[order(dat$celltype,dat$condition,abs(dat$p)),]
dat$Description <- factor(dat$Description,levels = unique(dat$Description))
ggplot(dat,aes(Description,p,fill=celltype))+
  geom_col()+coord_flip()+facet_wrap(~celltype)

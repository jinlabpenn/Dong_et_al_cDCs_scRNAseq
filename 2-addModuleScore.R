library(data.table)
library(Seurat)
library(RColorBrewer)
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
#FigureS3F markers from literature####
dat <- as.data.frame(fread('annotation-markerOfDC.txt'))
for (i in colnames(dat)){
  DC_features <- list(dat[,i][dat[,i]!=''])
  DC.proj <- AddModuleScore(
    object = DC.proj,
    features = DC_features,
    name = i)
  print(i)
}
head(DC.proj@meta.data)
colnames(DC.proj@meta.data)
ids <- grep('DC',colnames(DC.proj@meta.data))
colnames(DC.proj@meta.data)[ids] <- substring(colnames(DC.proj@meta.data)[ids], 1, nchar(colnames(DC.proj@meta.data)[ids]) - 1)
FeaturePlot(DC.proj, feature = colnames(DC.proj@meta.data)[ids])
#FigureS3H heatmap####
cdc1 <- subset(DC.proj,subset=celltype=='cDC1')
DimPlot(cdc1,label = T)
DefaultAssay(cdc1) <- 'RNA'
cdc1 <- FindVariableFeatures(cdc1)
cdc1 <- ScaleData(cdc1,features = rownames(cdc1))
marker.all.cdc1 <- FindAllMarkers(cdc1)
topG <- marker.all.cdc1[marker.all.cdc1$p_val_adj < 0.01& 
                          marker.all.cdc1$avg_log2FC>1&marker.all.cdc1$pct.1>0.1,] %>%
  group_by(cluster) %>% top_n(n = 15, avg_log2FC) %>% as.data.frame()
table(topG$cluster)
DoHeatmap(subset(cdc1,downsample=200),features = unique(topG[order(topG$cluster),]$gene))
DoHeatmap(cdc1,features = unique(topG[order(topG$cluster),]$gene),
          group.colors = c("#9ECAE1", "#6BAED6", "#08519C", "#FB9A99", "#E31A1C"))
#The results of GO and GSEA were produced in 1-geneFunctionalAnnotation.R
#Figure3E and S3I selected terms to show AMS####
dat <- as.data.frame(fread('gene_function_terms_AMS.txt'))
for (i in dat$Description){
  genes <- unlist(strsplit(unique(dat[dat$Description==i,]$geneID),'/'))
  print(length(genes))
  genes <- list(genes)
  cdc1 <- AddModuleScore(
    object = cdc1,
    features = genes,
    name = i)
  ids <- ncol(cdc1@meta.data)
  colnames(cdc1@meta.data)[ids] <- substring(colnames(cdc1@meta.data)[ids],
                                             1, nchar(colnames(cdc1@meta.data)[ids]) - 1)
  print(i)
}

cdc1@meta.data$V1 <- cdc1@reductions$umap@cell.embeddings[,1]
cdc1@meta.data$V2 <- cdc1@reductions$umap@cell.embeddings[,2]
ls <- list()
for(i in 1:length(terms)){
  ls[[i]] <- FeaturePlot(subset(cdc1,subset=V1>2.5 & V2<2),features = terms[i])+
    scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
}
patchwork::wrap_plots(ls[1:length(ls)])

#Figure3I cell cycle####
datMH <- as.data.frame(fread('mouse2Human.txt'))
datCC <- as.data.frame(fread('cellCycleGenes.txt'))
datCC <- merge(datCC,datMH,by.x='Gene',by.y='Human gene stable ID')
datCC$Phase <- factor(datCC$Phase,levels = c('G1/S','S','G2/M','M','M/G1'))
datCC <- datCC[order(datCC$Phase),]

for (i in as.character(unique(datCC$Phase))){
  genes <- unique(datCC[datCC$Phase==i,]$`Gene name`)
  print(length(genes))
  genes <- list(genes)
  cdc1 <- AddModuleScore(
    object = cdc1,
    features = genes,
    name = i)
  ids <- ncol(cdc1@meta.data)
  colnames(cdc1@meta.data)[ids] <- substring(colnames(cdc1@meta.data)[ids],
                                             1, nchar(colnames(cdc1@meta.data)[ids]) - 1)
  print(i)
}

colnames(cdc1@meta.data)
terms <- c('G1/S','S','G2/M','M','M/G1')
ls <- list()
for(i in 1:length(terms)){
  ls[[i]] <- FeaturePlot(subset(cdc1,subset=V1>2.5 & V2<2),features = terms[i])+
    scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
}
patchwork::wrap_plots(ls[1:length(ls)])

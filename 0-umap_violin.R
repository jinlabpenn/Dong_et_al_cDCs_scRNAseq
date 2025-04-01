library(Seurat)
library(ggplot2)
library(ggsignif)
library(data.table)
library(Seurat)
library(patchwork)
library(ggpubr)
theme_set(theme_bw())
#Figure3A UMAP plot for DC####
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
Idents(DC.proj) <- 'new_cluster'
colors <- c("#9ECAE1", "#6BAED6", "#3182BD", "#FB9A99", "#EF5A5A",
            "#A1D99B", "#41AB5D", "#238B45", "#FFED6F", "#B15928", "#CAB2D6", #"#6A3D9A",
            "#FF7F00","#FDBF6F", 'grey80')
DimPlot(DC.proj,label = T,group.by = 'new_cluster',pt.size = 1.1)+
  scale_color_manual(values = colors)+coord_fixed()

#Figure3B Violin plot####
my_comparisons <- list(c(12,13))
modify_vlnplot<- function(obj,feature, 
                          pt.size = 0,plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    stat_compare_means(comparisons = my_comparisons,vjust = 0.2,method = 'wilcox.test')+ # Add pairwise comparisons p-value
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin )+
    scale_fill_manual(values = colors)
  return(p)
}
## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}
## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm"),
                          ...) {
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
   plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
   ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}
genes <- c('Xcr1','Irf8','Clec9a','Sirpa','Irf4','Ccr7')
StackedVlnPlot(obj = subset(DC.proj,subset = new_cluster!=14), features = genes)

#Figure3C Dot plot####
dat <- DC.proj@meta.data
dat$ct_cluster <- paste0(dat$new_cluster)
dat$ct_cluster <- factor(dat$ct_cluster,levels = unique(dat[order(dat$new_cluster),]$ct_cluster))
head(dat)
tmp1 <- as.data.frame(table(dat$expCond,dat$ct_cluster))
tmp2 <- as.data.frame(table(dat$expCond))
tmp <- merge(tmp1,tmp2,by='Var1')
tmp$Percentage <- tmp$Freq.x/tmp$Freq.y
colnames(tmp)[1] <- 'Rep'
colnames(tmp)[2] <- 'Cluster'
tmp$Condition <- gsub('[0-9]','',tmp$Rep)
tmp$Condition <- factor(tmp$Condition,levels = c('SPF','GF'))
ggplot(tmp,aes(Cluster,Percentage,color=Condition))+
  geom_jitter()+
  scale_color_manual(values =  c("black","#ff0080"))

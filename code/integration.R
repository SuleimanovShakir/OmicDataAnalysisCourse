library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggpubr)
library(UpSetR)
library(ComplexUpset)
library(reshape)
library(ggplot2)
library(dplyr)
source('src.R')


dataset_diff_meth_all <- read.csv('../data/diff_meth_genes_cpg.csv', sep=' ')
dataset_meth_all <- read.csv('../data/meth_genes_cpg.csv', sep=' ')
dataset_meth_all_cpg_island <- read.csv('../data/meth_genes_cpg_islands.csv', sep=' ')
dataset_rnaseq <- read.csv('../data/results_significant_tb_without_LFC.csv')[2:9]

intersect(dataset_rnaseq$gene, dataset_diff_meth_all$nearest_gene)

#-------Scatterplot RNA-seq vs DiffMeth all
rnaseq_meth_intersection <- intersect(dataset_rnaseq$gene, dataset_meth_all$nearest_gene)

colnames(dataset_meth_all) <- c("chr","start","end","strand","pvalue","qvalue","meth.diff","nearest_tss","gene")
data_for_plot <- merge(dataset_rnaseq[(dataset_rnaseq$gene %in% rnaseq_meth_intersection), ], dataset_meth_all[(dataset_meth_all$gene %in% rnaseq_meth_intersection), ], by='gene')

ggplot(data_for_plot, aes(x=as.factor(gene_symbol), y=meth.diff, fill = log2FoldChange))+
  geom_violin(width=2.5) +
  theme_bw() +
  scale_fill_gradientn(colours = c('#fde725', '#f6b092', '#5ec962', '#21918c', '#3b528b')) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', linewidth=1)

ggplot(data_for_plot, aes(x=as.factor(gene_symbol), y=meth.diff, fill = log2FoldChange))+
  geom_boxplot() +
  theme_bw() +
  scale_fill_gradientn(colours = c('#fde725', '#f6b092', '#5ec962', '#21918c', '#3b528b')) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', linewidth=1)

data_for_scatterplot <- data_for_plot %>%
  group_by(gene) %>%
  summarise(rnaLFC = median(log2FoldChange), meth.diff_median = median(meth.diff))

ggplot(data_for_scatterplot, aes(x=rnaLFC, y=meth.diff_median))+
  geom_point() +
  theme_bw() +
  geom_smooth(method='lm')
  #scale_fill_gradientn(colours = c('#fde725', '#f6b092', '#5ec962', '#21918c', '#3b528b')) +
  #geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', linewidth=1)


#-------Scatterplot RNA-seq vs DiffMeth CpG islands
rnaseq_meth_cpg_islands_intersection <- intersect(dataset_rnaseq$gene, dataset_meth_all_cpg_island$nearest_gene)

colnames(dataset_meth_all_cpg_island) <- c("chr","start","end","strand","pvalue","qvalue","meth.diff","nearest_tss","gene")
data_for_plot_cpg_islands <- merge(dataset_rnaseq[(dataset_rnaseq$gene %in% rnaseq_meth_cpg_islands_intersection), ], dataset_meth_all_cpg_island[(dataset_meth_all_cpg_island$gene %in% rnaseq_meth_cpg_islands_intersection), ], by='gene')


ggplot(data_for_plot_cpg_islands, aes(x=as.factor(gene_symbol), y=meth.diff, fill = log2FoldChange))+
  geom_violin(width=1.5) +
  theme_bw() +
  scale_fill_gradientn(colours = c('#fde725', '#f6b092', '#5ec962', '#21918c', '#3b528b')) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', linewidth=1)

ggplot(data_for_plot_cpg_islands, aes(x=as.factor(gene_symbol), y=meth.diff, fill = log2FoldChange))+
  geom_boxplot() +
  theme_bw() +
  scale_fill_gradientn(colours = c('#fde725', '#f6b092', '#5ec962', '#21918c', '#3b528b')) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', linewidth=1)


data_for_scatterplot_cpg_islands <- data_for_plot_cpg_islands %>%
  group_by(gene) %>%
  summarise(rnaLFC = median(log2FoldChange), meth.diff_median = median(meth.diff))

ggplot(data_for_scatterplot_cpg_islands, aes(x=rnaLFC, y=meth.diff_median))+
  geom_point() +
  theme_bw() +
  geom_smooth(method='lm')
  #scale_fill_gradientn(colours = c('#fde725', '#f6b092', '#5ec962', '#21918c', '#3b528b')) +
  #geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', linewidth=1)


#------------------------------------------------Annotate bed
diffChIP_to_bed(path_to_file = '../data/bed/diff_H3K27ac_padj0.05.csv')
annotate_bed(in_path = '../data/bed/diff_H3K27ac_padj0.05.csv.bed', out_dir = '../data/annotation', filename = 'H3K27ac')


#------------------------------------------------Concatenate annotation with data from bed
ENSEMBLE <- read.csv("../data/raw_counts.csv", sep = ';')[1:2]

H3K4me1_ann <- read.delim('../data/annotation/H3K4me1')
H3K4me1_csv <- read.csv('../data/bed/diff_H3K4me1_padj0.05.csv')
H3K4me1_data <- left_join(H3K4me1_ann, H3K4me1_csv, by = c("start" = "Start", "end" = "End"), relationship = "many-to-many")
#write.csv(H3K4me1_data, file='../data/annotation/H3K4me1.annotated.csv')


H3K4me3_ann <- read.delim('../data/annotation/H3K4me3')
H3K4me3_csv <- read.csv('../data/bed/diff_H3K4me3_padj0.05.csv')
H3K4me3_data <- left_join(H3K4me3_ann, H3K4me3_csv, by = c("start" = "Start", "end" = "End"), relationship = "many-to-many")
#write.csv(H3K4me3_data, file='../data/annotation/H3K4me3.annotated.csv')


H3K27me3_ann <- read.delim('../data/annotation/H3K27me3')
H3K27me3_csv <- read.csv('../data/bed/diff_H3K27me3_padj0.05.csv')
H3K27me3_data <- left_join(H3K27me3_ann, H3K27me3_csv, by = c("start" = "Start", "end" = "End"), relationship = "many-to-many")
#write.csv(H3K27me3_data, file='../data/annotation/H3K27me3.annotated.csv')


H3K36me3_ann <- read.delim('../data/annotation/H3K36me3')
H3K36me3_csv <- read.csv('../data/bed/diff_H3K36me3_padj0.05.csv')
H3K36me3_data <- left_join(H3K36me3_ann, H3K36me3_csv, by = c("start" = "Start", "end" = "End"), relationship = "many-to-many")
#write.csv(H3K36me3_data, file='../data/annotation/H3K36me3.annotated.csv')


H3K27ac_ann <- read.delim('../data/annotation/H3K27ac')
H3K27ac_csv <- read.csv('../data/bed/diff_H3K27ac_padj0.05.csv')
H3K27ac_data <- left_join(H3K27ac_ann, H3K27ac_csv, by = c("start" = "Start", "end" = "End"), relationship = "many-to-many")
#write.csv(H3K27ac_data, file='../data/annotation/H3K27ac.annotated.csv')


#-------RNA-seq vs H3K4me3
rnaseq_H3K4me3_intersection <- intersect(dataset_rnaseq$gene_symbol, H3K4me3_data$GENENAME)

data_for_plot_rna_seq_H3K4me3 <- 
  merge(dataset_rnaseq[(dataset_rnaseq$gene_symbol %in% rnaseq_H3K4me3_intersection), ], 
        H3K4me3_data[(H3K4me3_data$GENENAME %in% rnaseq_H3K4me3_intersection), ], 
        by.x = 'gene_symbol', by.y = 'GENENAME')

ggplot(data_for_plot_rna_seq_H3K4me3, aes(x=as.factor(gene_symbol), y=log2FC, fill = log2FoldChange))+
  geom_violin(width=2.5) +
  theme_bw() +
  scale_fill_gradientn(colours = c('#fde725', '#f6b092', '#5ec962', '#21918c', '#3b528b')) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', linewidth=1)

ggplot(data_for_plot_rna_seq_H3K4me3, aes(x=as.factor(gene_symbol), y=log2FC, fill = log2FoldChange))+
  geom_boxplot() +
  theme_bw() +
  scale_fill_gradientn(colours = c('#fde725', '#f6b092', '#5ec962', '#21918c', '#3b528b')) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', linewidth=1)

data_for_scatterplot_rna_seq_H3K4me3 <- data_for_plot_rna_seq_H3K4me3 %>%
  group_by(gene_symbol) %>%
  summarise(rnaLFC = median(log2FoldChange), H3K4me3 = median(log2FC))

ggplot(data_for_scatterplot_rna_seq_H3K4me3, aes(x=rnaLFC, y=H3K4me3))+
  geom_point(size=2, color = 'darkgreen') +
  theme_bw() +
  geom_smooth(method='lm', alpha = 0.3) +
  xlab('RNA_LFC') +
  ylab('H3K4me3_LFC') +
  stat_cor(method="spearman", cor.coef.name = c("rho"), 
           label.y = max(data_for_scatterplot_rna_seq_H3K4me3$H3K4me3) * 1.5, 
           label.x = max(data_for_scatterplot_rna_seq_H3K4me3$rnaLFC) * 0.1)
  

ggplot(data_for_plot_rna_seq_H3K4me3, aes(x=log2FoldChange, y=log2FC))+
  geom_point(size=2, aes(color = as.factor(annotation))) +
  theme_bw() +
  geom_smooth(method='lm', alpha=0.3) +
  xlab('RNA_LFC') +
  ylab('H3K4me1_LFC')

#-------RNA-seq vs H3K4me1
rnaseq_H3K4me1_intersection <- intersect(dataset_rnaseq$gene_symbol, H3K4me1_data$GENENAME)

data_for_plot_rna_seq_H3K4me1 <- 
  merge(dataset_rnaseq[(dataset_rnaseq$gene_symbol %in% rnaseq_H3K4me1_intersection), ], 
        H3K4me1_data[(H3K4me1_data$GENENAME %in% rnaseq_H3K4me1_intersection), ], 
        by.x = 'gene_symbol', by.y = 'GENENAME')

ggplot(data_for_plot_rna_seq_H3K4me1, aes(x=as.factor(gene_symbol), y=log2FC, fill = log2FoldChange))+
  geom_violin(width=2.5) +
  theme_bw() +
  scale_fill_gradientn(colours = c('#fde725', '#f6b092', '#5ec962', '#21918c', '#3b528b')) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', linewidth=1)

ggplot(data_for_plot_rna_seq_H3K4me1, aes(x=as.factor(gene_symbol), y=log2FC, fill = log2FoldChange))+
  geom_boxplot() +
  theme_bw() +
  scale_fill_gradientn(colours = c('#fde725', '#f6b092', '#5ec962', '#21918c', '#3b528b')) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', linewidth=1)

data_for_scatterplot_rna_seq_H3K4me1 <- data_for_plot_rna_seq_H3K4me1 %>%
  group_by(gene_symbol) %>%
  summarise(rnaLFC = median(log2FoldChange), H3K4me1 = median(log2FC))

ggplot(data_for_scatterplot_rna_seq_H3K4me1, aes(x=rnaLFC, y=H3K4me1))+
  geom_point(size=2, color = 'darkgreen') +
  theme_bw() +
  geom_smooth(method='lm', alpha = 0.3) +
  xlab('RNA_LFC') +
  ylab('H3K4me1_LFC') +
  stat_cor(method="spearman", cor.coef.name = c("rho"), 
           label.y = max(data_for_scatterplot_rna_seq_H3K4me1$H3K4me1) * 1.5, 
           label.x = max(data_for_scatterplot_rna_seq_H3K4me1$rnaLFC) * 0.1)

ggplot(data_for_plot_rna_seq_H3K4me1, aes(x=log2FoldChange, y=log2FC))+
  geom_point(size=2, aes(color = as.factor(annotation))) +
  theme_bw() +
  geom_smooth(method='lm', alpha=0.3)

         
#-------RNA-seq vs H3K36me3
rnaseq_H3K36me3_intersection <- intersect(dataset_rnaseq$gene_symbol, H3K36me3_data$GENENAME)

data_for_plot_rna_seq_H3K36me3 <- 
  merge(dataset_rnaseq[(dataset_rnaseq$gene_symbol %in% rnaseq_H3K36me3_intersection), ], 
        H3K36me3_data[(H3K36me3_data$GENENAME %in% rnaseq_H3K36me3_intersection), ], 
        by.x = 'gene_symbol', by.y = 'GENENAME')

ggplot(data_for_plot_rna_seq_H3K36me3, aes(x=as.factor(gene_symbol), y=log2FC, fill = log2FoldChange))+
  geom_violin(width=2.5) +
  theme_bw() +
  scale_fill_gradientn(colours = c('#fde725', '#f6b092', '#5ec962', '#21918c', '#3b528b')) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', linewidth=1)

ggplot(data_for_plot_rna_seq_H3K36me3, aes(x=as.factor(gene_symbol), y=log2FC, fill = log2FoldChange))+
  geom_boxplot() +
  theme_bw() +
  scale_fill_gradientn(colours = c('#fde725', '#f6b092', '#5ec962', '#21918c', '#3b528b')) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', linewidth=1)

data_for_scatterplot_rna_seq_H3K36me3 <- data_for_plot_rna_seq_H3K36me3 %>%
  group_by(gene_symbol) %>%
  summarise(rnaLFC = median(log2FoldChange), H3K36me3 = median(log2FC))

ggplot(data_for_scatterplot_rna_seq_H3K36me3, aes(x=rnaLFC, y=H3K36me3))+
  geom_point(size=2, color = 'darkgreen') +
  theme_bw() +
  geom_smooth(method='lm') +
  xlab('RNA_LFC') +
  ylab('H3K36me3_LFC') +
  stat_cor(method="spearman", cor.coef.name = c("rho"), 
           label.y = max(data_for_scatterplot_rna_seq_H3K36me3$H3K36me3) * 1.5, 
           label.x = max(data_for_scatterplot_rna_seq_H3K36me3$rnaLFC) * 0.1)

ggplot(data_for_plot_rna_seq_H3K36me3, aes(x=log2FoldChange, y=log2FC))+
  geom_point(size=2, aes(color = as.factor(annotation))) +
  theme_bw() +
  geom_smooth(method='lm', alpha=0.3)


#-------RNA-seq vs H3K27me3
rnaseq_H3K27me3_intersection <- intersect(dataset_rnaseq$gene_symbol, H3K27me3_data$GENENAME)

data_for_plot_rna_seq_H3K27me3 <- 
  merge(dataset_rnaseq[(dataset_rnaseq$gene_symbol %in% rnaseq_H3K27me3_intersection), ], 
        H3K27me3_data[(H3K27me3_data$GENENAME %in% rnaseq_H3K27me3_intersection), ], 
        by.x = 'gene_symbol', by.y = 'GENENAME')

ggplot(data_for_plot_rna_seq_H3K27me3, aes(x=as.factor(gene_symbol), y=log2FC, fill = log2FoldChange))+
  geom_violin(width=2.5) +
  theme_bw() +
  scale_fill_gradientn(colours = c('#fde725', '#f6b092', '#5ec962', '#21918c', '#3b528b')) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', linewidth=1)

ggplot(data_for_plot_rna_seq_H3K27me3, aes(x=as.factor(gene_symbol), y=log2FC, fill = log2FoldChange))+
  geom_boxplot() +
  theme_bw() +
  scale_fill_gradientn(colours = c('#fde725', '#f6b092', '#5ec962', '#21918c', '#3b528b')) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', linewidth=1)

data_for_scatterplot_rna_seq_H3K27me3 <- data_for_plot_rna_seq_H3K27me3 %>%
  group_by(gene_symbol) %>%
  summarise(rnaLFC = median(log2FoldChange), H3K36me3 = median(log2FC))

ggplot(data_for_scatterplot_rna_seq_H3K27me3, aes(x=rnaLFC, y=H3K27me3))+
  geom_point(size=2, color = '#f6b092') +
  theme_bw() +
  geom_smooth(method='lm') +
  stat_cor(method="spearman", cor.coef.name = c("rho"), 
           label.y = max(data_for_scatterplot_rna_seq_H3K27me3$H3K27me3) * 1.5, 
           label.x = max(data_for_scatterplot_rna_seq_H3K27me3$rnaLFC) * 0.1)

ggplot(data_for_plot_rna_seq_H3K27me3, aes(x=log2FoldChange, y=log2FC))+
  geom_point(size=2, aes(color = as.factor(annotation))) +
  theme_bw() +
  geom_smooth(method='lm', alpha=0.3)


#-------RNA-seq vs H3K27ac
rnaseq_H3K27ac_intersection <- intersect(dataset_rnaseq$gene_symbol, H3K27ac_data$GENENAME)

data_for_plot_rna_seq_H3K27ac <- 
  merge(dataset_rnaseq[(dataset_rnaseq$gene_symbol %in% rnaseq_H3K27ac_intersection), ], 
        H3K27ac_data[(H3K27ac_data$GENENAME %in% rnaseq_H3K27ac_intersection), ], 
        by.x = 'gene_symbol', by.y = 'GENENAME')

ggplot(data_for_plot_rna_seq_H3K27ac, aes(x=as.factor(gene_symbol), y=log2FC, fill = log2FoldChange))+
  geom_violin(width=2.5) +
  theme_bw() +
  scale_fill_gradientn(colours = c('#fde725', '#f6b092', '#5ec962', '#21918c', '#3b528b')) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', linewidth=1)

ggplot(data_for_plot_rna_seq_H3K27ac, aes(x=as.factor(gene_symbol), y=log2FC, fill = log2FoldChange))+
  geom_boxplot() +
  theme_bw() +
  scale_fill_gradientn(colours = c('#fde725', '#f6b092', '#5ec962', '#21918c', '#3b528b')) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', linewidth=1)

data_for_scatterplot_rna_seq_H3K27ac <- data_for_plot_rna_seq_H3K27ac %>%
  group_by(gene_symbol) %>%
  summarise(rnaLFC = median(log2FoldChange), H3K36me3 = median(log2FC))

ggplot(data_for_scatterplot_rna_seq_H3K27ac, aes(x=rnaLFC, y=H3K27ac))+
  geom_point(size=2, color = '#f6b092') +
  theme_bw() +
  geom_smooth(method='lm')

ggplot(data_for_plot_rna_seq_H3K27ac, aes(x=log2FoldChange, y=log2FC))+
  geom_point(size=2, aes(color = as.factor(annotation))) +
  theme_bw() +
  geom_smooth(method='lm', alpha=0.3)


#-------Intersection plot H3K4me3
unique_diff_meth_genes <- 
  as.data.frame(ENSEMBLE[(ENSEMBLE$gene_id %in% unique(dataset_diff_meth_all$nearest_gene)),]$gene_symbol) %>%
  mutate(count = 1)
colnames(unique_diff_meth_genes) <- c('genes', 'METH')

unique_diff_rna_genes <- as.data.frame(unique(dataset_rnaseq$gene_symbol), col.names='rna') %>%
  mutate(count = 1)
colnames(unique_diff_rna_genes) <- c('genes', 'RNA')

unique_diff_H3K4me3_genes <- as.data.frame(unique(H3K4me3_data$GENENAME)) %>%
  mutate(count = 1)
colnames(unique_diff_H3K4me3_genes) <- c('genes', 'H3K4me3')

data_meth_rna_H3K4me3 <- merge(unique_diff_meth_genes, unique_diff_rna_genes, all=TRUE) %>%
  merge(unique_diff_H3K4me3_genes, all = TRUE) %>%
  replace(is.na(.), 0)

data_meth_rna_H3K4me3[2:4] = data_meth_rna_H3K4me3[2:4] == 1

intersection_plot_meth_rna_H3K4me3 <- upset(data_meth_rna_H3K4me3[2:4], colnames(data_meth_rna_H3K4me3[2:4]), name = "Intersection", width_ratio=0.1)

#-------Intersection plot H3K4me1
unique_diff_meth_genes <- 
  as.data.frame(ENSEMBLE[(ENSEMBLE$gene_id %in% unique(dataset_diff_meth_all$nearest_gene)),]$gene_symbol) %>%
  mutate(count = 1)
colnames(unique_diff_meth_genes) <- c('genes', 'METH')

unique_diff_rna_genes <- as.data.frame(unique(dataset_rnaseq$gene_symbol), col.names='rna') %>%
  mutate(count = 1)
colnames(unique_diff_rna_genes) <- c('genes', 'RNA')

unique_diff_H3K4me1_genes <- as.data.frame(unique(H3K4me1_data$GENENAME)) %>%
  mutate(count = 1)
colnames(unique_diff_H3K4me1_genes) <- c('genes', 'H3K4me1')

data_meth_rna_H3K4me1 <- merge(unique_diff_meth_genes, unique_diff_rna_genes, all=TRUE) %>%
  merge(unique_diff_H3K4me1_genes, all = TRUE) %>%
  replace(is.na(.), 0)

data_meth_rna_H3K4me1[2:4] = data_meth_rna_H3K4me1[2:4] == 1

intersection_plot_meth_rna_H3K4me1 <- upset(data_meth_rna_H3K4me1[2:4], colnames(data_meth_rna_H3K4me1[2:4]), name = "Intersection", width_ratio=0.1)

#-------Intersection plot H3K36me3
unique_diff_meth_genes <- 
  as.data.frame(ENSEMBLE[(ENSEMBLE$gene_id %in% unique(dataset_diff_meth_all$nearest_gene)),]$gene_symbol) %>%
  mutate(count = 1)
colnames(unique_diff_meth_genes) <- c('genes', 'METH')

unique_diff_rna_genes <- as.data.frame(unique(dataset_rnaseq$gene_symbol), col.names='rna') %>%
  mutate(count = 1)
colnames(unique_diff_rna_genes) <- c('genes', 'RNA')

unique_diff_H3K36me3_genes <- as.data.frame(unique(H3K36me3_data$GENENAME)) %>%
  mutate(count = 1)
colnames(unique_diff_H3K36me3_genes) <- c('genes', 'H3K36me3')

data_meth_rna_H3K36me3 <- merge(unique_diff_meth_genes, unique_diff_rna_genes, all=TRUE) %>%
  merge(unique_diff_H3K36me3_genes, all = TRUE) %>%
  replace(is.na(.), 0)

data_meth_rna_H3K36me3[2:4] = data_meth_rna_H3K36me3[2:4] == 1

intersection_plot_meth_rna_H3K36me3 <- upset(data_meth_rna_H3K36me3[2:4], colnames(data_meth_rna_H3K36me3[2:4]), name = "Intersection", width_ratio=0.1)


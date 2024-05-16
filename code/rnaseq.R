library(DESeq2)
library(DOSE)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)


#-----------------------------------------DESeq2_GLM_DGE_analysis
# read the data
data <- read.csv("../data/raw_counts.csv", sep = ';', row.names = 1)

#create metadata table
metadata <- colnames(data) %>% .[-c(1)] %>% strsplit(split = 'D') %>% as.data.frame() %>% t() %>% as.data.frame
rownames(metadata) <- colnames(data) %>% .[-c(1)]
colnames(metadata) <- c('Age', 'Sample')

metadata <- metadata %>%
  mutate(Age = ifelse(Age == 'Y', 'YD', 'OD'))

metadata$Age <- factor(metadata$Age, levels = c('YD', 'OD'))

# filter extremely low expressed genes - expression level in sum by all groups is lower than a threshold
data_filtered <- data[rowSums(data[,-c(1)]) > 10, ]

# create DESeq object
dds <- DESeqDataSetFromMatrix(countData = round(data_filtered[,-c(1)]), 
                              colData = metadata, 
                              design = ~ Age)

# estimate scaling factors
dds <- estimateSizeFactors(dds)

# visualize the difference between sizeFactors of samples
size_factors <- as.data.frame(dds$sizeFactor)
colnames(size_factors) <- c('sizeFactor')
ggplot(size_factors, aes(x = rownames(size_factors), y = sizeFactor, fill = as.factor(rownames(size_factors)))) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 1.0, linetype = 'dashed') +
  theme_bw() +
  theme(legend.position = 'None',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  xlab('Samples') +
  ylim(c(0.0, 1.4))

# obtain normalized counts
normalized_counts <- counts(dds, normalized=TRUE)


#-------------------------------------------------------------PCA

plotPCA.mystyle <- function (norm_mat = FALSE, object, ntop = 500, returnData = FALSE)
{
  font.size <- 18
  if (norm_mat == FALSE){
    rv <- rowVars(assay(object))
    r <- assay(object)
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
  } else {
    rv <- rowVars(object)
    r <- object
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
                                                       length(rv)))]
    pca <- prcomp(t(object[select, ]))
  }

  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  d1 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], 
                   Age = metadata$Age, 
                   name = colnames(object))
  
  ggplot(data = d1, aes_string(x = "PC1", y = "PC2")) +
    geom_point(aes_string(color = "Age", shape = "Age"), size = 6, alpha = 0.7) + 
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
    ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    theme_dose(font.size = font.size)+ #+ geom_label_repel(aes(label = colnames(data)), label.size = 0.1, box.padding = 0.2)
    theme(
      legend.key = element_rect(colour = NA, fill = NA), 
      legend.title= element_blank(), 
      legend.text=element_text(size=font.size-2),
      plot.title = element_text(size = 18, hjust = 0.5, face="bold")
    ) +
    ggtitle(glue::glue('{as.character(substitute(object))}')) +
    geom_text(aes(label = rownames(pca$x)), size=4)
} 

# calculate rld object for PCA
rld <- rlog(dds, blind = TRUE)

# PCA
plotPCA.mystyle(norm_mat = FALSE, object = rld)
plotPCA.mystyle(norm_mat = FALSE, object = dds)
plotPCA.mystyle(norm_mat = TRUE, object = normalized_counts)
# For both PCA we can observe 3 samples from YD group that are significantly far away from
# other samples. It might be interesting for further analysis

# Heatmap
rld_matrix <- assay(rld)
rld_cor_mat <- cor(rld_matrix)
pheatmap(rld_cor_mat, annotation = metadata[,c(1,2)])

# Show the design matrix. Check that design matrix is okay
model.matrix(design(dds), data = colData(dds))

# DESeq2 analysis
dds_analysis <- DESeq(dds)
plotDispEsts(dds_analysis)
results_before_shrinkage <- results(dds_analysis, contrast=c('Age', 'YD', 'OD'), alpha = 0.05)

results <- lfcShrink(dds_analysis, 
                     contrast=c('Age', 'YD', 'OD'), 
                     res=results_before_shrinkage, 
                     type='normal')

summary(results, alpha = 0.05)

# Set the cutoffs
padj.cutoff <- 0.05
lfc.cutoff <- 0.405

# Filter significant genes with high LFC
results_tb <- results %>% 
  data.frame() %>% 
  tibble::rownames_to_column(var="gene") %>% 
  as_tibble()

results_tb <- merge(results_tb, (data_filtered %>% tibble::rownames_to_column(., var = 'gene'))[, c(1,2)], by = 'gene')

results_significant_tb <- results_tb %>% 
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)

results_significant_tb_without_LFC <- results_tb %>%
  dplyr::filter(padj < padj.cutoff)


# Visualize volcano-plots
EnhancedVolcano(results_tb,
                lab = results_tb$gene_symbol, 
                x = 'log2FoldChange',
                y = 'padj', 
                title = 'KO vs Control, PBS',
                subtitle = NULL,
                pCutoff = padj.cutoff, 
                FCcutoff = lfc.cutoff,
                xlim = c(-2, 2),
                ylim = c(NA, 4))

EnhancedVolcano(results_tb,
                lab = results_tb$gene_symbol, 
                x = 'log2FoldChange',
                y = 'padj', 
                title = 'KO vs Control, PBS',
                subtitle = NULL,
                pCutoff = padj.cutoff, 
                FCcutoff = 0,
                xlim = c(-2, 2),
                ylim = c(NA, 4))

#-----------------------------------------DESeq2_MWtest_DGE_analysis

# Create DF of normalized counts
normalized_counts_df <- as.data.frame(normalized_counts)

# Perform Wilcoxon-Mann-Whitney and Student's T-test
for (k in 1:nrow(normalized_counts_df)){
  OD <- as.numeric(as.vector(select(normalized_counts_df, starts_with('OD'))[k,]))
  YD <- as.numeric(as.vector(select(normalized_counts_df, starts_with('YD'))[k,]))
  normalized_counts_df$p_value_ttest[k] <- t.test(OD, YD, paired = FALSE, alternative = "two.sided", exact = TRUE)$p.value
  normalized_counts_df$p_value_wilcoxon[k] <- wilcox.test(OD, YD, paired = FALSE, alternative = "two.sided", exact = TRUE)$p.value
}

# Perform multiple testing correction according to Benjamini-Hochberg method
normalized_counts_df$p_value_ttest_adj <- p.adjust(normalized_counts_df$p_value_ttest, method='BH', n = nrow(normalized_counts_df))
normalized_counts_df$p_value_wilcoxon_adj <- p.adjust(normalized_counts_df$p_value_wilcoxon, method='BH', n = nrow(normalized_counts_df))

# Calculate the fold change
normalized_counts_df$FC <- rowMeans(select(normalized_counts_df, starts_with('OD'))) / rowMeans(select(normalized_counts_df, starts_with('YD')))
normalized_counts_df$FC <- log(normalized_counts_df$FC)

simple_DEG_volcano_tb <- normalized_counts_df %>%
  tibble::rownames_to_column(var = 'gene') %>%
  select(!starts_with('OD') & !starts_with('YD'))

simple_DEG_volcano_tb <- merge(simple_DEG_volcano_tb, (data_filtered %>% tibble::rownames_to_column(., var = 'gene'))[, c(1,2)], by = 'gene')

#Unfortunately, the number of sample is not big enough to detect such small difference at statistically significant level



#---------------------------------GOenrichment
enrichmentGO <- enrichGO(gene = results_significant_tb$gene, 
                        universe = results_tb$gene, 
                        keyType = "ENSEMBL",
                        OrgDb = org.Hs.eg.db, 
                        ont = "BP",
                        pAdjustMethod = "BH", 
                        pvalueCutoff = 0.1)


cowplot::plot_grid(barplot(enrichmentGO))

enrichmentGO_no_LFC <- enrichGO(gene = results_significant_tb_without_LFC$gene, 
                         universe = results_tb$gene, 
                         keyType = "ENSEMBL",
                         OrgDb = org.Hs.eg.db, 
                         ont = "BP",
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 0.1)

cowplot::plot_grid(barplot(enrichmentGO_no_LFC))



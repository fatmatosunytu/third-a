library(DESeq2)

# Simulated expression data and sample metadata
count_data <- matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10)
col_data <- data.frame(
  condition = factor(rep(c("control", "treatment"), each = 5))
)
rownames(count_data) <- paste0("gene", 1:100)
rownames(col_data) <- colnames(count_data) <- paste0("sample", 1:10)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)

# Perform differential expression analysis
dds <- DESeq(dds)
results <- results(dds)

# Filter significant genes
significant_genes <- results[results$padj < 0.05 & abs(results$log2FoldChange) > 1, ]
head(significant_genes)


# Differential Expression Analysis in R

This document provides instructions for running the R script for alignment, gene counting, and differential expression analysis of _Arabidopsis thaliana_ RNA-seq data.

## Prerequisites

### Required Software

- R (version 4.0.0 or higher recommended)
- Bioconductor packages:
    - Rsubread
    - DESeq2
    - EnhancedVolcano
- Additional R packages:
    - dplyr
    - ggplot2

### Required Data Files

The script expects the following input files:

1. **Reference genome files**:
    
    - Reference genome in FASTA format (e.g., `pathForReference.fna`)
    - Genome annotation file in GTF format (e.g., `Arabidopsis_annotation.gtf`)
2. **RNA-seq data files**:
    
    - Treatment FASTQ files (e.g., `Treatment_file.fastq`)
    - Control FASTQ files (e.g., `Control_file.fastq.gz`)

## Installation of Required Packages

Before running the analysis, ensure all required packages are installed:

```R
# Install Bioconductor and required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Rsubread", "DESeq2", "EnhancedVolcano"))

# Install additional packages
install.packages(c("dplyr", "ggplot2"))
```

## Analysis Steps

### 1. Alignment to Reference Genome

```R
library(Rsubread)

# Build index for reference genome
buildindex(basename = "Arabidopsis_index", 
           reference = "pathForReference.fna")  # Update path to your reference FASTA

# Align treatment sample
align(index = "Arabidopsis_index", 
      readfile1 = "Treatment_file.fastq",       # Update path to your treatment FASTQ
      output_file = "Aligned_Treatment.bam", 
      nthreads = 4)

# Align control sample            
align(index = "Arabidopsis_index", 
      readfile1 = "Control_file.fastq.gz",      # Update path to your control FASTQ
      output_file = "Aligned_Control.bam", 
      nthreads = 4)
```

### 2. Generate Count Matrix

```R
# Define BAM files
bam_files <- c("Aligned_Treatment.bam", 
               "Aligned_Control.bam")            # Make sure these match your output BAM files

# Run featureCounts
counts <- featureCounts(files = bam_files,
                        annot.ext = "Arabidopsis_annotation.gtf",  # Update path to your GTF file
                        isGTFAnnotationFile = TRUE,
                        GTF.featureType = "exon",
                        GTF.attrType = "gene_id",
                        useMetaFeatures = TRUE,
                        allowMultiOverlap = TRUE,
                        nthreads = 4)

# Extract the count matrix
count_matrix <- counts$counts

# Save the count matrix
write.csv(count_matrix, "count_matrix.csv")
```

### 3. Differential Expression Analysis with DESeq2

```R
library(DESeq2)

# Create sample information
sample_info <- data.frame(
  condition = factor(c("treatment", "control")),
  row.names = colnames(count_matrix)
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ condition)

# Run differential expression analysis
dds <- DESeq(dds)

# Obtain results
res <- results(dds)

# Filter significant DEGs (adjusted p-value < 0.05)
res_sig <- res[which(res$padj < 0.05), ]

# Export results
write.csv(as.data.frame(res_sig), "DEGs.csv")
```

### 4. Create Visualization

#### Volcano Plot

```R
library(EnhancedVolcano)

# Create volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                xlab = 'Log2 Fold Change (Effect Size)',
                ylab = '-Log10 Adjusted P-value (Significance)',
                title = 'Volcano Plot of Differential Gene Expression',
                subtitle = 'Highlights of RNA-seq Analysis',
                pointSize = 2.0,
                labSize = 2.5,
                col = c('gray', 'orange', 'skyblue', 'purple'),
                selectLab = rownames(res)[1:10],  # Top 10 genes
                drawConnectors = TRUE,
                widthConnectors = 0.4,
                boxedLabels = TRUE,
                max.overlaps = Inf)
```

#### Top Differentially Expressed Genes Plot

```R
library(dplyr)
library(ggplot2)

# Convert DESeqResults to a data frame
significant_genes_df <- as.data.frame(res_sig)
significant_genes_df$Gene <- rownames(significant_genes_df)

# Add regulation status
significant_genes_df$Regulation <- ifelse(significant_genes_df$log2FoldChange > 0, 
                                         "Upregulated", "Downregulated")

# Get top 10 upregulated and downregulated genes by baseMean
top_genes <- significant_genes_df %>% 
  group_by(Regulation) %>% 
  slice_max(baseMean, n = 10)

# Create the plot
ggplot(top_genes, aes(x = reorder(Gene, -baseMean), y = baseMean, fill = Regulation)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Top 20 Upregulated & Downregulated Genes",
    x = "Genes",
    y = "BaseMean Expression") +
  scale_fill_manual(values = c("Upregulated" = "skyblue", "Downregulated" = "orange")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
```

## Notes

- This pipeline may require significant memory for the alignment step, particularly for large genomes.
- The analysis provided is for a simple two-condition comparison (treatment vs. control).
- Since this analysis uses samples without biological replicates, results should be interpreted with caution and considered exploratory.
- The script assumes that the treatment condition is the first column in the count matrix and control is the second.

## Customization

### For Multiple Treatments or Replicates

If you have multiple treatment groups or biological replicates, modify the `sample_info` data frame:

```R
sample_info <- data.frame(
  condition = factor(c("treatment1", "treatment1", "treatment2", "treatment2", "control", "control")),
  replicate = factor(c("rep1", "rep2", "rep1", "rep2", "rep1", "rep2")),
  row.names = colnames(count_matrix)
)

# Update the design formula to include replicates if needed
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ condition)
```

### For Different Contrasts

To compare specific conditions:

```R
# For example, to compare treatment1 vs control
res <- results(dds, contrast=c("condition", "treatment1", "control"))
```

## Troubleshooting

- If alignment fails, check that your reference genome path is correct and the file is not corrupted.
- Memory issues during alignment or counting can be addressed by reducing the number of threads.
- If DESeq2 produces errors about low counts, consider pre-filtering low-count genes:
    
    ```R
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    ```
    
- For issues with the GTF file, ensure it's properly formatted and contains the feature types and attribute types specified in the featureCounts command.
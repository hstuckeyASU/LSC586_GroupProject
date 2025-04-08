# LSC586_GroupProject
Final group project git repository for code used in differential expression analysis of Arabidopsis

# Arabidopsis thaliana Differential Expression Analysis

This repository contains the computational pipeline for performing differential expression analysis on _Arabidopsis thaliana_ (Thale Cress) RNA-seq data. The implementation offers two distinct approaches for sequence alignment and gene quantification:

1. **Python-based pipeline** - Utilizes STAR for alignment, optimized for performance with large datasets
2. **R-based workflow** - Alternative methodology for comparison and validation purposes

Both pipelines are fully documented with clear instructions for execution. The downstream differential expression analysis is implemented exclusively in R using the DESeq2 package.

## Key Features

- Dual implementation allowing flexibility based on dataset size and computational resources
- Optimized STAR alignment method (Python) for handling large-scale RNA-seq data
- Complete workflow from raw reads to differential expression results
- Comprehensive documentation and usage examples

## Important Disclaimer

This analysis is performed on a dataset without biological replicates. As such:

- Statistical robustness cannot be guaranteed
- Results should be interpreted as exploratory only
- No definitive biological inferences should be drawn regarding _A. thaliana_ gene expression
- The pipeline is provided primarily for methodological demonstration and educational purposes

## Getting Started

Refer to the instruction documentation in either the `Python/` or `R/` directories for installation requirements and execution instructions.

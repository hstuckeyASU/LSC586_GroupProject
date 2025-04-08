# RNA-Seq Analysis Pipeline for Arabidopsis thaliana

This document provides instructions for running the `LSC586_Align_and_Gene_Counts.ipynb` Jupyter notebook, which performs RNA-seq alignment and gene counting for _Arabidopsis thaliana_ samples.

## Prerequisites

### Required Software

- Python 3.x with the following packages:
    - gzip
    - shutil
    - os
    - subprocess
    - pandas (for the optional CSV export)
- STAR aligner (automatically installed by the notebook)
- Samtools (installed via apt-get in the notebook)
- Subread package (for featureCounts, installed via apt-get in the notebook)

### Required Data Files

The notebook expects the following input files:

1. **RNA-seq data files (gzipped FASTQ format)**:
    
    - Arabidopsis Lead Treatment.gz
    - Arabidopsis Untreated Control.gz
    - Arabidopsis Vanadium Treatment.gz
2. **Reference genome files**:
    
    - GCF_000001735.4_TAIR10.1_genomic.fna (FASTA format)
    - GCF_000001735.4_TAIR10.1_genomic.gtf.gz (GTF format, gzipped)

> It should be noted that in theory you can replace the RNA-seq data files with any of your chosing along with their appropriate genome files
## Google Colab Setup

This notebook is designed to run in Google Colab with data stored in Google Drive. If you're using Colab:

1. Upload all required data files to your Google Drive
2. Modify the file paths in the first cell to match your Google Drive directory structure:
    
    ```python
    Lead_path = '/content/drive/MyDrive/YOUR_DIRECTORY/Arabidopsis Lead Treatment.gz'Control_Path = '/content/drive/MyDrive/YOUR_DIRECTORY/Arabidopsis Untreated Control.gz'Vanadium_Path = '/content/drive/MyDrive/YOUR_DIRECTORY/Arabidopsis Vanadium Treatment.gz'fasta_path = '/content/drive/MyDrive/YOUR_DIRECTORY/GCF_000001735.4_TAIR10.1_genomic.fna'gtf_gz_path = '/content/drive/MyDrive/YOUR_DIRECTORY/GCF_000001735.4_TAIR10.1_genomic.gtf.gz'
    ```
    

## Local Environment Setup

To run this notebook locally instead of in Google Colab:

1. Install Jupyter Notebook or JupyterLab
2. Install all required Python packages
3. Ensure STAR, Samtools, and Subread are installed on your system
4. Modify the file paths to point to your local file locations
5. Remove or modify the Google Drive mounting code

## Pipeline Steps

The notebook follows these steps:

1. **Data Preparation**:
    
    - Mount Google Drive
    - Decompress gzipped RNA-seq files to FASTQ format
2. **STAR Installation and Setup**:
    
    - Download and compile STAR aligner
    - Install Samtools and Subread packages
3. **Genome Indexing**:
    
    - Decompress GTF annotation file
    - Build STAR index using the reference genome and annotation
4. **Read Alignment**:
    
    - Align each RNA-seq sample to the reference genome using STAR
    - This is the most time-consuming step (approximately 2 hours for the three samples)
5. **SAM to BAM Conversion**:
    
    - Convert SAM alignment files to BAM format
    - Sort BAM files for further processing
6. **Gene Expression Quantification**:
    
    - Use featureCounts to count reads mapped to genomic features
    - Generate a gene count matrix for differential expression analysis
    - Optionally export counts to CSV format

## Expected Output Files

After running the notebook, you should have:

1. STAR index files in the `star_index` directory
2. Alignment files in the `alignments_star` directory
3. Sorted BAM files in the `bam_dir` directory
4. Gene count matrix in `gene_counts.txt` and `gene_counts.csv`

## Notes

- This pipeline requires significant computational resources, especially for the STAR alignment step
- The total execution time is approximately 2-3 hours depending on your hardware
- Gene count data can be used for downstream differential expression analysis in R using DESeq2
- This analysis uses samples without biological replicates, so results should be interpreted as exploratory only

## Troubleshooting

- If file paths cause errors, ensure they match the actual locations of your files
- If you encounter memory issues during STAR indexing or alignment, try reducing the number of threads
- For any featureCounts errors, check the compatibility of your GTF annotation file

## Next Steps

After generating the gene count matrix, proceed to differential expression analysis:

1. Load the `gene_counts.csv` or `gene_counts.txt` file in R
2. Use DESeq2 to identify differentially expressed genes between treatments
3. Visualize and interpret the results
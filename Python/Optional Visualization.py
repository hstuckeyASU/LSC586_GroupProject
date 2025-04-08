def visualize_results(lead_results, vanadium_results, count_matrix):

    # Set up plotting style
    sns.set(style="whitegrid")
    plt.rcParams.update({'font.size': 12})
    
    # 1. MA Plot (Log2 Fold Change vs Mean Expression)
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.scatter(
        np.log2(lead_results['mean_expression']),
        lead_results['Lead_vs_Control_log2FC'],
        alpha=0.5,
        s=10,
        c='gray'
    )
    plt.scatter(
        np.log2(lead_results[lead_results['Lead_significant']]['mean_expression']),
        lead_results[lead_results['Lead_significant']]['Lead_vs_Control_log2FC'],
        alpha=0.7,
        s=20,
        c='red'
    )
    plt.axhline(y=0, color='blue', linestyle='-', alpha=0.3)
    plt.axhline(y=1, color='red', linestyle='--', alpha=0.3)
    plt.axhline(y=-1, color='red', linestyle='--', alpha=0.3)
    plt.xlabel('Log2 Mean Expression')
    plt.ylabel('Log2 Fold Change')
    plt.title('Lead vs Control')
    
    plt.subplot(1, 2, 2)
    plt.scatter(
        np.log2(vanadium_results['mean_expression']),
        vanadium_results['Vanadium_vs_Control_log2FC'],
        alpha=0.5,
        s=10,
        c='gray'
    )
    plt.scatter(
        np.log2(vanadium_results[vanadium_results['Vanadium_significant']]['mean_expression']),
        vanadium_results[vanadium_results['Vanadium_significant']]['Vanadium_vs_Control_log2FC'],
        alpha=0.7,
        s=20,
        c='red'
    )
    plt.axhline(y=0, color='blue', linestyle='-', alpha=0.3)
    plt.axhline(y=1, color='red', linestyle='--', alpha=0.3)
    plt.axhline(y=-1, color='red', linestyle='--', alpha=0.3)
    plt.xlabel('Log2 Mean Expression')
    plt.ylabel('Log2 Fold Change')
    plt.title('Vanadium vs Control')
    
    plt.tight_layout()
    plt.savefig("MA_plot.png", dpi=300)
    plt.close()
    
    # 2. Heatmap of top differentially expressed genes
    # Combine significant genes from both comparisons
    significant_genes = set(
        lead_results[lead_results['Lead_significant']].index
    ).union(set(
        vanadium_results[vanadium_results['Vanadium_significant']].index
    ))
    
    # Limit to top 50 genes by absolute fold change
    top_genes_lead = lead_results.loc[lead_results['Lead_significant']].nlargest(25, 'Lead_vs_Control_log2FC').index
    top_genes_van = vanadium_results.loc[vanadium_results['Vanadium_significant']].nlargest(25, 'Vanadium_vs_Control_log2FC').index
    top_genes = list(set(top_genes_lead).union(set(top_genes_van)))[:50]
    
    if len(top_genes) > 0:
        # Extract data for heatmap
        heatmap_data = count_matrix.loc[top_genes, sample_names]
        
        # Add gene names
        gene_labels = [f"{idx} ({count_matrix.loc[idx, 'gene_name']})" for idx in heatmap_data.index]
        
        # Log transform for better visualization
        log_data = np.log2(heatmap_data + 1)
        
        # Plot heatmap
        plt.figure(figsize=(10, 14))
        sns.heatmap(
            log_data,
            cmap='viridis',
            yticklabels=gene_labels,
            xticklabels=sample_names,
            cbar_kws={'label': 'Log2(Count + 1)'}
        )
        plt.title('Top Differentially Expressed Genes')
        plt.tight_layout()
        plt.savefig("top_genes_heatmap.png", dpi=300)
        plt.close()
    
   # 3. Expression comparison plot
    plt.figure(figsize=(10, 8))
    plt.scatter(
        np.log2(count_matrix['Control'] + 1),
        np.log2(count_matrix['Lead'] + 1),
        alpha=0.5,
        s=10,
        c='blue',
        label='All Genes'
    )
    
    # Highlight significant genes
    significant_indices = lead_results[lead_results['Lead_significant']].index
    plt.scatter(
        np.log2(count_matrix.loc[significant_indices, 'Control'] + 1),
        np.log2(count_matrix.loc[significant_indices, 'Lead'] + 1),
        alpha=0.7,
        s=30,
        c='red',
        label='Significant Genes'
    )
    
    # Add diagonal line
    max_val = max(
        np.log2(count_matrix['Control'] + 1).max(),
        np.log2(count_matrix['Lead'] + 1).max()
    )
    plt.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
    
    plt.xlabel('Log2(Control Count + 1)')
    plt.ylabel('Log2(Lead Count + 1)')
    plt.title('Expression Comparison: Lead vs Control')
    plt.legend()
    plt.tight_layout()
    plt.savefig("lead_vs_control_expression.png", dpi=300)
    plt.close()
    
    # Create a similar plot for Vanadium
    plt.figure(figsize=(10, 8))
    plt.scatter(
        np.log2(count_matrix['Control'] + 1),
        np.log2(count_matrix['Vanadium'] + 1),
        alpha=0.5,
        s=10,
        c='green',
        label='All Genes'
    )
    
    # Highlight significant genes
    significant_indices = vanadium_results[vanadium_results['Vanadium_significant']].index
    plt.scatter(
        np.log2(count_matrix.loc[significant_indices, 'Control'] + 1),
        np.log2(count_matrix.loc[significant_indices, 'Vanadium'] + 1),
        alpha=0.7,
        s=30,
        c='red',
        label='Significant Genes'
    )
    
    # Add diagonal line
    max_val = max(
        np.log2(count_matrix['Control'] + 1).max(),
        np.log2(count_matrix['Vanadium'] + 1).max()
    )
    plt.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
    
    plt.xlabel('Log2(Control Count + 1)')
    plt.ylabel('Log2(Vanadium Count + 1)')
    plt.title('Expression Comparison: Vanadium vs Control')
    plt.legend()
    plt.tight_layout()
    plt.savefig("vanadium_vs_control_expression.png", dpi=300)
    plt.close()
    
    print("Visualizations saved as PNG files")

# Generate visualizations
visualize_results(lead_results, vanadium_results, count_matrix)
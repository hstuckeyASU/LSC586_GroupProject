### This is a basic non robust differential analysis method. It is recommended to perform this step in R ###

def perform_differential_analysis(count_matrix, sample_names):

    analysis_data = count_matrix.copy()
    
    # Add count to avoid division by zero
    for col in sample_names:
        analysis_data[col] = analysis_data[col] + 1
    
    # log2 fold changes
    analysis_data['Lead_vs_Control_log2FC'] = np.log2(
        analysis_data['Lead'] / analysis_data['Control']
    )
    analysis_data['Vanadium_vs_Control_log2FC'] = np.log2(
        analysis_data['Vanadium'] / analysis_data['Control']
    )

    analysis_data['Lead_significant'] = abs(analysis_data['Lead_vs_Control_log2FC']) > 1
    analysis_data['Vanadium_significant'] = abs(analysis_data['Vanadium_vs_Control_log2FC']) > 1

    analysis_data['mean_expression'] = analysis_data[sample_names].mean(axis=1)

    lead_results = analysis_data.sort_values(
        by=['Lead_significant', 'Lead_vs_Control_log2FC'], 
        ascending=[False, False]
    )
    
    vanadium_results = analysis_data.sort_values(
        by=['Vanadium_significant', 'Vanadium_vs_Control_log2FC'],
        ascending=[False, False]
    )
    
    # Save results
    lead_results.to_csv("lead_vs_control_results.csv")
    vanadium_results.to_csv("vanadium_vs_control_results.csv")
    
    print("Differential analysis complete")
    print(f"Lead vs Control: {lead_results['Lead_significant'].sum()} significant genes")
    print(f"Vanadium vs Control: {vanadium_results['Vanadium_significant'].sum()} significant genes")
    
    return lead_results, vanadium_results

lead_results, vanadium_results = perform_differential_analysis(count_matrix, sample_names)

The function and parameters in Slim-TPCA package
=============================================================

List of function:
-------------------------------
The list of function can be found using:

The function will be explained in the following text.

preproc(table, ref_col=1): 
    Calculates soluble fraction at each temperature. 
    
    Parameters: 
        table : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures
        ref_col : Reference column index for conversion to soluble fraction
    
    Return: 
        A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures

dist(table, ref_col=1, method='cityblock'): 
    Calculates distance between every two proteins. 
    
    Parameters: 
        table : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures
        ref_col : Reference column index for conversion to soluble fraction
        method : Distance calculation method
    
    Return: 
        A two-dimensional matrix table where the data in columns i,j represent the distance between protein i and protein j
    
pair_found(table, pair_table, ref_col=1): 
    Based on the protein pair interaction Database, look for protein pairs where both proteins appear in the data. 
    
    Parameters: 
        table : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures
        pair_table : A table with the first protein and the second protein of each protein pair in the first two columns
        ref_col : Reference column index for conversion to soluble fraction
    
    Return: 
        A DataFrame table containing protein pairs where both proteins appear in the data
    
roc(table, pair_table, ref_col=1, method='cityblock'): 
    Calculate parameters of the ROC curve. 
    
    Parameters: 
        table : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures
        pair_table : A table with the first protein and the second protein of each protein pair in the first two columns
        ref_col : Reference column index for conversion to soluble fraction
        method : Distance calculation method
    
    Return: 
        Three parameters of the ROC curve
    
roc_plot(table, pair_table, ref_col=1, method='cityblock'): 
    Draw ROC plot based on parameters.
    
    Parameters: 
        table : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures
        pair_table : A table with the first protein and the second protein of each protein pair in the first two columns
        ref_col : Reference column index for conversion to soluble fraction
        method : Distance calculation method
    
    Return: 
        Return None
    
complex_found(table, complex_table, ref_col=1): 
    Look for complexes that meet the requirements of the analysis.
    
    Parameters: 
        table : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures
        complex_table : A DataFrame table with complex-related information, where the subunits (UniProt IDs) column contains the protein IDs in the complex (intervals by;)
        ref_col : Reference column index for conversion to soluble fraction
    
    Return: 
        A table containing the complexes that meet the analysis requirements, with the discovered and undiscovered proteins in columns Subunit_Found and No_Subunit_Found, respectively
    
complex_dist(table, complex_table, ref_col=1, method='cityblock'): 
    Calculate average distance between the subunit proteins of the complex.
    
    Parameters: 
        table : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures
        complex_table : A DataFrame table with complex-related information, where the subunits (UniProt IDs) column contains the protein IDs in the complex (intervals by;)
        ref_col : Reference column index for conversion to soluble fraction
        method : Distance calculation method
    
    Return: 
        A table including the average distance and z-score values for each protein complex
    
random_n(table, complex_table, ref_col=1, method='cityblock', samplesize=10000): 
    Sample virtual random complexes for calculation.
    
    Parameters: 
        table : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures
        complex_table : A DataFrame table with complex-related information, where the subunits (UniProt IDs) column contains the protein IDs in the complex (intervals by;)
        ref_col : Reference column index for conversion to soluble fraction
        method : Distance calculation method
        samplesize : Number of random samples
    Return: 
        One list contain virtual random complexes with the same size are sampled
    
complex_signature_sample(table, complex_table, ref_col=1, method='cityblock', samplesize=10000): 
    Calculate TPCA signatures of complexes by sampling.
    
    Parameters: 
        table : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures
        complex_table : A DataFrame table with complex-related information, where the subunits (UniProt IDs) column contains the protein IDs in the complex (intervals by;)
        ref_col : Reference column index for conversion to soluble fraction
        method : Distance calculation method
        samplesize : Number of random samples
    Return: 
        A table including p value and z-score values for each protein complex
    
complex_signature_beta(table, complex_table, ref_col=1, method='cityblock', samplesize=500): 
    Calculate TPCA signatures of complexes by fitting a beta distribution to random complexes.
    
    Parameters: 
        table : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures
        complex_table : A DataFrame table with complex-related information, where the subunits (UniProt IDs) column contains the protein IDs in the complex (intervals by;)
        ref_col : Reference column index for conversion to soluble fraction
        method : Distance calculation method
        samplesize : Number of random samples
    Return: 
        A table including p value and z-score values for each protein complex
    
align(table_1, table_2, ref_col=1): 
    Multiple sets of data may identify different proteins and align them here.
    
    Parameters: 
        table_1 : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures in one status
        table_2 : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures in other status
        ref_col : Reference column index for conversion to soluble fraction
    Return: 
        Two table for the protein ids of table1 and table2 after alignment
    
dynamic_complex_absolute_sample(table_1, table_2, complex_table, ref_col=1, method='cityblock', samplesize=10000): 
    Calculate TPCA dynamic modulation signatures of complexes by sampling and absolute distance.
    
    Parameters: 
        table_1 : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures in one status
        table_2 : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures in other status
        complex_table : A DataFrame table with complex-related information, where the subunits (UniProt IDs) column contains the protein IDs in the complex (intervals by;)
        ref_col : Reference column index for conversion to soluble fraction
        method : Distance calculation method
        samplesize : Number of random samples
    Return: 
        A DataFrame table contain average Manhattan distance between melting curves among all pairs of subunits of a protein complex in table1(col: Avg_Dist_1) and table2(col: Avg_Dist_2), z-scores value in table1(col: Avg_Dist_Derived_1) and table2(col: Avg_Dist_2), (col: Avg_Dist_Derived_2), Avg Dist relative change the dynamic p-values of the protein complex changes(col: Dynamic_P)

    
dynamic_complex_relative_sample(table_1, table_2, complex_table, ref_col=1, method='cityblock', samplesize=10000): 
    Calculate TPCA dynamic modulation signatures of complexes by sampling and relative distance.
    
    Parameters: 
        table_1 : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures in one status
        table_2 : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures in other status
        complex_table : A DataFrame table with complex-related information, where the subunits (UniProt IDs) column contains the protein IDs in the complex (intervals by;)
        ref_col : Reference column index for conversion to soluble fraction
        method : Distance calculation method
        samplesize : Number of random samples
    Return: 
        A DataFrame table contain average Manhattan distance between melting curves among all pairs of subunits of a protein complex in table1(col: Avg_Dist_1) and table2(col: Avg_Dist_2), z-scores value in table1(col: Avg_Dist_Derived_1) and table2(col: Avg_Dist_2), (col: Avg_Dist_Derived_2), Avg Dist relative change the dynamic p-values of the protein complex changes(col: Dynamic_P)

    
dynamic_complex_absolute_beta(table_1, table_2, complex_table, ref_col=1, method='cityblock', samplesize=500): 
    Calculate TPCA dynamic modulation signatures of complexes by Beta distribution fitting and absolute distance.
    
    Parameters: 
        table_1 : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures in one status
        table_2 : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures in other status
        complex_table : A DataFrame table with complex-related information, where the subunits (UniProt IDs) column contains the protein IDs in the complex (intervals by;)
        ref_col : Reference column index for conversion to soluble fraction
        method : Distance calculation method
        samplesize : Number of random samples
    Return: 
        A DataFrame table contain average Manhattan distance between melting curves among all pairs of subunits of a protein complex in table1(col: Avg_Dist_1) and table2(col: Avg_Dist_2), z-scores value in table1(col: Avg_Dist_Derived_1) and table2(col: Avg_Dist_2), (col: Avg_Dist_Derived_2), Avg Dist relative change the dynamic p-values of the protein complex changes(col: Dynamic_P)

    
dynamic_complex_relative_beta(table_1, table_2, complex_table, ref_col=1, method='cityblock', samplesize=500): 
    Calculate TPCA dynamic modulation signatures of complexes by Beta distribution fitting and relative distance.
    
    Parameters: 
        table_1 : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures in one status
        table_2 : A DataFrame table containing the soluble concentrations (or intensities) of all proteins at different temperatures in other status
        complex_table : A DataFrame table with complex-related information, where the subunits (UniProt IDs) column contains the protein IDs in the complex (intervals by;)
        ref_col : Reference column index for conversion to soluble fraction
        method : Distance calculation method
        samplesize : Number of random samples
    Return: 
        A DataFrame table contain average Manhattan distance between melting curves among all pairs of subunits of a protein complex in table1(col: Avg_Dist_1) and table2(col: Avg_Dist_2), z-scores value in table1(col: Avg_Dist_Derived_1) and table2(col: Avg_Dist_2), (col: Avg_Dist_Derived_2), Avg Dist relative change the dynamic p-values of the protein complex changes(col: Dynamic_P)

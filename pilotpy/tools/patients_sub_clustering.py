import numpy as np
import pandas as pd
from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from .Trajectory import *
from ..plot.ploting import *



def check_raw(data=None):
    '''
    Checks if input array contains raw gene counts or normalised counts.

    Parameters
    ----------
    data : np.ndarray or None, optional
        Array containing gene counts, with shape n_cells x n_genes.


    Returns
    -------
    is_int : bool
        If True, all values in array are integers.
    is_non_negative : bool
        If True, no negative values were found in the array.
    '''
    is_int = np.all(np.equal(np.floor(data), data))
    is_non_negative = np.all(data >= 0)
    return is_int, is_non_negative


def filter_genes(raw_counts, gene_names, threshold):
    '''
    Filters genes based on their expression frequency across cells.

    Removes genes that are not detected in a minimum fraction of the total 
    cell population. This helps reduce noise and improves the power of 
    downstream differential expression analysis.

    Parameters
    ----------
    raw_counts : np.ndarray
        A 2D array of gene expression counts with shape n_cells x n_genes.
    gene_names : pd.Index or np.ndarray
        An array or Index containing the names of the genes corresponding 
        to the columns of raw_counts.
    threshold : float
        The minimum fraction of cells (between 0 and 1) in which a gene 
        must be detected (count > 0) to be retained.

    Returns
    -------
    raw_counts : np.ndarray
        The filtered count matrix containing only the genes that passed 
        the threshold.
    gene_names : pd.Index or np.ndarray
        The names of the retained genes, matching the new column 
        dimensions of raw_counts.
    '''
    min_cells = int(threshold * raw_counts.shape[0])
    detected = np.count_nonzero(raw_counts > 0, axis=0)
    keep_genes = detected >= min_cells
    
    print("Total genes found: ", gene_names.shape)
    print("Total cell count: ", raw_counts.shape[0])
    print(f"Genes must be found in at least {min_cells} cells.")  # Should be 10
    print(f"detected min/max/mean: {detected.min()}/{detected.max()}/{detected.mean():.1f}")
    print(f"keep_genes sum: {keep_genes.sum()}/{len(keep_genes)}")
    
    raw_counts = raw_counts[:, keep_genes]
    gene_names = gene_names[keep_genes]
    return raw_counts, gene_names


def compute_diff_expressions(adata, 
                             cell_type: str = None, 
                             proportions: pd.DataFrame = None,
                             selected_genes: list = None, 
                             font_size:int=18,
                             group1: str = 'Tumor 1',
                             group2: str = 'Tumor 2',
                             fc_thr=0.5, 
                             pval_thr=0.05, 
                             exp_thr=0.1, 
                             sample_col='sampleID',
                             col_cell='cell_types', 
                             shrinkage=False, 
                             number_n=5, 
                             number_p=5, 
                             marker='o', 
                             color='w', 
                             markersize=8, 
                             font_weight_legend='normal', 
                             size_legend=12,
                             figsize=(15,15), 
                             dpi=100, 
                             **kwargs):
    '''
    Computes pseudobulk differential expression analysis for a specific cell type 
    using PyDESeq2 and generates a volcano plot.

    Parameters
    ----------
    adata : AnnData
        Annotated data object containing single-cell gene expression and metadata.
    cell_type : str, optional
        The specific cell type (found in col_cell) to perform DE analysis on.
    proportions : pd.DataFrame, optional
        Metadata table containing sample IDs and condition labels (Predicted_Labels).
    selected_genes : list, optional
        Subset of gene names to restrict the analysis to. If None, all filtered genes are used.
    font_size : int, optional
        Base font size for the plot text and labels. Default is 18.
    group1 : str, optional
        The name of the first experimental group (e.g., 'Tumor 1').
    group2 : str, optional
        The name of the second experimental group (e.g., 'Tumor 2').
    fc_thr : float, optional
        Log2 Fold Change threshold for significance in the volcano plot. Default is 0.5.
    pval_thr : float, optional
        Adjusted p-value threshold for significance. Default is 0.05.
    exp_thr : float, optional
        Expression threshold; filters genes expressed in fewer than this fraction of cells. Default is 0.1.
    sample_col : str, optional
        The column in adata.obs containing sample identifiers. Default is 'sampleID'.
    col_cell : str, optional
        The column in adata.obs containing cell type annotations. Default is 'cell_types'.
    shrinkage : bool, optional
        If True, applies apeGLM log2 fold change shrinkage. Default is False.
    number_n : int, optional
        Number of top down-regulated genes to label in the volcano plot. Default is 5.
    number_p : int, optional
        Number of top up-regulated genes to label in the volcano plot. Default is 5.
    marker : str, optional
        Marker style for the scatter plot points. Default is 'o'.
    color : str, optional
        Edge color for the markers in the plot. Default is 'w'.
    markersize : int, optional
        Size of the markers used in the legend. Default is 8.
    font_weight_legend : str, optional
        Font weight for the legend text. Default is 'normal'.
    size_legend : int, optional
        Font size for the legend text. Default is 12.
    figsize : tuple, optional
        Width and height of the resulting figure in inches. Default is (15, 15).
    dpi : int, optional
        Resolution of the saved figure in dots per inch. Default is 100.
    **kwargs : dict
        Additional keyword arguments passed to internal functions.

    Returns
    -------
    None
        The function saves a CSV of differential expression statistics and a 
        PDF volcano plot to the 'Results_PILOT' directory.
    '''
        
    result_path = Path('Results_PILOT')

    
    # Get gene counts for cell type
    cell_mask = adata.obs[col_cell] == cell_type
    raw_counts = adata.raw[cell_mask].X if adata.raw is not None else adata[cell_mask].X
    gene_names = adata.raw.var_names if adata.raw is not None else adata.var_names
    sample_ids = adata.obs.loc[cell_mask, sample_col]

    # Check for raw counts, filter for genes expressed in at least 10% of cells
    count_matrix = raw_counts.toarray() if hasattr(raw_counts, 'toarray') else raw_counts
    is_int, is_non_negative = check_raw(count_matrix)
    if not is_int or not is_non_negative:
        raise TypeError("Gene counts seem to be normalised, requires raw gene counts instead.")   
    
    
    count_matrix, gene_names = filter_genes(count_matrix, gene_names, threshold=exp_thr)
    
    pseudobulk = (pd.DataFrame(count_matrix.T, index=gene_names, columns=sample_ids)
                .groupby(level=0, axis=1).sum().T)
    
    # Free up RAM
    del count_matrix, raw_counts 
    
    relevant_samples = proportions[proportions['Predicted_Labels'].isin([group1, group2])]['sampleID']
    pseudobulk = pseudobulk.loc[pseudobulk.index.intersection(relevant_samples)].fillna(0).astype(int)
    metadata = proportions.loc[pseudobulk.index, 'Predicted_Labels'].to_frame()
    
    # Filter for selected genes
    if selected_genes is not None:
        selected_genes = list(set(selected_genes) & set(pseudobulk.columns))
        pseudobulk = pseudobulk[selected_genes]
    
    # Perform dispersion and log fold-change (LFC) estimation
    dds = DeseqDataSet(counts=pseudobulk, metadata=metadata, design_factors='Predicted_Labels')
    dds.deseq2()

    # Statistical Analysis: p-value computation using Wald test
    stat = DeseqStats(dds, contrast=['Predicted_Labels', group2, group1])  # Normal vs Tumor 1
    stat.summary()
    # LFC shrinkage with apeGLM
    if shrinkage:
        stat.lfc_shrink(coeff=f'Predicted_Labels[T.{group2}]')
    
    # Save to csv
    result_path = result_path / 'Diff_Expressions_Results' / cell_type 
    result_path.mkdir(parents=True, exist_ok=True)
    stat.results_df.to_csv(result_path / "Diff_expressions_stats.csv")

    pv_thr = -np.log10(pval_thr)  

    # Volcano Plot for viualisation
    volcano_plot(
        scores=pseudobulk.T,                    # Series/array mode: index = gene symbols
        foldchanges=stat.results_df['log2FoldChange'],
        p_values=stat.results_df['padj'],
        cell_type=cell_type,
        feature1=group1,
        feature2=group2,
        fc_thr=fc_thr,
        pv_thr=pv_thr,
        label_mode="topN",                     
        n_p=number_p,
        n_n=number_n,
        figsize=figsize,
        output_path=result_path,
        font_size=font_size,
        marker=marker,
        marker_edge_color=color,               
        markersize_legend=markersize,
        font_weight_legend=font_weight_legend,
        size_legend=size_legend,
        dpi=dpi
    )
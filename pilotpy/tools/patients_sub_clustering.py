from adjustText import adjust_text
from gprofiler import GProfiler
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.stats import ttest_ind
import seaborn as sns
from statsmodels.stats.multitest import multipletests
from .Trajectory import *
from ..plot.ploting import *


    
def filter_genes(raw_counts, gene_names, threshold):
    """
    Filters out genes that are expressed in fewer cells than a specified threshold.

    Parameters
    ----------
    raw_counts : np.ndarray
        Raw gene expression count matrix of shape (cells, genes).
    gene_names : np.ndarray or pd.Index
        Array of gene names corresponding to columns in raw_counts.
    threshold : float
        Minimum fraction of cells in which a gene must be detected (non-zero) to be retained.

    Returns
    -------
    filtered_counts : np.ndarray
        Subset of raw_counts containing only genes meeting the expression threshold.
    filtered_gene_names : np.ndarray or pd.Index
        Subset of gene_names corresponding to filtered_counts.

    
    Prints summary statistics about gene detection before filtering.
    Ensures that only genes expressed in at least `threshold * total_cells` are kept.
    """
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
    
    """
    Performs differential gene expression analysis using DESeq2 methodology (via PyDESeq2).
    Counts are aggregated into pseudobulk per sample filtered for lowly expressed genes,
    and analyzed using a negative binomial model. Log2 fold changes and adjusted p-values
    are computed using Wald tests. Optionally applies LFC shrinkage.
    Note: This function requires RAW gene counts. s

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.Requires storage of raw, non-normalised gene counts either in adata.X or adata.raw.X. 
    cell_type : str
        Specify cell type name to check its differential expression genes. The default is None.
    proportions : pd.DataFrame
        Cell types proportions in each sample. The default is None.
    selected_genes : list, optional
        Specify gene names to be considered for checking their differentiation. 
        If provided, restricts the analysis to the intersection of these genes and the dataset.
    font_size : int, optional
        Font size for plot labels and legends. The default is 18.
    group1 : str, optional
        Name of the first patient sub-group for comparison. The default is 'Tumor 1'.
    group2 : str, optional
        Name of the second patient sub-group for comparison. The default is 'Tumor 2'.
    fc_thr : float, optional
        Specify the fold change threshold. The default is 0.5.
    pval_thr : float, optional
        Specify the adjusted p-value threshold. The default is 0.05.
    exp_thr : float, optional
        Minimum fraction of cells in which a gene must be expressed to be retained. The defeult is 0.1=10%.
    sample_col : str, optional
        Name of the column containing sample IDs. The default is 'sampleID'.
    col_cell : str, optional
        Name of the column containing cell type annotations. The default is 'cell_types'.
    shrinkage: bool, optional
        Defines if LFC shrinkage is to be applied. The default is False.
    number_n : int, optional
        The number of labels that the user wants to show over the plot for negative thresholds. The default is 5.
    number_p : int, optional
        The number of labels that the user wants to show over the plot for positive thresholds. The default is 5.
    marker : str, optional
        Marker style for the labels in the volcano plot. The default is 'o'.
    color : str, optional
        Marker color for the labels in the volcano plot. The default is 'w'.
    markersize : int, optional
        Marker size for the labels in the volcano plot. The default is 8.
    font_weight_legend : str, optional
        Font weight for legend labels. The default is 'normal'.
    size_legend : int, optional
        Font size for legend labels. The default is 12.
    figsize: tuple, optional
        Figure size. The default is (15,15).
    dpi : int, optional
        Dots per inch for the saved plot image. Default is 100.

    Returns
    -------
    None

    Generates and displays a volcano plot of fold changes between two interested patient sub-groups.
    Saves a statistical table of each gene.
    Saves significantly differentiated genes in each group.
    """
    
    path_to_result='Results_PILOT'

    def check_raw(data):
        is_int = np.all(np.equal(np.floor(data), data))
        is_non_negative = np.all(data >= 0)
        return is_int, is_non_negative
    
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
    stat = DeseqStats(dds, contrast=['Predicted_Labels', group2, group1])  
    stat.summary()
    # LFC shrinkage with apeGLM
    if shrinkage:
        stat.lfc_shrink(coeff=f'Predicted_Labels[T.{group2}]')

    # Save statistical table of each gene
    path_to_result = path_to_result + f'/Diff_Expressions_Results/{cell_type}'
    os.makedirs(path_to_result, exist_ok=True)
    stat.results_df.to_csv(path_to_result + "/diff_expressions_stats_" + cell_type + ".csv")

    # Save up- and downregulated genes to csv
    sig = stat.results_df.dropna()
    sig_filtered = sig[
        (abs(sig['log2FoldChange']) >= fc_thr) &
        (sig['padj'] <= pval_thr)
    ]

    up = sig_filtered[sig_filtered['log2FoldChange'] > 0]
    down = sig_filtered[sig_filtered['log2FoldChange'] < 0]

    up.to_csv(path_to_result + "upregulated_genes.csv")
    down.to_csv(path_to_result + "downregulated_genes.csv")
    

    # Volcano Plot for viualisation
    pv_thr = -np.log10(pval_thr)  
    volcano_plot(
        pseudobulk.T,  
        stat.results_df['log2FoldChange'], 
        stat.results_df['padj'],
        cell_type, group1, group2, fc_thr, pv_thr,  
        figsize=figsize, output_path=path_to_result,
        n_p=number_p, n_n=number_n, font_size=font_size,
        marker=marker, color=color, markersize=markersize,
        font_weight_legend=font_weight_legend, size_legend=size_legend, dpi=dpi
    )
    


    
def install_r_packages():
    """
    Install R packages using rpy2.

    This function installs the "limma" R package using the "BiocManager" package manager.

    Parameters:
        None

    Returns:
        None
    """
    # Install R packages using rpy2
    import rpy2.robjects as robjects

    robjects.r('''
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    ''')

    robjects.r('''
    BiocManager::install("limma")
    ''')
   

   

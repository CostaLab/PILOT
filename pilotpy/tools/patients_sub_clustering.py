import numpy as np
import pandas as pd
from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from .Trajectory import *
from ..plot.ploting import *

## CHANGES
# replaced os with Pathlib
# replaced R's limma with pydeseq2
# added docstrings
# extract_cells... and install_R_... functions have become obsolete
# removed unnecessary imports
# added functions check_raw and filter_genes instead



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

    

#--------------------------------------------------------------------------------------------------------------------
#                                R Versions/old versions
#--------------------------------------------------------------------------------------------------------------------

def extract_cells_from_gene_expression_for_clustering_R(adata,sample_col,col_cell,cell_list,path_results=None,normalization=True,n_top_genes=2000,highly_variable_genes_=False):
    
    
    """
    Extract and save gene expression data for specific cells for clustering analysis.

    Parameters:
        adata : AnnData object
            An Annotated Data (AnnData) object containing gene expression data.
        sample_col : str
            The column name in the adata.obs DataFrame containing sample IDs.
        col_cell : str
            The column name in the adata.obs DataFrame containing cell type labels.
        cell_list : list
            A list of cell types for which gene expression data will be extracted and saved.
        path_results : str or None, optional (default=None)
            The path to the directory where the extracted data will be saved as CSV files.
            If None, the default path 'Results_PILOT/cells/' will be used.
        normalization : bool, optional (default=True)
            Whether to normalize the gene expression data by total count and apply log1p transformation.
        n_top_genes : int, optional (default=2000)
            The number of top highly variable genes to select. Only applicable if highly_variable_genes_ is True.
        highly_variable_genes_ : bool, optional (default=False)
            Whether to select highly variable genes for analysis.

    Returns:
        Gene exp. dataframe
    """
    for cell in cell_list:
        adata_new = adata[adata.obs[col_cell].isin([cell]),:]
        if normalization:
            sc.pp.normalize_total(adata_new, target_sum=1e4)
            sc.pp.log1p(adata_new)
        
        if highly_variable_genes_:
            
            sc.pp.highly_variable_genes(adata_new, n_top_genes=n_top_genes)
                # Access the list of highly variable genes
            highly_variable_genes = adata_new.var['highly_variable']
            df=adata_new[:,highly_variable_genes].X
            df=pd.DataFrame(df.toarray())
            highly_variable_gene_names = adata_new.var_names[np.array(adata_new.var['highly_variable'])]
            df.columns=list(highly_variable_gene_names)
        else:
            
            df=adata_new[:,adata_new.var_names].X
            df=pd.DataFrame(df.toarray())
            df.columns=adata_new.var_names
            
    
        df['sampleID']=list(adata_new.obs[sample_col])
        
        if path_results==None:
            if not os.path.exists('Results_PILOT/cells/'):
                os.makedirs('Results_PILOT/cells/')
            path_results='Results_PILOT/cells/'
        else:
            path_results=path_results

        
        df.to_csv(path_results+cell+'.csv')
        return df


def compute_diff_expressions_R(adata,cell_type: str = None,
                             proportions: pd.DataFrame = None,
                             selected_genes: list = None,
                             font_size:int=18,
                             group1: str = 'Tumor 1',
                             group2: str = 'Tumor 2',
                             label_name: str = 'Predicted_Labels',
                             fc_thr: float = 0.5,
                             pval_thr: float = 0.01,
                             sample_col:str='sampleID',
                             col_cell:str ='cell_types',
                             path=None,
                             normalization=False,
                             n_top_genes=2000,
                             highly_variable_genes_=True,
                             number_n=5,
                             number_p=5,
                             marker='o',
                             color='w',
                             markersize=8,
                             font_weight_legend='normal',
                             size_legend=12,
                             figsize=(15,15),dpi=100
                             ):

    """
    Using limma R package, lmFit fits a linear model using weighted least squares for each gene.
    Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models.
    Empirical Bayes smoothing of standard errors (shrinks standard errors
    that are much larger or smaller than those from other genes towards the average standard error).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    cell_type : str, optional
        Specify cell type name to check its differential expression genes. The default is None.
    proportions : pd.DataFrame, optional
        Cell types proportions in each sample. The default is None.
    selected_genes : list, optional
        Specify gene names to be considered for checking their differentiation.
    font_size : int, optional
        Font size for plot labels and legends. The default is 18.
    group1 : str, optional
        Name of the first patient sub-group for comparison. The default is 'Tumor 1'.
    group2 : str, optional
        Name of the second patient sub-group for comparison. The default is 'Tumor 2'.
    label_name : str, optional
        Name of the column containing the labels of patient sub-groups. The default is 'Predicted_Labels'.
    fc_thr : float, optional
        Specify the fold change threshold. The default is 0.5.
    pval_thr : float, optional
        Specify the adjusted p-value threshold. The default is 0.01.
    sample_col : str, optional
        Name of the column containing sample IDs. The default is 'sampleID'.
    col_cell : str, optional
        Name of the column containing cell type annotations. The default is 'cell_types'.
    path : str, optional
        Path to save the results. The default is None.
    normalization : bool, optional
        Perform gene expression normalization. The default is False.
    n_top_genes : int, optional
        Number of top variable genes to consider. The default is 2000.
    highly_variable_genes_ : bool, optional
        Determine highly variable genes. The default is True.
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
    
    if os.path.exists(path_to_result + "/cells/" + cell_type + ".csv"):

        cells = pd.read_csv(path_to_result + "/cells/" + cell_type + ".csv", index_col = 0)  
   
    elif cell_type not in adata.uns:
        cells=extract_cells_from_gene_expression_for_clustering(adata,sample_col=sample_col,col_cell=col_cell,cell_list=[cell_type],path_results=path,normalization=normalization,n_top_genes=n_top_genes,highly_variable_genes_=highly_variable_genes_)
        #cells = pd.read_csv(path_to_result + "/cells/" + cell_type + ".csv", index_col = 0)
    
    elif cell_type in adata.uns :
         cells=adata.uns[cell_type] 
    
    
    
    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    pandas2ri.activate()
    
    # prepare data for R
    proportions.index = proportions['sampleID']

    if selected_genes is None:
        selected_genes = cells.iloc[:,1:-1].columns
    data = cells[selected_genes]
    pred_labels = pd.DataFrame()
    pls = proportions.loc[cells['sampleID']]
    pred_labels['Predicted_Labels'] = pls[label_name]
    pred_labels['sampleID'] = pls['sampleID']
    
    # load R packages and data
    R=robjects.r 
    R('library(limma)') 
    R.assign('data',data) 
    R.assign('pred_labels', pred_labels) 
    R.assign('selected_groups', [group1, group2]) 
    R('selected_pred_labels <- pred_labels[which(pred_labels$Predicted_Labels %in% selected_groups),]')
    R('subresult <- data[row.names(selected_pred_labels),]')

    # delete for memory
    del data
    del pred_labels
    
    # run limma
    print('run limma lmFit')
    
    R('fit <- limma::lmFit(t(subresult), design = unclass(as.factor(selected_pred_labels$Predicted_Labels)))')
    print('run limma eBayes')
  
    R('fit <-  limma::eBayes(fit)')
    R('res <- limma::topTable(fit, n = 2000)')
    R('res <- res[colnames(data), ]')
    
    # get results
    res = R('''res''')
    
    if not os.path.exists(path_to_result+'/Diff_Expressions_Results/'+cell_type):
            os.makedirs(path_to_result+'/Diff_Expressions_Results/'+cell_type)
    path_to_result=path_to_result+'/Diff_Expressions_Results/'+cell_type+'/'
    res.to_csv(path_to_result + "/diff_expressions_stats_" + cell_type + ".csv")
      
    pv_thr = -np.log10(pval_thr)
    volcano_plot(cells[selected_genes].transpose(), res['logFC'], res['adj.P.Val'],
                 cell_type, group1, group2, fc_thr, pv_thr,
                 figsize = figsize, output_path = path_to_result,n_p=number_p,n_n=number_n,font_size=font_size, marker=marker,
                             color=color,
                             markersize=markersize,
                             font_weight_legend=font_weight_legend,
                             size_legend=size_legend,dpi=dpi)


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
   
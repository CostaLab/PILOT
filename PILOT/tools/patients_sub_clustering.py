import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from adjustText import adjust_text
from matplotlib.lines import Line2D
from gprofiler import GProfiler
from .Trajectory import *
from ..plot.ploting import *


    
def extract_cells_from_gene_expression_for_clustering(adata,sample_col,col_cell,cell_list,path_results=None,normalization=True,n_top_genes=2000,highly_variable_genes_=False):
    
    
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
def compute_diff_expressions(adata,cell_type: str = None,
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
    proportions.index = proportions['sampIeD']
   
    if selected_genes is None:
        selected_genes = cells.iloc[:,1:-1].columns
    data = cells[selected_genes]
    pred_labels = pd.DataFrame()
    pls = proportions.loc[cells['sampleID']]
    pred_labels['Predicted_Labels'] = pls[label_name]
    pred_labels['sampleID'] = pls['sampIeD']
    
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
   

   

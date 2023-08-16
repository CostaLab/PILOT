#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 08:49:06 2023

@authors: Mina&Mehdi
"""

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




    
def plot_hor_vs_vert(data, subplot, x, y, c, xlabel, ylabel, rotation,
                     tick_bottom, tick_left, title,fontsize=24):
    ax=plt.subplot(1,2,subplot)
    cols = ['tab:red' if x >= 0 else 'tab:blue' for x in data[x]]

    sns.barplot(x=x, y=y, data=data, ci=None, palette=cols)
    plt.title(title, fontsize=fontsize, fontweight = 'bold')
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.xticks(fontsize=fontsize, rotation=rotation)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    sns.despine(bottom=False, left=True)
    ax.grid(False)
    ax.tick_params(bottom=tick_bottom, left=tick_left)
    return None



    
def map_color(a, fc_thrr, pv_thrr):
    log2FoldChange, symbol, nlog10 = a
    if log2FoldChange >= fc_thrr and nlog10 >= pv_thrr:
        return 'very higher'
    elif log2FoldChange <= -fc_thrr and nlog10 >= pv_thrr:
        return 'very lower'
    elif log2FoldChange >= fc_thrr and nlog10 < pv_thrr:
        return 'higher'
    elif log2FoldChange <= -fc_thrr and nlog10 < pv_thrr:
        return 'lower'
    elif abs(log2FoldChange) < fc_thrr and nlog10 >= pv_thrr:
        return 'mix'
    else:
        return 'no'

def volcano_plot(scores, foldchanges, p_values, cell_type, feature1, feature2, fc_thr = 1, pv_thr = 1,
                 figsize = (20,20), output_path = None):
    df = pd.DataFrame(columns=['log2FoldChange', 'nlog10', 'symbol'])
    df['log2FoldChange'] = foldchanges
    df['nlog10'] = -np.log10(p_values.values)
    df['symbol'] = scores.index.values
    
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.dropna(subset=["nlog10"], how="all", inplace=True)
    

    selected_labels = df.loc[ (np.abs(df.log2FoldChange) >= fc_thr) & (df['nlog10'] >= pv_thr)]['symbol'].values
    group1_selected_labels = df.loc[ (df.log2FoldChange <= -fc_thr) & (df['nlog10'] >= pv_thr)]['symbol'].values
    pd.DataFrame(group1_selected_labels).to_csv(output_path + "/significant_genes_" + str(cell_type) + "_" + str(feature1) + ".csv")
    
    group2_selected_labels = df.loc[ (df.log2FoldChange >= fc_thr) & (df['nlog10'] >= pv_thr)]['symbol'].values
    pd.DataFrame(group2_selected_labels).to_csv(output_path + "/significant_genes_" + str(cell_type) + "_" + str(feature2) + ".csv")
    
    def map_shape(symbol):
        if symbol in selected_labels:
            return 'important'
        return 'not'
    
    df['color'] = df[['log2FoldChange', 'symbol', 'nlog10']].apply(map_color, fc_thrr = fc_thr, pv_thrr = pv_thr, axis = 1)
    df['shape'] = df.symbol.map(map_shape)
    df['baseMean'] = df.nlog10*10

    
    plt.figure(figsize = figsize, frameon=False, dpi=100)
    plt.style.use('default')

    ax = sns.scatterplot(data = df, x = 'log2FoldChange', y = 'nlog10', 
                         hue = 'color', hue_order = ['no', 'very higher','higher', 'mix', 'very lower', 'lower'],
                         palette = ['lightgrey', '#d62a2b', '#D62A2B7A',
                                    'lightgrey', '#1f77b4', '#1F77B47D'],
                         style = 'shape', style_order = ['not', 'important'],
                         markers = ['o', 'o'], 
                         size = 'baseMean', sizes = (40, 800)
                        )

    ax.axhline(pv_thr, zorder = 0, c = 'k', lw = 2, ls = '--')
    ax.axvline(fc_thr, zorder = 0, c = 'k', lw = 2, ls = '--')
    ax.axvline(-fc_thr, zorder = 0, c = 'k', lw = 2, ls = '--')

    texts = []
    for i in range(len(df)):
        if df.iloc[i].nlog10 >= pv_thr and (df.iloc[i].log2FoldChange >= fc_thr):
            texts.append(plt.text(x = df.iloc[i].log2FoldChange, y = df.iloc[i].nlog10, s = df.iloc[i].symbol,
                                 fontsize = 16, weight = 'bold', family = 'sans-serif'))
        if df.iloc[i].nlog10 >= pv_thr and ( df.iloc[i].log2FoldChange <= -fc_thr):
            texts.append(plt.text(x = df.iloc[i].log2FoldChange, y = df.iloc[i].nlog10, s = df.iloc[i].symbol,
                                 fontsize = 16, weight = 'bold', family = 'sans-serif'))
    adjust_text(texts)

    custom_lines = [Line2D([0], [0], marker='o', color='w', markerfacecolor='#d62a2b', markersize=8),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4', markersize=8)]

    plt.legend(custom_lines, ['Higher expressions in ' + feature2, 'Higher expressions in ' + feature1],loc = 1,
               bbox_to_anchor = (1,1.1), frameon = False, prop = {'weight': 'normal', 'size': 12})

    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(2)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.tick_params(width = 2)

    plt.title("Expression Score \n "+feature1+" - "+feature2, fontsize = 24)
    plt.xticks(size = 18, weight = 'bold')
    plt.yticks(size = 18, weight = 'bold')
    plt.xlabel("$log_{2}$ (Fold Change)", size = 18)
    plt.ylabel("-$log_{10}$ (P-value)", size = 18)

#     plt.savefig(filename, dpi = 100, bbox_inches = 'tight', facecolor = 'white')
    plt.savefig(output_path + "/volcano_" + str(feature1) + "-" + str(feature2) + "_FC.pdf",
                dpi = 100, bbox_inches = 'tight', facecolor = 'white')
    plt.show()

 


def extract_cells_from_gene_expression_for_clustering(adata,sample_col,col_cell,cell_list,path_results=None,normalization=True,n_top_genes=2000,highly_variable_genes_=False):
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
            df.columns=list(highly_variable_genes)
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
        

def compute_diff_expressions(adata,cell_type: str = None,
                             proportions: pd.DataFrame = None,
                             selected_genes: list = None,
                             group1: str = 'Tumor 1',
                             group2: str = 'Tumor 2',
                             label_name: str = 'Predicted_Labels',
                             fc_thr: float = 0.5,
                             pval_thr: float = 0.01,
                             sample_col:str='sampleID',
                             col_cell:str ='cell_types',
                             path=None,
                             normalization=True,
                             n_top_genes=2000,
                             highly_variable_genes_=False

                             ):
    """
    Using limma R package, lmFit fits a linear model using weighted least squares for each gene
    Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models
    Empirical Bayes smoothing of standard errors (shrinks standard errors 
    that are much larger or smaller than those from other genes towards the average standard error)

    Parameters
    ----------
    cell_type : str, optional
        Specify cell type name to check its differential expression genes. The default is None.
    proportions : pd.DataFrame, optional
        cell types proportion in each samples. The default is None.
    selected_genes : list, optional
        Specify genes names to be considered to checking their differentiation
    group1 : str, optional
        Name of first patients sub-group for comparison. The default is 'Tumor 1'.
    group2 : str, optional
        Name of second patients sub-group for comparison. The default is 'Tumor 2'.
    label_name : str, optional
        name of the column containing the labels of patient sub-groups. The default is 'Predicted_Labels'.
    fc_thr : float, optional
        Specify the fold change threshold. The default is 0.5.
    pval_thr : float, optional
        Specify the asj. p-value thrshold. The default is 0.01.
    
    Returns
    -------
    Plot volcano plot of fold changes between two intrested patients sub-groups.
    Save statistical table of each gene
    Save significantly differentiate genes in each group

    """
    path_to_result='Results_PILOT'
    
    if os.path.exists(path_to_result + "/cells/" + cell_type + ".csv"):

        cells = pd.read_csv(path_to_result + "/cells/" + cell_type + ".csv", index_col = 0)  
   
    elif cell_type not in adata.uns:
        extract_cells_from_gene_expression_for_clustering(adata,sample_col=sample_col,col_cell=col_cell,cell_list=[cell_type],path_results=path,normalization=normalization,n_top_genes=n_top_genes,highly_variable_genes_=highly_variable_genes_)
        cells = pd.read_csv(path_to_result + "/cells/" + cell_type + ".csv", index_col = 0)
    
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
                 figsize = (15,15), output_path = path_to_result)
    

    
def install_r_packages():
    # Install R packages using rpy2
    import rpy2.robjects as robjects

    robjects.r('''
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    ''')

    robjects.r('''
    BiocManager::install("limma")
    ''')
   

   
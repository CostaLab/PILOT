#!python
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



def plot_cell_types_distributions(proportions: pd.DataFrame = None,
                                  cell_types: list = None,
                                  labels:str = 'Predicted_Labels',
                                  file_path: str = None,
                                  figsize: tuple = (15, 7)):
    """
    

    Parameters
    ----------
    proportions : pd.DataFrame, optional
        cell types proportion in each samples. The default is None.
    cell_types : list, optional
        number of cell types to be considered. The default is None.
    labels : str, optional
        name of the column containing the labels of patient sub-groups. The default is 'Predicted_Labels'.
    file_path : str, optional
        determine the path to store the figure. The default is None.
    figsize : tuple, optional
        determine the figure size. The default is (15, 7).

    Returns
    -------
    None.
    Save of the plot cell types distribution in each patients sub-groups

    """
    col_names = list(cell_types)
    col_names.append(labels)
    df_long = pd.melt(proportions.loc[:,col_names], labels)
    fig = plt.figure(figsize = figsize)
    sns.boxplot(x = "variable", hue = labels, y = "value", data = df_long)
    plt.ylabel("Cell proportions", fontsize = 24)
    plt.xlabel("")
    plt.xticks(fontsize = 24, rotation = 45, ha = 'right', rotation_mode = 'anchor')
    plt.legend(fontsize = 24)
    fig.tight_layout()
    plt.savefig(file_path + "/Cell_types_distributions.png",
                          facecolor = 'white')
    plt.show()
    
def plot_hor_vs_vert(data, subplot, x, y, c, xlabel, ylabel, rotation,
                     tick_bottom, tick_left, title):
    ax=plt.subplot(1,2,subplot)
    cols = ['tab:blue' if x >= 0 else 'tab:red' for x in data[x]]

    sns.barplot(x=x, y=y, data=data, ci=None, palette=cols)
    plt.title(title, fontsize=24, fontweight = 'bold')
    plt.xlabel(xlabel, fontsize=24)
    plt.xticks(fontsize=24, rotation=rotation)
    plt.ylabel(ylabel, fontsize=24)
    plt.yticks(fontsize=24)
    sns.despine(bottom=False, left=True)
    ax.grid(False)
    ax.tick_params(bottom=tick_bottom, left=tick_left)
    return None


def cell_type_diff_two_sub_patient_groups(proportions: pd.DataFrame = None,
                                          cell_types: list = None,
                                          labels:str = 'Predicted_Labels',
                                          group1: str = 'Tumor 1',
                                          group2: str = 'Tumor 2',
                                          pval_thr: float = 0.05,
                                          figsize: tuple = (15, 4),
                                          file_path: str = None):
    """
    

    Parameters
    ----------
    proportions : pd.DataFrame, optional
        cell types proportion in each samples. The default is None.
    cell_types : list, optional
        number of cell types to be considered. The default is None.
    labels : str, optional
        name of the column containing the labels of patient sub-groups. The default is 'Predicted_Labels'.
    group1 : str, optional
        Name of the first patients sub-group to check the differentiations. The default is 'Tumor 1'.
    group2 : str, optional
        Name of the second patients sub-group to check the differentiations. The default is 'Tumor 2'.
    pval_thr : float, optional
        P-value threshold. The default is 0.05.
    figsize : tuple, optional
        Determine the figure size. The default is (15, 4).
    file_path : str, optional
        Determine the path to store the figure. The default is None.

    Returns
    -------
    None.
    Plot the statistical scores showing the how significant each cell type 
    differentiated in each patients sub-group (filter based on p-value threshold)
    Save table of statisical information of all cell types,
    sorted by their p-value and statistic score

    """
      
    group1_proportions = proportions[proportions['Predicted_Labels'] == group1]
    group2_proportions = proportions[proportions['Predicted_Labels'] == group2]
    
    stats_group12 = []
    scores_group12 = []
    # using Welch’s t-test to detect the statistically significant difference
    ## between group 1 and group 2
    for cell_type in cell_types:
        ttest_result = ttest_ind(group1_proportions[cell_type],
                                 group2_proportions[cell_type], equal_var = False)
        scores_group12.append(ttest_result[0])
        stats_group12.append(ttest_result[1])
        
    # adjust the p-values
    stats_group12 = multipletests(list(stats_group12), method='fdr_bh')[1]

    stats_bar_group12 = pd.DataFrame()
    stats_bar_group12['cell_type'] = cell_types
    stats_bar_group12['adjPval'] = stats_group12
    stats_bar_group12['-logPval'] = -np.log(stats_group12)
    stats_bar_group12['score'] = scores_group12
    
    # sort the data based on p-value and statistic score
    stats_bar_group12 = stats_bar_group12.sort_values(by = ['score', '-logPval'],
                                                      ascending = [False, False])
    
    # save the table of statistical differentiation of two groups
    ## as positive score shows cell types statistically significant on the first group
    ## and negative scrore shows cell types statistically significant on the second group
    stats_bar_group12.to_csv(file_path + "/Cell_type_diff_" + group1 + "_vs_" + group2 +".csv",
                             header = True, index = None)
    
    # filter data based on a p-value threshold
    stats_bar_group12 = stats_bar_group12[stats_bar_group12['adjPval'] < pval_thr]
    
    
    fig, ax = plt.subplots(figsize = figsize)
    plot_hor_vs_vert(stats_bar_group12, 1, x = 'score', y = 'cell_type', c = 'type',
                     xlabel = 'statistic score', ylabel = None,
                     rotation = None, tick_bottom = True, tick_left = False,
                     title = "Cell type rank " + group1 + " vs " + group2)
    fig.tight_layout()
    plt.savefig(file_path + "/Cell_type_diff_" + group1 + "_vs_" + group2 +".png",
                          facecolor = 'white')
    
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
                 figsize = (20,10), output_path = None):
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
                         size = 'baseMean', sizes = (40, 400)
                        )

    ax.axhline(pv_thr, zorder = 0, c = 'k', lw = 2, ls = '--')
    ax.axvline(fc_thr, zorder = 0, c = 'k', lw = 2, ls = '--')
    ax.axvline(-fc_thr, zorder = 0, c = 'k', lw = 2, ls = '--')

    texts = []
    for i in range(len(df)):
        if df.iloc[i].nlog10 >= pv_thr and (df.iloc[i].log2FoldChange >= fc_thr):
            texts.append(plt.text(x = df.iloc[i].log2FoldChange, y = df.iloc[i].nlog10, s = df.iloc[i].symbol,
                                 fontsize = 14, weight = 'bold', family = 'sans-serif'))
        if df.iloc[i].nlog10 >= pv_thr and ( df.iloc[i].log2FoldChange <= -fc_thr):
            texts.append(plt.text(x = df.iloc[i].log2FoldChange, y = df.iloc[i].nlog10, s = df.iloc[i].symbol,
                                 fontsize = 16, weight = 'bold', family = 'sans-serif'))
    adjust_text(texts)

    custom_lines = [Line2D([0], [0], marker='o', color='w', markerfacecolor='#d62a2b', markersize=15),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor='#1f77b4', markersize=15)]

    plt.legend(custom_lines, ['Higher expressions in ' + feature2, 'Higher expressions in ' + feature1],loc = 1,
               bbox_to_anchor = (1,1.1), frameon = False, prop = {'weight': 'normal', 'size': 16})

    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(2)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.tick_params(width = 2)

    plt.title("Expression Score \n "+feature1+" - "+feature2, fontsize = 24)
    plt.xticks(size = 15, weight = 'bold')
    plt.yticks(size = 15, weight = 'bold')

    plt.xlabel("$log_{2}$ (Fold Change)", size = 18)
    plt.ylabel("-$log_{10}$ (P-value)", size = 18)

#     plt.savefig(filename, dpi = 100, bbox_inches = 'tight', facecolor = 'white')
    plt.savefig(output_path + "/volcano_" + str(feature1) + "-" + str(feature2) + "_FC.pdf",
                dpi = 100, bbox_inches = 'tight', facecolor = 'white')
    plt.show()

 


    """

    Parameters
    ----------
    EMD : W distance,
    proportions : pd.DataFrame, optional
        cell types proportion in each samples. The default is None.
    annot : dataframe, annoation file.
    real_labels : list, status of samples.

    res : float, resulotion for leiden clsutering.
    metric: str, metric for leiden clustering and calculating sil.

    groupby_col: name of the column to do groupby by scanpy.pl.heatmap and show the grouping based on.

    Other parametrs are for scanpy.pl.heatmap function. https://scanpy.readthedocs.io/en/stable/generated/scanpy.pl.heatmap.html


    Returns
    -------
    Proportion matrix with predicted lables for each sample based on leiden clustering over EMD data.
    """



def clustering_EMD(EMD,proportions,annot,real_labels,res=0.3,metric='cosine',groupby_col='status',swap_axes=False,cmap="Blues_r",dendrogram=True,show_gene_labels=False,var_group_rotation=90,figsize=[12,12],save=False):
    import anndata
    EMD_df=pd.DataFrame(EMD,columns=proportions.keys())
    EMD_df['sampleID']=proportions.keys()
    EMD_df['status']=list(real_labels)
    
    
    adata_emd = sc.AnnData(EMD)
    sc.pp.neighbors(adata_emd, metric=metric)
    sc.tl.leiden(adata_emd, resolution = res)
    predicted_labels = np.array(adata_emd.obs.leiden)
    Silhouette = Sil_computing(EMD/EMD.max(), predicted_labels,metric='cosine')
    
    proportion_df=pd.DataFrame(proportions)
    proportion_df=proportion_df.T
    proportion_df.columns=annot.cell_type.unique()
    
    proportion_df['Predicted_Labels']=predicted_labels
    proportion_df['sampIeD']=list(proportions.keys())
    
    
    
    
    
    
    EMD_df['Leiden']=predicted_labels
    if groupby_col=='status':
        sorter=np.unique(real_labels)
        EMD_df['status'] = EMD_df.status.astype("category")
        EMD_df['status'] = EMD_df['status'].cat.set_categories(sorter)
        EMD_df=EMD_df.sort_values(["status"])
    elif groupby_col=='Leiden':
        sorter=EMD_df.Leiden.unique()
        EMD_df['Leiden'] = EMD_df.Leiden.astype("category")
        EMD_df['Leiden'] = EMD_df['Leiden'].cat.set_categories(sorter)
        EMD_df=EMD_df.sort_values(["Leiden"])
    obs = pd.DataFrame()
    obs['sampleID']=EMD_df.sampleID.astype(str)
    obs['status']=EMD_df.status.astype(str)
    obs['Leiden']=EMD_df.Leiden.astype(str)
    df_genes = pd.DataFrame(index = EMD_df.columns[0:EMD_df.shape[0]])
    adata_emd = anndata.AnnData(X = EMD_df[ EMD_df.columns[0:EMD_df.shape[0]]].values, var =df_genes, obs = obs )
    sc.pl.heatmap(adata_emd,adata_emd.obs.sampleID,groupby=[groupby_col],swap_axes=swap_axes,cmap=cmap,dendrogram=dendrogram,show_gene_labels=show_gene_labels,var_group_rotation=var_group_rotation,figsize=figsize,save=save)
    return proportion_df


def extract_cells_from_gene_expression_for_clustering(adata,sample_col,col_cell,cell_list,path_results,normalization=True):
    for cell in cell_list:
        adata_new = adata[adata.obs[col_cell].isin([cell]),:]
        if normalization:
            sc.pp.normalize_total(adata_new, target_sum=1e4)
            sc.pp.log1p(adata_new)
        df=adata_new[:,adata_new.var_names].X
        df=pd.DataFrame(df.toarray())
        df.columns=adata_new.var_names
        df['sampleID']=list(adata_new.obs[sample_col])

        if not os.path.exists(path_results+'/Cells/'):
            os.makedirs(path_results+'/Cells/') 
        df.to_csv(path_results+'/Cells/'+cell+'.csv')
        

def compute_diff_expressions(cell_type: str = None,
                             proportions: pd.DataFrame = None,
                             selected_genes: list = None,
                             group1: str = 'Tumor 1',
                             group2: str = 'Tumor 2',
                             label_name: str = 'Predicted_Labels',
                             fc_thr: float = 0.5,
                             pval_thr: float = 0.01,
                             path_to_result: str = None):
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
    path_to_result : str, optional
        Path to store the results. The default is None.

    Returns
    -------
    Plot volcano plot of fold changes between two intrested patients sub-groups.
    Save statistical table of each gene
    Save significantly differentiate genes in each group

    """
    
    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    pandas2ri.activate()
    
    # prepare data for R
    proportions.index = proportions['sampIeD']
    cells = pd.read_csv(path_to_result + "/Cells/" + cell_type + ".csv", index_col = 0)
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
    res.to_csv(path_to_result + "/diff_expressions_stats_" + cell_type + ".csv")
      
    pv_thr = -np.log10(pval_thr)
    volcano_plot(cells[selected_genes].transpose(), res['logFC'], res['adj.P.Val'],
                 cell_type, group1, group2, fc_thr, pv_thr,
                 figsize = (15,15), output_path = path_to_result)
    

    
def gene_annotation_cell_type_subgroup(cell_type: str = None,
                                       group: str = None,
                                       path_to_results: str = None):
    """
    

    Parameters
    ----------
    cell_type : str, optional
        Specify cell type name to check its differential expression genes. The default is None.
    group : str, optional
        Name of patients sub-group of interest . The default is None.
    path_to_results : str, optional
        Path to store the results. The default is None.

    Returns
    -------
    None.
    Plot to show the most relative GO terms for specifc cell-type of determind patient sub-group
    """


    group_genes = pd.read_csv(path_to_results + "/significant_genes_"+cell_type+"_"+group+".csv",
                               index_col=0)
    gp = GProfiler(return_dataframe=True)
    gprofiler_results = gp.profile(organism = 'hsapiens',
                query = list(group_genes['0'].values))
    num_gos = 19
    if(gprofiler_results.shape[0] < 19):
        num_gos = gprofiler_results.shape[0]
    selected_gps = gprofiler_results.loc[0:num_gos,['name', 'p_value']]
    selected_gps['nlog10'] = -np.log10(selected_gps['p_value'].values)

    figsize = (15,12)

    plt.figure(figsize = figsize, dpi = 100)
    plt.style.use('default')
    sns.scatterplot(data = selected_gps, x="nlog10", y="name", s = 200, color = 'tab:blue')

    plt.title('GO enrichment in ' + cell_type + ' associated with ' + group, fontsize = 32)

    plt.xticks(size = 24)
    plt.yticks(size = 24)

    plt.ylabel("GO Terms", size = 24)
    plt.xlabel("-$log_{10}$ (P-value)", size = 24)
    if not os.path.exists(path_to_results + "/Cells"+"/GO_analysis/"):
            os.makedirs(path_to_results + "/Cells"+"/GO_analysis/") 
    plt.savefig(path_to_results + "/Cells"+"/GO_analysis/" + cell_type + "_" + group + ".png", bbox_inches = 'tight', facecolor='white', transparent=False)
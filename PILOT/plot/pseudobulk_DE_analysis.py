#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 14:43:35 2024

@author: Mina Shaigan
"""

import os
import pandas as pd
import numpy as np
import anndata as ad
import itertools
from sklearn.feature_selection import VarianceThreshold
from sklearn.decomposition import PCA
from matplotlib.lines import Line2D

from adjustText import adjust_text
from gprofiler import GProfiler
import textwrap as tw
import seaborn as sns
import matplotlib.pyplot as plt

import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
from rpy2.robjects import pandas2ri

from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
import logging
rpy2_logger.setLevel(logging.ERROR)
pandas2ri.activate()


def plot_cell_numbers(adata, proportion_df,
                      cell_type: str = None,
                      cluster_col: str = "Predicted_Labels",
                      celltype_col: str = "cell_types",
                      sample_col: str = "sampleID",
                      my_pal = None):
    """
    

    Parameters
    ----------
    adata : TYPE
        DESCRIPTION.
    proportion_df : TYPE
        DESCRIPTION.
    cell_type : str, optional
        DESCRIPTION. The default is None.
    cluster_col : str, optional
        DESCRIPTION. The default is "Predicted_Labels".
    celltype_col : str, optional
        DESCRIPTION. The default is "cell_types".
    sample_col : str, optional
        DESCRIPTION. The default is "sampleID".
    my_pal : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.

    """
    
    copy_cells = adata.obs.copy()
    copy_cells= copy_cells[copy_cells[celltype_col] == cell_type]
    copy_cells['group'] = proportion_df.loc[copy_cells[sample_col]][cluster_col].values
    data = copy_cells.groupby(['group', sample_col])[sample_col].count()
    data = pd.DataFrame(data)
    data = data.loc[~(data==0).all(axis=1)]

    n_groups = np.unique(data.index.get_level_values("group").values)
    if my_pal is None:
        if len(n_groups) == 3:
            my_pal = dict(zip(n_groups, ["tab:red", "skyblue", "tab:blue"]))
        else:
            my_pal = dict(zip(n_groups, sns.color_palette("tab10", len(n_groups))))

    plt.figure(figsize=(20,5))
    x_values = data.index.get_level_values(sample_col).values
    plt.bar(range(data.shape[0]),data[sample_col].values,
                      color = [my_pal[key] for key in data.index.get_level_values("group").values],
                      tick_label = x_values)
    plt.xticks(fontsize=24, rotation = 45, ha= 'center')
    plt.yticks(fontsize=24)

    plt.title(cell_type, fontsize = 24)
    plt.ylabel('Number of cells', fontsize = 24)    
    colors = my_pal     
    labels = list(colors.keys())
    handles = [plt.Rectangle((0,0),1,1, color=colors[label]) for label in labels]
    plt.legend(handles, labels, fontsize = 24)
    plt.show()
    
def compute_pseudobulk_DE(
        cluster_counts: pd.DataFrame = None,
        cluster_metadata: pd.DataFrame = None,
        group1: str = None,
        group2: str = None,
        cluster_col: str = None):
    
    """
    Parameters
    ----------
    aggr_counts : pd.DataFrame, optional
        DESCRIPTION. The default is None.
    metadata : pd.DataFrame, optional
        DESCRIPTION. The default is None.
    cell_type : str, optional
        DESCRIPTION. The default is None.
    group1 : str, optional
        DESCRIPTION. The default is None.
    group2 : str, optional
        DESCRIPTION. The default is None.
    n_cpus : int, optional
        DESCRIPTION. The default is 8.

    Returns
    -------
    my_stat_res : TYPE
        DESCRIPTION.

    """

    # consider DE between two group of interest
    my_cluster_metadata = cluster_metadata[ (cluster_metadata[cluster_col] == group1 ) | (cluster_metadata[cluster_col] == group2)]
    my_cluster_counts = cluster_counts.loc[my_cluster_metadata.index]

    R = robjects.r
    R('library(SingleCellExperiment)')
    R('library(DESeq2)')
    R('library(apeglm)')
    R('library(tidyverse, verbose = FALSE)')
    R.assign('cluster_counts', my_cluster_counts)
    R.assign('cluster_metadata', my_cluster_metadata)

    try:

        R('dds <- DESeqDataSetFromMatrix(round(t(cluster_counts)), colData = cluster_metadata, design = ~ stage)')
        R('rld <- rlog(dds, blind = TRUE)')
        R('dds <- DESeq(dds)')
        R(' \
        mylist <- list(resultsNames(dds)); \
        for(coef in resultsNames(dds)){ \
            if(coef != "Intercept"){ \
                print(coef); \
                res <- results(dds, name = coef, alpha = 0.05); \
                res <- lfcShrink(dds, coef = coef, res = res, type = "apeglm"); \
                res_tbl <- res %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% arrange(padj); \
                mylist[[coef]] <- res_tbl; \
            } \
        } \
        ')
        
        res = R('''mylist''')
        return res
    except rpy2.rinterface_lib.embedded.RRuntimeError:
        return None
    
def compute_pseudobulk_PCA(
        cluster_counts: pd.DataFrame = None,
        cluster_metadata: pd.DataFrame = None):
    
    """
    Parameters
    ----------
    aggr_counts : pd.DataFrame, optional
        DESCRIPTION. The default is None.
    metadata : pd.DataFrame, optional
        DESCRIPTION. The default is None.
    cell_type : str, optional
        DESCRIPTION. The default is None.
    group1 : str, optional
        DESCRIPTION. The default is None.
    group2 : str, optional
        DESCRIPTION. The default is None.
    n_cpus : int, optional
        DESCRIPTION. The default is 8.

    Returns
    -------
    my_stat_res : TYPE
        DESCRIPTION.

    """

    
    # consider DE between two group of interest

    R = robjects.r
    R('library(SingleCellExperiment)')
    R('library(DESeq2)')
    R('library(apeglm)')
    R('library(tidyverse, verbose = FALSE)')
    R.assign('cluster_counts', cluster_counts)
    R.assign('cluster_metadata', cluster_metadata)

    try:

        R('dds <- DESeqDataSetFromMatrix(round(t(cluster_counts)), colData = cluster_metadata, design = ~ stage)')
        R('rld <- rlog(dds, blind = TRUE)')
        R('dds <- DESeq(dds)')
        rld = R('''as.data.frame(assay(rld))''')        
        return rld
    except rpy2.rinterface_lib.embedded.RRuntimeError:
        return None
    
def plotPCA_subgroups(proportions, deseq2_counts, cell_type, my_pal, cluster_col):
    # consider top variances features
    selector = VarianceThreshold(0.2)
    new_deseq2_counts = selector.fit_transform(deseq2_counts)
    new_deseq2_counts = pd.DataFrame(new_deseq2_counts, index = deseq2_counts.index)
    
    # reduce dimension by PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(new_deseq2_counts)
    
    # plot PCA
    color_map = [my_pal[val] for val in proportions.loc[new_deseq2_counts.index, cluster_col]]
    fig, ax = plt.subplots()
    ax.scatter(X_pca[:, 0], X_pca[:, 1],
                c = color_map,
                cmap='viridis', edgecolor='k', s = 200)
    plt.xlabel('PC1: ' + str(round(pca.explained_variance_ratio_[0]*100)) + "% variance", fontsize = 24)
    plt.ylabel('PC2: ' + str(round(pca.explained_variance_ratio_[1]*100)) + "% variance", fontsize = 24)
    plt.title('PCA ' + str(cell_type))
    legend_elements = []
    for k in my_pal.keys():
        legend_elements.append(Line2D([0], [0], marker='o', color='w', label=k,
                              markerfacecolor=my_pal[k], markersize=15))
    ax.legend(handles=legend_elements, loc=1)
    plt.show()
    
def map_color_ps(a, low_fc_thrr, high_fc_thrr, pv_thrr):
    log2FoldChange, symbol, nlog10 = a
    if log2FoldChange >= high_fc_thrr and nlog10 >= pv_thrr:
        return 'very higher'
    elif log2FoldChange <= -low_fc_thrr and nlog10 >= pv_thrr:
        return 'very lower'
    else:
        return 'no'

def volcano_plot_ps(data, symbol, foldchange, p_value,
                 cell_type,
                 feature1,
                 feature2,
                 low_fc_thr = 1,
                 high_fc_thr = 1,
                 pv_thr = 1,
                 figsize = (20,10),
                 output_path = None,
                 my_pal = None,
                 fontsize: int = 14
                ):
    """
    

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    symbol : TYPE
        DESCRIPTION.
    foldchange : TYPE
        DESCRIPTION.
    p_value : TYPE
        DESCRIPTION.
    cell_type : TYPE
        DESCRIPTION.
    feature1 : TYPE
        DESCRIPTION.
    feature2 : TYPE
        DESCRIPTION.
    low_fc_thr : TYPE, optional
        DESCRIPTION. The default is 1.
    high_fc_thr : TYPE, optional
        DESCRIPTION. The default is 1.
    pv_thr : TYPE, optional
        DESCRIPTION. The default is 1.
    figsize : TYPE, optional
        DESCRIPTION. The default is (20,10).
    output_path : TYPE, optional
        DESCRIPTION. The default is None.
    my_pal : TYPE, optional
        DESCRIPTION. The default is None.
    fontsize : int, optional
        DESCRIPTION. The default is 14.

    Returns
    -------
    str
        DESCRIPTION.

    """
    
    df = pd.DataFrame(columns=['log2FoldChange', 'nlog10', 'symbol'])
    df['log2FoldChange'] = data[foldchange]
    df['nlog10'] = -np.log10(data[p_value].values)
    df['symbol'] = data[symbol].values
    
    color1 = my_pal[feature1]
    color2 = my_pal[feature2]    
    
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.dropna(subset=["nlog10"], how="all", inplace=True)
    

    selected_labels = df.loc[ (df.log2FoldChange <= low_fc_thr) & (df.log2FoldChange >= high_fc_thr) & \
                             (df['nlog10'] >= pv_thr)]['symbol'].values
    
    def map_shape(symbol):
        if symbol in selected_labels:
            return 'important'
        return 'not'
    
    df['color'] = df[['log2FoldChange', 'symbol', 'nlog10']].apply(map_color_ps, low_fc_thrr = low_fc_thr, 
                                                                   high_fc_thrr = high_fc_thr,
                                                                   pv_thrr = pv_thr, axis = 1)
    df['shape'] = df.symbol.map(map_shape)
    df['baseMean'] = df.nlog10*10

    
    plt.figure(figsize = figsize, frameon=False, dpi=100)
    plt.style.use('default')

    ax = sns.scatterplot(data = df, x = 'log2FoldChange', y = 'nlog10', 
                         hue = 'color', hue_order = ['no', 'very higher', 'very lower'],
                         palette = ['lightgrey', color2, color1],
                         style = 'shape', style_order = ['not', 'important'],
                         markers = ['o', 'o'], 
                         size = 'baseMean', sizes = (40, 400)
                        )

    ax.axhline(pv_thr, zorder = 0, c = 'k', lw = 2, ls = '--')
    ax.axvline(high_fc_thr, zorder = 0, c = 'k', lw = 2, ls = '--')
    ax.axvline(-low_fc_thr, zorder = 0, c = 'k', lw = 2, ls = '--')

    texts = []
    for i in range(len(df)):
        if df.iloc[i].nlog10 >= pv_thr and (df.iloc[i].log2FoldChange >= high_fc_thr):
            texts.append(plt.text(x = df.iloc[i].log2FoldChange, y = df.iloc[i].nlog10, s = df.iloc[i].symbol,
                                 fontsize = fontsize, weight = 'bold', family = 'sans-serif'))
        if df.iloc[i].nlog10 >= pv_thr and ( df.iloc[i].log2FoldChange <= -low_fc_thr):
            texts.append(plt.text(x = df.iloc[i].log2FoldChange, y = df.iloc[i].nlog10, s = df.iloc[i].symbol,
                                 fontsize = fontsize + 2, weight = 'bold', family = 'sans-serif'))
    adjust_text(texts)

    custom_lines = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color2, markersize=fontsize),
                   Line2D([0], [0], marker='o', color='w', markerfacecolor=color1, markersize=fontsize)]

    plt.legend(custom_lines, ['Higher expressions in ' + feature2, 'Higher expressions in ' + feature1], loc = 1,
               bbox_to_anchor = (1,1.1), frameon = False, prop = {'weight': 'normal', 'size': fontsize})

    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(2)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.tick_params(width = 2)
    ax.set_ylim(bottom=0)
    plt.title("Expression Score \n " + feature1 + " - " + feature2, fontsize = fontsize + 4)
    plt.xticks(size = fontsize, weight = 'bold')
    plt.yticks(size = fontsize, weight = 'bold')

    plt.xlabel("$log_{2}$ (Fold Change)", size = fontsize + 2)
    plt.ylabel("-$log_{10}$ (P-value)", size = fontsize + 2)

    if output_path is not None:
        plt.savefig(output_path + "/volcano_" + str(feature1) + "-" + str(feature2) + "_FC.pdf",
                    dpi = 100, bbox_inches = 'tight', facecolor = 'white')
    plt.show()
    
def gene_annotation_cell_type_subgroup(data: pd.DataFrame = None,
                                       symbol: str = 'gene',
                                       sig_col: str = 'significant_gene',
                                       cell_type: str = None,
                                       group: str = None,
                                       sources: str = None,
                                       num_gos: int = 10,
                                       fig_h: int = 6,
                                       fig_w: int = 4,
                                       font_size: int = 14,
                                       max_length:int = 50,
                                       path_to_results: str = None,
                                       my_pal = None
                                     ):
    """
    Plot to show the most relative GO terms for specifc cell-type of determind patient sub-group

    Parameters
    ----------
    data : pd.DataFrame
        DESCRIPTION. The default is None.
    symbol : str, optional
        DESCRIPTION. The default is 'gene'.
    sig_col : str, optional
        DESCRIPTION. The default is 'significant_gene'.
    cell_type : str
        DESCRIPTION. The default is None.
    group : str
        DESCRIPTION. The default is None.
    sources : str, optional
        DESCRIPTION. The default is None.
    num_gos : int, optional
        DESCRIPTION. The default is 10.
    fig_h : int, optional
        DESCRIPTION. The default is 6.
    fig_w : int, optional
        DESCRIPTION. The default is 4.
    font_size : int, optional
        DESCRIPTION. The default is 14.
    max_length : int, optional
        DESCRIPTION. The default is 50.
    path_to_results : str, optional
        DESCRIPTION. The default is None.
    my_pal : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.

    """

#     path_to_results = 'Results_PILOT'

    color = my_pal[group]

#     group_genes = pd.read_csv(path_to_results + \
#                               "/significant_genes_" + cell_type + "_" + group + ".csv")

    group_genes = data.loc[data[sig_col] == group, symbol].values
    gp = GProfiler(return_dataframe = True)
    if list(group_genes):
        gprofiler_results = gp.profile(organism = 'hsapiens',
                                       query = list(group_genes),
                                       no_evidences = False,
                                       sources = sources)
    else:
        return "Genes list is empty!"
    
    if(gprofiler_results.shape[0] == 0):
        return "Not enough information!"
    elif(gprofiler_results.shape[0] < num_gos):
        num_gos = gprofiler_results.shape[0]

    all_gprofiler_results = gprofiler_results.copy()
    # display(all_gprofiler_results.head())
       
    # print(len(list(group_genes['symbol'].values)))
    # selected_gps = gprofiler_results.loc[0:num_gos,['name', 'p_value']]
    selected_gps = gprofiler_results.head(num_gos)[['name', 'p_value']]
    
    selected_gps['nlog10'] = -np.log10(selected_gps['p_value'].values)

    for i in selected_gps.index:
        split_name = "\n".join(tw.wrap(selected_gps.loc[i, 'name'], max_length))
        selected_gps.loc[i, 'name'] = split_name
    
    figsize = (fig_h, fig_w)

    plt.figure(figsize = figsize, dpi = 100)
    plt.style.use('default')
    sns.scatterplot(data = selected_gps, x = "nlog10", y = "name", s = 300, color = color)

    plt.title('GO enrichment in ' + cell_type + ' associated with ' + group + \
              '\n (number of genes: ' + str(len(list(group_genes))) + ")", fontsize = font_size + 2)

    plt.xticks(size = font_size)
    plt.yticks(size = font_size)

    plt.ylabel("GO Terms", size = font_size)
    plt.xlabel("-$log_{10}$ (P-value)", size = font_size)
    
    save_path = path_to_results + '/'
    if not os.path.exists(save_path):
            os.makedirs(save_path)
#     plt.savefig(save_path + group + ".pdf", bbox_inches = 'tight',
#                 facecolor = 'white', transparent = False)
    plt.show()
    
    all_gprofiler_results.to_csv(save_path + group + ".csv")
    
def get_sig_genes(data, symbol, foldchange, p_value, cell_type,
                 feature1, feature2,
                 low_fc_thr = 1, high_fc_thr = 1, pv_thr = 1):
    df = pd.DataFrame(columns=['log2FoldChange', 'nlog10', 'symbol'])
    df['log2FoldChange'] = data[foldchange]
    df['nlog10'] = -np.log10(data[p_value].values)
    df['symbol'] = data[symbol].values
    
    df.replace([np.inf, -np.inf], np.nan, inplace = True)
    df.dropna(subset = ["nlog10"], how = "all", inplace = True)
    
    data['significant_gene'] = ""
    group1_selected_labels = df.loc[ (df.log2FoldChange <= -low_fc_thr) & (df['nlog10'] >= pv_thr), 'symbol'].values
    data.loc[data[symbol].isin(group1_selected_labels), 'significant_gene'] = feature1
    
    group2_selected_labels = df.loc[ (df.log2FoldChange >= high_fc_thr) & (df['nlog10'] >= pv_thr), 'symbol'].values
    data.loc[data[symbol].isin(group2_selected_labels), 'significant_gene'] = feature2

    return data

def get_pseudobulk_DE(adata: ad.AnnData,
                      proportion_df: pd.DataFrame,
                      cell_type: str,
                      fc_thr: list,
                      pv_thr: float = 0.05,
                      celltype_col: str = "cell_types",
                      sample_col: str = "sampleID",
                      cluster_col: str = "Predicted_Labels",
                      remove_samples: list = [],
                      my_pal: dict = None,
                      path_to_results: str = 'Results_PILOT/',
                      figsize: tuple = (30, 15),
                      num_gos: int = 10,
                      fig_h: int = 6,
                      fig_w: int = 4,
                      sources: list = ['GO:CC', 'GO:PB', 'GO:MF'],
                      fontsize: int = 14,
                      load: bool = False
                     ):
    """
    

    Parameters
    ----------
    adata : ad.AnnData
        DESCRIPTION.
    proportion_df : pd.DataFrame
        DESCRIPTION.
    cell_type : str
        DESCRIPTION.
    fc_thr : list
        DESCRIPTION.
    pv_thr : float, optional
        DESCRIPTION. The default is 0.05.
    celltype_col : str, optional
        DESCRIPTION. The default is "cell_types".
    sample_col : str, optional
        DESCRIPTION. The default is "sampleID".
    cluster_col : str, optional
        DESCRIPTION. The default is "Predicted_Labels".
    remove_samples : list, optional
        DESCRIPTION. The default is [].
    my_pal : dict, optional
        DESCRIPTION. The default is None.
    path_to_results : str, optional
        DESCRIPTION. The default is 'Results_PILOT/'.
    figsize : tuple, optional
        DESCRIPTION. The default is (30, 15).
    num_gos : int, optional
        DESCRIPTION. The default is 10.
    fig_h : int, optional
        DESCRIPTION. The default is 6.
    fig_w : int, optional
        DESCRIPTION. The default is 4.
    sources : list, optional
        DESCRIPTION. The default is ['GO:CC', 'GO:PB', 'GO:MF'].
    fontsize : int, optional
        DESCRIPTION. The default is 14.
    load : bool, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    

    n_clusters = np.unique(proportion_df[cluster_col])
    if my_pal is None:
        if len(n_clusters) == 3:
            my_pal = dict(zip(n_clusters, ["tab:red", "skyblue", "tab:blue"]))
        else:
            my_pal = dict(zip(n_clusters, sns.color_palette("tab10", len(n_clusters))))

    save_path = path_to_results + "/Diff_Expressions_Results/" + str(cell_type) + "/pseudobulk/"
    log_pv_thr = -np.log10(pv_thr)

    print("Plot cells frequency for each sample... ")
    plot_cell_numbers(adata, proportion_df, cell_type = cell_type,
                  cluster_col = cluster_col, celltype_col = celltype_col,
                      sample_col = sample_col, my_pal= my_pal)
    
    if load == False:
        print("Aggregating the counts and metadata to the sample level...")
        counts_df = adata.to_df()
        counts_df[[celltype_col, sample_col]] = adata.obs[[celltype_col, sample_col]].values
        
        aggr_counts = counts_df.groupby([celltype_col, sample_col]).sum()
    
    
        cluster_counts = aggr_counts.loc[cell_type]
        cluster_metadata = proportion_df.loc[cluster_counts.index.values]
        cluster_metadata['stage'] = cluster_metadata[cluster_col].values
    
        # remove unwanted samples
        if not (remove_samples is None):
            for sample in remove_samples:
                if sample in cluster_metadata.index:
                    cluster_metadata = cluster_metadata.drop(index = sample)
                if sample in cluster_counts.index:
                    cluster_counts = cluster_counts.drop(index = sample)
    
        cluster_metadata = cluster_metadata.loc[cluster_counts.index]
        cluster_counts = cluster_counts.loc[:, (cluster_counts != 0).any(axis=0)]
    
        print("Use the median of ratios method for count normalization from DESeq2")
        print("Use regularized log transform (rlog) of the normalized counts from DESeq2")
        rld = compute_pseudobulk_PCA(cluster_counts, cluster_metadata)
    
        if rld is not None:
            
            if not os.path.exists(save_path):
                os.makedirs(save_path)
            rld.to_csv(save_path + "rld_PCA.csv")
    else:
        rld = pd.read_csv(save_path + "rld_PCA.csv", index_col = 0)
        
    deseq2_counts = rld.transpose()
    print("Plot the first two principal components... ")
    plotPCA_subgroups(proportion_df, deseq2_counts, cell_type, my_pal, cluster_col)

    print("Performing the DE analysis... ")
    j = 0
    for groups in itertools.combinations(n_clusters, 2):
        data = None
        if load == False:
            res = compute_pseudobulk_DE(cluster_counts, cluster_metadata,
                                        group1 = groups[0],
                                        group2 = groups[1],
                                        cluster_col = cluster_col)
            if res is not None:
                with (robjects.default_converter + pandas2ri.converter).context():
                    data = robjects.conversion.get_conversion().rpy2py(res[1])

                data = get_sig_genes(data, 'gene', 'log2FoldChange', 'padj', cell_type, 
                                     groups[0], groups[1], fc_thr[j], fc_thr[j], log_pv_thr)
                
                data.to_csv(save_path + "/" + str(groups[1]) + "vs" + str(groups[0]) + "_DE.csv")
        else:
            data = pd.read_csv(save_path + "/" + str(groups[1]) + "vs" + str(groups[0]) + "_DE.csv", index_col = 0)

        if data is not None:
            print("Plot volcano plot for " + str(groups[1]) + " vs " + str(groups[0]))
            volcano_plot_ps(data, 'gene', 'log2FoldChange', 'padj', cell_type, 
                         groups[0], groups[1], fc_thr[j], fc_thr[j], log_pv_thr, figsize = figsize,
                         output_path = save_path + "/",
                         my_pal = my_pal, fontsize = fontsize)

            print("Plot GO analysis for " + str(groups[1]) + " vs " + str(groups[0]))
            gene_annotation_cell_type_subgroup(data, cell_type = cell_type, group = groups[0],
                                               sources = sources, num_gos = num_gos,
                                               fig_h = fig_h, fig_w = fig_w, font_size = fontsize,
                                               path_to_results = save_path + "/" + str(groups[1]) + "vs" + str(groups[0]) + "/GOs/",
                                               my_pal = my_pal)
            gene_annotation_cell_type_subgroup(data, cell_type = cell_type, group = groups[1],
                                               sources = sources, num_gos = num_gos,
                                               fig_h = fig_h, fig_w = fig_w, font_size = fontsize,
                                               path_to_results = save_path + "/" + str(groups[1]) + "vs" + str(groups[0]) + "/GOs/",
                                               my_pal = my_pal)
        j += 1
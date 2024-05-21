#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:29:07 2024

@author: Mina Shaigan
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import anndata as ad
from sklearn.preprocessing import StandardScaler
import scipy.cluster.hierarchy as sch
import matplotlib.ticker as tick
from scipy.stats import zscore
from scipy.stats import norm
from gprofiler import GProfiler
import textwrap as tw
import matplotlib.colors as pltcolors
from decimal import Decimal
import json
import requests

from .curve_activity import _curvesnamic_network_char_terminal_logfc_, \
    _curvesnamic_network_char_transient_logfc_, \
        _curvesnamic_network_char_switching_time_, \
            _curvesnamic_network_char_area_


def generate_feature_list(func_type, X):
    """
    construct the features matrix
    """
    assert func_type in ['linear', 'quadratic', 'linear_quadratic'], 'func type parameter must be linear, quadratic, cubic, sqrt, reciprocal'
    
    if func_type == 'linear':
        X = np.column_stack((np.ones(len(X)), X))
        
    if func_type == 'linear_quadratic':
        X = np.column_stack((np.ones(len(X)), X, np.power(X,2) ))
    
    if func_type == 'quadratic':
        X = np.column_stack((np.ones(len(X)), np.power(X,2) ))

    return X

def make_curves(table, points):
    """
    define curves

    Parameters
    ----------
    table : TYPE
        DESCRIPTION.
    points : TYPE
        DESCRIPTION.

    Returns
    -------
    df_TFs_curves : TYPE
        DESCRIPTION.

    """
    TFs_coefs = np.array([table.iloc[:,]['Intercept'], 
                          table.iloc[:,]['Treat'], 
                          table.iloc[:,]['Treat2']])    

    TFs_curves = np.empty((0,len(points)))
    for TFi in range(TFs_coefs.shape[1]):
        funType = table.iloc[TFi,:]['Fitted function']
        polyline = generate_feature_list(funType, points)
        if(funType == "linear" or funType == "quadratic"):
            TF_coefs = np.array([TFs_coefs[:,TFi][0], TFs_coefs[:,TFi][1]], dtype = float)
        else:
            TF_coefs = np.array([TFs_coefs[:,TFi][0], TFs_coefs[:,TFi][1], TFs_coefs[:,TFi][2]], dtype = float)
        TFs_curves = np.append(TFs_curves, np.array([list(np.matmul(polyline, TF_coefs))]), axis = 0)

    df_TFs_curves = pd.DataFrame(TFs_curves, index = table['Gene ID'], columns = points)
    return df_TFs_curves

def get_noised_curves(adata: ad.AnnData = None,
                      cell_type: str = None,
                      filter_table_feature: str = 'R-squared',
                      filter_table_feature_pval: str = 'adjusted P-value',
                      table_filter_thr: float = 0.1,
                      table_filter_pval_thr: float = 0.05,
                      path_to_results: str = 'Results_PILOT/'):
    """
      Based on the fitted function we generate the curve again by considering the
    variance of the cells at each disease progression time point

    Parameters
    ----------
    adata : ad.AnnData, optional
        DESCRIPTION. The default is None.
    cell_type : str, optional
        DESCRIPTION. The default is None.
    filter_table_feature : str, optional
        DESCRIPTION. The default is 'R-squared'.
    filter_table_feature_pval : str, optional
        DESCRIPTION. The default is 'adjusted P-value'.
    table_filter_thr : float, optional
        DESCRIPTION. The default is 0.1.
    table_filter_pval_thr : float, optional
        DESCRIPTION. The default is 0.05.
    path_to_results : str, optional
        DESCRIPTION. The default is 'Results_PILOT/'.

    Returns
    -------
    scaled_noised_curves : TYPE
        DESCRIPTION.
    pseudotime_sample_names : TYPE
        DESCRIPTION.

    """

    # Get cells of one cell_type 
    cells = pd.read_csv(path_to_results + "/cells/" + str(cell_type) + ".csv", index_col = 0)

    # Get the pseudotime points
    pseudotime_sample_names = cells[['sampleID', 'Time_score']].groupby('Time_score').first()
    pseudotime_sample_names = pseudotime_sample_names.sort_index()
    
    # Read the table of fitted function for each gene
    table = pd.read_csv(path_to_results + "/Markers/" + str(cell_type) + "/Whole_expressions.csv", index_col = 0)
    table = table.fillna(0)

    # Filter the table based on the r-square and its p-value
    selected_table = table[ (np.abs(table[filter_table_feature]) >= table_filter_thr) & \
                            (table[filter_table_feature_pval] <= table_filter_pval_thr) ]

    # Get curves of each gene
    curves = make_curves(selected_table, pseudotime_sample_names.index.values)

    # Compute cells standard deviation in each time point
    cells_genes = cells.loc[:, ~cells.columns.isin(['sampleID']) ]
    sample_genes_var = cells_genes.groupby('Time_score').std() / 10

    # Get the sum of the covariates for each gene
    sum_covariates = np.sum(list(selected_table[['Treat', 'Treat2', 'Intercept']].values * [1, 1, -1]), axis = 1)
    
    noise_values = sample_genes_var[selected_table['Gene ID']].mul(sum_covariates)
                                  
    noised_curves = curves + noise_values.transpose()                              
    noised_curves.columns = pseudotime_sample_names.sampleID

    noised_curves = noised_curves.fillna(0)
    
    scaler = StandardScaler()
    scaled_noised_curves = pd.DataFrame(scaler.fit_transform(noised_curves.transpose()).transpose())
    scaled_noised_curves.columns = noised_curves.columns
    scaled_noised_curves.index = noised_curves.index
    
    scaler = StandardScaler()
    scaled_curves = pd.DataFrame(scaler.fit_transform(curves.transpose()).transpose())
    scaled_curves.columns = curves.columns
    scaled_curves.index = curves.index
    
    return scaled_curves, scaled_noised_curves, pseudotime_sample_names

def cluster_genes_curves(curves: pd.DataFrame = None,
                         cluster_method: str = 'complete',
                         cluster_metric: str = 'correlation',
                         scaler_value: float = 0.65):
    """
    

    Parameters
    ----------
    curves : pd.DataFrame, optional
        DESCRIPTION. The default is None.
    cluster_method : str, optional
        DESCRIPTION. The default is 'average'.
    cluster_metric : str, optional
        DESCRIPTION. The default is 'euclidean'.
    scaler_value : float, optional
        DESCRIPTION. The default is 0.65.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    try:
        # retrieve clusters using fcluster 
        d = sch.distance.pdist(curves)
        L = sch.linkage(d, method = cluster_method, metric = cluster_metric)
        
        # scaler_value can be modified to retrieve more stringent or relaxed clusters
        clusters = sch.fcluster(L, scaler_value * d.max(), 'distance')
    except ValueError:
        print("Please change scaler_value... you're getting too many clusters")
        clusters = [1] * len(curves.index)

    return pd.DataFrame({'Gene ID': curves.index, 'cluster': clusters})

def plot_heatmap_curves(curves: pd.DataFrame = None,
                        genes_clusters: pd.DataFrame = None,
                        cluster_method: str = 'complete',
                        cluster_metric: str = 'correlation',
                        cmap_color: str = 'RdBu_r',
                        figsize: tuple = (7, 9),
                        fontsize: int = 14
                        ):
    """
    

    Parameters
    ----------
    curves : pd.DataFrame, optional
        DESCRIPTION. The default is None.
    genes_clusters : pd.DataFrame, optional
        DESCRIPTION. The default is None.
    cluster_method : str, optional
        DESCRIPTION. The default is 'average'.
    cluster_metric : str, optional
        DESCRIPTION. The default is 'euclidean'.
    cmap_color : str, optional
        DESCRIPTION. The default is 'RdBu_r'.
    figsize : tuple, optional
        DESCRIPTION. The default is (7, 9).
    fontsize : int, optional
        DESCRIPTION. The default is 14.

    Returns
    -------
    None.

    """

    my_palette = dict(zip(genes_clusters['cluster'].unique(),
                          sns.color_palette("tab10", len(genes_clusters['cluster'].unique()))))
    row_colors = genes_clusters['cluster'].map(my_palette)
    row_colors.index = genes_clusters['Gene ID']
    
    g = sns.clustermap(
    curves,
    method = cluster_method,
    figsize = figsize,
    metric = cluster_metric,
    col_cluster = False,
    yticklabels = False,
    xticklabels = True,
    z_score = 0,
    cmap = cmap_color,
    dendrogram_ratio = (0.2, 0),
    cbar_pos = (-0.08, 0.8, .03, .1),
    row_colors = row_colors.loc[curves.index],
    annot_kws = {"size": fontsize + 2}
    )
    
    reordered_labels = curves.index[g.dendrogram_row.reordered_ind].tolist()
    genes_clusters = genes_clusters.iloc[g.dendrogram_row.reordered_ind]
    clusters_numbers = (genes_clusters.groupby(genes_clusters['cluster'])
                        .apply(lambda x: x.iloc[(len(x))//2])
                       )
    use_ticks = []
    for label in clusters_numbers['Gene ID']:
        use_ticks.append(reordered_labels.index(label))
    
    g.ax_heatmap.set(yticks = use_ticks, yticklabels = clusters_numbers['cluster'])

    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.set_xlabel("Samples", fontsize = fontsize + 2)
    g.ax_heatmap.tick_params(axis='x', labelsize = fontsize)
    g.ax_heatmap.tick_params(axis='y', labelsize = fontsize)
    g.ax_cbar.tick_params(labelsize = fontsize)
    g.ax_row_colors.tick_params(labelsize = fontsize + 2)
    
def plot_each_cluster_activities(curves: pd.DataFrame = None,
                                 genes_clusters: pd.DataFrame = None,
                                 pseudotime_sample_names: pd.DataFrame = None,
                                 fontsize: int = 14
                                ):
    """
    

    Parameters
    ----------
    curves : pd.DataFrame, optional
        DESCRIPTION. The default is None.
    genes_clusters : pd.DataFrame, optional
        DESCRIPTION. The default is None.
    pseudotime_sample_names : pd.DataFrame, optional
        DESCRIPTION. The default is None.
    fontsize : int, optional
        DESCRIPTION. The default is 14.

    Returns
    -------
    None.

    """
    
    n_clusters = np.unique(genes_clusters['cluster'])
    fig, axs = plt.subplots(nrows = 1, ncols = len(n_clusters), figsize = (len(n_clusters) * 3, 2))
    plt.subplots_adjust(hspace = 0.5)

    j = 0
    tickers = list(n_clusters)

    if len(n_clusters) == 1:
        my_axs = [axs]
    else:
        my_axs = axs.ravel()
    
    for ticker, ax in zip(tickers, my_axs):
        cluster_features = genes_clusters[genes_clusters['cluster'] == ticker]['Gene ID'].values
        
        curves_mean = curves.loc[cluster_features].mean()
        curves_std = curves.loc[cluster_features].std()
    
        ax.plot(pseudotime_sample_names.index, curves_mean, '-')
        
        ax.xaxis.set_major_locator(tick.NullLocator())
        ax.yaxis.set_major_locator(tick.NullLocator())
        
        if(len(cluster_features) < 40 and len(cluster_features) >= 30):
            alpha = 2.042
        elif(len(cluster_features) < 60 and len(cluster_features) >= 40):
            alpha = 2.021
        elif(len(cluster_features) < 80 and len(cluster_features) >= 60):
            alpha = 2.0
        elif(len(cluster_features) < 100 and len(cluster_features) >= 80):
            alpha = 1.99
        elif(len(cluster_features) < 1000 and len(cluster_features) >= 100):
            alpha = 1.984
        else:
            alpha = 1.962
             
        ax.fill_between(pseudotime_sample_names.index,
                        curves_mean - alpha * (curves_std),
                        curves_mean + alpha * (curves_std),
                        alpha = 0.2)
        ax.set_title('cluster ' + str(ticker), fontsize = fontsize)

        if j == 0:
            ax.set_ylabel('gene expression', fontsize = fontsize)
        ax.set_xlabel('disease progression', fontsize = fontsize)
        j += 1
    
    plt.show()
    
def adjust_p_values(p_values):
    p = np.asfarray(p_values)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

def compute_curves_activities(curves: pd.DataFrame = None,
                              genes_clusters: pd.DataFrame = None,
                              pseudotime_sample_names: pd.DataFrame = None,
                              cell_type: str = None,
                              path_to_results: str = 'Results_PILOT/'):
    """
    

    Parameters
    ----------
    curves : pd.DataFrame, optional
        DESCRIPTION. The default is None.
    genes_clusters : pd.DataFrame, optional
        DESCRIPTION. The default is None.
    pseudotime_sample_names : pd.DataFrame, optional
        DESCRIPTION. The default is None.
    cell_type : str, optional
        DESCRIPTION. The default is None.
    path_to_results : str, optional
        DESCRIPTION. The default is 'Results_PILOT/'.

    Returns
    -------
    curves_activities : TYPE
        DESCRIPTION.

    """

    times = pseudotime_sample_names.index
    curves_activities = pd.DataFrame(index = curves.index,
                                     columns = ('Terminal_logFC', 'Terminal_pvalue', 'Terminal_adjPvalue',
                                                'Transient_logFC', 'Switching_time', 'area'))
    terminal_logFC = _curvesnamic_network_char_terminal_logfc_(np.array((curves)), times)
    curves_activities['Terminal_logFC'] = np.round(terminal_logFC.ravel(), 2)
    
    z_score = zscore(curves_activities['Terminal_logFC'])
    p_values = norm.sf(abs(z_score)) * 2
    curves_activities['Terminal_pvalue'] = p_values
    curves_activities['Terminal_adjPvalue'] = adjust_p_values(p_values)
    
    transient_logFC = _curvesnamic_network_char_transient_logfc_(np.array((curves)), times)
    curves_activities['Transient_logFC'] = np.round(transient_logFC.ravel(), 2)
    switching_time = _curvesnamic_network_char_switching_time_(np.array((curves)), times)
    curves_activities['Switching_time'] = np.round(switching_time.ravel(), 2)
    areas = _curvesnamic_network_char_area_(np.array((curves)), times)
    curves_activities['area'] = np.round(areas, 2)

    genes_clusters.index = genes_clusters['Gene ID'].values
    curves_activities['cluster'] = genes_clusters.loc[curves_activities.index, 'cluster']

    save_path = path_to_results + "/Markers/" + str(cell_type) 
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    curves_activities.to_csv(save_path + "/curves_activities.csv")
    
    return curves_activities

def plot_rank_genes_cluster(curves_activities: pd.DataFrame = None,
                            fontsize: int = 12):
    """
    

    Parameters
    ----------
    curves_activities : pd.DataFrame, optional
        DESCRIPTION. The default is None.
    fontsize : int, optional
        DESCRIPTION. The default is 12.

    Returns
    -------
    None.

    """
    
    rank_genes = curves_activities.sort_values(['Terminal_logFC', 'Terminal_pvalue'],
                                               ascending=[False, True],
                                               key = abs).groupby('cluster').head(10)
    clusters = rank_genes['cluster']
    n_clusters = np.unique(clusters)
    
    tickers = list(n_clusters)
    fig, axs = plt.subplots(nrows = 1, ncols = len(n_clusters), figsize = (len(n_clusters) * 3, 2))
    plt.subplots_adjust(hspace = 0.5)

    if len(n_clusters) == 1:
        my_axs = [axs]
    else:
        my_axs = axs.ravel()
    
    j = 0
    for ticker, ax in zip(tickers, my_axs):
    
        my_data = rank_genes[rank_genes['cluster'] == ticker].copy()
        my_data.sort_values(['Terminal_logFC', 'Terminal_pvalue'],
                            ascending=[True, False], inplace = True, key = abs)
        my_data['log_pval'] = -np.log10(my_data['Terminal_pvalue'].values)
    
        x = my_data['log_pval'].values
        y = my_data['Terminal_logFC'].values
        txts = my_data.index.values
    
        start = 0
        step = 1
        num = len(y)
    
        y = start + np.arange(0, num) * step
    
        ax.scatter(x, y + 0.05, color = 'k', s = 10)
        for i, txt in enumerate(txts):
            ax.annotate(txt, (x[i] + 0.05, y[i]), ha = 'left',
                        fontsize = fontsize - 2)
    
        margin = 1
        ax.set_xlim(0, np.ceil(np.max(x)) + margin)
        ax.set_ylim(-0.1, num * step + 0.2)
    
        ax.set_xlabel('$-log_{10}$(p-value)', fontsize = fontsize - 2)
        if j == 0:
            ax.set_ylabel('rank', fontsize = fontsize - 2)
        ax.tick_params(axis='x', labelsize = fontsize - 2)
        ax.tick_params(axis='y', labelsize = fontsize - 2)
    
        ax.set_title('cluster ' + str(ticker), fontsize = fontsize)
        j += 1
    
    plt.show()
    return rank_genes
    

def gene_annotation_cell_type_genes(cell_type: str = None,
                                    genes: list = None,
                                    group: str = None,
                                    num_gos: int = 10,
                                    fig_h: int = 6,
                                    fig_w: int = 4,
                                    font_size: int = 16,
                                    max_length:int = 50,
                                    sources: list = ['GO:CC', 'GO:PB', 'GO:MF']
                                    ):
    """
    
    Parameters
    ----------
    cell_type : str, optional
        DESCRIPTION. The default is None.
    genes : list, optional
        DESCRIPTION. The default is None.
    group : str, optional
        DESCRIPTION. The default is None.
    num_gos : int, optional
        DESCRIPTION. The default is 10.
    fig_h : int, optional
        DESCRIPTION. The default is 6.
    fig_w : int, optional
        DESCRIPTION. The default is 4.
    font_size : int, optional
        DESCRIPTION. The default is 16.
    max_length : int, optional
        DESCRIPTION. The default is 50.
    sources : list, optional
        DESCRIPTION. The default is ['GO:CC', 'GO:PB', 'GO:MF'].

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    gp = GProfiler(return_dataframe = True)
    if list(genes):
        gprofiler_results = gp.profile(organism = 'hsapiens', sources = sources,
                                       query = list(genes), no_evidences = False)
    else:
        return "Genes list is empty!"
    
    if(gprofiler_results.shape[0] == 0):
        return "Not enough information!"
    elif(gprofiler_results.shape[0] < num_gos):
        num_gos = gprofiler_results.shape[0]
        
    all_gprofiler_results = gprofiler_results.copy()    
    
    selected_gps = gprofiler_results.head(num_gos)[['name', 'p_value']]
    
    selected_gps['nlog10'] = -np.log10(selected_gps['p_value'].values)

    for i in selected_gps.index:
        split_name = "\n".join(tw.wrap(selected_gps.loc[i, 'name'], max_length))
        selected_gps.loc[i, 'name'] = split_name

    figsize = (fig_h, fig_w)

    plt.figure(figsize = figsize, dpi = 100)
    plt.style.use('default')
    sns.scatterplot(data = selected_gps, x = "nlog10", y = "name", s = 200, color = 'tab:blue')

    plt.title('GO enrichment in ' + cell_type + ' associated with ' + str(group) + \
              '\n Number of genes: ' + str(len(genes))
              , fontsize = font_size + 2)

    plt.xticks(size = font_size)
    plt.yticks(size = font_size)

    plt.ylabel("GO Terms", size = font_size)
    plt.xlabel("-$log_{10}$ (P-value)", size = font_size)
    plt.show()
    
    return all_gprofiler_results

def annotation_cluster_genes_by_curves(curves_activities: pd.DataFrame = None,
                                       cell_type: str = None,
                                       num_gos: int = 10,
                                       fig_h: int = 6,
                                       fig_w: int = 4,
                                       max_length: int = 50,
                                       sources: list = ['GO:CC', 'GO:PB', 'GO:MF'],
                                       path_to_results: str = 'Results_PILOT/',
                                       fontsize: int = 14):
    """
    

    Parameters
    ----------
    curves_activities : pd.DataFrame, optional
        DESCRIPTION. The default is None.
    cell_type : str, optional
        DESCRIPTION. The default is None.
    num_gos : int, optional
        DESCRIPTION. The default is 10.
    fig_h : int, optional
        DESCRIPTION. The default is 6.
    fig_w : int, optional
        DESCRIPTION. The default is 4.
    max_length : int, optional
        DESCRIPTION. The default is 50.
    sources : list, optional
        DESCRIPTION. The default is ['GO:CC', 'GO:PB', 'GO:MF'].
    path_to_results : str, optional
        DESCRIPTION. The default is 'Results_PILOT/'.
    fontsize : int, optional
        DESCRIPTION. The default is 14.

    Returns
    -------
    None.

    """
    
    for c in list(np.unique(curves_activities['cluster'])):
        genes = curves_activities.loc[curves_activities['cluster'] == c].index.values
        gprofiler_results = gene_annotation_cell_type_genes(cell_type, genes, "cluster " + str(c),
                                                            num_gos, fig_h, fig_w,
                                                            fontsize, max_length,
                                                            sources)
        
        if type(gprofiler_results) != str:
            save_path = path_to_results + "/Markers/" + str(cell_type) + "/GOs/"
            if not os.path.exists(save_path):
                os.makedirs(save_path)
            gprofiler_results.to_csv(save_path + "/cluster_" + str(c) + ".csv")
        else:
            print("No information for cluster " + str(c) + "!")
    
def  plot_top_genes_patterns(rank_genes: pd.DataFrame,
                             pseudotime_sample_names: pd.DataFrame,
                             cell_type: str,
                             curves: pd.DataFrame,
                             sample_col: str = 'sampleID',
                             time_col: str = 'Time_score',
                             path_to_results: str = 'Results_PILOT/',
                             plot_color: str = 'tab:orange',
                             points_color: str = 'viridis',
                             fontsize: str = 14):
    """
    

    Parameters
    ----------
    rank_genes : pd.DataFrame
        DESCRIPTION.
    pseudotime_sample_names : pd.DataFrame
        DESCRIPTION.
    cell_type : str
        DESCRIPTION.
    curves : pd.DataFrame
        DESCRIPTION.
    sample_col : str, optional
        DESCRIPTION. The default is 'sampleID'.
    time_col : str, optional
        DESCRIPTION. The default is 'Time_score'.
    path_to_results : str, optional
        DESCRIPTION. The default is 'Results_PILOT/'.
    plot_color : TYPE, optional
        DESCRIPTION. The default is 'tab:orange'.
    points_color : TYPE, optional
        DESCRIPTION. The default is 'viridis'.
    fontsize : TYPE, optional
        DESCRIPTION. The default is 14.

    Returns
    -------
    None.

    """
    
    n_clusters = np.unique(rank_genes['cluster'])
    rank_genes.sort_values(['Terminal_logFC', 'Terminal_pvalue'],
                           ascending=[True, False], inplace = True, key = abs)
    rank_genes['rank'] = rank_genes.groupby('cluster').cumcount(ascending=True)
    
    cells = pd.read_csv(path_to_results + "/cells/" + str(cell_type) + ".csv",
                        usecols = np.concatenate((np.array([sample_col, time_col]), rank_genes.index.values)))
    
    scaler = StandardScaler()
    scaled_cells = pd.DataFrame(scaler.fit_transform(cells.iloc[:, 2:]))
    scaled_cells.columns = cells.iloc[:, 2:].columns
    scaled_cells.index = cells.iloc[:, 2:].index
    scaled_cells[[sample_col, time_col]] = cells[[sample_col, time_col]].values

    table = pd.read_csv(path_to_results + "/Markers/" + str(cell_type) + "/Whole_expressions.csv", index_col = 0)
    table.index = table['Gene ID']
    table = table.loc[rank_genes.index]
    
    n_px = 10

    fig, axes = plt.subplots(10, len(n_clusters), constrained_layout = True,
                             figsize = (len(n_clusters) * n_px * 2, 10 * n_px))
    
    for c in n_clusters:
        for g in range(10):
            if (rank_genes[['cluster','rank']].values == [c, g]).all(axis=1).any():
                gene_name = rank_genes[(rank_genes[['cluster','rank']].values == [c, g]).all(axis=1)].index.values[0]
                axes[g, c - 1].scatter(scaled_cells[time_col], scaled_cells[gene_name], c = scaled_cells[gene_name],
                                       alpha = 0.5, cmap = points_color, s = 100 * len(n_clusters),
                                       norm = pltcolors.CenteredNorm(np.mean(scaled_cells[gene_name])))
                
                
                axes[g, c - 1].plot(pseudotime_sample_names.index.values, curves.loc[gene_name],
                                    c = plot_color, linewidth = 6.0)
                axes[g, c - 1].set_xticklabels(axes[g, c - 1].get_xticks().astype(int))

                axes[g, c - 1].set_title(gene_name, size = fontsize * len(n_clusters), weight = 'bold')

                for item in (axes[g, c - 1].get_xticklabels() + axes[g, c - 1].get_yticklabels()):
                        item.set_fontsize( (fontsize - 2) * len(n_clusters))

                if(float(table.loc[gene_name, 'Slope']) > 0):
                    axes[g, c - 1].text(.01, .99, 'adj. p-value = %.2e \n$R^{2}$ = %.2f' % (Decimal(table.loc[gene_name, 'adjusted P-value']),
                                                                                      table.loc[gene_name, 'mod_rsquared_adj'] ),
                                  ha = 'left', va = 'top',
                                  transform=axes[g, c - 1].transAxes, size = (fontsize - 2) * len(n_clusters))
                else:
                    axes[g, c - 1].text(.99, .99, 'adj. p-value = %.2e\n$R^{2}$ = %.2f' % (Decimal(table.loc[gene_name, 'adjusted P-value']),
                                                                                      table.loc[gene_name, 'mod_rsquared_adj'] ),
                                  ha = 'right', va = 'top',
                                  transform=axes[g, c - 1].transAxes, size = (fontsize - 2) * len(n_clusters))
                
                if c == 1:
                    axes[g, c - 1].set_ylabel('Gene expression', size = (fontsize - 2) * len(n_clusters))
            else:
                axes[g, c - 1].set_axis_off()
    
    plt.show()
  
def plot_gene_list_pattern(gene_list: list,
                            cell_type: str,
                            sample_col: str = 'sampleID',
                            time_col: str = 'Time_score',
                            path_to_results: str = 'Results_PILOT/',
                            plot_color: str = 'tab:orange',
                            points_color: str = 'viridis',
                            fontsize = 14):
    """
    

    Parameters
    ----------
    gene_list : list
        DESCRIPTION.
    cell_type : str
        DESCRIPTION.
    sample_col : str, optional
        DESCRIPTION. The default is 'sampleID'.
    time_col : str, optional
        DESCRIPTION. The default is 'Time_score'.
    path_to_results : str, optional
        DESCRIPTION. The default is 'Results_PILOT/'.
    plot_color : str, optional
        DESCRIPTION. The default is 'tab:orange'.
    points_color : str, optional
        DESCRIPTION. The default is 'viridis'.
    fontsize : TYPE, optional
        DESCRIPTION. The default is 14.

    Returns
    -------
    Plot
        Plot genes pattern for specific cell type

    """
    
    if len(gene_list) > len(set(gene_list)):
        return "The gene list is not unique!"
    
    try:
        cells = pd.read_csv(path_to_results + "/cells/" + str(cell_type) + ".csv",
                                usecols = np.concatenate((np.array([sample_col, time_col]), gene_list)))
    except ValueError:
        return "Some of the genes are not exists in cell type " + str(cell_type) + "!"
        
    scaler = StandardScaler()
    scaled_cells = pd.DataFrame(scaler.fit_transform(cells.iloc[:, 2:]))
    scaled_cells.columns = cells.iloc[:, 2:].columns
    scaled_cells.index = cells.iloc[:, 2:].index
    scaled_cells[[sample_col, time_col]] = cells[[sample_col, time_col]].values
    
    try:
        table = pd.read_csv(path_to_results + "/Markers/" + str(cell_type) + "/Whole_expressions.csv", index_col = 0)
        table.index = table['Gene ID']
        table = table.loc[gene_list]
    except KeyError:
        return "There is no good fit for some of the genes!"
    
    # Get the pseudotime points
    pseudotime_sample_names = cells[['sampleID', 'Time_score']].groupby('Time_score').first()
    pseudotime_sample_names = pseudotime_sample_names.sort_index()
    
    # Get curves of each gene
    curves = make_curves(table, pseudotime_sample_names.index.values)
    
    scaler = StandardScaler()
    scaled_curves = pd.DataFrame(scaler.fit_transform(curves.transpose()).transpose())
    scaled_curves.columns = curves.columns
    scaled_curves.index = curves.index
    
    n_px = 5
    n_clusters = range(3)
    
    fig, axes = plt.subplots(int(np.ceil(len(gene_list) / 3)), len(n_clusters), constrained_layout = True,
                             figsize=(len(n_clusters) * n_px * 2, int(np.ceil(len(gene_list) / 3)) * len(n_clusters) * n_px / 2))
    axes = np.atleast_2d(axes)
    
    
    for c in n_clusters:
        for g in range(int(np.ceil(len(gene_list) / 3))):
            if (g * len(n_clusters) + c) < len(gene_list):
                gene_name = gene_list[g * len(n_clusters) + c]
                axes[g, c].scatter(scaled_cells[time_col],
                                   scaled_cells[gene_name],
                                   c = scaled_cells[gene_name],
                                   alpha = 0.5, cmap = points_color,
                                   s = 100 * len(n_clusters),
                                   norm = pltcolors.CenteredNorm(np.mean(scaled_cells[gene_name])))
                
                
                axes[g, c].plot(pseudotime_sample_names.index.values,
                                scaled_curves.loc[gene_name],
                                c = plot_color, linewidth = 6.0)
                axes[g, c].set_xticklabels(axes[g, c].get_xticks().astype(int))
    
                axes[g, c].set_title(gene_name,
                                     size = (fontsize - 2) * len(n_clusters),
                                     weight = 'bold')
    
                for item in (axes[g, c].get_xticklabels() + axes[g, c].get_yticklabels()):
                        item.set_fontsize( (fontsize - 4) * len(n_clusters))
    
                if(float(table.loc[gene_name, 'Slope']) > 0):
                    axes[g, c].text(.01, .99, 'adj. p-value = %.2e \n$R^{2}$ = %.2f' % (Decimal(table.loc[gene_name, 'adjusted P-value']),
                                                                                      table.loc[gene_name, 'mod_rsquared_adj'] ),
                                  ha = 'left', va = 'top',
                                  transform=axes[g, c].transAxes,
                                  size = (fontsize - 4) * len(n_clusters))
                else:
                    axes[g, c].text(.99, .99, 'adj. p-value = %.2e\n$R^{2}$ = %.2f' % (Decimal(table.loc[gene_name, 'adjusted P-value']),
                                                                                      table.loc[gene_name, 'mod_rsquared_adj'] ),
                                    ha = 'right', va = 'top',
                                    transform=axes[g, c].transAxes,
                                    size = (fontsize - 4) * len(n_clusters))
                
                if c == 0:
                    axes[g, c].set_ylabel('Gene expression',
                                          size = (fontsize - 4) * len(n_clusters))
            else:
                axes[g, c].set_axis_off()
    
    plt.show()  
  
def genes_selection_analysis(
        adata: ad.AnnData,
        cell_type: str,
        filter_table_feature: str = 'R-squared',
        filter_table_feature_pval: str = 'adjusted P-value',
        table_filter_thr: float = 0.05,
        table_filter_pval_thr: float = 0.05,
        cluster_method: str = 'complete',
        cluster_metric: str = 'correlation',
        scaler_value: float = 0.4,
        cmap_color: str = 'RdBu_r',
        figsize: tuple = (7, 9),
        num_gos: int = 10,
        fig_h: int = 6,
        fig_w: int = 4,
        max_length: int = 50,
        sources: list = ['GO:CC', 'GO:PB', 'GO:MF'],
        fontsize: int = 14,
        path_to_results: str = 'Results_PILOT/'
    ):
    """
    

    Parameters
    ----------
    adata : ad.AnnData
        DESCRIPTION.
    cell_type : str
        DESCRIPTION.
    filter_table_feature : str, optional
        DESCRIPTION. The default is 'R-squared'.
    filter_table_feature_pval : str, optional
        DESCRIPTION. The default is 'adjusted P-value'.
    table_filter_thr : float, optional
        DESCRIPTION. The default is 0.05.
    table_filter_pval_thr : float, optional
        DESCRIPTION. The default is 0.05.
    cluster_method : str, optional
        DESCRIPTION. The default is 'average'.
    cluster_metric : str, optional
        DESCRIPTION. The default is 'euclidean'.
    scaler_value : float, optional
        DESCRIPTION. The default is 0.4.
    cmap_color : str, optional
        DESCRIPTION. The default is 'RdBu_r'.
    figsize : tuple, optional
        DESCRIPTION. The default is (7, 9).
    num_gos : int, optional
        DESCRIPTION. The default is 10.
    fig_h : int, optional
        DESCRIPTION. The default is 7.
    fig_w : int, optional
        DESCRIPTION. The default is 5.
    max_length : int, optional
        DESCRIPTION. The default is 50.
    sources : list, optional
        DESCRIPTION. The default is ['GO:CC', 'GO:PB', 'GO:MF'].
    fontsize : int, optional
        DESCRIPTION. The default is 14.
    path_to_results : str, optional
        DESCRIPTION. The default is 'Results_PILOT/'.

    Returns
    -------
    None.

    """

    print("Filter genes based on R-square and p-value...")
    curves, noised_curves, pseudotime_sample_names = get_noised_curves(adata, cell_type,
                                                                       filter_table_feature,
                                                                       filter_table_feature_pval,
                                                                       table_filter_thr,
                                                                       table_filter_pval_thr,
                                                                       path_to_results)

    print("Cluster genes using hierarchical clustering... ")
    genes_clusters = cluster_genes_curves(noised_curves,
                                          cluster_method,
                                          cluster_metric,
                                          scaler_value)

    print("Plot the heatmap of genes clustered... ")
    plot_heatmap_curves(noised_curves, genes_clusters,
                        cluster_method, cluster_metric,
                        cmap_color, figsize, fontsize)

    print("Plot patterns of clusters... ")
    print("Compute curves activities... ")
    print("Save curves activities... ")
    print("Plot top 10 genes for each cluster")
    print("Plot GO analysis for each cluster")
    
    plot_each_cluster_activities(noised_curves, genes_clusters,
                                 pseudotime_sample_names)
    curves_activities = compute_curves_activities(noised_curves, genes_clusters,
                                                  pseudotime_sample_names,
                                                  cell_type, path_to_results)
    rank_genes = plot_rank_genes_cluster(curves_activities, fontsize)
    
    plot_top_genes_patterns(rank_genes, pseudotime_sample_names, cell_type, curves,
                            path_to_results = path_to_results, fontsize = fontsize)
    
    annotation_cluster_genes_by_curves(curves_activities, cell_type, num_gos,
                                       fig_h, fig_w, max_length, sources,
                                       path_to_results, fontsize)

def genes_selection_heatmap(
        adata: ad.AnnData,
        cell_type: str,
        filter_table_feature: str = 'R-squared',
        filter_table_feature_pval: str = 'adjusted P-value',
        table_filter_thr: float = 0.05,
        table_filter_pval_thr: float = 0.05,
        cluster_method: str = 'complete',
        cluster_metric: str = 'correlation',
        scaler_value: float = 0.4,
        cmap_color: str = 'RdBu_r',
        figsize: tuple = (7, 9),
        fontsize: int = 14,
        path_to_results: str = 'Results_PILOT/'
    ):
    """
    

    Parameters
    ----------
    adata : ad.AnnData
        DESCRIPTION.
    cell_type : str
        DESCRIPTION.
    filter_table_feature : str, optional
        DESCRIPTION. The default is 'R-squared'.
    filter_table_feature_pval : str, optional
        DESCRIPTION. The default is 'adjusted P-value'.
    table_filter_thr : float, optional
        DESCRIPTION. The default is 0.05.
    table_filter_pval_thr : float, optional
        DESCRIPTION. The default is 0.05.
    cluster_method : str, optional
        DESCRIPTION. The default is 'complete'.
    cluster_metric : str, optional
        DESCRIPTION. The default is 'correlation'.
    scaler_value : float, optional
        DESCRIPTION. The default is 0.4.
    cmap_color : str, optional
        DESCRIPTION. The default is 'RdBu_r'.
    figsize : tuple, optional
        DESCRIPTION. The default is (7, 9).
    fontsize : int, optional
        DESCRIPTION. The default is 14.
    path_to_results : str, optional
        DESCRIPTION. The default is 'Results_PILOT/'.

    Returns
    -------
    None.

    """

    print("Filter genes based on R-square and p-value...")
    curves, noised_curves, pseudotime_sample_names = get_noised_curves(adata, cell_type,
                                                                       filter_table_feature,
                                                                       filter_table_feature_pval,
                                                                       table_filter_thr,
                                                                       table_filter_pval_thr,
                                                                       path_to_results)

    print("Cluster genes using hierarchical clustering... ")
    genes_clusters = cluster_genes_curves(noised_curves,
                                          cluster_method,
                                          cluster_metric,
                                          scaler_value)

    print("Compute curves activities... ")
    print("Save curves activities... ")
    curves_activities = compute_curves_activities(noised_curves, genes_clusters,
                              pseudotime_sample_names,
                              cell_type, path_to_results)

    print("Plot the heatmap of genes clustered... ")
    plot_heatmap_curves(noised_curves, genes_clusters,
                        cluster_method, cluster_metric,
                        cmap_color, figsize, fontsize)
    
    adata.uns['gene_selection_heatmap'] = {'cell_type': cell_type,
                                           'curves': curves,
                                           'noised_curves': noised_curves,
                                           'pseudotime_sample_names': pseudotime_sample_names,
                                           'curves_activities': curves_activities}

    

    
def plot_rank_genes_cluster_specific(adata: ad.AnnData,
                                     cell_type: str,
                                     cluster: int,
                                     fontsize: int = 14):
    """
    

    Parameters
    ----------
    adata : ad.AnnData
        DESCRIPTION.
    cell_type : str
        DESCRIPTION.
    cluster : int
        DESCRIPTION.
    fontsize : int, optional
        DESCRIPTION. The default is 14.

    Returns
    -------
    str
        DESCRIPTION.

    """
    
    save_cell_type = adata.uns['gene_selection_heatmap']['cell_type']
    if save_cell_type == cell_type:
        curves_activities = adata.uns['gene_selection_heatmap']['curves_activities']
        if cluster in curves_activities['cluster'].values:
            
            fig, axs = plt.subplots(nrows = 1, ncols = 2, figsize = (2 * 3, 3))
            
            ### plot average cluster pattern
            cluster_features = curves_activities[curves_activities['cluster'] == cluster].index.values
            curves = adata.uns['gene_selection_heatmap']['noised_curves']
            pseudotime_sample_names = adata.uns['gene_selection_heatmap']['pseudotime_sample_names']
            
            curves_mean = curves.loc[cluster_features].mean()
            curves_std = curves.loc[cluster_features].std()
        
            axs[0].plot(pseudotime_sample_names.index, curves_mean, '-')
            
            axs[0].xaxis.set_major_locator(tick.NullLocator())
            axs[0].yaxis.set_major_locator(tick.NullLocator())
            
            if(len(cluster_features) < 40 and len(cluster_features) >= 30):
                alpha = 2.042
            elif(len(cluster_features) < 60 and len(cluster_features) >= 40):
                alpha = 2.021
            elif(len(cluster_features) < 80 and len(cluster_features) >= 60):
                alpha = 2.0
            elif(len(cluster_features) < 100 and len(cluster_features) >= 80):
                alpha = 1.99
            elif(len(cluster_features) < 1000 and len(cluster_features) >= 100):
                alpha = 1.984
            else:
                alpha = 1.962
                 
            axs[0].fill_between(pseudotime_sample_names.index,
                            curves_mean - alpha * (curves_std),
                            curves_mean + alpha * (curves_std),
                            alpha = 0.2)
            axs[0].set_title('cluster {}'.format(cluster), fontsize = fontsize)
        
            axs[0].set_ylabel('gene expression', fontsize = fontsize)
            axs[0].set_xlabel('disease progression', fontsize = fontsize)
            
            ### plot top 10 genes of cluster
            rank_genes = curves_activities[curves_activities['cluster'] == cluster].sort_values(['Terminal_logFC',
                                                                                                 'Terminal_pvalue'],
                                                   ascending=[False, True],
                                                   key = abs).head(10)
        
            
            my_data = rank_genes.copy()
            my_data.sort_values(['Terminal_logFC', 'Terminal_pvalue'],
                                ascending=[True, False], inplace = True, key = abs)
            my_data['log_pval'] = -np.log10(my_data['Terminal_pvalue'].values)
        
            x = my_data['log_pval'].values
            y = my_data['Terminal_logFC'].values
            txts = my_data.index.values
        
            start = 0
            step = 1
            num = len(y)
        
            y = start + np.arange(0, num) * step
        
            axs[1].scatter(x, y + 0.05, color = 'k', s = 10)
            for i, txt in enumerate(txts):
                axs[1].annotate(txt, (x[i] + 0.05, y[i]), ha = 'left',
                                fontsize = fontsize - 2)
        
            margin = 1
            axs[1].set_xlim(0, np.ceil(np.max(x)) + margin)
            axs[1].set_ylim(-0.1, num * step + 0.2)
        
            axs[1].set_xlabel('$-log_{10}$(p-value)', fontsize = fontsize - 2)
            axs[1].set_ylabel('rank', fontsize = fontsize - 2)
            axs[1].tick_params(axis='x', labelsize = fontsize - 2)
            axs[1].tick_params(axis='y', labelsize = fontsize - 2)
        
            axs[1].set_title('cluster {}'.format(cluster), fontsize = fontsize)
        
            plt.tight_layout()
            plt.show()
            
            
        else:
            return "The cluster does not exist!"
    else:
        return "Please run the funtion genes_selection_heatmap first!"
    
    
    
def plot_top_genes_patterns_cluster_specific(adata: ad.AnnData,
                                             cell_type: str,
                                             cluster: int,
                                             sample_col: str = 'sampleID',
                                             time_col: str = 'Time_score',
                                             plot_color: str = 'tab:orange',
                                             points_color: str = 'viridis',
                                             fontsize: int = 14,
                                             path_to_results: str = 'Results_PILOT/'):
    """
    

    Parameters
    ----------
    adata : ad.AnnData
        DESCRIPTION.
    cell_type : str
        DESCRIPTION.
    cluster : int
        DESCRIPTION.
    sample_col : str, optional
        DESCRIPTION. The default is 'sampleID'.
    time_col : str, optional
        DESCRIPTION. The default is 'Time_score'.
    plot_color : str, optional
        DESCRIPTION. The default is 'tab:orange'.
    points_color : str, optional
        DESCRIPTION. The default is 'viridis'.
    fontsize : int, optional
        DESCRIPTION. The default is 14.
    path_to_results : str, optional
        DESCRIPTION. The default is 'Results_PILOT/'.

    Returns
    -------
    str
        DESCRIPTION.

    """
    save_cell_type = adata.uns['gene_selection_heatmap']['cell_type']
    if save_cell_type == cell_type:
    
        curves_activities = adata.uns['gene_selection_heatmap']['curves_activities']
        if cluster in curves_activities['cluster'].values:
            rank_genes = curves_activities[curves_activities['cluster'] == cluster]
            rank_genes = rank_genes.sort_values(['Terminal_logFC', 'Terminal_pvalue'], 
                                                ascending = [False, True], key = abs).head(10)
            
                
            my_data = rank_genes.copy()
            my_data.sort_values(['Terminal_logFC', 'Terminal_pvalue'],
                                ascending=[True, False], inplace = True, key = abs)
            plot_gene_list_pattern(my_data.index, cell_type,
                                   sample_col, time_col, path_to_results,
                                   plot_color, points_color, fontsize)
        else:
            return "The cluster does not exist!"
    else:
        return "Please run the funtion genes_selection_heatmap first!"
    
def annotation_genes_cluster_specific(adata: ad.AnnData,
                                      cell_type: str,
                                      cluster: int,
                                      num_gos: int = 10,
                                      fig_h: int = 6,
                                      fig_w: int = 4,
                                      max_length: int = 50,
                                      sources: list = ['GO:CC', 'GO:PB', 'GO:MF'],
                                      fontsize: int = 14,
                                      path_to_results: str = 'Results_PILOT/'):
    """
    

    Parameters
    ----------
    adata : ad.AnnData
        DESCRIPTION.
    cell_type : str
        DESCRIPTION.
    cluster : int
        DESCRIPTION.
    num_gos : int, optional
        DESCRIPTION. The default is 10.
    fig_h : int, optional
        DESCRIPTION. The default is 6.
    fig_w : int, optional
        DESCRIPTION. The default is 4.
    max_length : int, optional
        DESCRIPTION. The default is 50.
    sources : list, optional
        DESCRIPTION. The default is ['GO:CC', 'GO:PB', 'GO:MF'].
    fontsize : int, optional
        DESCRIPTION. The default is 14.
    path_to_results : str, optional
        DESCRIPTION. The default is 'Results_PILOT/'.

    Returns
    -------
    str
        DESCRIPTION.

    """
    save_cell_type = adata.uns['gene_selection_heatmap']['cell_type']
    if save_cell_type == cell_type:
    
        curves_activities = adata.uns['gene_selection_heatmap']['curves_activities']
        genes = curves_activities.loc[curves_activities['cluster'] == cluster].index.values
        gprofiler_results = gene_annotation_cell_type_genes(cell_type, genes, "cluster " + str(cluster),
                                                            num_gos, fig_h, fig_w,
                                                            fontsize, max_length,
                                                            sources)
        
        if type(gprofiler_results) != str:
            save_path = path_to_results + "/Markers/" + str(cell_type) + "/GOs/"
            if not os.path.exists(save_path):
                os.makedirs(save_path)
            gprofiler_results.to_csv(save_path + "/cluster_" + str(cluster) + ".csv")
        else:
            print("No information for cluster " + str(cluster) + "!")
    else:
        return "Please run the funtion genes_selection_heatmap first!"

def star_sig(x):
    if x < 0.001:
        return "**"
    elif x < 0.01:
        return "*"
    else:
        return 
        
def plot_hallmark_genes_clusters(adata: ad.AnnData,
                                 cell_type: str,
                                 gene_set_library: str = 'MSigDB_Hallmark_2020',
                                 cmap: str = 'coolwarm',
                                 font_size: int = 14):
    """
    

    Parameters
    ----------
    adata : ad.AnnData
        DESCRIPTION.
    cell_type : str
        DESCRIPTION.
    gene_set_library : str, optional
        DESCRIPTION. The default is 'MSigDB_Hallmark_2020'.
    cmap : str, optional
        DESCRIPTION. The default is 'coolwarm'.
    font_size : int, optional
        DESCRIPTION. The default is 14.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    None.

    """
    save_cell_type = adata.uns['gene_selection_heatmap']['cell_type']
    if save_cell_type == cell_type:
        curves_activities = adata.uns['gene_selection_heatmap']['curves_activities']
        col_names = ['rank', 'term', 'p-value', 'odds ratio', 'combined score',
                     'evidence genes', 'q-value', 'unknown1', 'unknown2']
        
        GO_terms = pd.DataFrame(columns = col_names)
        for cluster in np.unique(curves_activities['cluster']):
            ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
            genes_str = '\n'.join(curves_activities[curves_activities['cluster'] == cluster].index.values)
            description = 'Example gene list'
            payload = {
                'list': (None, genes_str),
                'description': (None, description)
            }
            
            response = requests.post(ENRICHR_URL, files=payload)
            if not response.ok:
                raise Exception('Error analyzing gene list')
            
            data = json.loads(response.text)
            
            ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
            query_string = '?userListId=%s&backgroundType=%s'
            user_list_id = data['userListId']
            gene_set_library = gene_set_library#'MSigDB_Hallmark_2020'#'KEGG_2021_Human'
            response = requests.get(
                ENRICHR_URL + query_string % (user_list_id, gene_set_library)
             )
            if not response.ok:
                raise Exception('Error fetching enrichment results')
            
            data = json.loads(response.text)
            data_df = pd.DataFrame(data[gene_set_library], columns = col_names)
            data_df['cluster'] = cluster
            GO_terms = pd.concat([GO_terms, data_df])
            
        GO_terms['cluster'] = GO_terms['cluster'].astype(int)
        GO_terms_clusters = pd.DataFrame(0, index = np.unique(GO_terms['term']), columns = np.unique(curves_activities['cluster']))
        for cluster in np.unique(curves_activities['cluster']):
            data = GO_terms[GO_terms['cluster'] == cluster]
            data.index = data['term'].values
            GO_terms_clusters.loc[data.index, cluster] = data['combined score'].values
            
        
        scaler = StandardScaler()
        scaled_GO_terms_clusters = pd.DataFrame(scaler.fit_transform(GO_terms_clusters.transpose()).transpose())
        scaled_GO_terms_clusters.columns = GO_terms_clusters.columns
        scaled_GO_terms_clusters.index = GO_terms_clusters.index
        
        GO_terms_annot = pd.DataFrame(0, index = np.unique(GO_terms['term']), columns = np.unique(curves_activities['cluster']))
        for cluster in np.unique(curves_activities['cluster']):
            data = GO_terms[GO_terms['cluster'] == cluster]
            data.index = data['term'].values
            GO_terms_annot.loc[data.index, cluster] = data['q-value'].values
            
        
        plt.figure(figsize=(scaled_GO_terms_clusters.shape[1],
                            int(np.round(scaled_GO_terms_clusters.shape[0]/3))))
        plt.title("{}".format(gene_set_library))
        g = sns.heatmap(scaled_GO_terms_clusters.sort_values(list(scaled_GO_terms_clusters.columns[::-1])),
                        cmap = cmap,
                        cbar_kws={"shrink": 0.3, 'label': 'Z-scaled score'},
                        yticklabels = True, center = 0, vmin = -1.5, vmax = 1.5,
                        annot = GO_terms_annot.map(star_sig).fillna(""), fmt="")
        g.set_xticklabels(g.get_xmajorticklabels(), fontsize = font_size);
        g.set_yticklabels(g.get_ymajorticklabels(), fontsize = font_size);    
    else:
        return "Please run the funtion genes_selection_heatmap first!"
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
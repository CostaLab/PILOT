#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:29:07 2024

@author: mina
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

import curve_activity as ca


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
    # define curves
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
    
    return noised_curves, pseudotime_sample_names

def cluster_genes_curves(curves: pd.DataFrame = None,
                         cluster_method: str = 'average',
                         cluster_metric: str = 'euclidean',
                         scaler_value: float = 0.65):
    
    scaler = StandardScaler()
    df = pd.DataFrame(scaler.fit_transform(curves.transpose()).transpose())
    df.columns = curves.columns
    df.index = curves.index

    try:
        # retrieve clusters using fcluster 
        d = sch.distance.pdist(df)
        L = sch.linkage(d, method = cluster_method, metric = cluster_metric)
        
        # scaler_value can be modified to retrieve more stringent or relaxed clusters
        clusters = sch.fcluster(L, scaler_value * d.max(), 'distance')
    except ValueError:
        print("Please change scaler_value... you're getting too many clusters")
        clusters = [1] * len(curves.index)

    return pd.DataFrame({'Gene ID': curves.index, 'cluster': clusters})

def plot_heatmap_curves(curves: pd.DataFrame = None,
                        genes_clusters: pd.DataFrame = None,
                        cluster_method: str = 'average',
                        cluster_metric: str = 'euclidean',
                        cmap_color: str = 'RdBu_r',
                        figsize: tuple = (7, 9),
                        fontsize: int = 14
                        ):

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

    times = pseudotime_sample_names.index
    curves_activities = pd.DataFrame(index = curves.index,
                                     columns = ('Terminal_logFC', 'Terminal_pvalue', 'Terminal_adjPvalue',
                                                'Transient_logFC', 'Switching_time', 'area'))
    terminal_logFC = ca._curvesnamic_network_char_terminal_logfc_(np.array((curves)), times)
    curves_activities['Terminal_logFC'] = np.round(terminal_logFC.ravel(), 2)
    
    z_score = zscore(curves_activities['Terminal_logFC'])
    p_values = norm.sf(abs(z_score)) * 2
    curves_activities['Terminal_pvalue'] = p_values
    curves_activities['Terminal_adjPvalue'] = adjust_p_values(p_values)
    
    transient_logFC = ca._curvesnamic_network_char_transient_logfc_(np.array((curves)), times)
    curves_activities['Transient_logFC'] = np.round(transient_logFC.ravel(), 2)
    switching_time = ca._curvesnamic_network_char_switching_time_(np.array((curves)), times)
    curves_activities['Switching_time'] = np.round(switching_time.ravel(), 2)
    areas = ca._curvesnamic_network_char_area_(np.array((curves)), times)
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
    
    rank_genes = curves_activities.sort_values('Terminal_pvalue').groupby('cluster').head(10)
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
        my_data.sort_values(['Terminal_logFC', 'Terminal_pvalue'], ascending=[True, False],
                                inplace = True, key = abs)
        my_data['log_pval'] = -np.log10(my_data['Terminal_pvalue'].values)
    
        x = my_data['log_pval'].values
        y = my_data['Terminal_logFC'].values
        txts = my_data.index.values
    
        start = 0
        step = 0.5
        num = len(y)
    
        y = start + np.arange(0, num) * step
    
        ax.scatter(x, y, color = 'white')
        for i, txt in enumerate(txts):
            ax.annotate(txt, (x[i], y[i]), ha='left', fontsize = fontsize - 2)
    
        margin = 0.5
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
    
def genes_selection_analysis(
    adata: ad.AnnData = None,
    cell_type: str = None,
    filter_table_feature: str = 'R-squared',
    filter_table_feature_pval: str = 'adjusted P-value',
    table_filter_thr: float = 0.05,
    table_filter_pval_thr: float = 0.05,
    cluster_method: str = 'average',
    cluster_metric: str = 'euclidean',
    scaler_value: float = 0.4,
    cmap_color: str = 'RdBu_r',
    figsize: tuple = (7, 9),
    fontsize: int = 14,
    path_to_results: str = 'Results_PILOT/'
    ):

    print("Filter genes based on R-square and p-value...")
    noised_curves, pseudotime_sample_names = get_noised_curves(adata, cell_type,
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
    print("Plot top 10 genes for cluster")
    
    plot_each_cluster_activities(noised_curves, genes_clusters,
                                 pseudotime_sample_names)
    curves_activities = compute_curves_activities(noised_curves, genes_clusters,
                                                  pseudotime_sample_names, cell_type,
                                                  path_to_results)
    plot_rank_genes_cluster(curves_activities, fontsize)
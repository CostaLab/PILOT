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
from pathlib import Path

from adjustText import adjust_text
from gprofiler import GProfiler
import textwrap as tw
import seaborn as sns
import matplotlib.pyplot as plt

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from .ploting import gene_annotation_cell_type_subgroup, volcano_plot


def plot_cell_numbers(adata, 
                      proportion_df,
                      cell_type: str = None,
                      cluster_col: str = "Predicted_Labels",
                      celltype_col: str = "cell_types",
                      sample_col: str = "sampleID",
                      my_pal = None):
    """   

    Parameters
    ----------
    adata : AnnData
        Annotated data object containing single-cell gene expression and metadata.
    proportion_df : pandas.DataFrame
        A dataframe indexed by sample IDs that contains group or cluster assignments for each sample.
    cell_type : str, optional
        The specific cell type to filter for from the `celltype_col`. 
        The default is None.
    cluster_col : str, optional
        The column name in `proportion_df` representing the group/cluster 
        labels used for coloring the bars. The default is "Predicted_Labels".
    celltype_col : str, optional
        The column name in `adata.obs` that contains cell type annotations. 
        The default is "cell_types".
    sample_col : str, optional
        The column name in `adata.obs` representing unique sample identifiers. 
        The default is "sampleID".
    my_pal : TYPE, optional
        A custom color palette mapping groups to colors. If None, it defaults 
        to a red/blue scheme for 3 groups or 'tab10' otherwise. 
        The default is None.

    Returns
    -------
    None.
        The function generates and displays a matplotlib bar plot.
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
        group1: str = "Tumor1",
        group2: str = "Tumor1",
        cluster_col: str = None,
        n_cpus: int = 8):
    
    """
    Computes Differential Expression using PyDESeq2, mimicking the R DESeq2 pipeline with apeglm LFC shrinkage.
    Performs Wald test for the contrast stage_group1_vs_group2 group2 as reference level).

    Parameters
    ----------
    aggr_counts : pd.DataFrame
        Pseudobulk counts matrix with rows = samples, columns = genes, non-negative integers. The default is None.
    metadata : pd.DataFrame
        Sample metadata with 'stage' column containing group1/group2 labels. The default is None.
    group1 : str, optional
        Name of experimental group (numerator of fold change). The default is "Tumor1".
    group2 : str, optional
        Name of reference/control group (denominator of fold change). The default is "Tumor2".
    n_cpus : int, optional
        Number of CPU cores for dispersion estimation and Wald tests. The default is 8.

    Returns
    -------
    dict
        Single key 'stage_{group1}_vs_{group2}' mapping to pandas DataFrame 
        with columns: 'gene', 'log2FoldChange', 'pvalue', 'padj', etc.
        Sorted by adjusted p-value (padj).


    Note
    ----
    Expects raw integer counts (DESeq2-style normalization applied internally).
    Sets explicit categorical reference: group2 < group1 for consistent 
    contrast direction matching DESeq2 coefficient naming.
    """

    # Filter metadata and counts for two groups of interest
    mask = cluster_metadata[cluster_col].isin([group1, group2])
    my_cluster_metadata = cluster_metadata[mask].copy()
    
    # Ensure counts match filtered metadata
    my_cluster_counts = cluster_counts.loc[my_cluster_metadata.index].copy()

    my_cluster_metadata["stage"] = pd.Categorical(
        #my_cluster_metadata["stage"],
        my_cluster_metadata[cluster_col],
        categories=[group2, group1],  # reference first
        ordered=True,
    )

    try:
        dds = DeseqDataSet(
            counts=my_cluster_counts.astype(int),
            metadata=my_cluster_metadata,
            design_factors="stage", 
            refit_cooks=True,
            n_cpus=n_cpus
        )

        # Run DESeq2 pipeline
        dds.deseq2()

         # Compute results, apply LFC Shrinkage (apeglm)
        stat_res = DeseqStats(dds, 
                              contrast=['stage', group2, group1], 
                              n_cpus=n_cpus)
        
        # Executes Wald test and applies shrinkage
        stat_res.summary()

        res_df = stat_res.results_df
        res_df = res_df.reset_index().rename(columns={'index': 'gene'})
        res_df = res_df.sort_values('padj')

        return {f"stage_{group1}_vs_{group2}": res_df}

    except Exception as e:
        print(f"Error during DE analysis: {e}")
        return None


def compute_pseudobulk_PCA(
        cluster_counts: pd.DataFrame = None,
        cluster_metadata: pd.DataFrame = None):
    """
    Computes regularized log-transformed pseudobulk counts for PCA using PyDESeq2.
    Replaces DESeq2's rlog() with PyDESeq2's equivalent regularized transformation.
    
    Parameters
    ----------
    cluster_counts : pd.DataFrame
        Pseudobulk counts matrix (rows = samples, columns = genes). The default is None.
    cluster_metadata : pd.DataFrame  
        Sample metadata with 'stage' column. The default is None.

    Returns
    -------
    pd.DataFrame
        Regularized log-transformed counts (genes x samples), suitable for PCA
        Column names = sample names, row names = gene names
    """
    
    try:
        # Initialize DeseqDataSet (PyDESeq2 expects samples x genes)
        dds = DeseqDataSet(
            counts=cluster_counts.astype(int),
            metadata=cluster_metadata,
            design_factors="stage",
            refit_cooks=True
        )
        dds.deseq2()

        # Get regularized log-transformed counts 
        dds.vst()

        rld_df = pd.DataFrame(
            dds.layers["vst_counts"].T,
            index=cluster_counts.columns,    # genes as index (from input)
            columns=cluster_counts.index     # samples as columns (from input)
        )
        return rld_df
        
    except Exception as e:
        print(f"Error during rlog computation: {e}")
        return None
    

def plotPCA_subgroups(proportions, deseq2_counts, cell_type, my_pal, cluster_col):

    """
    Performs feature selection via variance thresholding and visualizes 
    samples using Principal Component Analysis (PCA).

    Parameters
    ----------
    proportions : pandas.DataFrame
        A dataframe containing metadata or cluster assignments for the samples. 
        Must be indexed by sample IDs that match the index of `deseq2_counts`.
    deseq2_counts : pandas.DataFrame
        A dataframe of normalized expression counts (or similar features) 
        where rows are samples and columns are features (e.g., genes).
    cell_type : str
        The name of the cell type being analyzed, used primarily for the 
        plot title.
    my_pal : dict
        A dictionary mapping cluster/group labels to specific hex colors 
        or matplotlib color names.
    cluster_col : str
        The column name in the `proportions` dataframe used to color 
        the samples in the PCA plot.

    Returns
    -------
    None
        Displays a PCA scatter plot with explained variance ratios on 
        the axes and a custom legend.
    """

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


def get_sig_genes(data, 
                  symbol: str, 
                  foldchange: str, 
                  p_value: str, 
                  feature1: str, 
                  feature2: str,
                  low_fc_thr: float = 1, 
                  high_fc_thr: float = 1, 
                  pv_thr: float = 1):
    """
    Identifies and labels significantly differentially expressed genes based on 
    log2 fold-change and p-value thresholds.

    Parameters
    ----------
    data : AnnData
        Annotated data object containing single-cell gene expression and metadata..
    symbol : str
        The column name in `data` containing gene symbols or identifiers.
    foldchange : str
        The column name in `data` containing log2 fold-change values.
    p_value : str
        The column name in `data` containing p-values (unadjusted or adjusted).
    feature1 : str
        The label to assign to significantly down-regulated genes (e.g., "Downregulated").
    feature2 : str
        The label to assign to significantly up-regulated genes (e.g., "Upregulated").
    low_fc_thr : float, optional
        The absolute log2 fold-change threshold for down-regulation. 
        Genes must be <= -low_fc_thr. The default is 1.
    high_fc_thr : float, optional
        The absolute log2 fold-change threshold for up-regulation. 
        Genes must be >= high_fc_thr. The default is 1.
    pv_thr : float, optional
        The threshold for significance on the -log10(p-value) scale. 
        The default is 1 (which corresponds to p < 0.1).

    Returns
    -------
    pandas.DataFrame
        The original dataframe with an additional 'significant_gene' column containing 
        the labels provided in `feature1` and `feature2` for genes passing thresholds.
    """

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
                      path_to_results: str = 'Results_PILOT',
                      figsize: tuple = (15, 15),
                      num_gos: int = 10,
                      fig_h: int = 6,
                      fig_w: int = 4,
                      sources: list = ['GO:CC', 'GO:PB', 'GO:MF'],
                      fontsize: int = 14,
                      load: bool = False,
                      label_mode: str = "all",
                      n_p: int = 10,
                      n_n: int = 10
                     ):
    """
    Pseudobulk differential expression analysis pipeline using PyDESeq2.
    Aggregates single-cell counts to sample-level pseudobulk, performs PCA QC,
    runs pairwise DE between clusters, volcano plots, and GO enrichment.
    
    Parameters
    ----------
    adata : ad.AnnData
        Single-cell count matrix and metadata.
    proportion_df : pd.DataFrame
        Cell type proportions per sample (index = sampleID).
    cell_type : str
        Cell type/cluster for pseudobulk aggregation and DE.
    fc_thr : list
        Fold-change thresholds for volcano plots (one per cluster pair).
    pv_thr : float, optional
        Adjusted p-value threshold. The default is 0.05.
    celltype_col : str, optional
        Column name for cell types in adata.obs. The default is "cell_types".
    sample_col : str, optional
        Column name for sample IDs in adata.obs. The default is "sampleID".
    cluster_col : str, optional
        Column name for cluster labels in proportion_df. The default is "Predicted_Labels".
    remove_samples : list, optional
        Sample IDs to exclude from analysis. The default is [].
    my_pal : dict, optional
        Color palette for clusters. The default is None.
    path_to_results : str, optional
        Base directory for saving results. The default is 'Results_PILOT/'.
    figsize : tuple, optional
        Volcano plot figure size. The default is (30, 15).
    num_gos : int, optional
        Number of top GO terms to plot. The default is 10.
    fig_h : int, optional
        GO plot height. The default is 6.
    fig_w : int, optional
        GO plot width. The default is 4.
    sources : list, optional
        GO term sources. The default is ['GO:CC', 'GO:PB', 'GO:MF'].
    fontsize : int, optional
        Plot font size. The default is 14.
    load : bool, optional
        Load precomputed results instead of recomputing. The default is False.
    label_mode : {'all', 'topN'}, optional
        How to select labels for display.
        'all': label all points that pass the thresholds.
        'topN': label only the top n_p and n_n most significant points per side.
    n_p : int, optional
        Number of highest‑significance points on the positive log2FC side to label (used when label_mode == 'topN').
    n_n : int, optional
        Number of highest‑significance points on the negative log2FC side to label (used when label_mode == 'topN').

    Returns
    -------
    None
        Saves PCA, DE results, volcano plots, and GO analyses to disk.

    Notes
    -----
    - Uses PyDESeq2's rlog_norm() for PCA-ready normalized counts
    - Performs all pairwise DE between unique clusters in proportion_df[cluster_col]
    - Assumes raw counts in adata.X/layers (no normalization applied before aggregation)
    - Saves results to: Results_PILOT/Diff_Expressions_Results/{cell_type}/pseudobulk/
    """

    # Color palette for clusters
    n_clusters = np.unique(proportion_df[cluster_col])
    if my_pal is None:
        if len(n_clusters) == 3:
            my_pal = dict(zip(n_clusters, ["tab:red", "skyblue", "tab:blue"]))
        else:
            my_pal = dict(zip(n_clusters, sns.color_palette("tab10", len(n_clusters))))

    # Output dir
    save_path = Path(path_to_results) / "Diff_Expressions_Results" / str(cell_type) / "pseudobulk"
    log_pv_thr = -np.log10(pv_thr)

    # Plot cell frequency QC
    print("Plot cells frequency for each sample... ")
    plot_cell_numbers(adata, proportion_df, cell_type=cell_type,
                     cluster_col=cluster_col, celltype_col=celltype_col,
                     sample_col=sample_col, my_pal=my_pal)
    
    # Aggregate to pseudobulk and compute rlog-normalized counts
    if load == False:
        print("Aggregating the counts and metadata to the sample level...")
        counts_df = adata.to_df()
        counts_df[[celltype_col, sample_col]] = adata.obs[[celltype_col, sample_col]].values
        
        aggr_counts = counts_df.groupby([celltype_col, sample_col]).sum()
        cluster_counts = aggr_counts.loc[cell_type]
        cluster_metadata = proportion_df.loc[cluster_counts.index.values]
        cluster_metadata['stage'] = cluster_metadata[cluster_col].values
    
        # Remove unwanted samples
        if remove_samples:
            for sample in remove_samples:
                cluster_metadata = cluster_metadata.drop(index=sample, errors='ignore')
                cluster_counts = cluster_counts.drop(index=sample, errors='ignore')
    
        cluster_metadata = cluster_metadata.loc[cluster_counts.index]
        cluster_counts = cluster_counts.loc[:, (cluster_counts != 0).any(axis=0)]
    
        print("Computing rlog-normalized counts using PyDESeq2...")
        rld = compute_pseudobulk_PCA(cluster_counts, cluster_metadata)
    
        if rld is not None:
            save_path.mkdir(parents=True, exist_ok=True)
            rld.to_csv(save_path / "rld_PCA.csv")
    else:
        rld = pd.read_csv(save_path / "rld_PCA.csv", index_col=0)
        
    # PCA plot QC
    deseq2_counts = rld.transpose()
    print("Plot the first two principal components... ")
    plotPCA_subgroups(proportion_df, deseq2_counts, cell_type, my_pal, cluster_col)

    # Pairwise DE analysis
    print("Performing the DE analysis... ")
    j = 0
    for groups in itertools.combinations(n_clusters, 2):
        data = None
        if load == False:
            # Use adapted PyDESeq2 function1_adapted_to_function2
            res = compute_pseudobulk_DE(cluster_counts, cluster_metadata,
                                        group1=groups[0], group2=groups[1],
                                        cluster_col=cluster_col)
            
            if res is not None:
                # Extract single DataFrame 
                data = list(res.values())[0]  

                data = get_sig_genes(data, 'gene', 'log2FoldChange', 'padj', 
                                   groups[0], groups[1], fc_thr[j], fc_thr[j], log_pv_thr)
                
                data.to_csv(save_path / f"{str(groups[1])}_VS_{str(groups[0])}_DE.csv")
        else:
            data = pd.read_csv(save_path / f"{str(groups[1])}_VS_{str(groups[0])}_DE.csv", index_col=0)

        # Plot results
        if data is not None:
            print(f"Plot volcano plot for {str(groups[1])} vs {str(groups[0])}")
            volcano_plot(data=data, 
                         symbol_col='gene', 
                         fc_col='log2FoldChange', 
                         pval_col='padj', 
                         cell_type=cell_type, 
                         feature1=groups[0], 
                         feature2=groups[1], 
                         fc_thr=fc_thr[j], 
                         pv_thr=log_pv_thr, 
                         label_mode=label_mode,
                         n_p=n_p,
                         n_n=n_n,
                         my_pal=my_pal, 
                         figsize=figsize, 
                         font_size=fontsize,
                         output_path=save_path
                        )

            # GO analysis for both groups
            print(f"Plot GO analysis for {str(groups[1])} vs {str(groups[0])}")
            Path(save_path / f"{str(groups[1])}vs{str(groups[0])}" / "GOs").mkdir(parents=True, exist_ok=True)
            for group in [groups[0], groups[1]]:
                print("About to run gene annotation...")
                print(cell_type)
                gene_annotation_cell_type_subgroup(data=data,
                                                   cell_type=cell_type, 
                                                   group=group,
                                                   sources=sources, 
                                                   num_gos=num_gos,
                                                   figsize=(fig_h,fig_w),
                                                   font_size=fontsize,
                                                   path_to_results=Path(save_path) / f"{str(groups[1])}vs{str(groups[0])}" / "GOs",
                                                   my_pal=my_pal)
        j += 1

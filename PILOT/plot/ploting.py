import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import os
import seaborn as sns
import scipy
from matplotlib import pyplot as plt
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from adjustText import adjust_text
from matplotlib.lines import Line2D
from gprofiler import GProfiler
import pydiffmap
from sklearn import metrics
from pydiffmap import diffusion_map
from sklearn.neighbors import NearestCentroid
from sklearn.metrics.cluster import rand_score
from scipy.spatial import distance
from sknetwork.clustering import Louvain
from sklearn.preprocessing import label_binarize
import time
import sys 
from genericpath import isfile
from matplotlib.image import imread
import warnings
import ot
from logging import info, warn
from cycler import cycler
warnings.filterwarnings('ignore')
import elpigraph 
from scipy.stats import ttest_ind
from matplotlib.font_manager import FontProperties
from ..tools.Gene_cluster_specific_functions import *



def trajectory(adata,n_evecs = 2, epsilon =1, alpha = 0.5,knn= 64, sample_col=1, clusters = 'status',label_act = False,colors=['#377eb8','#ff7f00','#e41a1c'],location_labels='center', figsize=(12,12),font_size=24,axes_line_width=1,axes_color='black',facecolor='white',point_size=100,cmap='viridis',fontsize_legend=24,alpha_trans=1,plot_titel = "Trajectory of the disease progression"):
    
    
    
    
    """
    Find trajectories using Diffusion Maps and visualize the trajectory plot.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing EMD, annotations, and other necessary data.
    n_evecs : int, optional
        Number of embedding vectors, by default 2.
    epsilon : float or str, optional
        Method for choosing the epsilon in diffusion maps, by default 1.
    alpha : float, optional
        Normalization parameter in the bandwidth function, by default 0.5.
    knn : int, optional
        Number of nearest neighbors for constructing the kernel, by default 64.
    sample_col : int, optional
        Index of the column representing sample IDs in annotation data, by default 1.
    clusters : str, optional
        Name of the column representing clusters/categories in annotation data, by default 'status'.
    label_act : bool, optional
        Whether to label data points with sample IDs, by default False.
    colors : list, optional
        List of colors for different categories, by default ['#377eb8', '#ff7f00', '#e41a1c'].
    location_labels : str, optional
        Location of the legend labels, by default 'center'.
    figsize : tuple, optional
        figsize, by default (12,12).
    font_size : int, optional
        Font size for labels and annotations, by default 24.
    axes_line_width : float, optional
        Line width of the axes, by default 1.
    axes_color : str, optional
        Color of the axes lines, by default 'black'.
    facecolor : str, optional
        Background color of the figure, by default 'white'.
    point_size : int, optional
        Size of the data points in the plot, by default 100.
    cmap : str, optional
        Colormap for scatter plot, by default 'viridis'.
    fontsize_legend : int, optional
        Font size of the legend, by default 24.
    alpha_trans : float, optional
        Transparency level for data points, by default 1.
    plot_title : str, optional
        Title of the plot, by default "Trajectory of the disease progression".

    Returns
    -------
    None
        Visualizes and saves the trajectory plot.
    """
    
    EMD=adata.uns['EMD']/adata.uns['EMD'].max()
    df=adata.uns['annot']
    
    path='Results_PILOT/plots'
    
    if not os.path.exists(path):
        os.makedirs(path)
        
        
    custom_cycler = (cycler(color=colors))
    
    plot_titel = plot_titel

   
    mydmap = diffusion_map.DiffusionMap.from_sklearn(n_evecs = n_evecs, epsilon =epsilon, alpha = alpha, k=knn)
    embedding = mydmap.fit_transform(EMD)
        
    plt.rcParams.update({'font.size': font_size})
    plt.rcParams["axes.edgecolor"] = axes_color
    plt.rcParams["axes.linewidth"] = axes_line_width

    df = df.drop_duplicates(subset =[df.columns[sample_col]])
    fig = plt.figure(figsize=figsize)
    
    ax = plt.gca()
    
    ax.set(facecolor = facecolor)
    ax.set_prop_cycle(custom_cycler)
   
    for category in df[clusters].unique():
        aux = np.array(df[clusters] == category)
        group_data = embedding[aux]
        scatter = ax.scatter(group_data[:,0],group_data[:,1], alpha=alpha_trans,label=category, cmap=cmap,s=point_size) 
        
        if label_act == True:
            names = np.array(df[df[clusters] == category].sampleID)
            k = 0
            for txt in names:
                ax.annotate(txt, (group_data[k,0], group_data[k,1]),fontsize=font_size)
                k=k+1

    ax.legend(loc=location_labels, fontsize=fontsize_legend)
    plt.title(plot_titel)
    plt.savefig(path+"/"+plot_titel+'.pdf')
    plt.show()         
    plt.close(fig)
    
    
    adata.uns['embedding']=embedding
    

    
def heatmaps(adata,figsize=(12,12),col_cluster=True,row_cluster=True,cmap='Blues_r',font_scale=2):
    
    """
    Plot heatmaps of cost matrix and Wasserstein distances.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing annotations, cost matrix, and Wasserstein distances.
    
    figsize : tuple, optional
        Figure size for  heatmap, by default (12,12).
    col_cluster : bool, optional
        Whether to cluster the columns, by default True.
    row_cluster : bool, optional
        Whether to cluster the rows, by default True.
    cmap : str, optional
        Colormap for the heatmaps, by default 'Blues_r'.
    font_scale : int, optional
        Font scale for labels and annotations, by default 2.

    Returns
    -------
    None
        Plots and saves heatmaps of the cost matrix and Wasserstein distances.
    """

    annot=adata.uns['annot']
    cost=adata.uns['cost']
    path='Results_PILOT/plots'
    
    if not os.path.exists(path):
        os.makedirs(path)
        
    fig = plt.figure()
   # sns.set(font_scale=font_scale)
    sns.clustermap(cost[annot.cell_type.unique()],cmap=cmap,figsize=figsize,col_cluster=col_cluster,row_cluster=row_cluster);
    plt.title('Cost Matrix',loc='center')
    plt.savefig(path+'/Cost_matrix.pdf') 
    plt.close(fig)
    
    fig = plt.figure()
    #sns.set(font_scale=font_scale)
    emd=adata.uns['EMD_df']
    sns.clustermap(emd,cmap=cmap,figsize=figsize,col_cluster=col_cluster,row_cluster=row_cluster)
    plt.title('Wasserstein distance',loc='center')
    plt.savefig(path+'/Wasserstein distance.pdf') 
    plt.close(fig)
 

 
def heatmaps_df(df, figsize=(12, 12), col_cluster=True, row_cluster=True, cmap='Blues_r'):
    """
    Plot heatmaps of cost matrix and Wasserstein distances.

    Parameters
    ----------
    df : DataFrame
        Input DataFrame containing data for heatmap plotting.
    figsize : tuple, optional
        Figure size for the heatmap, by default (12, 12).
    col_cluster : bool, optional
        Whether to cluster the columns, by default True.
    row_cluster : bool, optional
        Whether to cluster the rows, by default True.
    cmap : str, optional
        Colormap for the heatmaps, by default 'Blues_r'.

    Returns
    -------
    None
        Plots and saves heatmaps based on the input DataFrame.
    """

    path='Results_PILOT/plots'
    sns.clustermap(df,
                                row_cluster=row_cluster,col_cluster=col_cluster,annot=False,cmap=cmap,figsize=figsize,xticklabels=True);
    plt.savefig(path+"/"+'Proportions_of_cell_types_for_samples_over_trajectory.pdf')
    
      
    

def fit_pricipla_graph(adata,NumNodes=20,source_node=0,show_text=True,Do_PCA=False,figsize=(12,12),X_color='r', Node_color='k', DimToPlot=[0, 1],facecolor='white',title='Principal graph'):
    
    """
    Fit an Elastic Principal Graph to the data and extract pseudotime information.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing embeddings and other necessary data.
    NumNodes : int, optional
        Number of nodes for building the backbone, by default 20.
    source_node : int, optional
        Index of the source node to start pseudotime estimation, by default 0.
    show_text : bool, optional
        Whether to show numbers in the backbone plot, by default True.
    Do_PCA : bool, optional
        Whether to perform PCA on the nodes, by default False.
    figsize : tuple, optional
        Figsize, by default (12,12).
    X_color : str, optional
        Color of the X-axis in the plot, by default 'r'.
    Node_color : str, optional
        Color of the backbone's nodes, by default 'k'.
    DimToPlot : list, optional
        List of integers specifying the PCs or dimensions to plot, by default [0, 1].
    facecolor : str, optional
        Background color of the figure, by default 'white'.
    title : str, optional
        Title of the plot, by default 'Principal graph'.

    Returns
    -------
    None
        Fits an Elastic Principal Graph, plots it, and extracts pseudotime information.
    """
    
    path='Results_PILOT/plots'
    
    if not os.path.exists(path):
        os.makedirs(path)
        
    emb=adata.uns['embedding']
    pg_curve = elpigraph.computeElasticPrincipalTree(emb,NumNodes=NumNodes)[0]
    fig = plt.figure(figsize=figsize)
    ax = plt.gca()
    ax.set(facecolor = facecolor)
    elpigraph.plot.PlotPG(emb,pg_curve,Do_PCA=Do_PCA,show_text=show_text,DimToPlot=DimToPlot,Node_color=Node_color,X_color=X_color)
    plt.title(title)
    plt.savefig(path+"/"+'Principal graph'+'.pdf')
    plt.show()         
    plt.close(fig)
    elpigraph.utils.getPseudotime(emb,pg_curve,source=source_node,target=None)
    pseudotime = pg_curve['pseudotime']
    adata.uns['pseudotime']=pseudotime
    


 
def clustering_emd(adata,res=0.3,metric='cosine',groupby_col='Leiden',swap_axes=False,cmap="Blues_r",dendrogram=True,show_gene_labels=True,var_group_rotation=45,figsize=[12,12],save=False,sorter_leiden=None):
    
    """
    Perform clustering and visualization of EMD (Earth Mover's Distance) data in AnnData object.

    Parameters:
    adata (AnnData): Input AnnData object containing EMD data.
    res (float): Resolution parameter for Leiden clustering. Default is 0.3.
    metric (str): Distance metric for clustering. Default is 'cosine'.
    groupby_col (str): Grouping variable for plotting. 'Leiden' groups by predicted clusters, 'status' groups by real labels. Default is 'Leiden'.
    swap_axes (bool): Swap the axes in the heatmap. Default is False.
    cmap (str): Colormap for the heatmap. Default is "Blues_r".
    dendrogram (bool): Display dendrograms in the heatmap. Default is True.
    show_gene_labels (bool): Show gene labels in the heatmap. Default is True.
    var_group_rotation (int): Rotation angle for gene labels. Default is 45 degrees.
    figsize (list): Size of the heatmap figure. Default is [12, 12].
    save (bool): Save the heatmap figure. Default is False.
    sorter_leiden (list or None): Custom order for Leiden clusters. If not provided, the default order is used.

    Returns:
    proportion_df (DataFrame): DataFrame containing proportions of sub-clusters in each sample.
    """
    
    EMD=adata.uns['EMD']
    proportions=adata.uns['proportions']
    annot=adata.uns['annot']
    real_labels=adata.uns['real_labels']
    
    EMD_df=pd.DataFrame(EMD,columns=proportions.keys())
    EMD_df['sampleID']=proportions.keys()
    EMD_df['status']=list(real_labels)
    
    
    adata_emd = sc.AnnData(EMD)
    sc.pp.neighbors(adata_emd, metric=metric)
    sc.tl.leiden(adata_emd, resolution = res)
    predicted_labels = np.array(adata_emd.obs.leiden)
    Silhouette = Sil_computing(EMD/EMD.max(), predicted_labels,metric=metric)
    
    proportion_df=pd.DataFrame(proportions)
    proportion_df=proportion_df.T
    proportion_df.columns=annot.cell_type.unique()
    
    proportion_df['Predicted_Labels']=predicted_labels
    proportion_df['sampleID']=list(proportions.keys())
    
    EMD_df['Leiden']=predicted_labels
    if groupby_col=='status':
        sorter=np.unique(real_labels)
        EMD_df['status'] = EMD_df.status.astype("category")
        EMD_df['status'] = EMD_df['status'].cat.set_categories(sorter)
        EMD_df=EMD_df.sort_values(["status"])
    elif groupby_col=='Leiden':
        if sorter_leiden==None:
            sorter=EMD_df.Leiden.unique()
        else:
            sorter=sorter_leiden
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
    

    
    


def Sil_computing(EMD, real_labels, metric='cosine'):
    """
    Compute the Silhouette score based on Wasserstein distances.

    Parameters
    ----------
    EMD : numpy.ndarray
        Wasserstein distances matrix.
    real_labels : list or numpy.ndarray
        True labels or ground truth.
    metric : str, optional
        Metric for calculating pairwise distances, by default 'cosine'.

    Returns
    -------
    float
        Silhouette score indicating cluster quality.
    """
    Silhouette = metrics.silhouette_score(EMD, real_labels, metric =metric)
    #print("Silhouette score: ", Silhouette) 
    return Silhouette


def select_best_sil(adata,resolutions=[],marker='o',figsize=(10,10),facecolor="white",metric='cosine',path=None,start=0.2,step=0.1,end=2):
    """
    Parameters
    ----------
    adata : adata,
    EMD : W distance,
    path_to_results:str, path to save the plot
    resolutions: list,a list of your desire resulotions
    marker : str, optional (default='o')
            Marker style for data points in the plot.
    figsize : tuple, optional
        Figure size (width, height) in inches. Default is (10, 10).
    facecolor : str, optional
        Background color of the saved plot image. Default is 'white'.

    metric: str, metric for leiden clustering and calculating sil. 

    Returns
    -------
    None,
    plot Silhouette Score vs. Resolution to figure out the best sil.
    plot Silhouette Score vs. Number of clusters.
    """
    
    resolutions = [start + step * i for i in range(int((end - start) / step) + 1)]
    # Create a list to store the Silhouette Scores
    best_res=start
    best_sil=0
    if path==None:
        if not os.path.exists('Results_PILOT/plots'):
            os.makedirs('Results_PILOT/plots')
        path_to_results='Results_PILOT/plots'
    else:
        path_to_results=path
    sil_scores = []
    number_of_clusters=[]
    EMD=adata.uns['EMD']
    # Define the range of resolutions/number of clusters to test
    # Loop through the resolutions and calculate the Silhouette Score for each
    for resolution in resolutions:
        # Cluster the data using Louvain clustering at the specified resolution
        adata_emd = sc.AnnData(EMD)
        sc.pp.neighbors(adata_emd, metric=metric)
        sc.tl.leiden(adata_emd, resolution = resolution)
        # Calculate the Silhouette Score
        predicted_labels = np.array(adata_emd.obs.leiden)
        
        Silhouette = metrics.silhouette_score(EMD, predicted_labels, metric =metric)
        if float(Silhouette) > best_sil:
            best_sil=Silhouette
            best_res=resolution
            

        # Append the Silhouette Score to the list
        sil_scores.append(Silhouette)
        number_of_clusters.append(len(np.unique(predicted_labels)))
    adata.uns['best_res']=best_res
    # Plot the Silhouette Scores against the resolutions
    plt.figure(figsize = figsize,facecolor=facecolor)
    plt.plot(resolutions, sil_scores, marker=marker)
    plt.xlabel('Resolution')
    plt.ylabel('Silhouette Score')
    plt.title('Silhouette Score vs. Resolution')
    plt.savefig(path_to_results+'/silhouette_score_vs_resolution.pdf')
    plt.show()
    
    
    # Plot the Silhouette Scores against the resolutions
    plt.figure(figsize = figsize,facecolor=facecolor)
    plt.plot(resolutions, number_of_clusters, marker=marker)
    plt.xlabel('Resolution')
    plt.ylabel('Number of Clusters')
    plt.title('Number of Clusters vs. Resolution')
    plt.savefig(path_to_results+'/number_of_clusters_vs_resolution.pdf')
    plt.show()
    
def cell_type_diff_two_sub_patient_groups(proportions: pd.DataFrame = None,
                                          cell_types: list = None,
                                          labels:str = 'Predicted_Labels',
                                          group1: str = 'Tumor 1',
                                          group2: str = 'Tumor 2',
                                          pval_thr: float = 0.05,
                                          figsize: tuple = (15, 4),
                                          file_path: str = None,fontsize:int = 24):
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
    fontsize:int, the size of the font
    Returns
    -------
    None.
    Plot the statistical scores showing the how significant each cell type 
    differentiated in each patients sub-group (filter based on p-value threshold)
    Save table of statisical information of all cell types,
    sorted by their p-value and statistic score

    """
    if file_path==None:
        if not os.path.exists('Results_PILOT/plots'):
            os.makedirs('Results_PILOT/plots')
        file_path='Results_PILOT/plots'
    else:
        file_path=file_path
    
        
     
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
    
    
    if not os.path.exists('Results_PILOT/Diff_Expressions_Results'):
                os.makedirs('Results_PILOT/Diff_Expressions_Results')
    stats_bar_group12.to_csv('Results_PILOT/Diff_Expressions_Results' + "/Cell_type_diff_" + group1 + "_vs_" + group2 +".csv",
                             header = True, index = None)
    
    # filter data based on a p-value threshold
    stats_bar_group12 = stats_bar_group12[stats_bar_group12['adjPval'] < pval_thr]
    
    
    fig, ax = plt.subplots(figsize = figsize)
    plot_hor_vs_vert(stats_bar_group12, 1, x = 'score', y = 'cell_type', c = 'type',
                     xlabel = 'statistic score', ylabel = None,
                     rotation = None, tick_bottom = True, tick_left = False,
                     title = "Cell type rank " + group1 + " vs " + group2,fontsize=fontsize)
    fig.tight_layout()
    plt.savefig(file_path + "/Cell_type_diff_" + group1 + "_vs_" + group2 +".pdf",
                          facecolor = 'white')
    
    
    
def plot_cell_types_distributions(proportions: pd.DataFrame = None,
                                  cell_types: list = None,
                                  labels:str = 'Predicted_Labels',
                                  file_path: str = None,
                                  figsize: tuple = (15, 7),label_order=None,label_colors=None,fontsize=24,rotation = 45):
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
    label_order: list, order of lables
    label_colors: list, color of lables based
    fontsize:int, size of fonts
    rotation:int, rotation
    Returns
    -------
    None.
    Save of the plot cell types distribution in each patients sub-groups

    """
    
    if file_path==None:
        if not os.path.exists('Results_PILOT/plots'):
            os.makedirs('Results_PILOT/plots')
        file_path='Results_PILOT/plots'
    else:
        file_path=file_path
    
    col_names = list(cell_types)
    col_names.append(labels)
    df_long = pd.melt(proportions.loc[:,col_names], labels)
    
    if label_order:
        df_long['Predicted_Labels'] = df_long.Predicted_Labels.astype("category")
        sorter=label_order
        df_long['Predicted_Labels'] = df_long['Predicted_Labels'].cat.set_categories(sorter)
        df_long=df_long.sort_values(["Predicted_Labels"])

    
    fig = plt.figure(figsize = figsize)
    sns.boxplot(x = "variable", hue = labels, y = "value", data = df_long,palette=label_colors)
    plt.ylabel("Cell proportions", fontsize = fontsize)
    plt.xlabel("")
    plt.xticks(fontsize = fontsize, rotation = rotation, ha = 'right', rotation_mode = 'anchor')
    plt.legend(fontsize = fontsize)
    fig.tight_layout()
    plt.savefig(file_path + "/Cell_types_distributions.pdf",
                          facecolor = 'white')
    plt.show()
 


    
def gene_annotation_cell_type_subgroup(cell_type: str = None,
                                   group: str = None,
                                   source: str = None,
                                   num_gos: int = 15,
                                   figsize=(12,12),
                                   font_size: int = 24,
                                   bbox_inches: str = 'tight',
                                   facecolor: str = 'white',
                                   transparent: bool = False,
                                   organism: str = 'hsapiens',
                                   dpi: int = 100,
                                   s: int = 200,
                                   color: str = 'tab:blue'):
    """
    Perform Gene Ontology (GO) enrichment analysis and create a scatterplot of enriched terms.

    Parameters:
    ----------
    cell_type : str, optional
        Specify cell type name to check its differential expression genes. The default is None.
    group : str, optional
        Name of patients sub-group of interest. The default is None.
    source : str, optional
        Specify the source of GO terms. The default is None.
    num_gos: int, optional
        Number of GO terms to plot. Default is 5.
    figsize: tuple, optional
        figsize. Default is (12,12).
    font_size: int, optional
        Font size for labels. Default is 24.
    bbox_inches: str, optional
        Bounding box for saving the plot. Default is 'tight'.
    facecolor: str, optional
        Background color of the figure. Default is 'white'.
    transparent: bool, optional
        Set to True for a transparent figure. Default is False.
    organism: str, optional
        The organism for GO analysis. Default is 'hsapiens'.
    dpi: int, optional
        Dots per inch for the saved plot image. Default is 100.
    s: int, optional
        Marker size for the scatterplot. Default is 200.
    color: str, optional
        Color of the scatterplot markers. Default is 'tab:blue'.

    Returns:
    --------
    None
        Saves the scatterplot of enriched GO terms as a PDF file.
    """

    path_to_results='Results_PILOT'
    group_genes = pd.read_csv(path_to_results +'/Diff_Expressions_Results/'+cell_type+"/significant_genes_"+cell_type+"_"+group+".csv",
                               index_col=0)
    
    
    
    gp = GProfiler(return_dataframe=True)
    if list(group_genes['0'].values):
        gprofiler_results = gp.profile(organism = organism,
                                       query = list(group_genes['0'].values))
    else:
        return "Genes list is empty!"
    
    
    if(gprofiler_results.shape[0] == 0):
        return "Not enough information!"

    
    if(gprofiler_results.shape[0] < num_gos):
        num_gos = gprofiler_results.shape[0]
    
    if source: 
        gprofiler_results = gprofiler_results[gprofiler_results['source']==source]
       
    
  
    selected_gps = gprofiler_results.head(num_gos)[['name', 'p_value']]
    
    selected_gps['nlog10'] = -np.log10(selected_gps['p_value'].values)

    plt.figure(figsize = figsize, dpi = dpi)
    plt.style.use('default')
    sns.scatterplot(data = selected_gps, x= "nlog10", y= "name", s = s, color = color)

    plt.title('GO enrichment in ' + cell_type + ' associated with ' + group, fontsize = font_size)

    plt.xticks(size = font_size)
    plt.yticks(size = font_size)

    plt.ylabel("GO Terms", size = font_size)
    plt.xlabel("-$log_{10}$ (P-value)", size = font_size)
    
    if not os.path.exists(path_to_results+'/Diff_Expressions_Results/'+cell_type+'/GO_analysis/'):
            os.makedirs(path_to_results+'/Diff_Expressions_Results/'+cell_type+'/GO_analysis/')
    path_to_results=path_to_results+'/Diff_Expressions_Results/'+cell_type+'/GO_analysis/'
    plt.savefig(path_to_results+ group + ".pdf", bbox_inches = bbox_inches, facecolor=facecolor, transparent=transparent)
    gprofiler_results.to_csv(path_to_results+group+'_'+cell_type+"_all_gprofiler_results.csv")



def exploring_specific_genes(cluster_name='cell_type',font_size=24,gene_list=[],fig_size=(64, 56),p_value=0.01,create_new_plot_folder=True,fc_ther=0.5):
    """
    Explore specific genes within a cluster to analyze their patterns in comparison to other cell types.

    Parameters:
        cluster_name : str
            The name of the cluster you're interested in exploring.
        sort : list
            List of column names for sorting the results.
        number_genes : int
            Number of top genes to display for each pattern (linear, quadratic, etc.).
        cluster_names : list
            List of cluster names for exploration (if empty, all clusters will be considered).
        font_size : int
            Font size for text and labels in plots.
        gene_list : list
            List of specific genes to explore within the cluster.
        fig_size: tuple,optional
              Size of the plot.

        fc_ther: float, optional
            threshold for FC.

    Returns:
        Show the genes for the interested cell types
    """
    path='Results_PILOT/'
    file_name = "/Whole_expressions.csv"
    cluster_names = [os.path.splitext(f)[0] for f in listdir(path + '/cells/') \
                         if isfile(join(path + '/cells/', f))]
    
    all_stats_extend = pd.read_csv(path + "/gene_clusters_stats_extend.csv", sep = ",")
    
    with open(path + '/genes_dictionary.pkl', 'rb') as handle:
        gene_dict = pickle.load(handle)
    
    pline = np.linspace(1, 20, 20)
    filtered_all_stats_extend=all_stats_extend[all_stats_extend['gene'].isin(gene_list)]
    filtered_all_stats_extend=filtered_all_stats_extend[filtered_all_stats_extend['cluster'].isin([cluster_name])]
    
 
    plot_stats_by_pattern(cluster_names, filtered_all_stats_extend, gene_dict, pline, path, file_name,font_size=font_size,p_value=p_value,create_new_plot_folder=create_new_plot_folder,fc_ther=fc_ther)
    
    # Load the PNG image file
    if create_new_plot_folder:
        image_path =path+'/plot_genes_for_'+str(cluster_name)+'/'+str(cluster_name) + ".png"  # Replace with the actual path to your PNG image
        image = imread(image_path)
    else:
        
        image_path =path+'plots_gene_cluster_differentiation/'+cluster_name+'.png'  # Replace with the actual path to your PNG image
        image = imread(image_path)
    
    # Set the size of the figure
    fig, ax = plt.subplots(figsize=fig_size)  # Adjust the width and height as needed

    # Display the PNG image
    ax.imshow(image)
    ax.axis('off')  # Turn off axis labels and ticks
    plt.show()
    
   
    
    
def go_enrichment(df,num_gos=20,source=None,cell_type='cell_type',fontsize=32,s=200, figsize = (15,12),color = 'tab:blue',dpi=100,bbox_inches = 'tight', facecolor='white', transparent=False,organism='hsapiens'):
    
    """
    Perform Gene Ontology (GO) enrichment analysis and create a scatterplot of enriched terms.

    Parameters:
    -----------
    df : DataFrame
        Input DataFrame containing gene information.
    num_gos : int, optional
        Number of top enriched GO terms to visualize. Default is 20.
    source : str, optional
        Specify the source of GO terms. The default is None.
    cell_type : str, optional
        Name of the cell type for which the GO enrichment is performed. Default is 'cell_type'.
    fontsize : int, optional
        Font size for the plot title and labels. Default is 32.
    s : int, optional
        Marker size for the scatterplot. Default is 200.
    figsize : tuple, optional
        Figure size (width, height) in inches. Default is (15, 12).
    color : str, optional
        Color of the scatterplot markers. Default is 'tab:blue'.
    dpi : int, optional
        Dots per inch for the saved plot image. Default is 100.
    bbox_inches : str, optional
        Bounding box options for saving the plot. Default is 'tight'.
    facecolor : str, optional
        Background color of the saved plot image. Default is 'white'.
    transparent : bool, optional
        Whether the saved plot image has a transparent background. Default is False.
    organism : str, optional
        Organism for GO enrichment analysis. Default is 'hsapiens'.

    Returns:
    --------
    None
        Saves the scatterplot of enriched GO terms as a PDF file.
    """

    path='Results_PILOT/'
    df_sorted =df
    gp = GProfiler(return_dataframe=True)
    gprofiler_results = gp.profile(organism = organism,
                query = list(df_sorted['gene'].values))
    
    if not os.path.exists(path+'GO/'+cell_type+'/'):
        os.makedirs(path+'GO/'+cell_type+'/')
    gprofiler_results.to_csv(path+'GO/'+cell_type+'/'+cell_type+"_all_gprofiler_results.csv")
    if(gprofiler_results.shape[0] < num_gos):
        num_gos = gprofiler_results.shape[0]

    #selected_gps = gprofiler_results.loc[0:num_gos,['name', 'p_value']]
    if source:   
        gprofiler_results = gprofiler_results[gprofiler_results['source']==source]
    selected_gps = gprofiler_results.head(num_gos)[['name', 'p_value']]

    selected_gps['nlog10'] = -np.log10(selected_gps['p_value'].values)

    figsize = figsize

    plt.figure(figsize = figsize, dpi = dpi)
    plt.style.use('default')
    sns.scatterplot(data = selected_gps, x="nlog10", y="name", s = s, color =color)

    plt.title('GO enrichment in ' + cell_type, fontsize = fontsize)

    plt.xticks(size = fontsize)
    plt.yticks(size = fontsize)

    plt.ylabel("GO Terms", size = fontsize)
    plt.xlabel("-$log_{10}$ (P-value)", size = fontsize)
    #if not os.path.exists(path+'GO/'):
       # os.makedirs(path+'GO/')
    #plt.savefig(path+'GO/'+cell_type+".pdf", bbox_inches = 'tight', facecolor='white', transparent=False)
    if not os.path.exists(path+'GO/'+cell_type+'/'):
        os.makedirs(path+'GO/'+cell_type+'/')
    plt.savefig(path+'GO/'+cell_type+'/'+cell_type+".pdf", bbox_inches = bbox_inches, facecolor=facecolor, transparent=transparent)
    
def plt_gene_cluster_differentiation(cellnames=['healthy_CM','Myofib'],font_size=22,p_value=0.01,fc_ther=0.5):
    """
    Generate and save plots showcasing gene expression patterns for selected cell clusters.

    Parameters:
    -----------
    cellnames : list, optional
        List of cluster names for which gene expression patterns will be visualized. Default is ['healthy_CM', 'Myofib'].
    font_size : int, optional
        Font size for labels and titles in the plots. Default is 16.
    p_value : float, optional
        P-value threshold for significance. Default is 0.05.
    fc_ther: float, optional
            threshold for FC.
    Returns:
    --------
    None
        Generates and saves scatterplots of gene expression patterns for selected clusters.
    """
    
    file_name = "/Whole_expressions.csv"
    path='Results_PILOT/'
    all_stats_extend=pd.read_csv(path+'gene_clusters_stats_extend.csv')
    with open(path + '/genes_dictionary.pkl', 'rb') as handle:
        gene_dict = pickle.load(handle)
    
    pline = np.linspace(1, 20, 20)
    plot_stats_by_pattern(cluster_names=cellnames, all_stats_extend=all_stats_extend, gene_dict=gene_dict, 
                          pline=pline, path_to_results=path, file_name=file_name,font_size=font_size,
                          p_value=p_value,fc_ther=fc_ther)

def qq_plot_gene(target, data, sorted_best, gene_name):
    """
    Generate a QQ plot for a specific gene's performance.

    Parameters:
        target (pd.DataFrame): The target DataFrame containing gene labels.
        data (pd.DataFrame): The data DataFrame containing gene data.
        sorted_best (dict): A dictionary containing sorted gene results.
        gene_name (str): The name of the gene to create the QQ plot for.

    Returns:
        None: Displays the QQ plot for the specified gene's performance.
    """

    x = np.array(list(data['label']))

    best_tf = list(target.loc[[gene_name]].values.ravel())

    best_results = sorted_best[gene_name][1]
    best_func = sorted_best[gene_name][0]
    # pattern = sorted_best[gene_name][2]
    data = pd.DataFrame({'x': x, 'y':best_tf})

    X = generate_feature(best_func, x)
    new_X = np.append(np.ones((len(X),1)), X, axis=1)
    # print(X.shape)
    # print(best_results['params'])
    predictions = np.matmul(new_X, best_results['params'])
    
    plt.figure(figsize = (8,6))
    stats.probplot( (best_tf - predictions), dist="norm", plot=pylab)
    pylab.show()

def plot_best_matches_cell_types(target, data,df,sorted_best, scale_name, min_target=0, max_target=35,num=11,width=25,height=25,xlim=4,point_size=100,color_back=None,fontsize=28,alpha=1,cmap='viridis'):
    """
    Plot the best-fitted models for cell types.

    Parameters:
        target (pd.DataFrame): The target data for gene activity.
        data (pd.DataFrame): The data containing cell type labels.
        df (pd.DataFrame): The data frame containing sample information.
        sorted_best (dict): A dictionary containing the best-fitted model results for each factor, sorted by R-squared or modified R-squared.
        scale_name (str): The name of the scale.
        min_target (int): The minimum target value.
        max_target (int): The maximum target value.
        num (int): Number of models to plot.
        width (int): Width of the figure.
        height (int): Height of the figure.
        xlim (int): X-axis limit.
        point_size (int): Size of data points.
        color_back (str): Background color of the plot.
        fontsize (int): Font size for the titles and labels.
        alpha (float): Alpha value for data points.
        cmap (str): Colormap for data points.

    Returns:
        None: This function generates the plot but does not return any value.
    """

    
    sorted_best = dict(sorted(sorted_best.items(), key=lambda x: x[1][-1]))

    x = np.array(list(data['label']))
    min_x = min(x)
    max_x = max(x)
    start=df[df['Time_score']==min_x]['sampleID'].unique()
    end=df[df['Time_score']==max_x]['sampleID'].unique()
    xlabel='From  '+ start+'  to  '+end
    
    plt.figure(figsize=((width, height)))
    plt.subplots_adjust(wspace = 0.5, hspace = 1 )
   # plt.suptitle(xlabel, fontsize=20, y=0.95)
    plt.tight_layout()
    
    j = 1
    for i in range(num):
        
      
        
        
        best_tf_name = list(sorted_best)[i]
        best_tf = target.loc[[best_tf_name]].values
        best_results = list(sorted_best.items())[i][1][1]
        best_func = list(sorted_best.items())[i][1][0]
        p_vals_best_model = sorted_best[best_tf_name][1]['pvalues']
        polyline = np.linspace(min_target,max_target)
        polyline = generate_feature(best_func, polyline)
        
        ax = plt.subplot(math.ceil(num/4), 4, j)
        plt.xticks(np.arange(min_x,max_x,xlim))
        if color_back!=None:
            ax.set_facecolor(color_back)
        ax.scatter(x, best_tf, c =best_tf ,cmap=cmap,alpha=alpha,s=point_size)
        new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
        curve = np.matmul(new_polyline, best_results['params'])
        ax.plot(np.linspace(min_target,max_target), curve)
        
        formatted_number= "{:.3e}".format(sorted_best[best_tf_name][-1])
        
        line1=str(best_tf_name)+"\n" +"P-val: "+ "\n"
        bold_word=formatted_number                     

        if float(sorted_best[best_tf_name][-1]) <= 0.05:
            title =line1+"{}".format(bold_word)
            font_props = FontProperties(weight='bold')
            ax.set_title(title, fontproperties=font_props,fontsize=fontsize)
        else:
            combined_line = "{} {}".format(line1,bold_word)
            ax.set_title(combined_line,fontsize=fontsize)
            
            
      #  ax.set_title(str(best_tf_name)+ "\n" +"P-val: "+ "\n" +str("{:.3e}".format(sorted_best[best_tf_name][-1])),fontsize=fontsize)
            
       # best_tf_name=str(best_tf_name)+"\n Modified adj R2 = " + str("{:.2f}".format(best_results['mod_rsquared_adj']))
       # ax.set_title(best_tf_name,fontsize=fontsize)
       
        ax.axis(xmin = min_target,xmax = max_target)
    
        if((j+4) % 4 == 1):
            #ax.set(ylabel = scale_name)
            ax.set_ylabel(scale_name, fontsize=fontsize)
            
       # if show_R:
        #    if x_lab==None:
         #       x_lab=np.mean(x)
          #  if y_lab==None:
           #     y_lab=np.mean(best_tf)
            #ax.annotate("adj R2 = " + str("{:.2f}".format(best_results['rsquared_adj'])), 
             #        (x_lab,y_lab),
              #       color=color_label,
               #      size=fontsize)
     
    
        j += 1



  
def plot_best_matches(target, data,df, sorted_best, scale_name, plot_color='tab:orange',num=16,width=25,height=25,x_lim=4,fontsize=24,alpha=0.5,cmap='viridis',color_back=None):
    """
    Plot the best-fitted models for different patterns.

    Parameters:
        target (pd.DataFrame): The target data for gene activity.
        data (pd.DataFrame): The data containing cell type labels.
        df (pd.DataFrame): The data frame containing sample information.
        sorted_best (dict): A dictionary containing the best-fitted model results for each factor, sorted by R-squared or modified R-squared.
        scale_name (str): The name of the scale.
        plot_color (str): Color of the plots.
        num (int): Number of models to plot.
        width (int): Width of the figure.
        height (int): Height of the figure.
        x_lim (int): X-axis limit.
        fontsize (int): Font size for titles and labels.
        alpha (float): Alpha value for data points.
        cmap (str): Colormap for data points.
        color_back (str): Background color of the plot.

    Returns:
        None: This function generates the plot but does not return any value.
    """

    
    plt.rcParams.update({'font.size': fontsize})

    
    all_motifs = np.array(list(sorted_best.keys()))
    all_motifs = all_motifs.reshape(all_motifs.shape[0],1)
    all_patterns = np.array([i[2] for i in list(sorted_best.values())])
    all_patterns = all_patterns.reshape(all_patterns.shape[0],1)
    mat = np.append(all_motifs, all_patterns, axis = 1)
    patterns=np.unique(mat[:,1])

    x = np.array(list(data['label']))
    min_x = min(x)
    max_x = max(x)
    start=df[df['Time_score']==min_x]['sampleID'].unique()
    end=df[df['Time_score']==max_x]['sampleID'].unique()

    plt_count=0
    plt.figure(figsize=((width, height)))
    xlabel='From  '+ start+'  to  '+end
   # plt.suptitle(xlabel, fontsize=20, y=0.05)
    #plt.tight_layout()
    plt.subplots_adjust(


                    wspace=0.5,
                    hspace=1)

    genes_expressions=patterns
    pltt=0
    for pattern in genes_expressions:
        counter=0
        flag=False
        j=0

        if (pltt+4)%4==3:
            pltt=pltt+1
        elif (pltt+4)%4==2:
            pltt=pltt+2
        elif (pltt+4)%4==1:
            pltt=pltt+3

        for key,values in sorted_best.items():

            if list(sorted_best.values())[counter][2]==pattern:
                pltt+=1
                ax = plt.subplot(len(patterns), 4, pltt)
                if color_back!=None:
                    ax.set_facecolor(color_back)
                plt.xticks(np.arange(min_x,max_x,x_lim))
                #flag=True
                best_tf_name = list(sorted_best)[counter]
                best_tf = target.loc[[best_tf_name]].values
                best_results = list(sorted_best.items())[counter][1][1]
                best_func = list(sorted_best.items())[counter][1][0]

                data = pd.DataFrame({'x': x.tolist(), 'y':best_tf.T.tolist()})

                pline = np.linspace(min_x, max_x)
                polyline = generate_feature(best_func, pline)
                new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
                curve = np.matmul(new_polyline, best_results['params'])


               
                    
                ax.scatter(x, best_tf, c = best_tf,cmap=cmap, alpha = alpha)
                 
                
                    
                formatted_number= "{:.3e}".format(sorted_best[best_tf_name][-1])

                ax.plot(pline, curve, c = plot_color)
                
                line1=str(best_tf_name)+"\n" +"P-val: "+ "\n"
                bold_word=formatted_number 
                
          
                
                
                
                
                if float(sorted_best[best_tf_name][-1]) <= 0.05:
                    
                    title =line1+"{}".format(bold_word)
                    font_props = FontProperties(weight='bold')
                    ax.set_title(title, fontproperties=font_props,fontsize=fontsize)
                    
                   
                else:
                    combined_line = "{} {}".format(line1,bold_word)
                    ax.set_title(combined_line,fontsize=fontsize)
                    
                
                
               # best_tf_name=best_tf_name+"\n Modified adj R2 = " + str("{:.2f}".format(best_results['mod_rsquared_adj']))
               # ax.set_title(best_tf_name,fontsize=fontsize)

               

                ax.axis(xmin = min(x)-0.01,xmax = max(x)+0.01)

                if((j+4) % 4 == 0):
                    ax.set_ylabel(pattern, fontsize=fontsize)
                    #ax[plt_count,j].set(ylabel = pattern,fontsize=8)
                    #ax[plt_count,j].set_xlabel(str(xlabel),fontsize=12)
            #    if show_R:
             #       if x_lab==None:
              #          x_lab=np.mean(x)
               #     if y_lab==None:
                #        y_lab=np.mean(best_tf)
                       
                 #   ax.annotate("Modified adj R2 = " + str("{:.2f}".format(best_results['mod_rsquared_adj'])),
                              
                  #           (x_lab,y_lab),
                   #          color=color_label,
                    #         size=fontsize)

                j += 1
                if j%4==0:
                    break




            counter=counter+1  

    
  

def plot_two_genes(adata, sorted_best_WT, sorted_best_KO, gene_name, scale_name, plot_color1 = 'tab:blue', plot_color2 = 'tab:red'):
    
    """
    Plot the gene activity for two different conditions (e.g., WT and KO) along with their best-fitted models.

    Parameters:
        adata (AnnData): An AnnData object containing the data.
        sorted_best_WT (dict): A dictionary containing the best-fitted model results for the wild-type (WT) condition.
        sorted_best_KO (dict): A dictionary containing the best-fitted model results for the knock-out (KO) condition.
        gene_name (str): The name of the gene for which the activity is plotted.
        scale_name (str): The name of the scale.
        plot_color1 (str): Color for the WT condition plot (default is 'tab:blue').
        plot_color2 (str): Color for the KO condition plot (default is 'tab:red').

    Returns:
        None: This function generates the plot but does not return any value.
    """

    WT_ATAC_data = pd.DataFrame()
    WT_ATAC_target = adata[adata.obs.group == 'WT'].to_df().transpose()
    WT_ATAC_data['label'] = list(adata[adata.obs.group == 'WT'].obs['label'].values)
    
    WT_x = WT_ATAC_data['label']
    WT_tf = list(WT_ATAC_target.loc[[gene_name]].values.ravel())
    WT_results = sorted_best_WT[gene_name][1]
    WT_func = sorted_best_WT[gene_name][0]
    # WT_pattern = sorted_best_WT[gene_name][2]
    
    KO_ATAC_data = pd.DataFrame()
    KO_ATAC_target = adata[adata.obs.group == 'KO'].to_df().transpose()
    KO_ATAC_data['label'] = list(adata[adata.obs.group == 'KO'].obs['label'].values)
    
    KO_x = KO_ATAC_data['label']
    KO_tf = list(KO_ATAC_target.loc[[gene_name]].values.ravel())
    KO_results = sorted_best_KO[gene_name][1]
    KO_func = sorted_best_KO[gene_name][0]
    # KO_pattern = sorted_best_KO[gene_name][2]
    
    pline = np.linspace(min(KO_x), max(KO_x))
    
    polyline = generate_feature(WT_func, pline)
    new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
    WT_curve = np.matmul(new_polyline, WT_results['params'])
    
    polyline = generate_feature(KO_func, pline)
    new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
    KO_curve = np.matmul(new_polyline, KO_results['params'])
    
    ax = plt.subplot(1, 1, 1)
    ax.scatter(WT_x, WT_tf, c = WT_tf, alpha = 0.5, cmap='Blues', norm=colors.CenteredNorm(-1),)
    ax.scatter(KO_x, KO_tf, c = KO_tf, alpha = 0.5, cmap='Reds', norm=colors.CenteredNorm(-1),)
    
    
    
    ax.axis(xmin=min(WT_x),xmax=max(WT_x))

    ax.plot(np.linspace(min(WT_x),max(WT_x)), WT_curve, c = plot_color1)
    ax.plot(np.linspace(min(WT_x),max(WT_x)), KO_curve, c = plot_color2)
    ax.set_title(re.sub('.*?:', '', gene_name))
    ax.set_ylabel(scale_name)
    

def plot_one_gene(target, data, sorted_best, gene_name, scale_name, plot_color):
    
    """
    Plot the gene activity and its best-fitted model for a single gene.

    Parameters:
        target (pd.DataFrame): A Pandas DataFrame containing the target gene expression data.
        data (pd.DataFrame): A Pandas DataFrame containing the data.
        sorted_best (dict): A dictionary containing the best-fitted model results for multiple genes.
        gene_name (str): The name of the gene for which the activity is plotted.
        scale_name (str): The name of the scale.
        plot_color (str): Color for the plot.

    Returns:
        None: This function generates the plot but does not return any value.
    """

    x = data['label']
    min_x = min(x)
    max_x = max(x)
    
    best_tf = list(target.loc[[gene_name]].values.ravel())

    best_results = sorted_best[gene_name][1]
    best_func = sorted_best[gene_name][0]
    pattern = sorted_best[gene_name][2]
    # data = pd.DataFrame({'x': x, 'y':best_tf})


    
    pline = np.linspace(min_x, max_x)
    polyline = generate_feature(best_func, pline)
    new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
    curve = np.matmul(new_polyline, best_results['params'])

    ax = plt.subplot(1, 1, 1)
    ax.scatter(x, best_tf, c = best_tf, alpha = 0.5, cmap='plasma')
    ax.axis(xmin=min(x),xmax=max(x))

    ax.plot(np.linspace(min(x),max(x)), curve, c = plot_color)
    ax.set_title(re.sub('.*?:', '', gene_name))
    ax.set_ylabel(scale_name)
    
    
    ax.annotate("adj R2 = " + str("{:.2f}".format(best_results['rsquared_adj'])), 
                  (0.2,np.max(best_tf)-0.1),
                  color='black',
                  size=14)
    ax.set_xlabel(pattern)

def plot_gene(target, data, sorted_best, gene_name, scale_name, plot_color):
    """
    Plot gene expression data along with the best-fitted curve and statistical information.

    Parameters:
        target (pd.DataFrame): A Pandas DataFrame containing gene expression data.
        data (pd.DataFrame): A Pandas DataFrame containing additional data such as labels.
        sorted_best (dict): A dictionary containing the best-fitted models for different genes.
        gene_name (str): The name of the gene to plot.
        scale_name (str): The name of the scale for the y-axis.
        plot_color (str): The color to use for plotting.

    Returns:
        None

    Note:
        This function plots gene expression data for a specific gene along with the best-fitted curve,
        statistical information, and other visual elements using the ggplot library.
        It includes the gene name, adjusted R-squared value, p-value, and the fitted function type in the plot.
        The appearance of the plot may vary based on the fitted function type (linear, quadratic, or linear-quadratic).
    """

    x = data['label']
    min_x = min(x)
    max_x = max(x)
    
    best_tf = list(target.loc[[gene_name]].values.ravel())

    best_results = sorted_best[gene_name][1]
    best_func = sorted_best[gene_name][0]
    pattern = sorted_best[gene_name][2]
    data = pd.DataFrame({'x': x, 'y':best_tf})

    pline = np.linspace(min_x, max_x)
    polyline = generate_feature(best_func, pline)
    new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
    curve = np.matmul(new_polyline, best_results['params'])
    
    p = ggplot(aes(x='x', y='y'), data,)
    p = p + geom_pointdensity(size = 3, show_legend=False)
    
    dt = pd.DataFrame({'x': pline, 'y': curve})
    p = p + geom_line(aes(x='x', y='y'), dt, colour = plot_color, size = 1)
    
    p = p + theme(panel_background=element_rect(fill='white'),              
          axis_line_x=element_line(color='black'),
          axis_line_y=element_line(color='black'),
          panel_grid=element_blank(),
          panel_border=element_blank())
    p = p + labs(title = gene_name + "\n\n" + \
                 "adj $R^{2}$ = " + str("{:.2f}".format(best_results['rsquared_adj'])) \
                            + "\n p-value = " + str("{:.2f}".format(sorted_best[gene_name][4]))
                 , y = scale_name, x = pattern)
    # p = p + annotate(geom = "text", 
    #                   label = "adj R2 = " + str("{:.2f}".format(best_results['rsquared_adj'])) \
    #                         + "\n p-value = " + str("{:.2f}".format(sorted_best[gene_name][4])),
    #                   x = 0, 
    #                   y = np.max(best_tf))
    if(best_func == 'linear'):
        p = p + annotate(geom = "text", 
                          label = "ge = " +
                              str(round(best_results['params'][0],2)) \
                              + ["", "+"][best_results['params'][1] > 0]+ \
                                  str(round(best_results['params'][1],2)) + "t" ,
                          x = pline[10], 
                          y = curve[10])
    elif(best_func == 'quadratic'):
        p = p + annotate(geom = "text", 
                          label = "ge = " +
                              str(round(best_results['params'][0],2)) \
                               + ["", "+"][best_results['params'][1] > 0]+ \
                                   str(round(best_results['params'][1],2)) + "$t^2$" ,
                          x = pline[10], 
                          y = curve[10])
    elif(best_func == 'linear_quadratic'):
        p = p + annotate(geom = "text", 
                          label = "ge = " + 
                              str(round(best_results['params'][0],2)) \
                               + ["", "+"][best_results['params'][1] > 0]+ \
                                   str(round(best_results['params'][1],2)) + "t" \
                                   + ["", "+"][best_results['params'][2] > 0]+\
                                       str(round(best_results['params'][2],2)) + "$t^2$" ,
                          x = pline[10], 
                          y = curve[10])
    p.draw()
    

def plot_gene_specific(target, data, sorted_best, gene_name, scale_name, plot_color):
    
    """
    Plot gene expression and its best-fitted model for a single gene.

    Parameters:
        target (pd.DataFrame): A Pandas DataFrame containing the target gene expression data.
        data (pd.DataFrame): A Pandas DataFrame containing the data.
        sorted_best (dict): A dictionary containing the best-fitted model results for multiple genes.
        gene_name (str): The name of the gene for which the expression is plotted.
        scale_name (str): The name of the scale.
        plot_color (str): Color for the plot.

    Returns:
        None: This function generates the plot but does not return any value.
    """

    x = data['label']
    min_x = min(x)
    max_x = max(x)
    
    best_tf = list(target.loc[[gene_name]].values.ravel())
    best_results = sorted_best[gene_name][1]
    best_func = sorted_best[gene_name][0]
    pattern = sorted_best[gene_name][2]
    
    new_x = generate_feature(best_func, np.array(list(data['label'])))
    all_main_polyline = np.append(np.ones((len(new_x),1)), new_x, axis=1)
    new_curve = np.matmul(all_main_polyline, best_results['params'])
    
    my_data = pd.DataFrame({'x': x, 'y':best_tf})
    pline = np.linspace(min_x, max_x)
    polyline = generate_feature(best_func, pline)
    new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
    curve = np.matmul(new_polyline, best_results['params'])
    
    
    # stdev = np.sqrt(sum((best_tf - target.loc[gene_name])**2) / (len(target.loc[gene_name]) - 2))
    
    
    print("mean value: " + str(np.mean(target.loc[gene_name])))
    ss_res = np.sum( (target.loc[gene_name] - new_curve )**2)
    print("SSres: " + str(ss_res) )
    
    delata = 1.35
    modified_ss_res = 0
    res_e = target.loc[gene_name] - new_curve
    for e in res_e:
        if abs(e) < delata:
            modified_ss_res += 1/2*np.sum(e**2)
        else:
            modified_ss_res += delata *(abs(e) - (1/2)*delata)
    print("Modified SSres: " + str(modified_ss_res))
    
    ss_tot = np.sum( (target.loc[gene_name] - np.mean(target.loc[gene_name]))**2)
    print("SStot: " + str(ss_tot) )
    r2 = 1 - ( ss_res / ss_tot )
    print("R2: " + str(r2) )
    
    r2 = 1 - ( modified_ss_res / ss_tot )
    r2 = 1 - (1-r2)*(len(best_tf)-1)/(len(best_tf)-new_x.shape[1]-1)
    print("Modified R2: " + str(r2) )
    
    mse = mean_absolute_error(new_curve, target.loc[gene_name])
    print("MAE: " + str(mse))
    
    p = ggplot(aes(x='x', y='y'), my_data,)
    p = p + geom_pointdensity(size = 3, show_legend=False)
    
    dt = pd.DataFrame({'x': pline, 'y': curve})
    p = p + geom_line(aes(x='x', y='y'), dt, colour=plot_color, size = 1)
    
    mean_dt = pd.DataFrame({'x': pline, 'y': np.mean(target.loc[gene_name])*np.ones(len(pline))})
    p = p + geom_line(aes(x='x', y='y'), mean_dt, colour='tab:red', size = 1)
    
    p = p + theme(panel_background=element_rect(fill='white'),              
          axis_line_x=element_line(color='black'),
          axis_line_y=element_line(color='black'),
          panel_grid=element_blank(),
          panel_border=element_blank())
    p = p + labs(title = gene_name, y = scale_name, x = pattern)
    p = p + annotate(geom = "text", 
                      label = "\n\n adj $R^2$ = " + str("{:.2f}".format(best_results['rsquared_adj'])) \
                            + "\n Modified adj $R^2$ = " + str("{:.2f}".format(r2))\
                            + "\n p-value = " + str("{:.2f}".format(sorted_best[gene_name][4])),
                      x = np.mean(pline), 
                      y = np.max(best_tf)+0.15)
    if(best_func == 'linear'):
        p = p + annotate(geom = "text", 
                          label = "ge = " + str(round(best_results['params'][0],2)) \
                              + ["", "+"][best_results['params'][1] > 0]+ \
                                  str(round(best_results['params'][1],2)) + "t" ,
                          color = "tab:orange",
                          size = 16,
                          x = np.max(pline)-15, 
                          y = np.max(curve)+0.01)
    elif(best_func == 'quadratic'):
        p = p + annotate(geom = "text", 
                          label = "ge = " + str(round(best_results['params'][0],2)) \
                               + ["", "+"][best_results['params'][1] > 0]+ \
                                   str(round(best_results['params'][1],2)) + "$t^2$" ,
                          x = pline[10], 
                          y = curve[10])
    elif(best_func == 'linear_quadratic'):
        p = p + annotate(geom = "text", 
                          label = "ge = " + str(round(best_results['params'][0],2)) \
                               + ["", "+"][best_results['params'][1] > 0]+ \
                                   str(round(best_results['params'][1],2)) + "t" \
                                   + ["", "+"][best_results['params'][2] > 0]+\
                                       str(round(best_results['params'][2],2)) + "$t^2$" ,
                          x = pline[10], 
                          y = curve[10])
    p = p + annotate(geom = "text", 
                    label = "Mean Model" ,
                        x = np.max(mean_dt['x'])-8, 
                        y = np.max(mean_dt['y'])-0.03,
                    color = 'tab:red',
                     size = 16
                        )
    p.draw()
   
  
def plot_gene_distribtion(target, gene_name):
    """
    Plot the distribution of gene expression for a specific gene.

    Parameters:
        target (pd.DataFrame): A Pandas DataFrame containing the target gene expression data.
        gene_name (str): The name of the gene for which the expression distribution is plotted.

    Returns:
        None: This function generates the plot but does not return any value.
    """

    data = pd.DataFrame({'gene': target.loc[gene_name]})
    plt.figure(figsize = (8,6))
    sns.histplot(data=data, x='gene', edgecolor='k')
    plt.title(gene_name + " expression distribtion")
    plt.xlabel('gene expression', size = 16)
    plt.ylabel('Frequency', size = 16)
    plt.yticks(size = 12, weight = 'bold')
    plt.xticks(size = 12, weight = 'bold')
    plt.show()
   
def plot_gene_density(target, data, sorted_best, gene_name, scale_name, plot_color):
    """
    Plot the density distribution of gene expression for a specific gene.

    Parameters:
        target (pd.DataFrame): A Pandas DataFrame containing the target gene expression data.
        data (pd.DataFrame): A Pandas DataFrame containing the data with 'label' and 'y' columns.
        sorted_best (dict): A dictionary containing sorted best-fitted models for genes.
        gene_name (str): The name of the gene for which the density distribution is plotted.
        scale_name (str): The name of the scale or measurement associated with the gene expression.
        plot_color (str): The color for the density plot.

    Returns:
        None: This function generates the density plot but does not return any value.
    """

    x = data['label']

    best_tf = list(target.loc[[gene_name]].values.ravel())
    data = pd.DataFrame({'x': x, 'y':best_tf})
    
    labels=np.sort([y for y in list(data.x.unique())])
    fig, axes = joypy.joyplot(data, by="x", column="y", labels=labels, 
                          grid="y", linewidth=1, legend=False, overlap=0.5, figsize=(6,5),kind="kde", bins=80,
                          title=gene_name, ylabels=False,
                          colormap=cm.autumn_r)

def plot_pval_rsq_correlation(table, feature1, feature2, show_fit = True, log_transform = False):
    
    """
    Plot the correlation between two features, optionally showing a linear fit.

    Parameters:
        table (pd.DataFrame): A Pandas DataFrame containing the data.
        feature1 (str): The name of the first feature to be plotted on the x-axis.
        feature2 (str): The name of the second feature to be plotted on the y-axis.
        show_fit (bool, optional): If True, a linear fit line will be plotted. Default is True.
        log_transform (bool, optional): If True, the second feature will be log-transformed (-log10) before plotting. Default is False.

    Returns:
        None: This function generates the scatter plot but does not return any value.
    """

    
    x = [float(x) for x in table[feature1].values]
    
    if(log_transform):
        y = [float(x) for x in table[feature2].values]
        y = -np.log10(y)
        data = {"x": x, "y": y}
    
    plt.figure(figsize=(6,6))
    sns.scatterplot(x, y, color='tab:blue')
    
    if(show_fit):
        results = ols("y ~ x", data=data).fit()
        print('ratio: ' + str(results.params[1]))
        pline = np.linspace(min(x),max(x))
        curve = np.matmul(np.array(results.params), [np.ones(len(pline)),pline])
        # y_pred = np.matmul(np.array(results.params), [np.ones(len(x)),x])
    
        sns.lineplot(pline, curve, color = 'tab:red')
    
        # textstr = '\n'.join((
        #         r'RMSE=%.2f' % (compute_metrics(y, y_pred)[1], ),
        #         r'MAE=%.2f' % (compute_metrics(y, y_pred)[2], ),
        #         r'r2=%.2f' % (compute_metrics(y, y_pred)[0], ),
        #         r'Cor=%.2f' % (compute_metrics(y, y_pred)[3], )))
        # plt.text(min(pline), max(y), textstr, fontsize=16,
        #         verticalalignment='top')
        
    plt.xlabel(feature1, fontsize=16)
    if(log_transform):
        plt.ylabel("$-log_{10}($" + str(feature2) + ")", fontsize=16)
    else:
        plt.ylabel(feature2, fontsize=16)
    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.show()

def plot_condition(target, data, sorted_best, condition_type, scale_name, plot_color):
    """
    Plot the best-fitted models based on up or downregulation patterns in gene expression.

    Parameters:
        target (pd.DataFrame): A Pandas DataFrame containing the target gene expressions.
        data (pd.DataFrame): A Pandas DataFrame containing the feature data.
        sorted_best (dict): A dictionary containing the best-fitted models for each gene.
        condition_type (str): The type of condition to plot, either 'increasing' or 'decreasing'.
        scale_name (str): The name of the scale for the y-axis.
        plot_color (str): The color to use for the plot.

    Returns:
        None: This function generates the plot but does not return any value.
    """

    assert condition_type in ['increasing', 'decreasing'], 'condition type parameter must be increasing, decreasing'
    
    x = np.array(list(data['label']))
    j = 1
    i = 0
    
    if(condition_type == 'increasing'):
        plt.figure(figsize=(15, 12))
        plt.subplots_adjust(hspace=0.5, wspace=0.3)
        plt.suptitle("Increasing Expression", fontsize=18, y=0.95)
        plt.tight_layout()
   
        while(j<=9):
            
            best_tf_name = list(sorted_best)[i]
            best_tf = target.loc[[best_tf_name]].values
            best_results = list(sorted_best.items())[i][1][1]
            best_func = list(sorted_best.items())[i][1][0]
            
            pline = np.linspace(min(x), max(x))
            polyline = generate_feature(best_func, pline)
            new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
            curve = np.matmul(new_polyline, best_results['params'])
            mit = list(sorted_best.items())[i][1][3]
            if(mit > 0):
                
                
                ax = plt.subplot(3, 3, j)

                ax.scatter(x, best_tf, c = best_tf, alpha = 0.5, cmap='RdBu', norm=colors.CenteredNorm(), )
                ax.axis(xmin = min(x),xmax = max(x))
                ax.set(ylabel = scale_name)
                ax.plot(pline, curve, c = plot_color)
                ax.set_title(best_tf_name)
                
                ax.annotate("adj R2 = " + str("{:.2f}".format(best_results['rsquared_adj'])), 
                      (0,np.max(best_tf)-0.2),
                      color='black',
                      size=14)
        
                j += 1
            i += 1
    else:
        plt.figure(figsize=(15, 12))
        plt.subplots_adjust(hspace=0.5, wspace=0.3)
        plt.suptitle("Decreasing Expression", fontsize=18, y=0.95)
        plt.tight_layout()
        
        while(j<=9):
            
            best_tf_name = list(sorted_best)[i]
            best_tf = target.loc[[best_tf_name]].values
            best_results = list(sorted_best.items())[i][1][1]
            best_func = list(sorted_best.items())[i][1][0]
            
            pline = np.linspace(min(x), max(x))
            polyline = generate_feature(best_func, pline)
            new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
            curve = np.matmul(new_polyline, best_results['params'])
            mit = list(sorted_best.items())[i][1][3]
            if(mit < 0):
                
                
                ax = plt.subplot(3, 3, j)
                ax.scatter(x, best_tf, c = best_tf, alpha = 0.5, cmap='RdBu', norm=colors.CenteredNorm(), )
                ax.axis(xmin = min(x),xmax = max(x))
                ax.set(ylabel = scale_name)
                ax.plot(pline, curve, c = plot_color)
                ax.set_title(best_tf_name)
                
                ax.annotate("adj R2 = " + str("{:.2f}".format(best_results['rsquared_adj'])), 
                      (1.8,np.max(best_tf)-0.1),
                      color='black',
                      size=14)
        
                j += 1
            i += 1
            
def plot_6_best(target, data, sorted_best, scale_name, plot_color):
    """
    Plot the 6 best-fitted models for each pattern type (linear, quadratic, linear_quadratic) and regulation direction (up, down) in gene expression.

    Parameters:
        target (pd.DataFrame): A Pandas DataFrame containing the target gene expressions.
        data (pd.DataFrame): A Pandas DataFrame containing the feature data.
        sorted_best (dict): A dictionary containing the best-fitted models for each gene.
        scale_name (str): The name of the scale for the y-axis.
        plot_color (str): The color to use for the plot.

    Returns:
        None: This function generates the plot but does not return any value.
    """

    x = np.array(list(data['label']))
    
    plt.figure(figsize=(15, 12))
    plt.subplots_adjust(hspace=0.5)
    
    pattern_types = ['linear', 'quadratic', 'linear_quadratic']
    mit_types = ['up', 'down']
    
    k = 1
    for mit in mit_types:
        for pattern in pattern_types:
            motifs = data_interst(sorted_best, pattern, mit)
            tf_name = list(motifs.keys())[0]
            tf = target.loc[[tf_name]].values
            results = list(motifs.items())[0][1][1]
            func = list(motifs.items())[0][1][0]
            patt = list(motifs.items())[0][1][2]
            pline = np.linspace(min(x), max(x))
            polyline = generate_feature(func, pline)
            new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
            curve = np.matmul(new_polyline, results['params'])
            
            ax = plt.subplot(3, 3, k)
            ax.scatter(x, tf, c = tf, alpha = 0.5, cmap='RdBu', norm=colors.CenteredNorm(), )
            ax.axis(xmin=min(x),xmax=max(x))

            ax.plot(np.linspace(min(x),max(x)), curve, c = plot_color)
            ax.set_title(re.sub('.*?:', '', tf_name))
            if(k % 3 == 1):
                ax.set(ylabel = scale_name)
            
            ax.annotate("adj R2 = " + str("{:.2f}".format(results['rsquared_adj'])), 
                          (0.2,np.max(tf)-0.1),
                          color='black',
                          size=14)
            ax.set_xlabel(patt)
            k += 1
            
    
def plot_hor_vs_vert(data, subplot, x, y, c, xlabel, ylabel, rotation,
                     tick_bottom, tick_left, title,fontsize=24):
    
    
    
    '''
    Plot horizontal and vertical bar charts using Seaborn.

    Parameters:
        data : pandas.DataFrame
            The input DataFrame containing the data to be plotted.
        subplot : int
            The subplot number (1 or 2) for the placement of the chart.
        x : str
            The name of the column to be plotted on the x-axis.
        y : str
            The name of the column to be plotted on the y-axis.
        c : str
            The color column used to determine the color of bars.
        xlabel : str
            The label for the x-axis.
        ylabel : str
            The label for the y-axis.
        rotation : float
            The rotation angle for x-axis labels in degrees.
        tick_bottom : bool
            If True, display ticks on the bottom of the plot.
        tick_left : bool
            If True, display ticks on the left side of the plot.
        title : str
            The title of the plot.
        fontsize : int, optional
            The font size for title, labels, and ticks (default is 24).

    Returns:
        None
    '''
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

def volcano_plot(scores, foldchanges, p_values, cell_type, feature1, feature2, fc_thr = 1, pv_thr = 1,
                 figsize = (20,20), output_path = None,n_p=5,n_n=5,font_size=18, marker='o',
                             color='w',
                             markersize=8,
                             font_weight_legend='normal',
                             size_legend=12,dpi=100
                 
                 
                             ):
   
    """
    Generate a volcano plot to visualize gene expression significance.

    Parameters:
        scores : pandas Series
            A pandas Series containing the expression scores for genes.
        foldchanges : array-like
            An array-like containing the fold changes for genes.
        p_values : pandas Series
            A pandas Series containing the p-values for genes.
        cell_type : str
            The name of the cell type being analyzed.
        feature1 : str
            The name of the first feature being compared.
        feature2 : str
            The name of the second feature being compared.
        fc_thr : float, optional (default=1)
            The threshold for log2FoldChange to determine significance.
        pv_thr : float, optional (default=1)
            The threshold for negative log10 of p-value to determine significance.
        figsize : tuple, optional (default=(15, 15))
            The size of the plot figure.
        output_path : str, optional (default=None)
            The path to save the output plot. If None, the plot will be displayed.
        n_p : int, optional (default=5)
            The number of labels that the user wants to show over the plot for positive threshold.
        n_n : int, optional (default=5)
            The number of labels that the user wants to show over the plot for negative threshold.
        font_size : int, optional (default=18)
            Font size for the plot.
        marker : str, optional (default='o')
            Marker style for data points in the plot.
        color : str, optional (default='w')
            Color for data points in the plot.
        markersize : int, optional (default=8)
            Marker size for data points in the plot.
        font_weight_legend : str, optional (default='normal')
            Font weight for legend text.
        size_legend : int, optional (default=12)
            Font size for legend text.
        dpi : int, optional
            Dots per inch for the saved plot image. Default is 100.

    Returns:
        None
    """

    
    
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
    
    
    #plt.xlim(-xlim, xlim)
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
    filtered_df = df.loc[df['nlog10'] >= pv_thr]
    subset_labels_fold_change_pos = filtered_df.loc[filtered_df['log2FoldChange'] >= fc_thr]
    subset_labels_fold_change_pos = subset_labels_fold_change_pos.sort_values(by='nlog10', ascending=False)
    subset_labels_fold_change_pos = subset_labels_fold_change_pos.head(n_p)['symbol'].values

    subset_labels_fold_change_neg = filtered_df.loc[filtered_df['log2FoldChange'] <= -fc_thr]
    subset_labels_fold_change_neg = subset_labels_fold_change_neg.sort_values(by='nlog10', ascending=False)
    subset_labels_fold_change_neg = subset_labels_fold_change_neg.head(n_n)['symbol'].values
    # Combine the subsets of genes
    subset_labels = np.concatenate([subset_labels_fold_change_pos, subset_labels_fold_change_neg])
    for i in range(len(df)):
        if df.iloc[i].symbol in subset_labels:
            if df.iloc[i].nlog10 >= pv_thr and (df.iloc[i].log2FoldChange >= fc_thr):
                texts.append(plt.text(x = df.iloc[i].log2FoldChange, y = df.iloc[i].nlog10, s = df.iloc[i].symbol,
                                     fontsize = font_size, weight = 'bold', family = 'sans-serif'))
            if df.iloc[i].nlog10 >= pv_thr and ( df.iloc[i].log2FoldChange <= -fc_thr):
                texts.append(plt.text(x = df.iloc[i].log2FoldChange, y = df.iloc[i].nlog10, s = df.iloc[i].symbol,
                                     fontsize = font_size, weight = 'bold', family = 'sans-serif'))
    adjust_text(texts)
   # for i in range(len(df)):
    #    if df.iloc[i].symbol in subset_labels:
     #       if df.iloc[i].nlog10 >= pv_thr and (df.iloc[i].log2FoldChange >= fc_thr):
      #          texts.append(plt.text(x = df.iloc[i].log2FoldChange, y = df.iloc[i].nlog10, s = df.iloc[i].symbol,
          #                           fontsize = 16, weight = 'bold', family = 'sans-serif'))
       #     if df.iloc[i].nlog10 >= pv_thr and ( df.iloc[i].log2FoldChange <= -fc_thr):
        #        texts.append(plt.text(x = df.iloc[i].log2FoldChange, y = df.iloc[i].nlog10, s = df.iloc[i].symbol,
         #                            fontsize = 16, weight = 'bold', family = 'sans-serif'))
    #adjust_text(texts)

    custom_lines = [Line2D([0], [0], marker=marker, color=color, markerfacecolor='#d62a2b', markersize=markersize),
                   Line2D([0], [0], marker=marker, color=color, markerfacecolor='#1f77b4', markersize=markersize)]

    plt.legend(custom_lines, ['Higher expressions in ' + feature2, 'Higher expressions in ' + feature1],loc = 1,
               bbox_to_anchor = (1,1.1), frameon = False, prop = {'weight': font_weight_legend, 'size': size_legend})

    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(2)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.tick_params(width = 2)

    plt.title("Expression Score \n "+feature1+" - "+feature2, fontsize = font_size)
    plt.xticks(size = font_size, weight = 'bold')
    plt.yticks(size = font_size, weight = 'bold')
    plt.xlabel("$log_{2}$ (Fold Change)", size = font_size)
    plt.ylabel("-$log_{10}$ (P-value)", size = font_size)

#     plt.savefig(filename, dpi = 100, bbox_inches = 'tight', facecolor = 'white')
    plt.savefig(output_path + "/volcano_" + str(feature1) + "-" + str(feature2) + "_FC.pdf",
                dpi = dpi, bbox_inches = 'tight', facecolor = 'white')
    plt.show()

 

                
def map_color(a, fc_thrr, pv_thrr):
    """
    Map colors based on specified thresholds for Fold Change and p-value.

    Parameters:
        a : tuple
            A tuple containing log2FoldChange, symbol, and negative log10 of p-value.
        fc_thrr : float
            The threshold for log2FoldChange to determine different color mappings.
        pv_thrr : float
            The threshold for negative log10 of p-value to determine different color mappings.

    Returns:
        str
            A string indicating the color mapping based on the provided thresholds and input values.
    """
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

def plot_stats_by_pattern(cluster_names: list = None,
                          all_stats_extend: pd.DataFrame = None,
                          gene_dict: dict = None,
                          pline: list = None,
                          path_to_results: str = None,
                          file_name: str = "/Whole_expressions.csv",
                          font_size: int = 24,p_value=0.01,create_new_plot_folder=False,fc_ther=0.5):
    """
    Parameters
    ----------
    cluster_names : list, optional
        List of cluster names to be considered.
    all_stats_extend : DataFrame, optional
        DataFrame containing stats for each gene and the cluster its belong.
    gene_dict : dict, optional
        Dictionary map gene to cluster its belong.
    pline : list, optional
        Values for axis x to plot curve. 
    path_to_results : str, optional
        path to save result contained Markers and cells.
    file_name : str, optional
        File name same for all clusters have fitting information for genes.
        The default is "/Whole_expressions.csv".
    font_size : int, optional
        Specify font size for labels and names. The default is 12.
    fc_ther: float, optional
        Threshold for FC.

    Returns
    -------
    return figures for each cluster ploting top 4 genes for each pattern.

    """
    plt.rcParams.update({'font.size': font_size})
    # plot results
    for cluster in cluster_names:
        my_data = all_stats_extend[ (all_stats_extend['cluster'] == str(cluster)) & (all_stats_extend['pvalue'] < p_value)]
        my_data['FC']=my_data['FC'].astype(float)
        my_data=my_data[my_data['FC'] > fc_ther]
        sort_my_data = my_data.sort_values(['pvalue'],
                ascending=[True]).groupby('Expression pattern').head(4)
        expression_patterns = np.unique(sort_my_data['Expression pattern'])
        if(len(expression_patterns)):

            n_col = 4
            n_row = len(expression_patterns)
            n_px = 10

            plt.rcParams["figure.facecolor"] = 'w'
            
                
            #plt.figure(figsize=(80, 80))
            #f, axs = plt.subplots(n_row, n_col, figsize=(16, 8))
            plt.figure(figsize=(80, 80))  # Set the overall figure size

            # Adjust the size of individual subplots
            subplot_width = 8  # Choose an appropriate value
            subplot_height = 6  # Choose an appropriate value

            f, axs = plt.subplots(n_row, n_col, figsize=(n_col * subplot_width, n_row * subplot_height))
            axs = np.atleast_2d(axs)
           
          

            p = 0
            for pattern in expression_patterns:
                pattern_stats = sort_my_data[sort_my_data['Expression pattern'] == pattern]
                k = 0
                for i, row in pattern_stats.iterrows():
                    cluster_n = row['cluster']
                    gene_name = row['gene']
                    WT_x, WT_tf, WT_curve, KO_x, KO_tf, KO_curve = get_two_tables(gene_name, cluster_n, gene_dict,
                                                                                  file_name, pline, path_to_results)

                    if(KO_x is not None):
                        axs[p, k].scatter([x+0.2 for x in KO_x], KO_tf, alpha = 0.2, c = "tab:gray")
                    axs[p, k].scatter(WT_x, WT_tf, alpha = 0.5, c = WT_tf, cmap = 'viridis',
                              norm = colors.CenteredNorm(np.mean(WT_tf)))

                    axs[p, k].axis(xmin = min(WT_x), xmax = max(WT_x))

                    if(KO_x is not None):
                        axs[p, k].plot(np.linspace(min(WT_x),max(WT_x)), KO_curve,
                                       c = "dimgray", linewidth = 4.0)
                    axs[p, k].plot(np.linspace(min(WT_x),max(WT_x)), WT_curve,
                                   c = "tab:orange", linewidth = 4.0)

                    axs[p, k].set_title(gene_name, size = font_size, weight = 'bold')
                    if( k == 0):
                        plt.rcParams.update({'font.size': font_size-4})
                        axs[p, k].set_ylabel(pattern, size = font_size-4)
                    #else:
                        
                     #   axs[p, k].set_ylabel('Gene expression', size = 16)
                    plt.rcParams.update({'font.size': font_size})
                    for item in (axs[p, k].get_xticklabels() + axs[p, k].get_yticklabels()):
                        item.set_fontsize(font_size)

                    plt.text(.01, .99, 'p-value = %.2e ' % + Decimal(str(row['pvalue'])),
                             ha = 'left', va = 'top',
                             transform=axs[p, k].transAxes, size = font_size)
                    axs[p, k].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    k += 1
                while(k != 4):
                    axs[p, k].set_axis_off()
                    k += 1
                p += 1
            plt.subplots_adjust(wspace = 0.5, hspace = 0.7)
            
            if create_new_plot_folder:
                if not os.path.exists(path_to_results+'/plot_genes_for_'+str(cluster)+'/'):  
                    os.makedirs(path_to_results+'/plot_genes_for_'+str(cluster)+'/')
            
                save_path = path_to_results+'/plot_genes_for_'+str(cluster)+'/'+str(cluster) + ".png"
                plt.savefig(save_path)
                plt.close()
            else:   
                if not os.path.exists(path_to_results+'/plots_gene_cluster_differentiation/'):  
                        os.makedirs(path_to_results+'/plots_gene_cluster_differentiation/')

                save_path = path_to_results+'/plots_gene_cluster_differentiation/'+str(cluster) + ".png"
                plt.savefig(save_path)
                print('Plot for '+str(cluster))
                plt.show()
                plt.close()


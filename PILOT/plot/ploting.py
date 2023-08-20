import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import os
import seaborn as sns
import scipy
import matplotlib.pyplot as plt
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
import warnings
import ot
from logging import info, warn
from cycler import cycler
warnings.filterwarnings('ignore')
import elpigraph
from ..tools.Trajectory import *
from ..tools.patients_sub_clustering import *
from ..tools.Gene_cluster_specific import *


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
fig_h : int, optional
    Height of the figure, by default 12.
fig_w : int, optional
    Width of the figure, by default 12.
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

def trajectory(adata,n_evecs = 2, epsilon =1, alpha = 0.5,knn= 64, sample_col=1, clusters = 'status',label_act = False,colors=['#377eb8','#ff7f00','#e41a1c'],location_labels='center', fig_h=12,fig_w=12,font_size=24,axes_line_width=1,axes_color='black',facecolor='white',point_size=100,cmap='viridis',fontsize_legend=24,alpha_trans=1,plot_titel = "Trajectory of the disease progression"):
    
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
    fig = plt.figure()
    fig.set_size_inches(fig_h, fig_w)
    
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
    
"""
Plot heatmaps of cost matrix and Wasserstein distances.

Parameters
----------
adata : AnnData
    Annotated data matrix containing annotations, cost matrix, and Wasserstein distances.
figsize_h : int, optional
    Height of the heatmap figure, by default 12.
figsize_w : int, optional
    Width of the heatmap figure, by default 12.
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

    
def heatmaps(adata,figsize_h=12,figsize_w=12,col_cluster=True,row_cluster=True,cmap='Blues_r',font_scale=2):
    annot=adata.uns['annot']
    cost=adata.uns['cost']
    path='Results_PILOT/plots'
    
    if not os.path.exists(path):
        os.makedirs(path)
        
    fig = plt.figure()
   # sns.set(font_scale=font_scale)
    sns.clustermap(cost[annot.cell_type.unique()],cmap=cmap,figsize=(figsize_h,figsize_w),col_cluster=col_cluster,row_cluster=row_cluster);
    plt.title('Cost Matrix',loc='center')
    plt.savefig(path+'/Cost_matrix.pdf') 
    plt.close(fig)
    
    fig = plt.figure()
    #sns.set(font_scale=font_scale)
    emd=adata.uns['EMD_df']
    sns.clustermap(emd,cmap='Blues_r',figsize=(figsize_h,figsize_w),col_cluster=True,row_cluster=True)
    plt.title('Wasserstein distance',loc='center')
    plt.savefig(path+'/Wasserstein distance.pdf') 
    plt.close(fig)
    
    
    
    

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
fig_x_size : int, optional
    Width of the figure, by default 12.
fig_y_size : int, optional
    Height of the figure, by default 12.
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

def fit_pricipla_graph(adata,NumNodes=20,source_node=0,show_text=True,Do_PCA=False,fig_x_size=12,fig_y_size=12,X_color='r', Node_color='k', DimToPlot=[0, 1],facecolor='white',title='Principal graph'):
    
    path='Results_PILOT/plots'
    
    if not os.path.exists(path):
        os.makedirs(path)
        
    emb=adata.uns['embedding']
    pg_curve = elpigraph.computeElasticPrincipalTree(emb,NumNodes=NumNodes)[0]
    fig = plt.figure()
    fig.set_size_inches(fig_x_size, fig_y_size)
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


    
    
def clustering_emd(adata,res=0.3,metric='cosine',groupby_col='Leiden',swap_axes=False,cmap="Blues_r",dendrogram=True,show_gene_labels=True,var_group_rotation=45,figsize=[12,12],save=False):
    
    
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


    """
    Parameters
    ----------
    EMD : W distance,
    path_to_results:str, path to save the plot
    resolutions: list,a list of your desire resulotions

    metric: str, metric for leiden clustering and calculating sil. 

    Returns
    -------
    None,
    plot Silhouette Score vs. Resolution to figure out the best sil.
    """


def select_best_sil(adata,resolutions=[],metric='cosine',path=None,start=0.2,step=0.1,end=2):
    
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
    plt.figure(facecolor="white")
    plt.plot(resolutions, sil_scores, marker='o')
    plt.xlabel('Resolution')
    plt.ylabel('Silhouette Score')
    plt.title('Silhouette Score vs. Resolution')
    plt.savefig(path_to_results+'/silhouette_score_vs_resolution.pdf')
    plt.show()
    
    
     # Plot the Silhouette Scores against the resolutions
    plt.figure(facecolor="white")
    plt.plot(resolutions, number_of_clusters, marker='o')
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
                                      source:str=None,num_gos:int=5,fig_h:int=15,fig_w:int=12,font_size:int=24):
    """
    

    Parameters
    ----------
    cell_type : str, optional
        Specify cell type name to check its differential expression genes. The default is None.
    group : str, optional
        Name of patients sub-group of interest . The default is None.
    path_to_results : str, optional
        Path to store the results. The default is None.
    num_gos: float, number of GO for ploting.
    fig_h:int, height of figure,
    fig_w:int, width of figure
    font_size:int, size of labels
    Returns
    -------
    None.
    Plot to show the most relative GO terms for specifc cell-type of determind patient sub-group
    """

    path_to_results='Results_PILOT'
    group_genes = pd.read_csv(path_to_results +'/Diff_Expressions_Results/'+cell_type+"/significant_genes_"+cell_type+"_"+group+".csv",
                               index_col=0)
    gp = GProfiler(return_dataframe=True)
    gprofiler_results = gp.profile(organism = 'hsapiens',
                query = list(group_genes['0'].values))
    num_gos = num_gos
    if(gprofiler_results.shape[0] < num_gos):
        num_gos = gprofiler_results.shape[0]
    
    if source: 
        gprofiler_results = gprofiler_results[gprofiler_results['source']==source]
       
    
    #selected_gps = gprofiler_results.loc[0:num_gos,['name', 'p_value']]
    selected_gps = gprofiler_results.head(num_gos)[['name', 'p_value']]
    
    selected_gps['nlog10'] = -np.log10(selected_gps['p_value'].values)

    figsize = (fig_h,fig_w)

    plt.figure(figsize = figsize, dpi = 100)
    plt.style.use('default')
    sns.scatterplot(data = selected_gps, x="nlog10", y="name", s = 200, color = 'tab:blue')

    plt.title('GO enrichment in ' + cell_type + ' associated with ' + group, fontsize = 32)

    plt.xticks(size = 24)
    plt.yticks(size = 24)

    plt.ylabel("GO Terms", size = 24)
    plt.xlabel("-$log_{10}$ (P-value)", size = 24)
    
    if not os.path.exists(path_to_results+'/Diff_Expressions_Results/'+cell_type+'/GO_analysis/'):
            os.makedirs(path_to_results+'/Diff_Expressions_Results/'+cell_type+'/GO_analysis/')
    path_to_results=path_to_results+'/Diff_Expressions_Results/'+cell_type+'/GO_analysis/'
    plt.savefig(path_to_results+ group + ".pdf", bbox_inches = 'tight', facecolor='white', transparent=False)
    

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

Returns:
    Show the genes for the interested cell types
"""

def exploring_specific_genes(cluster_name='cell_type',font_size=24,gene_list=[],fig_size=(64, 56),p_value=0.01,create_new_plot_folder=True):
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
    
 
    plot_stats_by_pattern(cluster_names, filtered_all_stats_extend, gene_dict, pline, path, file_name,font_size=font_size,p_value=p_value,create_new_plot_folder=create_new_plot_folder)
    
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
    
    
    
def go_enrichment(df,num_gos=20,cell_type='cell_type',fontsize=32,s=200, figsize = (15,12),color = 'tab:blue',dpi=100):
    
    
    """
    Perform Gene Ontology (GO) enrichment analysis and create a scatterplot of enriched terms.

    Parameters:
    -----------
    df : DataFrame
        Input DataFrame containing gene information.
    num_gos : int, optional
        Number of top enriched GO terms to visualize. Default is 20.
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

    Returns:
    --------
    None
        Saves the scatterplot of enriched GO terms as a PDF file.
    """
    path='Results_PILOT/'
    df_sorted =df
    gp = GProfiler(return_dataframe=True)
    gprofiler_results = gp.profile(organism = 'hsapiens',
                query = list(df_sorted['gene'].values))

    if(gprofiler_results.shape[0] < num_gos):
        num_gos = gprofiler_results.shape[0]

    #selected_gps = gprofiler_results.loc[0:num_gos,['name', 'p_value']]
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
    if not os.path.exists(path+'GO/'):
        os.makedirs(path+'GO/')
    plt.savefig(path+'GO/'+cell_type+".pdf", bbox_inches = 'tight', facecolor='white', transparent=False)
    
def plt_gene_cluster_differentiation(cellnames=['healthy_CM','Myofib'],font_size=22,p_value=0.01):
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
                          p_value=p_value)
    

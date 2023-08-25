import pandas as pd
import numpy as np
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
from .Cell_gene_selection import *
from .Gene_cluster_specific import *
import warnings
import ot
from logging import info, warn
from cycler import cycler
from matplotlib.image import imread
warnings.filterwarnings('ignore')







def wasserstein_distance(adata,emb_matrix='X_PCA',
clusters_col='cell_types',sample_col='sampleID',status='status',
                              metric='cosine',
                               regulizer=0.2,normalization=True,
                               regularized='unreg',reg=0.1,
                               res=0.01,steper=0.01,data_type='scRNA',return_sil_ari=False):
    
    """
    Calculate the Wasserstein (W) distance among samples using PCA representation and clustering information.

    Parameters
    ----------
    adata : AnnData
        Loaded AnnData object containing the data.
    emb_matrix : numpy.ndarray
        PCA representation of data (variable).
    clusters_col : str
        Column name in the observation level of 'adata' that represents cell types or clustering.
    sample_col : str
        Column name in the observation level of 'adata' that represents samples or patients.
    status : str
        Column name in the observation level of 'adata' that represents status or disease, e.g., control/case.
    regulizer : float, optional
        Hyper-parameter of a Dirichlet distribution for regularization, by default 0.1.
    metric : str, optional
        Metric for calculating the cost matrix, by default 'cosine'.
    regularized : bool, optional
        Whether to use regularized optimal transport, by default True.
    reg : float, optional
        Regularization parameter if 'regularized' is True, by default 0.1.
    res : float, optional
        Resolution for Leiden clustering to achieve desired cluster count, by default 0.1.
    steper : float, optional
        Stepper value for finding the best Leiden resolution, by default 0.01.
    data_type : str, optional
        Type of your data, e.g., 'scRNA' or 'pathomics', by default 'scRNA'.
    return_sil_ari : bool, optional
        Whether to return ARI (Adjusted Rand Index) or Silhouette score for assessing W distance effects, by default False.

    Returns
    -------
    None
        Calculates and stores the W distance among samples in the adata object.
    """
    
    
    if data_type=='scRNA':
        data,annot=extract_data_anno_scRNA_from_h5ad(adata,emb_matrix=emb_matrix,
clusters_col=clusters_col,sample_col=sample_col,status=status)
    else:
         data,annot=extract_data_anno_pathomics_from_h5ad(adata,var_names=list(adata.var_names),clusters_col=clusters_col,sample_col=sample_col,status=status)
           
        
       
        
        
    adata.uns['data'] =data
    adata.uns['annot'] =annot
    proportions = Cluster_Representations(annot,regulizer=regulizer,
    normalization=normalization)
    adata.uns['proportions'] = proportions
    
    cost,cost_df= cost_matrix(annot,data,metric=metric)
    adata.uns['cost'] = cost_df
    
    EMD,emd_df = wasserstein_d(proportions,cost/cost.max(),regularized=regularized,
    reg=reg)
    
    adata.uns['EMD_df']=emd_df
    adata.uns['EMD'] =EMD
    #Computing clusters and then ARI
    if return_sil_ari:
        predicted_labels, ARI, real_labels = Clustering(EMD/EMD.max(), annot,metric=metric,res=res,steper=steper)
        adata.uns['real_labels'] =real_labels
        #Computing Sil
        Silhouette = Sil_computing(EMD/EMD.max(), real_labels,metric=metric)
        adata.uns['Sil']=Silhouette
        adata.uns['ARI']=ARI
    else:
        adata.uns['real_labels']=return_real_labels(annot)
      
        



def load_h5ad(path):
    """
    Load h5ad data from the specified path.

    Parameters
    ----------
    path : str
        Path to the h5ad data file.

    Returns
    -------
    AnnData
        Loaded AnnData object containing the data from the specified path.
    """
    if isfile(path):
        adata=sc.read_h5ad(path)
        return adata
    else:
        print('There is no such data, check the path or name')
    
    
    



def set_path_for_results():
    """
    Create a folder named 'Results_PILOT' to save the generated results by PILOT.

    Parameters
    ----------
    dataset : str
        Name of your data.

    Returns
    -------
    Created path.
    """
    if not os.path.exists('Results_PILOT/plots'):
        os.makedirs('Results_PILOT/plots')
        return 'Results_PILOT/plots'
    else:
        
        return 'Results_PILOT/plots'

 




def extract_annot_expression(adata,columns=['cell_type_original','patient_region','region','X_pca'],reclustering=False,reduction=False,resu=0.1,max_value=10,target_sum=1e4): 
    """
    Extract annotations and expression data from a scRNA dataset.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    columns : list, optional
        List of column names to extract from observation annotations, by default ['cell_type_original', 'patient_region', 'region', 'X_pca'].
    reclustering : bool, optional
        Whether to perform reclustering based on annotations, by default False.
    reduction : bool, optional
        Whether to apply dimensionality reduction, by default False.
    resu : float, optional
        Resolution for Leiden clustering during reclustering, by default 0.1.
    max_value : int, optional
        Maximum value to clip expression data, by default 10.
    target_sum : float, optional
        Target sum for normalization, by default 1e4.

    Returns
    -------
    AnnData
        Annotated data matrix with extracted annotations and processed expression data.

    """
    
    if reduction:
                sc.pp.normalize_total(adata, target_sum=target_sum)
                sc.pp.log1p(adata)
                sc.pp.scale(adata, max_value=max_value)
                sc.tl.pca(adata, svd_solver='arpack')
                data=adata.obsm['X_pca']  
                cols=[]
                for i in range(1,data.shape[1]+1):
                        cols.append('PCA_'+str(i))
                data=pd.DataFrame(data,columns=cols)
    if not reduction:
                data=adata.obsm[columns[3]]

                cols=[]
                for i in range(1,data.shape[1]+1):
                        cols.append('PCA_'+str(i))
                data=pd.DataFrame(data,columns=cols)

    if reclustering:
                sc.pp.neighbors(adata)
                louvain = Louvain(resolution = resu)
                labels = louvain.fit_transform(adata.obsp['distances'])
                annot = pd.DataFrame()
                annot['cell_types']=labels 
                annot['sampleID']=list(adata.obs[columns[1]])
                annot['status']=list(adata.obs[columns[2]])
    if not reclustering:
                annot=adata.obs[columns[0:3]]
                annot.columns=['cell_types','sampleID','status']                          
    return data,annot    





def extract_data_anno_scRNA_from_h5ad(adata,emb_matrix='PCA',clusters_col='cell_type',sample_col='sampleID',status='status'):
    """
    This function extracts the needed inputs for PILOT(pathomics)
    Parameters
    ----------
    adata : AnnData
        Loaded AnnData object containing the data.
    emb_matrix : numpy.ndarray
        PCA representation of data (variable).
    clusters_col : str
        Column name in the observation level of 'adata' that represents cell types or clustering.
    sample_col : str
        Column name in the observation level of 'adata' that represents samples or patients.
    status : str
        Column name in the observation level of 'adata' that represents status or disease, e.g., control/case.
    Returns
    -------
    AnnData
        Annotated data matrix with extracted annotations and processed expression data.
    """
    global path_to_results
    data=adata.obsm[emb_matrix]  
    col_add=[]
    for i in range(1,adata.obsm[emb_matrix].shape[1]+1):
        col_add.append('PCA_'+str(i))
    data=pd.DataFrame(data,columns=col_add) 
    data = data.reset_index(drop=True)
    annot=adata.obs[[clusters_col,sample_col,status]]
    annot.columns=['cell_type','sampleID','status']
    annot = annot.reset_index(drop=True) 
    path_to_results=set_path_for_results()
    
    return data,annot



def extract_data_anno_pathomics_from_h5ad(adata,var_names=[],clusters_col='Cell_type',sample_col='sampleID',status='status'):
    """
    This function extracts the needed inputs for PILOT (scRNA)
    Parameters
    ----------
    adata : AnnData
        Loaded AnnData object containing the data.
    emb_matrix : numpy.ndarray
        PCA representation of data (variable).
    clusters_col : str
        Column name in the observation level of 'adata' that represents cell types or clustering.
    sample_col : str
        Column name in the observation level of 'adata' that represents samples or patients.
    status : str
        Column name in the observation level of 'adata' that represents status or disease, e.g., control/case.
    Returns
    -------
    AnnData
        Annotated data matrix with extracted annotations and processed expression data.
    """
    global path_to_results
    data=adata[:,var_names].X
    data=pd.DataFrame(data,columns=var_names)
    data = data.reset_index(drop=True)
    annot=adata.obs[[clusters_col,sample_col,status]]
    annot.columns=['cell_type','sampleID','status']
    annot = annot.reset_index(drop=True) 
    path_to_results=set_path_for_results()
    
    return data,annot

              

       
    
def load_annot(path,name):
    
    """
    Load the annotation matrix from a specified path.

    Parameters
    ----------
    path : str
    Path to the annotation matrix file.
    name : str
    Name of the annotation dataset.

    Returns
    -------
    AnnotationMatrix
    Loaded annotation matrix from the specified path.
    """
    
   
    for format in [".csv"]:
        p = os.path.join(path,name + format)
        if isfile(p):
            info("loading " + p)
            try: 
                if format == ".csv":
                    annot = pd.read_csv(p)
                annot = annot.loc[:, ~annot.columns.str.contains('^Unnamed')]  
                return annot
            except (Exception):
                warn("loading " + name + "failed, check the path/name of this data")
                
                

                     
def load_expression(path,name):
    """
    Load the gene expression matrix from a specified path.

    Parameters
    ----------
    path : str
        Path to the gene expression matrix file.
    name : str
        Name of the expression dataset.

    Returns
    -------
    ExpressionMatrix
        Loaded gene expression matrix from the specified path.
    """

    
    for format in [".csv"]:
                  
        p = os.path.join(path, name + format)
        if isfile(p):
            info("loading " + p)
            try: 
                if format == ".csv":
                    exper = pd.read_csv(p)
                
                exper = exper.loc[:, ~exper.columns.str.contains('^Unnamed')]
                return exper
            except (Exception):
                warn("loading " + name + " failed, check the path/name of this data")
        
   


              
    

def Cluster_Representations(df, cell_col = 0, sample_col = 1,regulizer=0.2,normalization=True): 
    """
    Calculate proportions of clusters per sample.

    Parameters
    ----------
    df : DataFrame
        Annotation data containing sample, cell/cluster, and status information.
    cell_col : int or str, optional
        Cell-type/cluster column index or name in the annotation data, by default 0 (the 0th column).
    sample_col : str, optional
        Sample column name in the annotation data, by default 1.
    regulizer : float, optional
        Regularizer for smoothing proportions, by default 0.2.
    normalization : bool, optional
        Whether to normalize the proportions, by default True.

    Returns
    -------
    Dictionary
        Cluster proportions per sample.
    """

    dict ={}

    cells = df[df.columns[cell_col]].unique()  #An array of cell_types
    cell_col = df.columns[cell_col]
    sample_col = df.columns[sample_col]
    prior=np.ones(len(df['cell_type'].unique()))
    for cell in cells:
        prior[cells ==cell] =len(df[df['cell_type']==cell])/(len(df)-1)

    prior=prior*(regulizer)  # prior= Nk/N*c
    
    ## Dictionary of Samples: Each sample stores the normalized counts of CellTypes.
    for sample in df[sample_col].unique():
        
        aux = df[df[sample_col]==sample]
        dict[sample] = aux[cell_col].value_counts()
        
    ## Vector: Stores the distribution over the whole cell-space --> makes zeros for non-appering CellTypes in a sample.
    for x in df[sample_col].unique():
        vector = np.zeros(len(df[cell_col].unique()))
        z = dict[x].keys()  #z includes cell_types
        for i in z:  
            vector[cells == i] = dict[x][i] #dict[x][i] number of counts
            
        
        dict[x] = vector  #set proportions for each sample
        
    
    if normalization == True:
        for x in df['sampleID'].unique():
            dict[x] = (dict[x]+prior)/(sum(dict[x])+sum(prior))
          
            
            
            
    
    return dict




def cost_matrix(annot,data,metric='cosine'):
    """
    Compute the cost matrix to find distances between clusters/cell-types.

    First, it calculates the median of each cluster, and then the dissimilarity between clusters is computed.

    Parameters
    ----------
    annot : DataFrame
        Annotation data containing cluster/cell-type information.
    data : numpy.ndarray
        PCA representations.
    path : str
        Path for saving the cost matrix.
    metric : str
        Metric for calculating the distances between centroids.

    Returns
    -------
    matrix of distances and dataframe    
    """
    cells = annot[annot.columns[0]].unique() # gets cells 
    centroids = []

    for i in cells:
        centroids.append(list(data[annot[annot.columns[0]] == i].median(axis=0))) #get the centriod of each cluster

    dis_t = scipy.spatial.distance.pdist(centroids, metric = metric) 
    dis = scipy.spatial.distance.squareform(dis_t, force ='no', checks = True)
    cost = pd.DataFrame.from_dict(dis).T
    cost.columns=annot.cell_type.unique()
    cost['cell_types']=annot.cell_type.unique()
    cost=cost.set_index('cell_types')
         
    return dis,cost



def wasserstein_d(Clu_rep, cost,regularized = "unreg", reg = 0.1):
    """
    Compute Wasserstein distances among samples.

    By default, the method is classical OT (unregularized). You can change the method to regularized OT.

    Parameters
    ----------
    Clu_rep : dict
        Dictionary containing proportions matrix for each sample.
    cost : numpy.ndarray
        Cost matrix.
    regularized : str, optional
        Type of OT method, by default "unreg" (classical OT without regularization). Use "reg" for regularized OT.
    reg : float, optional
        Regularization parameter for Sinkhorn regularization, by default 0.1.

    Returns
    -------
    numpy.ndarray, DataFrame
        Array of Wasserstein distances and a DataFrame of distances indexed by sample IDs.
    """
    
   
    samples_id = list(Clu_rep.keys())
    n_samples = len(samples_id)
    EMD = np.zeros((n_samples,n_samples))
    
    if regularized == "unreg":
        for i in range(n_samples):
            for j in range(n_samples):
              
                EMD[i,j] = ot.emd2(Clu_rep[samples_id[i]],Clu_rep[samples_id[j]], cost)
    else:
        for i in range(n_samples):
            for j in range(n_samples):
                EMD[i,j] = ot.sinkhorn2(Clu_rep[samples_id[i]],Clu_rep[samples_id[j]], cost, reg, method = "sinkhorn_stabilized")
                
   
    emd = pd.DataFrame.from_dict(EMD).T
    emd.columns=samples_id 
    emd['sampleID']=samples_id 
    emd=emd.set_index('sampleID')
    
    return EMD,emd



def Clustering(EMD, df, category = 'status', sample_col=1,res = 0.01,metric ='cosine',steper=0.01):
    """
    Perform clustering based on Wasserstein distances.

    Parameters
    ----------
    EMD : numpy.ndarray
        Wasserstein distances matrix.
    df : DataFrame
        Annotation data.
    category : str, optional
        Name of the column representing disease/status in annotation data, by default 'status'.
    sample_col : int, optional
        Index of the column representing sample IDs in annotation data, by default 1.
    res : float, optional
        Initial resolution for Leiden clustering, by default 0.01.
    metric : str, optional
        Metric for neighborhood graph construction in knn space, by default 'cosine'.
    steper : float, optional
        Stepper value for adjusting resolution to match the number of true labels, by default 0.01.

    Returns
    -------
    numpy.ndarray, float, list
        Cluster labels, Adjusted Rand Index (ARI), and true labels.
    """

    samples = df[df.columns[sample_col]].unique()
    condition = df[category]
    n_clust = len(df[category].unique())
    
    flag = True
    inc = steper
    while flag == True:
        adata = sc.AnnData(EMD)
        sc.pp.neighbors(adata, metric=metric)
        sc.tl.leiden(adata, resolution = res)
        labels = np.array(adata.obs.leiden)
        if len(df.status.unique()) > len(np.unique(labels)):

            res = res + inc

        elif len(df.status.unique()) < len(np.unique(labels)):

            flag = True
            res = res - 0.001

        else: 
            flag = False
            labels = labels.astype(int)

    
    #print("Cluster labels: ", df[category].unique())
    
    true_labels = []
    for i in range(len(samples)):
        a = condition[df[df.columns[sample_col]] == samples[i]].unique()
        true_labels.append(a[0])
    
    S = rand_score(true_labels, labels)
    #print("ARI: ", S)
    return labels, S, true_labels;




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






def plt_trajectory(adata,n_evecs = 2, epsilon =1, alpha = 0.5,knn= 64, sample_col=1, clusters = 'status',label_act = False,colors=['#377eb8','#ff7f00','#e41a1c'],location_labels='center', fig_h=12,fig_w=12,font_size=24,axes_line_width=1,axes_color='black',facecolor='white',point_size=100,cmap='viridis',fontsize_legend=24,alpha_trans=1,plot_titel = "Trajectory of the disease progression"):
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

    
    EMD=adata.uns['EMD']/adata.uns['EMD'].max()
    df=adata.uns['annot']
    path=path_to_results
    
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

    
    
    
    
      
        
def return_real_labels(df, category = 'status', sample_col=1):
    """
    Load Annotaion of the data.

    Parameters
    ----------
    df : dataframe
        the annotaion (samples/cells/status).

    Returns
    -------
    Real lables of samples

    """  

    samples = df[df.columns[sample_col]].unique()
    condition = df[category]
    n_clust = len(df[category].unique())

    true_labels = []
    for i in range(len(samples)):
        a = condition[df[df.columns[sample_col]] == samples[i]].unique()
        true_labels.append(a[0])
    
    
    return true_labels;
        


    
def plt_heatmaps(adata,figsize_h=12,figsize_w=12,col_cluster=True,row_cluster=True,cmap='Blues_r',font_scale=2):
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
    annot=adata.uns['annot']
    cost=adata.uns['cost']
    path=path_to_results
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
    
    
    
    




def fit_pricipla_graph(adata,NumNodes=20,source_node=0,show_text=True,Do_PCA=False,fig_x_size=12,fig_y_size=12,X_color='r', Node_color='k', DimToPlot=[0, 1],facecolor='white',title='Principal graph'):
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
    path=path_to_results
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
    





def cell_importance(adata,heatmap_h=20,heatmap_w=12,width=40,height=35,xlim=5,p_val=1,plot_cell=True,point_size=100,color_back='white',fontsize=20,alpha=1,cmap='viridis',save_as_pdf=True):
    """
    Order cells based on estimated time and visualize cell type importance.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing necessary data.
    heatmap_h : int, optional
        Height of the heatmap figure, by default 20.
    heatmap_w : int, optional
        Width of the heatmap figure, by default 12.
    width : int, optional
        Width of the plot, by default 40.
    height : int, optional
        Height of the plot, by default 35.
    xlim : int, optional
        Limit for x-axis in the plot, by default 5.
    p_val : float, optional
        P-value for filtering the fitting models, by default 1.
    plot_cell : bool, optional
        Whether to plot the cell type importance, by default True.
    point_size : int, optional
        Size of points in the plot, by default 100.
    color_back : str, optional
        Background color of the plot, by default 'white'.
    fontsize : int, optional
        Font size for labels and annotations, by default 20.
    alpha : float, optional
        Transparency level for plotting, by default 1.
    cmap : str, optional
        Colormap for plotting, by default 'viridis'.
    save_as_pdf : bool, optional
        Whether to save the plot as PDF, by default False.

    Returns
    -------
    None
        Visualizes and saves the cell type importance plot.
    """

    
    path=path_to_results
    real_labels=adata.uns['real_labels']
    pseudotime=adata.uns['pseudotime']
    annot=adata.uns['annot']
    bins=adata.uns['proportions']
    embedding_diff=adata.uns['embedding']
    cell_types_propo=bins
    patients_id=bins.keys()
    cell_types = annot['cell_type'].unique()
    emd=embedding_diff
    labels=real_labels
    emd_dataframe = pd.DataFrame({'sampleID':list(patients_id), 'pseudotime':pseudotime,'lables':list(labels)},dtype=object)
    emd_dataframe_sort = emd_dataframe.sort_values('pseudotime', ascending=True) 
    emd_dataframe_sort['pseudotime']=np.arange(1, len(emd_dataframe)+1, 1).tolist()
    pathies_cell_proportions = pd.DataFrame.from_dict(bins).T
    pathies_cell_proportions.columns=cell_types
    pathies_cell_proportions.index.name='sampleID'
    df_join = pd.merge(emd_dataframe_sort['sampleID'], pathies_cell_proportions, how='inner', on = 'sampleID')
    df_join=df_join.set_index('sampleID')
    #Normalizing the proportions for heat map
    normalized_df=(df_join-df_join.min())/(df_join.max()-df_join.min())
 
    #Saving Heat map based on sorte pseuduscores of the Trajectory 
    
    sns.clustermap(normalized_df[cell_types],row_cluster=False,annot=False,cmap='Blues',figsize=(heatmap_h,heatmap_w),xticklabels=True);
    plt.savefig(path+"/"+'Samples_over_trajectory.pdf')
    
    #Building a model based on Regression and pseuduscores 
    pathies_cell_proportions['Time_score']=list(emd_dataframe['pseudotime'])
    pathies_cell_proportions = pathies_cell_proportions.sort_values('Time_score', ascending=True)
    pathies_cell_proportions['Time_score']=np.arange(1, len(emd_dataframe)+1, 1).tolist()
    pathies_cell_proportions=pathies_cell_proportions.reset_index()
    RNA_data = pd.DataFrame()
    RNA_data['label'] = pathies_cell_proportions['Time_score']
    RNA_target = np.transpose(pathies_cell_proportions.iloc[:,1:len(annot.cell_type.unique())+1])
    min_target=min(RNA_data['label'])
    max_target=max(RNA_data['label'])
    #sorted_best = fit_best_model_cell_types(RNA_target, RNA_data,min_target=min_target, max_target=max_target)
    
    
    sorted_best=fit_best_model(RNA_target, RNA_data, model_type='LinearRegression',max_iter_huber=1,epsilon_huber=1,pval_thr=p_val,modify_r2=False)
    
    if save_as_pdf:
        suffix='Cell_types_importance.pdf'
    else:
        suffix='Cell_types_importance.png'
    
    if plot_cell:
    
        with plt.rc_context():
                plot_best_matches_cell_types(RNA_target, RNA_data,pathies_cell_proportions, sorted_best, "Cell Proportion",
        min_target=min_target, max_target=max_target,num=len(sorted_best.keys()),width=width,height=height,xlim=xlim,point_size=point_size,color_back=color_back,fontsize=fontsize,alpha=alpha,cmap=cmap)
                plt.savefig(path+"/"+suffix)
    
    
    cellnames=list(sorted_best.keys())
    #return orders of samples based on the trejctory and selected cell-types
    
    
    if not os.path.exists(path+'/Cell_type_Report'):
        os.makedirs(path+'/Cell_type_Report')
  
    
    save_data(sorted_best, 
                   ['Cell name', 'Expression pattern', 'Slope', 'Fitted function', 'Intercept', 'Treat', 'Treat2', 'adjusted P-value', 'R-squared','mod_rsquared_adj'],
                       path+'/Cell_type_Report/','Cells_Importance',p_val=p_val,pro=None)
    
    
    
    adata.uns['cellnames']=cellnames
    adata.uns['orders']=pathies_cell_proportions[['sampleID','Time_score']]
 



def extract_cells_from_gene_expression(adata,sample_col,col_cell,cell_list=[],normalize=True):
    """
    Extract gene expression data for specific cells and associate them with pseudotime.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing gene expression data.
    sample_col : str
        Name of the column representing sample IDs in annotation data.
    col_cell : str
        Name of the column representing cell types in annotation data.
    cell_list : list, optional
        List of cell types to extract, by default an empty list (extract all cell types).
    normalize : bool, optional
        Whether to normalize and log-transform the gene expression data, by default True.

    Returns
    -------
    None
        Extracts and saves gene expression data for specified cells associated with pseudotime.
    """

    
    path_results='Results_PILOT'
    orders=adata.uns['orders']

    if len(cell_list)==0:
        cell_types=adata.uns['cellnames']
    else:
        cell_types=cell_list


    for cell in cell_types:

        adata_new = adata[adata.obs[col_cell].isin([cell]),:]
        if normalize:
            sc.pp.normalize_total(adata_new, target_sum=1e4)
            sc.pp.log1p(adata_new)
        df=adata_new[:,adata_new.var_names].X
        df=pd.DataFrame(df.toarray())
        df.columns=adata_new.var_names
        df['sampleID']=list(adata_new.obs[sample_col])
        joint=orders.merge(df,on='sampleID')
        if not os.path.exists(path_results+'/cells/'):
            os.makedirs(path_results+'/cells/') 
        joint.to_csv(path_results+'/cells/'+cell+'.csv')

        return joint





            
def genes_importance(adata,name_cell,col='Time_score',genes_index=[],p_value=0.05,max_iter_huber=100,epsilon_huber=1.35,x_lim=4,width=20,height=30,store_data=True,genes_interesting=[],modify_r2 = False,model_type = 'HuberRegressor',fontsize=8,alpha=0.5,cmap='viridis',color_back=None,save_as_pdf=False,plot_genes=True,colnames=[],sample_col='sampleID',col_cell='cell_types',normalize=True):
    
    """
    Order genes based on estimated time and visualize gene importance.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing gene expression data.
    name_cell : str
        Name of the cell type for which to analyze gene importance.
    col : str, optional
        Name of the time column, by default 'Time_score'.
    genes_index : list, optional
        Indices of genes in the data file, by default an empty list (use all genes).
    p_value : float, optional
        P-value for filtering the fitting models, by default 0.05.
    max_iter_huber : int, optional
        Number of iterations for Huber model, by default 100.
    epsilon_huber : float, optional
        Epsilon parameter for Huber model, by default 1.35.
    x_lim : int, optional
        Limit for x-axis in the plot, by default 4.
    width : int, optional
        Width of the plot, by default 20.
    height : int, optional
        Height of the plot, by default 30.
    store_data : bool, optional
        Whether to save the data, by default True.
    genes_interesting : list, optional
        List of interesting gene names, by default an empty list.
    modify_r2 : bool, optional
        Use modified R squared, by default False.
    model_type : str, optional
        Model type ('HuberRegressor' or 'LinearRegression'), by default 'HuberRegressor'.
    fontsize : int, optional
        Font size for labels and annotations, by default 8.
    alpha : float, optional
        Transparency level for plotting, by default 0.5.
    cmap : str, optional
        Colormap for plotting, by default 'viridis'.
    color_back : str, optional
        Background color of the plot, by default None.
    save_as_pdf : bool, optional
        Whether to save the plot as PDF, by default False.
    plot_genes : bool, optional
        Whether to plot the gene importance, by default True.
    colnames : list, optional
        List of column names (for pathmics data), by default an empty list.
    sample_col : str, optional
        Name of the sample ID column, by default 'sampleID'.
    col_cell : str, optional
        Name of the cell type column, by default 'cell_types'.
    normalize: bool, optional
        Whether to normalize gene exp, by default True.

    Returns
    -------
    None
        Visualizes and saves gene importance plots.
    """

    
    path='Results_PILOT'
    if not os.path.exists(path+'/cells/'+name_cell+'.csv'):
        data=extract_cells_from_gene_expression(adata,sample_col=sample_col,col_cell=col_cell,cell_list=[name_cell],normalize=normalize)
    
    elif os.path.exists(path+'/cells/'+name_cell+'.csv'):
        
        data =loadTarget(path+'/cells/', name_cell)
        
        
    if len(colnames)!=0:  #for pathmics data
        for col in colnames:
            data[col]=np.log(data[col]+1)
        pro=cal_proportions(data)
    else:
        pro=cal_proportions(data)  #Calculate the proportion of Zero for genes
        fiternames_pro=list(pro[pro['proportion']> 0.95]['Gene ID']) #Filter genes with high probability zero fraction
        data=data[data.columns[~data.columns.isin(fiternames_pro)]]
    
    if len(genes_index)==0:
        genes_index=list(range(2, data.shape[1]))
        
    RNA_data = pd.DataFrame()
    RNA_data['label'] = data[col]
    print(f"Name of Cell type : \033[1m{name_cell}\033[0m")

    RNA_target = np.transpose(data.iloc[:,genes_index])
    print("sparsity:" + str(calculate_sparsity(RNA_target)))
    #%%%
    min_target=min(RNA_data['label'])
    max_target=max(RNA_data['label'])
    #%%% fit best models for RNA data

    sorted_best =fit_best_model(RNA_target, RNA_data, model_type,max_iter_huber,epsilon_huber,p_value,modify_r2)
    
 
        
    
    
    #%%% plot RNA data
    if not os.path.exists(path+'/Markers'):
        os.makedirs(path+'/Markers')
    if not os.path.exists(path+'/Markers/'+name_cell):
            os.makedirs(path+'/Markers/'+name_cell)
    
    if save_as_pdf:
        suffix='.pdf'
    else:
        suffix='.png'
    

    if store_data:
    
        
        save_data(sorted_best, 
                   ['Gene ID', 'Expression pattern', 'Slope', 'Fitted function', 'Intercept', 'Treat', 'Treat2', 'adjusted P-value', 'R-squared','mod_rsquared_adj'],
                       path+'/Markers/',name_cell,p_val=p_value,pro=pro)
        
        
        
    if plot_genes:
        with plt.rc_context():
            plot_best_matches(RNA_target, RNA_data,data, sorted_best, "Gene expression",plot_color='tab:orange',num=len(sorted_best.keys()),x_lim=x_lim,width=width,height=height,fontsize=fontsize,alpha=alpha,cmap=cmap,color_back=color_back)
            plt.savefig(path+'/Markers/'+name_cell+'/'+'genes_ranking for cell type '+name_cell+suffix)


        if len(genes_interesting):
            print(" Plots for interesting genes: ")
            filtered_dict = {k:v for (k,v) in sorted_best.items() if k in genes_interesting}
            with plt.rc_context():


                plot_best_matches(RNA_target, RNA_data,data, filtered_dict, "Gene expression",             plot_color='tab:orange',num=len(filtered_dict.keys()),width=width,height=height,x_lim=x_lim,fontsize=fontsize,alpha=alpha,cmap=cmap,color_back=color_back)
                plt.savefig(path+'/Markers/'+name_cell+'/'+'Interesting genes_ranking for cell type '+name_cell+'.png')



    
def reclustering_data(adata,resu=0.01,normalization=False,target_sum=1e6,n_neighbor=15,method_='umap', metric_t = 'cosine',mode='distances',origine_scr_rna=False,dimension_rect=False,n_component=25):
    """
    Recluster data using Louvain clustering algorithm.

    Parameters
    ----------
    adata : AnnData or ndarray
        Annotated data matrix or ndarray containing gene expression data.
    resu : float, optional
        Resolution parameter for Louvain clustering, by default 0.01.
    normalization : bool, optional
        Whether to perform data normalization, by default False.
    target_sum : float, optional
        Target sum for data normalization, by default 1e6.
    n_neighbor : int, optional
        Number of neighbors for neighborhood graph construction, by default 15.
    method_ : str, optional
        Method for constructing the neighborhood graph ('umap' or other), by default 'umap'.
    metric_t : str, optional
        Metric for calculating distances between cells, by default 'cosine'.
    mode : str, optional
        Mode for clustering ('distances' or 'connectivities'), by default 'distances'.
    origine_scr_rna : bool, optional
        Whether the data comes from original single-cell RNA-seq data, by default False.
    dimension_rect : bool, optional
        Whether to use dimensionality reduction (PCA), by default False.
    n_component : int, optional
        Number of components for PCA if dimension_rect is True, by default 25.

    Returns
    -------
    labels : ndarray
        Cluster labels assigned by Louvain clustering algorithm.
    """
    
    if origine_scr_rna:
        #Normalization
        if normalization:
            sc.pp.normalize_total(adata,target_sum)
            sc.pp.log1p(adata)
            adata.raw = adata
            sc.pp.scale(adata, max_value=10) 
        #Use dimention reduction (PCA)
        if  dimension_rect:
            sc.tl.pca(adata)
            sc.pp.neighbors(adata,metric=metric_t,n_neighbors=n_neighbor,method=method_,n_pcs=n_component)
        else:
            sc.pp.neighbors(adata,metric=metric_t,n_neighbors=n_neighbor,method=method_)

        #Using Louvain clustering 
        louvain = Louvain(resolution = resu)
    else:
        adata = sc.AnnData(adata)
        sc.pp.neighbors(adata)
        louvain = Louvain(resolution = resu)
        
        
    if mode=='distances':
        labels = louvain.fit_transform(adata.obsp['distances'])
    else:
        labels = louvain.fit_transform(adata.obsp['connectivities'])

    return labels


            


def cal_proportions(data):
    """
    Compute the proportion of sparsity for genes in the given data.

    Parameters
    ----------
    data : DataFrame
        DataFrame containing gene expression data.

    Returns
    -------
    pro : DataFrame
        DataFrame with columns 'Gene ID', 'proportion', and 'mean'.
        'Gene ID': Gene IDs from the data.
        'proportion': Proportion of sparsity (genes with value 0) for each gene.
        'mean': Mean expression value for each gene.
    """
    pro=pd.DataFrame()
    pro['Gene ID']=list(data.columns)
    pro['proportion']=0
    pro['mean']=0
    for name_genes in list(data.columns):
        if name_genes not in ['sampleID','Time_score']:
            means=data[name_genes].mean()
            count = (data[name_genes] == 0).sum()
            pro.loc[pro['Gene ID'] == name_genes, 'proportion'] = count/len(data)
            pro.loc[pro['Gene ID'] == name_genes, 'mean'] = means
    return pro            
            


def extract_cells_from_pathomics(adata,path=None):
    """
    Extract clusters along with their features and pseudotime.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing features and pseudotime information.
    path : str, optional
        Path to the directory where the extracted data will be saved, by default 'Results_PILOT/'.

    Returns
    -------
    None.
    Saves the extracted data, including features and pseudotime, to the specified path.
    """
    
    if path==None:
        path='Results_PILOT/'
    annot=adata.uns['annot']
    data=adata.uns['data']
    orders=adata.uns['orders']
    data['sampleID']=list(annot['sampleID'])
    joint=pd.merge(data,orders,on='sampleID')
    if not os.path.exists(path+'/cells/'):  
                    os.makedirs(path+'/cells/')
    joint.to_csv(path+'/cells/'+'All.csv')
 


def gene_cluster_differentiation(adata,cellnames=[],sort=['Expression pattern', 'adjusted P-value', 'R-squared'],number_genes=10,cluster_names=[],font_size=14,gene_list=[]):
    
    """
    Perform gene cluster differentiation analysis.

    Parameters
    ----------
    cellnames : list, optional
        List of cell names for which you want to perform gene cluster differentiation analysis, by default [].
    sort : list, optional
        List of criteria for sorting gene expressions, by default ['Expression pattern', 'adjusted P-value', 'R-squared'].
    number_genes : int, optional
        Number of top genes to consider for each expression pattern, by default 10.
    cluster_names : list, optional
        List of cluster names to consider for gene cluster differentiation, by default [].
    font_size : int, optional
        Font size for plots, by default 12.
        gene_list=list, optional
        Your interested genes. If you have any intereted genes!

    Returns
    -------
    None.
    Performs gene cluster differentiation analysis based on specified parameters and saves the results.
    """
    path='Results_PILOT/'
    start=min(adata.uns['orders']['Time_score'])
    end=max(adata.uns['orders']['Time_score'])
    if len(gene_list)==0:
        gene_list=[]
        for cell in cellnames: #Your interested cell, this gets the genes of these cells and compares with others
            data = pd.read_csv(path + '/Markers/' + cell + '/Whole_expressions.csv', index_col = 0)
            specific_data = data.sort_values(sort,
                       ascending=[True, True, False]).groupby('Expression pattern').head(number_genes)
            gene_list.extend(specific_data['Gene ID'].tolist())
        
        gene_list = np.unique(gene_list)
        infer_gene_cluster_differentiation(gene_list,path_to_results = path,font_size=font_size,start=start,
                                       end=end)
    else:
        gene_list = np.unique(gene_list)
        infer_gene_cluster_differentiation(gene_list,path_to_results = path,font_size=font_size,start=start,
                                       end=end)


    
    
def morphological_features_importance(data,name_cell='All',col='Time_score',genes_index=[],p_value=0.05,max_iter_huber=100,epsilon_huber=1.35,x_lim=135,width=20,height=8,store_data=True,genes_interesting=[],modify_r2 = False,model_type = 
'HuberRegressor',fontsize=10,alpha=0.5,cmap='viridis',color_back=None,save_as_pdf=False,plot_genes=True,path=None):
    
    
    """
    Perform importance analysis for morphological features.

    Parameters
    ----------
    data : pandas DataFrame
        Data containing morphological features and time information.
    name_cell : str, optional
        Name of the cell or cluster for which to perform the analysis, by default 'All'.
    col : str, optional
        Name of the column representing time information, by default 'Time_score'.
    genes_index : list, optional
        List of indices of features to include in the analysis, by default [].
    p_value : float, optional
        P-value threshold for significance, by default 0.05.
    max_iter_huber : int, optional
        Maximum number of iterations for the HuberRegressor model, by default 100.
    epsilon_huber : float, optional
        Epsilon parameter for the HuberRegressor model, by default 1.35.
    x_lim : int, optional
        Limit for the x-axis in the plots, by default 135.
    width : int, optional
        Width of the plots, by default 20.
    height : int, optional
        Height of the plots, by default 8.
    store_data : bool, optional
        Whether to store the analysis results, by default True.
    genes_interesting : list, optional
        List of interesting genes to include in the plots, by default [].
    modify_r2 : bool, optional
        Whether to use modified R squared, by default False.
    model_type : str, optional
        Type of regression model, by default 'HuberRegressor'.
    fontsize : int, optional
        Font size for the plots, by default 10.
    alpha : float, optional
        Transparency for the plots, by default 0.5.
    cmap : str, optional
        Color map for the plots, by default 'viridis'.
    color_back : str, optional
        Background color for the plots, by default None.
    save_as_pdf : bool, optional
        Whether to save the plots in PDF format, by default False.
    plot_genes : bool, optional
        Whether to plot the gene importance analysis, by default True.
    path : str, optional
        Path to store the analysis results and plots, by default None.

    Returns
    -------
    None.
    Performs gene importance analysis for morphological features based on specified parameters and saves the results.
    """
    
   
   
    if path==None:
        path='Results_PILOT/'
        
       
    RNA_data = pd.DataFrame()
    RNA_data['label'] = data[col]
    print(f"Name of Cluster : \033[1m{name_cell}\033[0m")
   
    if len(genes_index)!=0:
         RNA_target = np.transpose(data.iloc[:,genes_index])
    else:
        RNA_target = np.transpose(data.iloc[:,0:data.shape[1]-2])
   
    print("Sparsity:" + str(calculate_sparsity(RNA_target)))
    #%%%
    min_target=min(RNA_data['label'])
    max_target=max(RNA_data['label'])
    #%%% fit best models for RNA data

    sorted_best =fit_best_model(RNA_target, RNA_data, model_type,max_iter_huber,epsilon_huber,p_value,modify_r2)
    
 
        
    
    
    #%%% plot RNA data
    if not os.path.exists(path+'/Markers'):
        os.makedirs(path+'/Markers')
    if not os.path.exists(path+'/Markers/'+name_cell):
            os.makedirs(path+'/Markers/'+name_cell)
    
    if save_as_pdf:
        suffix='.pdf'
    else:
        suffix='.png'
    
    pro=cal_proportions(data) #calculate sparsity  
    if store_data:
    
        
        save_data(sorted_best, 
                   ['Gene ID', 'Expression pattern', 'Slope', 'Fitted function', 'Intercept', 'Treat', 'Treat2', 'adjusted P-value', 'R-squared','mod_rsquared_adj'],
                       path+'/Markers/',name_cell,p_val=p_value,pro=pro)
        
        
       
    if plot_genes:
        with plt.rc_context():
            plot_best_matches(RNA_target, RNA_data,data, sorted_best, "Gene expression",plot_color='tab:orange',num=len(sorted_best.keys()),x_lim=x_lim,width=width,height=height,fontsize=fontsize,alpha=alpha,cmap=cmap,color_back=color_back)
            plt.savefig(path+'/Markers/'+name_cell+'/'+'Morphological_ranking for cell type '+name_cell+suffix)

        
        if len(genes_interesting):
            print(" Plots for interesting Morphological features: ")
            filtered_dict = {k:v for (k,v) in sorted_best.items() if k in genes_interesting}
            with plt.rc_context():


                plot_best_matches(RNA_target, RNA_data,data, filtered_dict, "Gene expression",             plot_color='tab:orange',num=len(filtered_dict.keys()),width=width,height=height,x_lim=x_lim,fontsize=fontsize,alpha=alpha,cmap=cmap,color_back=color_back)
                plt.savefig(path+'/Markers/'+name_cell+'/'+'Interesting Morphological_features_ranking for cell type '+name_cell+'.png')

                

                


def norm_morphological_features(path=None,column_names=[],name_cell=None):
    """
    Normalize morphological features.

    Parameters
    ----------
    path : str, optional
        Path to the directory containing the cell data, by default None.
    column_names : list, optional
        List of column names representing morphological features, by default [].
    name_cell : str, optional
        Name of the cell or cluster for which to normalize features, by default None.

    Returns
    -------
    data : pandas DataFrame
        Data with normalized morphological features.
    """
    if path==None:
        path='Results_PILOT/'
    data =loadTarget(path+'/cells/', name_cell)
    for col in column_names:
            data[col]=np.log(data[col]+1)
            
    return data
                
             
def results_gene_cluster_differentiation(cluster_name=None,sort_columns=['pvalue'],ascending=[True],threshold=0.5,p_value=0.01):
    """
    Retrieve and sort gene cluster statistics based on specified criteria.

    Parameters:
        cluster_name : str, optional
            The name of the gene cluster to retrieve statistics for. Default is None.
        sort_columns : list of str, optional
            List of column names to sort the data by. Default is ['pvalue'].
        ascending : list of bool, optional
            List indicating sorting order for each corresponding column. Default is [True].
        threshold: float, optional
             Select genes with fold changes higher than the defined threshold. Default is 0.5.
        p_value: float, optional
               Select genes with the Wald test p-value less than the defined one. Default is 0.01.

    Returns:
        pandas.DataFrame
            A sorted dataframe containing gene cluster statistics based on the specified criteria.
    """
    path='Results_PILOT/'
    statistics=pd.read_csv(path+'/gene_clusters_stats_extend.csv')
    df=statistics[statistics['cluster']==cluster_name]
    df['FC']=df['FC'].astype(float)
    df=df[df['FC'] > threshold]
    df['pvalue']=df['pvalue'].astype(float)
    df=df[df['pvalue'] < p_value]
    df_sorted = df.sort_values(by=sort_columns, ascending=ascending)
 
    return df_sorted[['gene','cluster','waldStat','pvalue','FC','Expression pattern','fit-pvalue','fit-mod-rsquared']]  





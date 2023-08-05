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
warnings.filterwarnings('ignore')


'''
This function loads the h5ad data, set the path.   
Indicate the @path for reading the data
'''


def load_h5ad(path):
    if isfile(path):
        adata=sc.read_h5ad(path)
        return adata
    else:
        print('There is no such data, check the path or name')
    
    
    
'''
This function creates a folder named Results_PILOT for saving the whole generated results by PILOT 
 @param:name dataset: Name of your data
'''

def set_path_for_results(name_dataset):
    if not os.path.exists('Results_PILOT/'+name_dataset):
        os.makedirs('Results_PILOT/'+name_dataset)
        return 'Results_PILOT/'+name_dataset
    else:
        print('The path has already been set.')
        return 'Results_PILOT/'+name_dataset
        
'''
For extracting annotation from scRNA
'''
def extract_annot_expression(adata,columns=['cell_type_original','patient_region','region','X_pca'],reclustering=False,reduction=False,resu=0.1,max_value=10,target_sum=1e4):
    
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

 
'''
This function extracts the needed inputs for PILOT
 @param:adata: loaded Anndata.
 @param:emb_matrix: PCA representation of data (variable)
 @param:clusters_col: cell_type/clustering column name in observation level of your  Anndata
 @param:sample_col: samples/patients column name in observation level of your Anndata
 @param:status: status/disease column name, e.g. control/case 
 @param:name_dataset: name of your data ,PILOT creates a folder with this name for saving the whole generated results 
'''

def extract_data_anno_scRNA_from_h5ad(adata,emb_matrix='PCA',clusters_col='cell_type',sample_col='sampleID',status='status' ,name_dataset='unnamed_dataset'):
    data=adata.obsm[emb_matrix]  
    col_add=[]
    for i in range(1,adata.obsm[emb_matrix].shape[1]+1):
        col_add.append('PCA_'+str(i))
    data=pd.DataFrame(data,columns=col_add) 
    data = data.reset_index(drop=True)
    annot=adata.obs[[clusters_col,sample_col,status]]
    annot.columns=['cell_type','sampleID','status']
    annot = annot.reset_index(drop=True) 
    path_to_results=set_path_for_results(name_dataset)
    
    return data,annot,path_to_results


'''
For extracting annotation from pathomics data
'''
def extract_data_anno_pathomics_from_h5ad(adata,var_names=[],clusters_col='Cell_type',sample_col='sampleID',status='status' ,name_dataset='unnamed_dataset'):
    data=adata[:,var_names].X
    data=pd.DataFrame(data,columns=var_names)
    data = data.reset_index(drop=True)
    annot=adata.obs[[clusters_col,sample_col,status]]
    annot.columns=['cell_type','sampleID','status']
    annot = annot.reset_index(drop=True) 
    path_to_results=set_path_for_results(name_dataset)
    
    return data,annot,path_to_results

              
       
    
def load_annot(path,name):
    """
    loads the annotation matrix
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
    loads the gene expression matrix
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
        
                  
      
'''
Calculating preportions of culsters per sample.
 @param:df is annot data (includes sample,cells/clusters and status)
 @param: cell_col is cell-type/clusters column of your, defualt is the 0th column of your annotation data
 @param: sample  is sample column 
 @param: regulizer is regulizers 
'''                  
    

def Cluster_Representations(df, cell_col = 0, sample_col = 1,regulizer=0.2,normalization=True):    

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

'''
Computing cost matrix, it finds distances between clusters/cell-types,
First calculates the median of each cluster and then dissimilarity between clusters is computed.  
 @param: annot is annotation data
 @param: data are PCAs 
 @param: path for saving matrix
 @param: cell_col is cells/cluster column
 @metric: metric, the measurement for calculating the distances between centroids 
'''

def cost_matrix(annot,data,path,cell_col = 0,metric='cosine',figsize_h=12,figsize_w=12,col_cluster=True,row_cluster=True,cmap='Blues_r',font_scale=2):
    

    cells = annot[annot.columns[cell_col]].unique() # gets cells 
    centroids = []

    for i in cells:
        centroids.append(list(data[annot[annot.columns[cell_col]] == i].median(axis=0))) #get the centriod of each cluster

    dis_t = scipy.spatial.distance.pdist(centroids, metric = metric) 
    dis = scipy.spatial.distance.squareform(dis_t, force ='no', checks = True)
    
    fig = plt.figure()
    
    cost = pd.DataFrame.from_dict(dis).T
    cost.columns=annot.cell_type.unique()
    cost['cell_types']=annot.cell_type.unique()
    cost=cost.set_index('cell_types')
    sns.set(font_scale=font_scale)
    sns.clustermap(cost[annot.cell_type.unique()],cmap=cmap,figsize=(figsize_h,figsize_w),col_cluster=col_cluster,row_cluster=row_cluster, tree_kws={"linewidths": 0.});
    
    plt.title('Cost Matrix',loc='center')
    plt.savefig(path+'/Cost_matrix.pdf') 
    plt.close(fig)     
    return dis

'''
Computes Wassertein disatcnes among samples, the defualt method if callasical OT (unreg),
you can change OT with Regularation one.
@param: Clu_rep is proportions matrix
@param: cost is cost matrix
@param: regularized  is type of OT, by default uses calssical OT without any regulariztion, otherwise applys sinkhorn_stabilized
'''
def wasserstein_d(Clu_rep, cost,regularized = "unreg", reg = 0.1, path=None,figsize_h=12,figsize_w=12):
    
   
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
                
    fig = plt.figure()
    emd = pd.DataFrame.from_dict(EMD).T
    emd.columns=samples_id 
    emd['sampleID']=samples_id 
    emd=emd.set_index('sampleID')
    sns.set(font_scale=2)
    sns.clustermap(emd,cmap='Blues_r',figsize=(figsize_h,figsize_w),col_cluster=True,row_cluster=True, tree_kws={"linewidths": 0.})
    plt.title('Wasserstein distance',loc='center')
    plt.savefig(path+'/Wasserstein distance.pdf') 
    plt.close(fig)


    return EMD

'''
Doing the clustering based on the OT distance
@param: EMD, wasserstein distances
@param: df , annotation data
@param: category name of the column as disease/status in annot data
@param: sample_col, is the sample column in annot data
@parma: res, resolution od Leiden clustering
@param: metric, is the measurement for neighborhood graph in the knn space
@parma: steper, is the steper for increasin/decreasing resolution till getting the number of real lables
'''
def Clustering(EMD, df, category = 'status', sample_col=1,res = 0.01,metric ='cosine',steper=0.01):

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

    
    print("Cluster labels: ", df[category].unique())
    
    true_labels = []
    for i in range(len(samples)):
        a = condition[df[df.columns[sample_col]] == samples[i]].unique()
        true_labels.append(a[0])
    
    S = rand_score(true_labels, labels)
    print("ARI: ", S)
    return labels, S, true_labels;



def Sil_computing(EMD, real_labels, metric='cosine'):
    Silhouette = metrics.silhouette_score(EMD, real_labels, metric =metric)
    print("Silhouette score: ", Silhouette) 
    return Silhouette



'''
Finding the trajecories by applying DM
@param: EMD, is the output of OT,
@param: n_evecs, number of embeddings
@param: k : int, optional Number of nearest neighbors over which to construct the kernel.
       
@param: epsilon, string or scalar, optional
            Method for choosing the epsilon.  Currently, the only options are to provide a scalar (epsilon is set to the provided scalar) 'bgh' (Berry, Giannakis and Harlim), and 'bgh_generous' ('bgh' method, with answer multiplied by 2.

@param: alpha, normalization in bandwith function

Reference:
 .. [1] T. Berry, and J. Harlim, Applied and Computational Harmonic Analysis 40, 68-96
           (2016).
'''
def trajectory(EMD, predicted_labels,df, n_evecs = 2, epsilon =1, alpha = 0.5,knn= 64, sample_col=1, clusters = 'status',label_act = False, path=None,colors=['#377eb8','#ff7f00','#e41a1c'],location_labels='center', fig_h=12,fig_w=12,font_size=24,axes_line_width=1,axes_color='black',facecolor='white',point_size=100,cmap='viridis',fontsize_legend=24,alpha_trans=1,plot_titel = "Trajectory of the disease progression"):
    
    custom_cycler = (cycler(color=colors))
    
    plot_titel = plot_titel

   
    mydmap = diffusion_map.DiffusionMap.from_sklearn(n_evecs = n_evecs, epsilon =epsilon, alpha = alpha, k=knn)
    embedding = mydmap.fit_transform(EMD)
        
    plt.rcParams.update({'font.size': font_size})
    plt.rcParams["axes.edgecolor"] = axes_color
    plt.rcParams["axes.linewidth"] = axes_line_width

    fig = plt.figure()
    fig.set_size_inches(fig_h, fig_w)
    
    ax = plt.gca()
    
    ax.set(facecolor = facecolor)
    ax.set_prop_cycle(custom_cycler)
    scatter = ax.scatter(embedding[:,0],embedding[:,1], c = predicted_labels, label = predicted_labels,s=point_size)
    legend= ax.legend(*scatter.legend_elements(), loc="lower left", title="Classes")
    ax.add_artist(legend)
   
    
    plt.title(plot_titel+" (predicted_lables) ")
    plt.savefig(path+"/"+plot_titel+"(predicted_lables)"+'.pdf')    
    plt.close(fig)

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
    frame = legend.get_frame()
    frame.set_color(facecolor)
    

    plt.title(plot_titel)
    plt.savefig(path+"/"+plot_titel+'.pdf')
    plt.show()         
    plt.close(fig)
    

    
    return embedding


'''
Getting the pseudotime based on the pricipla graph

@param: emb, embeddings of DM.
@param: NumNodes, number of node for buiding the backbone.
@param: source_node, start from this node
@param: show_text, showing the numbers in the backbone
@param: Do_PCA Perform PCA on the nodes
@param: DimToPlot a integer vector specifing the PCs (if Do_PCA=TRUE) or dimension (if Do_PCA=FALSE) to plot. 
@param: Node_color, colors of backbon's  node.
@return orders of samples based on the trajectory and selected cell-types
'''



def fit_pricipla_graph(emb,NumNodes=10,source_node=0,path=None,show_text=True,Do_PCA=False,fig_x_size=12,fig_y_size=12,X_color='r', Node_color='k', DimToPlot=[0, 1],facecolor='white',title='Principal graph'):
 
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
    return pseudotime


'''
Ordering cells based on the estimated time by PILOT.
@param: bins, proportions.
@param: annot, annotation data.
@param: embedding_diff, emds of DM. 
@param: pseudotime, calculated pseudotime with fit_pricipla_graph.
@param: heatmap_h,heatmap_w, height and width of heatmap.
@param: p_val, p-value for filtering the found fitting models.
@param: plot_cell: ploting the plots
@param witdth and height, size of plots
@param xlim, for grouping the disease progression(pseudotime) and showing the groups over x-axis
@param color_back, for backgroundcolor
@param save_as_pdf, format for saving
@param alpha, adjust the transparency 
'''
def Cell_importance(bins,annot,embedding_diff,real_labels,path,pseudotime,heatmap_h=12,heatmap_w=12,width=25,height=25,xlim=5,p_val=1,plot_cell=True,point_size=100,color_back=None,fontsize=20,alpha=1,cmap='viridis',save_as_pdf=False):
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
    
    
    return pathies_cell_proportions[['sampleID','Time_score']],cellnames
   

'''
Ordering genes based on the estimated time by PILOT
@param: pro, the percentage of the sparsity of genes.
@param: data, returned data for each cell by extract_cells_from_gene_expression function.
@param: name_cell, name of the cell
@param:col, name of the time column
@param: genes_index, the indices of genes in the data file
@param: store_data, if save the plots and output
@param:genes_interesting, the name of your interested genes
@param:modify_r2, use modified R squared
@param:model_type, choose traditional lose function or HuberRegressor
@param:p_value, default is 0.05
@param: max_iter_huber, number of iteration for huber model
@param: epsilon_huber: epsilon parameter for huber model
@param:x_lim  for grouping the disease progression(pseudotime) and showing the groups over x-axis
@param color_back, for backgroundcolor
@param save_as_pdf, format for saving
@param alpha, adjust the transparency

'''
            
def genes_importance(pro,data,path,name_cell,col='Time_score',genes_index=[],p_value=0.05,max_iter_huber=500,epsilon_huber=1.5,x_lim=4,width=12,height=12,store_data=True,genes_interesting=[],modify_r2 = True,model_type = 'HuberRegressor',fontsize=20,alpha=0.5,cmap='viridis',color_back=None,save_as_pdf=False,plot_genes=False):
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

'''
Extracts genes expressed. with from object, this function after extracting genes joins them with pseudotime
@param: adata is the Anndata
@param: order, the samples with their pseudotime (output of Cell_importance function)
@cell_list, cells which are changing significantly over the trajectory(pseudotime)
Return each cell beside genes and their time associated with the trajectory of samples, they are saved in the cell folder.
'''

def extract_cells_from_gene_expression(adata,orders,sample_col,col_cell,cell_list,path_results):
            
        cell_types=cell_list
        for cell in cell_types:
            
            adata_new = adata[adata.obs[col_cell].isin([cell]),:]
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



            
'''
Computing the proportions of sparsity for genes
@param: data, the cell with its genes
'''
def cal_proportions(data):
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
            
'''
Extracting clusters with their features and time 
@param:orders, the order by PILOT,
@param:annot, annotation data
@param:data, features data
'''
def extract_cells_from_pathomics(orders,annot,data,path):
    
    data['sampleID']=list(annot['sampleID'])
    joint=pd.merge(data,orders,on='sampleID')
    if not os.path.exists(path+'/cells/'):  
                    os.makedirs(path+'/cells/')
    joint.to_csv(path+'/cells/'+'All.csv')
    
    

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
from sklearn.decomposition import PCA
from sklearn.cluster import SpectralClustering
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import NearestCentroid
from sklearn.metrics.cluster import rand_score
from sklearn.cluster import DBSCAN
from scipy.spatial import distance
from sklearn.manifold import MDS
from sknetwork.clustering import Louvain
import shap
from sklearn.ensemble import RandomForestClassifier,RandomForestRegressor
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import train_test_split
import time
import sys 
from genericpath import isfile
from .Cell_gene_selection import *
import warnings
import ot
from logging import info, warn
from cycler import cycler
#custom_cycler = (cycler(color=['#1589e8','#f0c3c3','#ec0c0e']))


def load_h5ad(path):
    if isfile(path):
        adata=sc.read_h5ad(path)
        return adata
    else:
        print('There is no such data, check the path or name')
    
    
    
    
def set_path_for_results(name_dataset):
    if not os.path.exists('Results_PILOT/'+name_dataset):
        os.makedirs('Results_PILOT/'+name_dataset)
        return 'Results_PILOT/'+name_dataset
    else:
        print('The path has already been set.')
        return 'Results_PILOT/'+name_dataset
        

def extract_annot_expression(adata,columns=['cell_type_original','patient_region','region','X_pca'],reclustering=False,reduction=False,resu=0.1):
    
        if reduction:
                    sc.pp.normalize_total(adata, target_sum=1e4)
                    sc.pp.log1p(adata)
                    sc.pp.scale(adata, max_value=10)
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
                warn("loading " + name + " failed, failed, check the path/name of this data")
       
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
        
                  
      
                  
    

def Cluster_Representations(df, cell_col = 0, sample_col = 1, regularization = None,regulizer=10):    

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
        
    
    if regularization == True:
        for x in df['sampleID'].unique():
            dict[x] = (dict[x]+prior)/(sum(dict[x])+sum(prior))
            #vector = np.ones(len(df['cell_type'].unique()))*max(bins[x])*0.1
           # vector = np.ones(len(df['cell_type'].unique()))*regulizer #len(df['cell_type'].unique())
           # dict[x] = (dict[x]+vector)/(sum(dict[x])+sum(vector))
            
            
            
    
    return dict


def bins_agreg(bins):
    collector = np.zeros(len(bins[list(bins.keys())[0]]))
    for i in bins.keys():
        collector = collector + bins[i]
    
    return collector

def find_outliers(arrayMatrix):
    q75,q25 = np.percentile(arrayMatrix,[75,25])
    intr_qr = q75-q25
 
    max_v = q75+(1.5*intr_qr)
    min_v = q25-(1.5*intr_qr)
    
    z_scores = [(y > max_v or y < min_v) for y in arrayMatrix]
    return np.where(np.array(z_scores) == True)


def filter_outliers(outliers,annot,data):
    
    n_cell = len(annot.cell_type.unique())

    to_remove = annot.sampleID.unique()[outliers]
    
    for i in to_remove:
        data = data.drop(annot[annot.sampleID == i].index)
        annot = annot.drop(annot[annot.sampleID == i].index)

    return data, annot


def data_regularization(annot, data, collector, min_count = 10, max_count = 100):
    
    significant_cells =[]
    name_cells = []
    boolean_list=[]
    names = annot.cell_type.unique()  
    k=0
    
    for i in collector:
        if min_count < i < max_count:
            significant_cells.append(i)
            name_cells.append(names[k])
        k=k+1 
    
    for i in annot.index.values:
        boolean_list.append(annot['cell_type'].loc[i] in name_cells)
        
    return annot[boolean_list],data[boolean_list]

    

def cost_matrix(annot, data,path,cell_col = 0):
    

    cells = annot[annot.columns[cell_col]].unique()
    centroids = []

    for i in cells:
        centroids.append(list(data[annot[annot.columns[cell_col]] == i].median(axis=0)))

    dis_t = scipy.spatial.distance.pdist(centroids, metric = 'cosine')
    dis = scipy.spatial.distance.squareform(dis_t, force ='no', checks = True)
    
    fig = plt.figure()
    
    cost = pd.DataFrame.from_dict(dis).T
    cost.columns=annot.cell_type.unique()
    cost['cell_types']=annot.cell_type.unique()
    cost=cost.set_index('cell_types')
    sns.set(font_scale=2)
    sns.clustermap(cost[annot.cell_type.unique()],cmap='Blues_r',figsize=(12,12),col_cluster=True,row_cluster=True, tree_kws={"linewidths": 0.});
    
    plt.title('Cost Matrix',loc='center')
    plt.savefig(path+'/Cost_matrix.pdf') 
    plt.close(fig)     
    return dis


def wasserstein_d(bins, cost,regularized = "unreg", reg = 0.1, path=None):
    
   
    samples_id = list(bins.keys())
    n_samples = len(samples_id)
    EMD = np.zeros((n_samples,n_samples))
    
    if regularized == "unreg":
        for i in range(n_samples):
            for j in range(n_samples):
                #EMD[i,j] = ot.unbalanced.mm_unbalanced2(bins[samples_id[i]],bins[samples_id[j]], cost, 0.8, 'kl')
                EMD[i,j] = ot.emd2(bins[samples_id[i]],bins[samples_id[j]], cost)
    else:
        for i in range(n_samples):
            for j in range(n_samples):
                EMD[i,j] = ot.sinkhorn2(bins[samples_id[i]],bins[samples_id[j]], cost, reg, method = "sinkhorn_stabilized")
                
    fig = plt.figure()
    emd = pd.DataFrame.from_dict(EMD).T
    emd.columns=samples_id 
    emd['sampleID']=samples_id 
    emd=emd.set_index('sampleID')
    sns.set(font_scale=2)
    sns.clustermap(emd,cmap='Blues_r',figsize=(12,12),col_cluster=True,row_cluster=True, tree_kws={"linewidths": 0.})
    plt.title('Wasserstein distance',loc='center')
    plt.savefig(path+'/Wasserstein distance.pdf') 
    plt.close(fig)


    return EMD


def Clustering(EMD, df, category = 'status', sample_col=1,res = 0.01):

    samples = df[df.columns[sample_col]].unique()
    condition = df[category]
    n_clust = len(df[category].unique())
    
    flag = True
    inc = 0.01
    while flag == True:
        adata = sc.AnnData(EMD)
        sc.pp.neighbors(adata, metric ='cosine')
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



def Sil_computing(EMD, real_labels, metric, space = 'diffusion'):

    n_comp = EMD.shape[0]
    
    if space == 'MDS':
        mydmap = MDS(n_components=n_comp)
        embedding = mydmap.fit_transform(EMD)
    
    elif space == 'PCA':
        mydmap = PCA(n_components=n_comp)
        embedding = mydmap.fit_transform(EMD)
        
    elif space == 'diffusion':
        embedding = EMD
    
    Silhouette = metrics.silhouette_score(embedding, real_labels, metric ='cosine')
    print("Silhouette score: ", Silhouette)
        
    return Silhouette


def trajectory(EMD, predicted_labels,df, knn= 64, sample_col=1, clusters = 'status', embed_coord = 'diffusion',label_act = False, path=None,colors=['#377eb8','#ff7f00','#e41a1c'],location_labels='center'):
    
    custom_cycler = (cycler(color=colors))
    
    plot_titel = "Trajectory of the disease progression"

    if embed_coord == 'diffusion':
        mydmap = diffusion_map.DiffusionMap.from_sklearn(n_evecs = 2, epsilon =1, alpha = 0.5, k=64)
        embedding = mydmap.fit_transform(EMD)
        
    elif embed_coord == 'MDS':
        mydmap = MDS(n_components=2)
        embedding = mydmap.fit_transform(EMD)
        
    elif embed_coord == 'PCA':
        mydmap = PCA(n_components=2)
        embedding = mydmap.fit_transform(EMD)
    
    else:
        print('Choose one of the embeddings options including diffusion map, MDs or PCA.')
      
    
    
    
    plt.rcParams.update({'font.size': 16})
    plt.rcParams["axes.edgecolor"] = "black"
    plt.rcParams["axes.linewidth"] = 1

    fig = plt.figure()
    fig.set_size_inches(12, 12)
    
    ax = plt.gca()
    
    ax.set(facecolor = "white")
    ax.set_prop_cycle(custom_cycler)
    scatter = ax.scatter(embedding[:,0],embedding[:,1], c = predicted_labels, label = predicted_labels,s=100)
    legend= ax.legend(*scatter.legend_elements(), loc="lower left", title="Classes")
    ax.add_artist(legend)
   
    
    plt.title(plot_titel+" (predicted_lables) ")
    plt.savefig(path+"/"+plot_titel+"(predicted_lables)"+'.pdf')    
    plt.close(fig)

    
    

    df = df.drop_duplicates(subset =[df.columns[sample_col]])
    fig = plt.figure()
    fig.set_size_inches(12, 12)
    
    ax = plt.gca()
    
    ax.set(facecolor = "white")
    ax.set_prop_cycle(custom_cycler)
    #plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['#1589e8','#f0c3c3','#ec0c0e'])
    for category in df[clusters].unique():
        aux = np.array(df[clusters] == category)
        group_data = embedding[aux]
        scatter = ax.scatter(group_data[:,0],group_data[:,1], alpha=1,label=category, cmap='viridis',s=100) 
        
        if label_act == True:
            names = np.array(df[df[clusters] == category].sampleID)
            k = 0
            for txt in names:
                ax.annotate(txt, (group_data[k,0], group_data[k,1]),fontsize=12)
                k=k+1

    ax.legend(loc=location_labels, fontsize=24)
    frame = legend.get_frame()
    frame.set_color('white')
    

    plt.title(plot_titel)
    plt.savefig(path+"/"+plot_titel+'.pdf')
    if embed_coord =='diffusion':
        plt.show()         
    plt.close(fig)
    

    
    return embedding

def plot_by_pheno(data, predicted_labels, class_id = 2):
    boolean = []
    k=0
    for i in predicted_labels:
        if i == class_id:
            boolean.append(k)
        else:
            pass       
        k=k+1
        
    category = data.loc[boolean]
    umap_object = umap.UMAP(n_components=2, random_state=42)
    embedding_umap = umap_object.fit_transform(category)

    fig,ax = plt.subplots()
    scatter = ax.scatter(embedding_umap[:,0],embedding_umap[:,1])
    ax.set_title('Cluster projection by phenotypic state: ' + str(class_id))
    plt.show() 

def contigency_mat(true_labels,predicted_labels, normalize = 'index'):
    contigency = {'predicted_labels': predicted_labels, 'true_labels': true_labels}
    cont_matrix = pd.DataFrame(data=contigency)
    data_crosstab = pd.crosstab(cont_matrix['predicted_labels'], cont_matrix['true_labels'], margins = False,normalize= normalize)
    print(data_crosstab)
    return data_crosstab




def Cell_importance(bins,annot, embedding_diff, real_labels, path, sort_axis='emb_x', p_val=0.05):
    
    cell_types_propo = bins
    patients_id = bins.keys()
    cell_types = annot['cell_type'].unique()
    emd = embedding_diff
    labels = real_labels
    #creat a dataframe of samples and their pseuduscores
    emd_dataframe = pd.DataFrame({'sampleID':list(patients_id), 'emb_x':emd[:,0],'emb_y':emd[:,1],'lables':list(labels)},dtype=object)
    #sort samples based on defined axis, you should choose correct one!
    emd_dataframe_sort = emd_dataframe.sort_values(sort_axis, ascending=True) 
    orders=list(emd_dataframe_sort['sampleID'])
    times=list(range(1, len(orders)+1))
    pathies_cell_proportions = pd.DataFrame.from_dict(bins).T
    pathies_cell_proportions.columns=cell_types
    pathies_cell_proportions.index.name='sampleID'
    df_join = pd.merge(emd_dataframe_sort['sampleID'], pathies_cell_proportions, how='inner', on = 'sampleID')
    df_join=df_join.set_index('sampleID')
    #Normalizing the proportions for heat map
    normalized_df=(df_join-df_join.min())/(df_join.max()-df_join.min())
 
    #Saving Heat map based on sorte pseuduscores of the Trajectory 
    
    sns.clustermap(normalized_df[cell_types],row_cluster=False,annot=False,cmap='Blues',figsize=(12,12),xticklabels=True);
    plt.savefig(path+"/"+'Samples_over_trajectory.pdf')
    
    #Building a model based on Regression and pseuduscores 
    pathies_cell_proportions['Time_score']=list(emd_dataframe[sort_axis])
    pathies_cell_proportions = pathies_cell_proportions.sort_values('Time_score', ascending=True)
    pathies_cell_proportions['Time_score']=list(times)
    pathies_cell_proportions=pathies_cell_proportions.reset_index()
    RNA_data = pd.DataFrame()
    RNA_data['label'] = pathies_cell_proportions['Time_score']
    RNA_target = np.transpose(pathies_cell_proportions.iloc[:,1:len(annot.cell_type.unique())+1])
    sorted_best=fit_best_model(RNA_target, RNA_data, model_type = 'LinearRegression', pval_thr = p_val, modify_r2 = False)
        
    with plt.rc_context():
            min_x = min(np.array(list(RNA_data['label'])))
            max_x = min(np.array(list(RNA_data['label'])))
            start = pathies_cell_proportions[pathies_cell_proportions['Time_score'] == min_x]['sampleID'].unique()
            end = pathies_cell_proportions[pathies_cell_proportions['Time_score'] == max_x]['sampleID'].unique()
            plot_name = 'From  ' + start + '  to  ' + end
            plot_best_matches_cell_types(RNA_target, RNA_data, sorted_best, plot_name, "Cell Proportion")
            plt.savefig(path+"/"+'Cell_types_importance.pdf')
    
    
    cellnames=list(sorted_best.keys())
    return pathies_cell_proportions[['sampleID','Time_score']],cellnames
   

def feature_importance_shap(EMD,bins,annot,embedding_diff,real_labels,path,sort_axis='emb_x'):
    Similarities=EMD
    cell_types_propo=bins
    patients_id=bins.keys()
    cell_types = annot['cell_type'].unique()
    emd=embedding_diff
    labels=real_labels
    #creat a dataframe of samples and their pseuduscores
    emd_dataframe = pd.DataFrame({'sampleID':list(patients_id), 'emb_x':emd[:,0],'emb_y':emd[:,1],'lables':list(labels)},dtype=object)
    #sort samples based on defined axis, you should choose correct one!
    emd_dataframe_sort = emd_dataframe.sort_values(sort_axis, ascending=False)
    orders=list(emd_dataframe_sort['sampleID'])
    
    #get the smaples proportions 
    pathies_cell_proportions = pd.DataFrame.from_dict(bins).T
    pathies_cell_proportions.columns=cell_types
    pathies_cell_proportions.index.name='sampleID'
    df_join = pd.merge(emd_dataframe_sort['sampleID'], pathies_cell_proportions, how='inner', on = 'sampleID')
    df_join=df_join.set_index('sampleID')
    #Normalizing the proportions for heat map
    normalized_df=(df_join-df_join.min())/(df_join.max()-df_join.min())
 
    #Saving Heat map based on sorte pseuduscores of the Trajectory 
    
    sns.clustermap(normalized_df[cell_types],row_cluster=False,annot=False,cmap='Blues',figsize=(12,12),xticklabels=True);
    plt.savefig(path+"/"+'heat_map_trjectory.png')
    
    #Building a model based on Regression and pseuduscores 
    pathies_cell_proportions['labels']=list(emd_dataframe[sort_axis])
    Y = pathies_cell_proportions['labels']

    X = pathies_cell_proportions[cell_types]
    # Split the data into train and test data:
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size = 0.3)
    # Build the model with the random forest regression algorithm:
    model = RandomForestRegressor(max_depth=10, random_state=0, n_estimators=10)
    #Training the model
    model.fit(X_train, Y_train)
    explainer = shap.TreeExplainer(model)
    
    fig2 = plt.figure()
  
    
    shap.summary_plot(explainer.shap_values(X_train),X_train,plot_type='dot',show=False,plot_size=[12,12]);
    
    
        
    fig2.savefig(path+"/"+'feature_importance_dot.pdf', bbox_inches='tight')   
    fig3 = plt.figure()
 
    shap.summary_plot(explainer.shap_values(X_train), X_train, plot_type="bar",show=False,plot_size=[12,12]);
   
    fig2.savefig(path+"/"+'feature_importance_dot.png')
    fig3.savefig(path+"/"+'feature_importance_bar.png')
 
    return orders

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

            
def genes_importance(pro,data,path,name_cell,col='Time_score',genes_index=[],p_value=0.05,max_iter_huber=500,epsilon_huber=1.5,x_lim=4,store_data=1,genes_interesting=[],modify_r2 = True,model_type = 'HuberRegressor'):
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
    
    
    

    if store_data:
    
        
        save_data(sorted_best, 
                   ['Gene ID', 'Expression pattern', 'Slope', 'Fitted function', 'Intercept', 'Treat', 'Treat2', 'adjusted P-value', 'R-squared','mod_rsquared_adj'],
                       path+'/Markers/',name_cell,p_val=p_value,pro=pro)

    with plt.rc_context():
        plot_best_matches(RNA_target, RNA_data,data, sorted_best, "Gene expression", plot_color='tab:orange',num=len(sorted_best.keys()),x_lim=x_lim)
        plt.savefig(path+'/Markers/'+name_cell+'/'+'genes_ranking for cell type '+name_cell+'.pdf')
    
    
    if len(genes_interesting):
        print(" Plots for interesting genes: ")
        filtered_dict = {k:v for (k,v) in sorted_best.items() if k in genes_interesting}
        with plt.rc_context():
            
            
            plot_best_matches(RNA_target, RNA_data,data, filtered_dict, "Gene expression",             plot_color='tab:orange',num=len(filtered_dict.keys()),x_lim=x_lim)
            plt.savefig(path+'/Markers/'+name_cell+'/'+'Interesting genes_ranking for cell type '+name_cell+'.pdf')

def proportions(data):
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
            
            
  

### Trajectory Analysis of Kidney IgAN Data with PILOT

<div class="alert alert-block alert-info">
<b>PILOT</b>

Welcome to the PILOT Package Tutorial for pathomics Data!
 
You can find the pathomics data [here](https://github.com/CostaLab/PILOT/tree/main/Tutorial/Datasets).

</div>


```python
import PILOT as pl
import scanpy as sc
```

#### Kidney_IgAN Tubuli

##### Reading Anndata
<div class="alert alert-block alert-info">
The Anndata (h5ad) file is in the Datasets folder.
</div>


```python
adata_T=sc.read_h5ad('Datasets/Kidney_IgAN_T.h5ad') 
```

#### Loading the required information and computing the Wasserstein distance:
<div class="alert alert-block alert-info"> In order to work with PILOT, ensure that your Anndata object is loaded and contains the required information.
    
Use the following parameters to configure PILOT for your analysis (Setting Parameters):
    
adata: Pass your loaded Anndata object to PILOT.
    
emb_matrix: Provide the name of the variable in the obsm level that holds the dimension reduction ( PCA representation).
    
clusters_col: Specify the name of the column in the observation level of your Anndata that corresponds to cell types or clusters.
    
sample_col: Indicate the column name in the observation level of your Anndata that contains information about samples or patients.
    
status: Provide the column name that represents the status or disease (e.g., "control" or "case").
  
</div>


```python
pl.tl.wasserstein_distance(adata_T,clusters_col='Cell_type',sample_col='sampleID',status='status' 
                           ,data_type='Pathomics')
```

##### plotting the Cost matrix and the Wasserstein distance:
<div class="alert alert-block alert-info"> 
 Here we show the heatmaps of the Cost matrix (cells) and Wasserstein distance (samples).      
</div>


```python
pl.pl.heatmaps(adata_T)
```


    
![png](Combination_Kidney_IgAN_files/Combination_Kidney_IgAN_9_0.png)
    



    
![png](Combination_Kidney_IgAN_files/Combination_Kidney_IgAN_9_1.png)
    


##### Trajectory:
<div class="alert alert-block alert-info"> 
 Here we show the Diffusion map of Wasserstein distance. In the upcoming trajectory analysis, the labels "<30" signify an eGFR level below 30 (considered as Low), "30-60" denotes a reduced eGFR range (marked as Reduced), and ">60" indicates a normal eGFR level
</div>


```python
pl.pl.trajectory(adata_T,colors=['red','blue','orange'])
```


    
![png](Combination_Kidney_IgAN_files/Combination_Kidney_IgAN_11_0.png)
    


#### Kidney_IgAN Glomeruli

##### Reading Anndata
<div class="alert alert-block alert-info">
The Anndata (h5ad) file is in the Datasets folder.
   </div>


```python
adata_G=sc.read_h5ad('Datasets/Kidney_IgAN_G.h5ad') #First read the object
```


```python
pl.tl.wasserstein_distance(adata_G,clusters_col='Cell_type',sample_col='sampleID',status='status'
                          ,data_type='Pathomics')
```


```python
pl.pl.heatmaps(adata_G)
```


    
![png](Combination_Kidney_IgAN_files/Combination_Kidney_IgAN_16_0.png)
    



    
![png](Combination_Kidney_IgAN_files/Combination_Kidney_IgAN_16_1.png)
    



```python
pl.pl.trajectory(adata_G,colors=['red','blue','orange'])
```


    
![png](Combination_Kidney_IgAN_files/Combination_Kidney_IgAN_17_0.png)
    


#### Combination:
<div class="alert alert-block alert-info"> 
Here, we combine the distances of samples. We get the sum of distances of samples based on Tubuli and Glomeruli distances.   
</div>


```python
adata_Com=adata_G
adata_Com.uns['EMD']=adata_G.uns['EMD']+adata_T.uns['EMD']
```


```python
pl.pl.trajectory(adata_Com,colors=['red','blue','orange'])
```


    
![png](Combination_Kidney_IgAN_files/Combination_Kidney_IgAN_20_0.png)
    


####  Fit a principal graph:
<div class="alert alert-block alert-info"> 
The difussion map creates an embedding that potentially reveals a trajectory in the data. Next, PILOT explores EIPLGraph to find the structure of the trajectory. An important parameter is the source_node, which indicates the start of the trajectory. Here, we selected a normal sample (node by id 2). This method returns a rank samples, which we define as a disease progression score (t = t1, ..., tn), where tl represents the ranking of the nth sample.
</div>


```python
pl.pl.fit_pricipla_graph(adata_Com,source_node=2)
```


    
![png](Combination_Kidney_IgAN_files/Combination_Kidney_IgAN_22_0.png)
    


##### Cell-type Importance Glomeruli:


```python
pl.tl.cell_importance(adata_Com,xlim=125)
```


    
![png](Combination_Kidney_IgAN_files/Combination_Kidney_IgAN_24_0.png)
    



    
![png](Combination_Kidney_IgAN_files/Combination_Kidney_IgAN_24_1.png)
    


#### Feature selection for Glomeruli based on Combination:

##### Saving morphological features and maps them with the obtained order by PILOT  (for Glomeruli):

<div class="alert alert-block alert-info"> 
This step extracts features associated with all clusters and map them with the obtained time by PILOT (based on the Trjaectory order of Samples).
    
    
* The function "extract_cells_from_pathomics"  automatically creates a cells folder and puts the extracted features associated with cells and obtained time by PILOT(orders).
</div>


```python
pl.tl.extract_cells_from_pathomics(adata_Com)
```

##### Getting the log scale of features 


```python
data=pl.tl.norm_morphological_features(column_names=['glom_sizes',
 'glom_distance_to_closest_glom','glom_diameters','glom_tuft_sizes','glom_bowman_sizes'],name_cell='All')
```

##### Morphological features changes for Glomeruli:

<div class="alert alert-block alert-info">
      
 Firstly, we should note that for pathmocis data, instead of using cluster-specific changes for features/genes (please see MI tutorial analysis), we use the fit models to find the structural changes. 
    
We apply the morphological_features_importance function to catch the features/structures that are changing over the trajectory(combination) for Glomeruli.
</div>


```python
pl.tl.morphological_features_importance(data,height=10,x_lim=160,width=20)    
```

    Name of Cluster : [1mAll[0m
    Sparsity:6.486269746268921e-05
    For this cell_type, p-value of  14 genes are statistically significant.
      Expression pattern  count
    0        linear down      6
    1          linear up      4
    2     quadratic down      2
    3       quadratic up      2
    data saved successfully



    
![png](Combination_Kidney_IgAN_files/Combination_Kidney_IgAN_31_1.png)
    


#### Feature slection for Tubuli based on Combination

#### Saving morphological features and map them with the obtained order by PILOT  (for Tubuli):



```python
adata_T.uns['orders']=adata_Com.uns['orders']
pl.tl.extract_cells_from_pathomics(adata_T)
```

#####  Getting the log scale of features 


```python
data=pl.tl.norm_morphological_features(column_names=['tubule_diameters',
 'tubule_sizes',
 'tubule_distance_to_closest_instance'],name_cell='All')
```

##### Morphological features changes for Tubuli


```python
pl.tl.morphological_features_importance(data,x_lim=160,height=5)    
```

    Name of Cluster : [1mAll[0m
    Sparsity:0.0
    For this cell_type, p-value of  3 genes are statistically significant.
      Expression pattern  count
    0        linear down      2
    1          linear up      1
    data saved successfully



    
![png](Combination_Kidney_IgAN_files/Combination_Kidney_IgAN_38_1.png)
    



```python

```

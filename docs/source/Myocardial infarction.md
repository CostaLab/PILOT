# Exploring Myocardial Infarction scRNA Data using PILOT: Unraveling Disease Trajectory and Cellular/Molecular Changes 

<div class="alert alert-block alert-info">
<b>PILOT</b>

Welcome to the PILOT Package Tutorial for scRNA Data!

Here we show the whole process for applying PILOT to scRNA data using Myocardial Infarction scRNA Data, you can download it from [here](https://costalab.ukaachen.de/open_data/PILOT/myocardial_infarction.h5ad).

</div>


```python
import PILOT as pl
import scanpy as sc
```

### Reading Anndata


```python
adata=sc.read_h5ad('Datasets/myocardial_infarction.h5ad')
```

### Loading the required information and computing the Wasserstein distance:
<div class="alert alert-block alert-info"> In order to work with PILOT, ensure that your Anndata object is loaded and contains the required information.
    
Use the following parameters to configure PILOT for your analysis (Setting Parameters):
    
adata: Pass your loaded Anndata object to PILOT.
    
emb_matrix: Provide the name of the variable in the obsm level that holds the PCA representation.
    
clusters_col: Specify the name of the column in the observation level of your Anndata that corresponds to cell types or clusters.
    
sample_col: Indicate the column name in the observation level of your Anndata that contains information about samples or patients.
    
status: Provide the column name that represents the status or disease (e.g., "control" or "case").
       
</div>


```python
pl.tl.wasserstein_distance(adata,emb_matrix='PCA',
clusters_col='cell_subtype',sample_col='sampleID',status='Status')
```

### Ploting the Cost matrix and the Wasserstein distance:
<div class="alert alert-block alert-info"> 
 Here we show the heatmaps of Cost matrix (cells) and Wasserstein distance (samples).      
</div>


```python
pl.pl.heatmaps(adata)
```


    
![png](Myocardial%20infarction_files/Myocardial%20infarction_8_0.png)
    



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_8_1.png)
    


### Trajectory:
<div class="alert alert-block alert-info"> 
 Here we show the Diffusion map of Wasserstein distance.
</div>


```python
pl.pl.trajectory(adata,colors=['Blue','red'])
```


    
![png](Myocardial%20infarction_files/Myocardial%20infarction_10_0.png)
    


###  Fit a principal graph:
<div class="alert alert-block alert-info"> 
Drawing the backbone of the trajectory with EIPLGraph.
Here the source_node is important to start ranking samples, simply you can choose the start point from control samples.  It also allows
us to rank samples with a disease progression score t = t1, ...,tn, where tl is the ranking of the sample n
</div>



```python
pl.pl.fit_pricipla_graph(adata,source_node=7)
```


    
![png](Myocardial%20infarction_files/Myocardial%20infarction_12_0.png)
    


###  Cell-type importance:
<div class="alert alert-block alert-info"> 
Here we get the critical cells that are changing over the disease progression(sorted samples based on the trajectory of PILOT with EIPLGraph).
</div>


```python
pl.tl.cell_importance(adata)
```


    
![png](Myocardial%20infarction_files/Myocardial%20infarction_14_0.png)
    



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_14_1.png)
    


# Applyin PILOT for finding Markers

### Gene selection:
<div class="alert alert-block alert-info"> 
In this step, we find genes that are changed specifically over the disease progression (order of Trajectory from Control to IZ) per specific cell. In other words, we uncover genes with different patterns.You need to reproduce whole markers for cells by running the following code. After running the code, you can  see a folder named 'Markers' that for each cell there is a folder inside that includes 'Whole_expressions.csv'. 
Whole_expressions file covers the found genes and their statistics.    
Please be patient for this part, it takes time for whole cells. Here we find the genes for 'Healthy_CM' cell type.
    
* You need to set names of columns that show cell_types/clusters and Samples/Patinets in your object.
</div>


```python
pl.tl.genes_importance(adata,name_cell=adata.uns['cellnames'][1],sample_col='sampleID',col_cell='cell_subtype')
```

    Name of Cell type : [1mhealthy_CM[0m
    sparsity:0.6932279564364608
    For this cell_type, p-value of  3206 genes are statistically significant.
               Expression pattern  count
    4    linear up quadratic down   1134
    2    linear down quadratic up    758
    0                 linear down    581
    3                   linear up    347
    6              quadratic down    216
    7                quadratic up    164
    1  linear down quadratic down      5
    5      linear up quadratic up      1
    data saved successfully



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_17_1.png)
    


### Gene Cluster Differentiation:
<div class="alert alert-block alert-info"> 
After finding genes per cell, in the next step, we use the Gene_Cluster_Differentiation function to uncover the pattern of the found genes for each cell compared to other cells. You need to set your interested cells, e.g., we use 'healthy_CM' in the following part. In other hands, we  want to see the found genes of just this cell. You can adjust it with your interest or select whole cells.

Next, PILOT picks the genes of selected cell and then finds the pattern of each gene in a distinct cell compared to other cells. 
 
* number_genes is the interest number of genes for each pattern.

* Please not you should do the previous step for whole cells types to run this part! We just do that for Healthy_CM!
 
 
</div>


```python
pl.tl.gene_cluster_differentiation(cellnames=['healthy_CM'],number_genes=1)
```


    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_1.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_3.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_5.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_7.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_9.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_11.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_13.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_15.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_17.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_19.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_21.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_23.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_25.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_27.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_29.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_31.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_33.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_35.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_37.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_39.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_41.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_43.png)
    



    <Figure size 640x480 with 0 Axes>



    
![png](Myocardial%20infarction_files/Myocardial%20infarction_19_45.png)
    



```python

```

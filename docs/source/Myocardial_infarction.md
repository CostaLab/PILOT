### Trajectory analysis of Myocardial Infarction using PILOT

<div class="alert alert-block alert-info">
<b>PILOT</b>

Welcome to the PILOT Package Tutorial for scRNA Data!

Here we show the whole process for applying PILOT to scRNA data using Myocardial Infarction scRNA Data, you can download the Anndata (h5ad) file from [here](https://costalab.ukaachen.de/open_data/PILOT/myocardial_infarction.h5ad).

</div>


```python
import PILOT as pl
import scanpy as sc
```

##### Reading Anndata


```python
adata=sc.read_h5ad('Datasets/myocardial_infarction.h5ad')
```

###### Loading the required information and computing the Wasserstein distance:
<div class="alert alert-block alert-info"> In order to work with PILOT, ensure that your Anndata object is loaded and contains the required information.
    
Use the following parameters to configure PILOT for your analysis (Setting Parameters):
    
adata: Pass your loaded Anndata object to PILOT.
    
emb_matrix: Provide the name of the variable in the obsm level that holds the dimension reduction (PCA representation).
    
clusters_col: Specify the name of the column in the observation level of your Anndata that corresponds to cell types or clusters.
    
sample_col: Indicate the column name in the observation level of your Anndata that contains information about samples or patients.
    
status: Provide the column name that represents the status or disease (e.g., "control" or "case").
       
</div>


```python
pl.tl.wasserstein_distance(adata,emb_matrix='PCA',
clusters_col='cell_subtype',sample_col='sampleID',status='Status')
```

##### Ploting the Cost matrix and the Wasserstein distance:
<div class="alert alert-block alert-info"> 
 Here we show the heatmaps of Cost matrix (cells) and Wasserstein distance (samples).      
</div>


```python
pl.pl.heatmaps(adata)
```


    
![png](Myocardial_infarction_files/Myocardial_infarction_8_0.png)
    



    
![png](Myocardial_infarction_files/Myocardial_infarction_8_1.png)
    


##### Trajectory:
<div class="alert alert-block alert-info"> 
 Here we show the Diffusion map of Wasserstein distance.
</div>


```python
pl.pl.trajectory(adata,colors=['Blue','red'])
```


    
![png](Myocardial_infarction_files/Myocardial_infarction_10_0.png)
    


#####  Fit a principal graph:
<div class="alert alert-block alert-info"> 
The difussion map creates an embeding that potentially reveals a trajectory in the data. Next, PILOT explores EIPLGraph to find the structure of the trajectory. An important parameter is the source_node, which indicate the start of the trajectory. Here, we selected a control sample. This method returns a rank samples, which we define as a disease progression score (t = t1, ..., tn), where tl represents the ranking of the nth sample.
</div>



```python
pl.pl.fit_pricipla_graph(adata,source_node=7)
```


    
![png](Myocardial_infarction_files/Myocardial_infarction_12_0.png)
    


#####  Cell-type importance:
<div class="alert alert-block alert-info"> 
Next, we can use the robust regression model to find cells whose proportions change linearly or non-linearly with disease progression. As indicated in the paper, major halmarks of MI progression are detected, i.e., a decrease of cardiomyocyte cells (CM) and an increase of fibroblasts and myeloid cells.
</div>


```python
pl.tl.cell_importance(adata,height=45,width=38,fontsize=28)
```


    
![png](Myocardial_infarction_files/Myocardial_infarction_14_0.png)
    



    
![png](Myocardial_infarction_files/Myocardial_infarction_14_1.png)
    


##### Applyin PILOT for finding Markers

##### Gene selection:
<div class="alert alert-block alert-info"> 
Given that we found interesting cell types, we would like next to investigate genes associated with these trajectories, i.e. genes, which expression changes linear or quadratically with the disease progression. After running the command, you can find a folder named ‘Markers’. There, we will have a folder for each cell type. The file ‘Whole_expressions.csv’ contains all statistics associated with genes for that cell type. Here, we run the genes_importance function for whole cell types.
    
* You need to set names of columns that show cell_types/clusters and Samples/Patinets in your object.
</div>


```python
for cell in adata.uns['cellnames']:
    pl.tl.genes_importance(adata,name_cell=cell,sample_col='sampleID',col_cell='cell_subtype',plot_genes=False)
```

    


##### Cluster Specific Marker Changes:
<div class="alert alert-block alert-info"> 
The previous test, only finds genes with significant changes over time for a given cell type. However, it does not consider if a similar pattern and expression values are found in other clusters. To further select genes, we use a Wald test that compares the fit of the gene in the cluster vs. the fit of the gene in other clusters.
In the code below, we consider top genes (regarding the regression fit) for two interesting cell types discussed in the manuscript (‘healthy CM’ and ‘Myofib’).
</div>


```python
pl.tl.gene_cluster_differentiation(cellnames=['healthy_CM','Myofib'],number_genes=70)
```

    



```python
pl.pl.plt_gene_cluster_differentiation(cellnames=['healthy_CM','Myofib'],font_size=20)
```


    <Figure size 8000x8000 with 0 Axes>



    <Figure size 8000x8000 with 0 Axes>



<div class="alert alert-block alert-info"> 
Test results are saved in ‘gene_clusters_stats_extend.csv’ and plots are saved at “plots_gene_cluster_differentiation”. To find a final list of genes, we only consider genes with a fold change higher than 0.5, i.e. genes which expression is increased in the cluster at hand; and we sort the genes based on the Wald test p-value. These can be seen bellow.
</div>


```python
df=pl.tl.results_gene_cluster_differentiation(cluster_name='Myofib').head(50)
df.head(15)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gene</th>
      <th>cluster</th>
      <th>waldStat</th>
      <th>pvalue</th>
      <th>FC</th>
      <th>Expression pattern</th>
      <th>fit-pvalue</th>
      <th>fit-mod-rsquared</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>5845</th>
      <td>RORA</td>
      <td>Myofib</td>
      <td>726.167859</td>
      <td>4.443328e-157</td>
      <td>0.940741</td>
      <td>quadratic down</td>
      <td>7.232834e-174</td>
      <td>0.587234</td>
    </tr>
    <tr>
      <th>2673</th>
      <td>GAS7</td>
      <td>Myofib</td>
      <td>405.944093</td>
      <td>1.141824e-87</td>
      <td>0.927879</td>
      <td>linear up quadratic down</td>
      <td>1.873033e-107</td>
      <td>0.570704</td>
    </tr>
    <tr>
      <th>6418</th>
      <td>SRSF11</td>
      <td>Myofib</td>
      <td>97.018608</td>
      <td>6.799134e-21</td>
      <td>0.865734</td>
      <td>linear down quadratic up</td>
      <td>2.897170e-83</td>
      <td>0.540767</td>
    </tr>
    <tr>
      <th>5027</th>
      <td>PKNOX2</td>
      <td>Myofib</td>
      <td>78.876085</td>
      <td>5.346793e-17</td>
      <td>0.855504</td>
      <td>quadratic down</td>
      <td>1.039404e-117</td>
      <td>0.544122</td>
    </tr>
    <tr>
      <th>2876</th>
      <td>GSN</td>
      <td>Myofib</td>
      <td>38.588189</td>
      <td>2.121753e-08</td>
      <td>0.633831</td>
      <td>linear up quadratic down</td>
      <td>2.942472e-279</td>
      <td>0.601684</td>
    </tr>
    <tr>
      <th>1282</th>
      <td>CHD9</td>
      <td>Myofib</td>
      <td>33.359494</td>
      <td>2.704617e-07</td>
      <td>0.566664</td>
      <td>linear up quadratic down</td>
      <td>7.658862e-77</td>
      <td>0.559604</td>
    </tr>
    <tr>
      <th>2518</th>
      <td>FN1</td>
      <td>Myofib</td>
      <td>27.507285</td>
      <td>4.608276e-06</td>
      <td>1.573680</td>
      <td>linear down quadratic up</td>
      <td>2.947389e-188</td>
      <td>0.633774</td>
    </tr>
    <tr>
      <th>1488</th>
      <td>COL6A3</td>
      <td>Myofib</td>
      <td>23.179833</td>
      <td>3.704328e-05</td>
      <td>1.069156</td>
      <td>linear down quadratic up</td>
      <td>3.514298e-172</td>
      <td>0.608543</td>
    </tr>
    <tr>
      <th>5983</th>
      <td>SEC24D</td>
      <td>Myofib</td>
      <td>18.968725</td>
      <td>2.775001e-04</td>
      <td>0.812640</td>
      <td>linear down quadratic up</td>
      <td>6.604860e-64</td>
      <td>0.522700</td>
    </tr>
    <tr>
      <th>1881</th>
      <td>DST</td>
      <td>Myofib</td>
      <td>16.720179</td>
      <td>8.068366e-04</td>
      <td>0.509683</td>
      <td>linear down quadratic up</td>
      <td>5.705082e-10</td>
      <td>0.559501</td>
    </tr>
    <tr>
      <th>3869</th>
      <td>MGP</td>
      <td>Myofib</td>
      <td>14.382736</td>
      <td>2.427875e-03</td>
      <td>0.838889</td>
      <td>quadratic down</td>
      <td>1.327779e-225</td>
      <td>0.571374</td>
    </tr>
    <tr>
      <th>1443</th>
      <td>COL3A1</td>
      <td>Myofib</td>
      <td>13.488091</td>
      <td>3.691628e-03</td>
      <td>1.240454</td>
      <td>linear down quadratic up</td>
      <td>0.000000e+00</td>
      <td>0.665616</td>
    </tr>
    <tr>
      <th>1423</th>
      <td>COL1A2</td>
      <td>Myofib</td>
      <td>13.174701</td>
      <td>4.273640e-03</td>
      <td>1.327753</td>
      <td>linear down quadratic up</td>
      <td>0.000000e+00</td>
      <td>0.655032</td>
    </tr>
    <tr>
      <th>2138</th>
      <td>EXT1</td>
      <td>Myofib</td>
      <td>13.100268</td>
      <td>4.424713e-03</td>
      <td>0.570081</td>
      <td>linear up quadratic down</td>
      <td>3.159831e-35</td>
      <td>0.555757</td>
    </tr>
    <tr>
      <th>1712</th>
      <td>DCN</td>
      <td>Myofib</td>
      <td>9.286334</td>
      <td>2.571648e-02</td>
      <td>1.137701</td>
      <td>linear up quadratic down</td>
      <td>1.866152e-284</td>
      <td>0.588602</td>
    </tr>
  </tbody>
</table>
</div>



<div class="alert alert-block alert-info"> 
Plots of genes are saved at 'plot_genes_for_Myofib' folder. We can also vizualise specfici genes, for example the ones discussed in PILOT manuscript (COL1A2, DCN and EXT1). In the plot, the orange line indicates the fit in the target cell type (shown as orange lines) compared to other cell types (represented by grey lines).
</div>


```python
pl.pl.exploring_specific_genes(cluster_name='Myofib',font_size=20,gene_list=['COL1A2','DCN','EXT1'])
```


    <Figure size 8000x8000 with 0 Axes>



    
![png](Myocardial_infarction_files/Myocardial_infarction_24_1.png)
    


<div class="alert alert-block alert-info"> 
Here is the GO enrichment for 50 the first top genes for Myofib.
</div>


```python
pl.pl.go_enrichment(df,cell_type='Myofib')
```


    
![png](Myocardial_infarction_files/Myocardial_infarction_26_0.png)
    


<div class="alert alert-block alert-info"> 
We can repeate the same analysis for healthy_CM cell type by using the following commands.
</div>


```python
df=pl.tl.results_gene_cluster_differentiation(cluster_name='healthy_CM').head(50)
df.head(15)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gene</th>
      <th>cluster</th>
      <th>waldStat</th>
      <th>pvalue</th>
      <th>FC</th>
      <th>Expression pattern</th>
      <th>fit-pvalue</th>
      <th>fit-mod-rsquared</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>7527</th>
      <td>ZSWIM6</td>
      <td>healthy_CM</td>
      <td>10732.437845</td>
      <td>0.0</td>
      <td>2.439606</td>
      <td>linear down quadratic up</td>
      <td>3.203533e-07</td>
      <td>0.539056</td>
    </tr>
    <tr>
      <th>755</th>
      <td>BRAF</td>
      <td>healthy_CM</td>
      <td>11070.052527</td>
      <td>0.0</td>
      <td>2.784995</td>
      <td>linear down quadratic up</td>
      <td>1.651067e-154</td>
      <td>0.572363</td>
    </tr>
    <tr>
      <th>6357</th>
      <td>SRPK2</td>
      <td>healthy_CM</td>
      <td>2899.114236</td>
      <td>0.0</td>
      <td>2.343306</td>
      <td>linear down quadratic up</td>
      <td>4.902958e-02</td>
      <td>0.529040</td>
    </tr>
    <tr>
      <th>6244</th>
      <td>SORBS1</td>
      <td>healthy_CM</td>
      <td>2038.349722</td>
      <td>0.0</td>
      <td>1.828163</td>
      <td>linear down quadratic up</td>
      <td>5.178918e-179</td>
      <td>0.560146</td>
    </tr>
    <tr>
      <th>1444</th>
      <td>COL4A2</td>
      <td>healthy_CM</td>
      <td>2697.352511</td>
      <td>0.0</td>
      <td>1.928972</td>
      <td>linear down quadratic up</td>
      <td>5.469598e-47</td>
      <td>0.547734</td>
    </tr>
    <tr>
      <th>1564</th>
      <td>CPNE5</td>
      <td>healthy_CM</td>
      <td>7898.219297</td>
      <td>0.0</td>
      <td>2.928602</td>
      <td>linear down quadratic up</td>
      <td>0.000000e+00</td>
      <td>0.765984</td>
    </tr>
    <tr>
      <th>1582</th>
      <td>CREB5</td>
      <td>healthy_CM</td>
      <td>14355.325817</td>
      <td>0.0</td>
      <td>3.577095</td>
      <td>linear down quadratic up</td>
      <td>2.412752e-295</td>
      <td>0.669214</td>
    </tr>
    <tr>
      <th>1659</th>
      <td>CUX1</td>
      <td>healthy_CM</td>
      <td>11836.217139</td>
      <td>0.0</td>
      <td>2.296052</td>
      <td>linear down quadratic up</td>
      <td>2.555804e-16</td>
      <td>0.584346</td>
    </tr>
    <tr>
      <th>1688</th>
      <td>DAB1</td>
      <td>healthy_CM</td>
      <td>2333.703271</td>
      <td>0.0</td>
      <td>3.148011</td>
      <td>linear down quadratic up</td>
      <td>0.000000e+00</td>
      <td>0.645854</td>
    </tr>
    <tr>
      <th>2520</th>
      <td>FNDC3A</td>
      <td>healthy_CM</td>
      <td>2999.657669</td>
      <td>0.0</td>
      <td>2.638358</td>
      <td>linear down quadratic up</td>
      <td>1.411204e-18</td>
      <td>0.541994</td>
    </tr>
    <tr>
      <th>2613</th>
      <td>FRMD5</td>
      <td>healthy_CM</td>
      <td>6497.378414</td>
      <td>0.0</td>
      <td>1.279514</td>
      <td>linear up quadratic up</td>
      <td>3.445313e-286</td>
      <td>0.543240</td>
    </tr>
    <tr>
      <th>3185</th>
      <td>KDM2A</td>
      <td>healthy_CM</td>
      <td>3498.453036</td>
      <td>0.0</td>
      <td>2.095851</td>
      <td>linear down quadratic up</td>
      <td>7.834409e-04</td>
      <td>0.519640</td>
    </tr>
    <tr>
      <th>3269</th>
      <td>L3MBTL4</td>
      <td>healthy_CM</td>
      <td>6616.669429</td>
      <td>0.0</td>
      <td>2.849918</td>
      <td>linear down quadratic up</td>
      <td>0.000000e+00</td>
      <td>0.643816</td>
    </tr>
    <tr>
      <th>3334</th>
      <td>LDLRAD4</td>
      <td>healthy_CM</td>
      <td>1950.046140</td>
      <td>0.0</td>
      <td>1.231544</td>
      <td>linear down quadratic up</td>
      <td>4.568236e-276</td>
      <td>0.567330</td>
    </tr>
    <tr>
      <th>4942</th>
      <td>PDLIM3</td>
      <td>healthy_CM</td>
      <td>3489.299461</td>
      <td>0.0</td>
      <td>3.290811</td>
      <td>linear down quadratic up</td>
      <td>2.326237e-257</td>
      <td>0.589730</td>
    </tr>
  </tbody>
</table>
</div>




```python
pl.pl.exploring_specific_genes(cluster_name='healthy_CM',gene_list=list(df['gene'][0:20]))
```


    <Figure size 8000x8000 with 0 Axes>



    
![png](Myocardial_infarction_files/Myocardial_infarction_29_1.png)
    



```python
pl.pl.go_enrichment(df,cell_type='healthy_CM')
```


    
![png](Myocardial_infarction_files/Myocardial_infarction_30_0.png)
    



```python

```

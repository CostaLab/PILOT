### Patients sub-group detection by clustering of the W disatnce of Pancreas data and ranking important cells/genes based on found sub-group



<div class="alert alert-block alert-info">
<b>PILOT</b>

In this tutorial, we are using single-cell data from pancreatic ductal adenocarcinomas to find sub-groups of patients using clustering of W distances. And then we rank the cells/genes based on clustering results.

* You can download data for this tutorial from [here](https://costalab.ukaachen.de/open_data/PILOT/PDAC.h5ad).

</div>

We first load our packages to use its functions:


```python
import PILOT as pl
import scanpy as sc
```

Load the data,


```python
adata=sc.read_h5ad('Datasets/PDAC_C.h5ad')
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
pl.tl.wasserstein_distance(adata,emb_matrix='X_pca',
clusters_col='cell_types',sample_col='sampleID',status='status')
```

### Ploting the Cost matrix and the Wasserstein distance:
<div class="alert alert-block alert-info"> 
 Here we show the heatmaps of Cost matrix (cells) and Wasserstein distance (samples).      
</div>


```python
pl.pl.heatmaps(adata)
```


    
![png](Patients_sub_group_detection_files/Patients_sub_group_detection_9_0.png)
    



    
![png](Patients_sub_group_detection_files/Patients_sub_group_detection_9_1.png)
    


### Trajectory:
<div class="alert alert-block alert-info"> 
 Here we show the Diffusion map of Wasserstein distance.
</div>


```python
pl.pl.trajectory(adata,colors=['red','Blue'])
```


    
![png](Patients_sub_group_detection_files/Patients_sub_group_detection_11_0.png)
    


### In this section, we should find the optimal number of clusters. 
<div class="alert alert-block alert-info"> 
The Silhouette Score Curve is used to find the optimal number of clusters by plotting the average Silhouette Score for different numbers of clusters. The number of clusters corresponding to the highest average Silhouette Score is considered the optimal number of clusters.
</div>



```python
pl.pl.select_best_sil(adata)
```


    
![png](Patients_sub_group_detection_files/Patients_sub_group_detection_13_0.png)
    



    
![png](Patients_sub_group_detection_files/Patients_sub_group_detection_13_1.png)
    


## Patients sub-group detection by clustering EMD. 
<div class="alert alert-block alert-info"> 
Using the Silhouette scores of the previous step, we can find the optimal number of cluster of patients to detect different stage of disease. 
</div>


```python
proportion_df=pl.pl.clustering_emd(adata,res=adata.uns['best_res'])
```

    WARNING: dendrogram data not found (using key=dendrogram_Leiden). Running `sc.tl.dendrogram` with default parameters. For fine tuning it is recommended to run `sc.tl.dendrogram` independently.



    
![png](Patients_sub_group_detection_files/Patients_sub_group_detection_15_1.png)
    



```python
proportion_df
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
      <th>Fibroblast cell</th>
      <th>Stellate cell</th>
      <th>Macrophage cell</th>
      <th>Endothelial cell</th>
      <th>T cell</th>
      <th>B cell</th>
      <th>Ductal cell type 2</th>
      <th>Endocrine cell</th>
      <th>Ductal cell type 1</th>
      <th>Acinar cell</th>
      <th>Predicted_Labels</th>
      <th>sampIeD</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>T1</th>
      <td>0.037589</td>
      <td>0.113576</td>
      <td>0.105891</td>
      <td>0.274959</td>
      <td>0.046117</td>
      <td>0.015377</td>
      <td>0.367178</td>
      <td>0.039278</td>
      <td>0.000031</td>
      <td>0.000006</td>
      <td>2</td>
      <td>T1</td>
    </tr>
    <tr>
      <th>T2</th>
      <td>0.023683</td>
      <td>0.050644</td>
      <td>0.080895</td>
      <td>0.079913</td>
      <td>0.219654</td>
      <td>0.366306</td>
      <td>0.053610</td>
      <td>0.000658</td>
      <td>0.103918</td>
      <td>0.020718</td>
      <td>0</td>
      <td>T2</td>
    </tr>
    <tr>
      <th>T3</th>
      <td>0.111619</td>
      <td>0.223975</td>
      <td>0.153370</td>
      <td>0.133641</td>
      <td>0.069095</td>
      <td>0.012913</td>
      <td>0.290039</td>
      <td>0.004557</td>
      <td>0.000786</td>
      <td>0.000006</td>
      <td>0</td>
      <td>T3</td>
    </tr>
    <tr>
      <th>T4</th>
      <td>0.075958</td>
      <td>0.093477</td>
      <td>0.100291</td>
      <td>0.050654</td>
      <td>0.184007</td>
      <td>0.073023</td>
      <td>0.402102</td>
      <td>0.009737</td>
      <td>0.007823</td>
      <td>0.002928</td>
      <td>2</td>
      <td>T4</td>
    </tr>
    <tr>
      <th>T5</th>
      <td>0.050237</td>
      <td>0.194602</td>
      <td>0.181150</td>
      <td>0.265452</td>
      <td>0.052019</td>
      <td>0.000905</td>
      <td>0.087912</td>
      <td>0.020626</td>
      <td>0.127364</td>
      <td>0.019734</td>
      <td>0</td>
      <td>T5</td>
    </tr>
    <tr>
      <th>T6</th>
      <td>0.028871</td>
      <td>0.115979</td>
      <td>0.132011</td>
      <td>0.214852</td>
      <td>0.036881</td>
      <td>0.005349</td>
      <td>0.155536</td>
      <td>0.017102</td>
      <td>0.259211</td>
      <td>0.034207</td>
      <td>0</td>
      <td>T6</td>
    </tr>
    <tr>
      <th>T7</th>
      <td>0.042859</td>
      <td>0.203452</td>
      <td>0.042852</td>
      <td>0.155289</td>
      <td>0.024106</td>
      <td>0.002689</td>
      <td>0.183404</td>
      <td>0.001341</td>
      <td>0.309202</td>
      <td>0.034807</td>
      <td>0</td>
      <td>T7</td>
    </tr>
    <tr>
      <th>T8</th>
      <td>0.012943</td>
      <td>0.054532</td>
      <td>0.121944</td>
      <td>0.030166</td>
      <td>0.077470</td>
      <td>0.017225</td>
      <td>0.674181</td>
      <td>0.000003</td>
      <td>0.010092</td>
      <td>0.001445</td>
      <td>2</td>
      <td>T8</td>
    </tr>
    <tr>
      <th>T9</th>
      <td>0.270590</td>
      <td>0.036114</td>
      <td>0.058202</td>
      <td>0.002137</td>
      <td>0.005952</td>
      <td>0.000004</td>
      <td>0.596398</td>
      <td>0.002550</td>
      <td>0.026351</td>
      <td>0.001702</td>
      <td>2</td>
      <td>T9</td>
    </tr>
    <tr>
      <th>T10</th>
      <td>0.025385</td>
      <td>0.311543</td>
      <td>0.071262</td>
      <td>0.153383</td>
      <td>0.169056</td>
      <td>0.152148</td>
      <td>0.111132</td>
      <td>0.002418</td>
      <td>0.002458</td>
      <td>0.001216</td>
      <td>0</td>
      <td>T10</td>
    </tr>
    <tr>
      <th>T11</th>
      <td>0.548986</td>
      <td>0.131443</td>
      <td>0.150856</td>
      <td>0.014013</td>
      <td>0.034693</td>
      <td>0.001276</td>
      <td>0.067799</td>
      <td>0.005093</td>
      <td>0.044884</td>
      <td>0.000957</td>
      <td>0</td>
      <td>T11</td>
    </tr>
    <tr>
      <th>T12</th>
      <td>0.032607</td>
      <td>0.135680</td>
      <td>0.204836</td>
      <td>0.114541</td>
      <td>0.173558</td>
      <td>0.265179</td>
      <td>0.004863</td>
      <td>0.000441</td>
      <td>0.059922</td>
      <td>0.008373</td>
      <td>0</td>
      <td>T12</td>
    </tr>
    <tr>
      <th>T13</th>
      <td>0.096212</td>
      <td>0.098639</td>
      <td>0.129248</td>
      <td>0.203105</td>
      <td>0.022355</td>
      <td>0.006321</td>
      <td>0.113711</td>
      <td>0.004860</td>
      <td>0.194362</td>
      <td>0.131186</td>
      <td>0</td>
      <td>T13</td>
    </tr>
    <tr>
      <th>T14</th>
      <td>0.031040</td>
      <td>0.047052</td>
      <td>0.026033</td>
      <td>0.040552</td>
      <td>0.004010</td>
      <td>0.001005</td>
      <td>0.848283</td>
      <td>0.001502</td>
      <td>0.000518</td>
      <td>0.000004</td>
      <td>2</td>
      <td>T14</td>
    </tr>
    <tr>
      <th>T15</th>
      <td>0.247942</td>
      <td>0.216757</td>
      <td>0.050618</td>
      <td>0.154397</td>
      <td>0.059305</td>
      <td>0.089464</td>
      <td>0.062386</td>
      <td>0.007158</td>
      <td>0.089478</td>
      <td>0.022496</td>
      <td>0</td>
      <td>T15</td>
    </tr>
    <tr>
      <th>T16</th>
      <td>0.208068</td>
      <td>0.215408</td>
      <td>0.163394</td>
      <td>0.178700</td>
      <td>0.020201</td>
      <td>0.010408</td>
      <td>0.093648</td>
      <td>0.007344</td>
      <td>0.066109</td>
      <td>0.036720</td>
      <td>0</td>
      <td>T16</td>
    </tr>
    <tr>
      <th>T17</th>
      <td>0.064274</td>
      <td>0.028304</td>
      <td>0.117983</td>
      <td>0.057564</td>
      <td>0.012475</td>
      <td>0.005280</td>
      <td>0.709303</td>
      <td>0.001919</td>
      <td>0.001936</td>
      <td>0.000963</td>
      <td>2</td>
      <td>T17</td>
    </tr>
    <tr>
      <th>T18</th>
      <td>0.104355</td>
      <td>0.106913</td>
      <td>0.014735</td>
      <td>0.106281</td>
      <td>0.004489</td>
      <td>0.000646</td>
      <td>0.658712</td>
      <td>0.001922</td>
      <td>0.001943</td>
      <td>0.000005</td>
      <td>2</td>
      <td>T18</td>
    </tr>
    <tr>
      <th>T19</th>
      <td>0.040320</td>
      <td>0.192341</td>
      <td>0.205322</td>
      <td>0.290049</td>
      <td>0.072087</td>
      <td>0.005811</td>
      <td>0.133247</td>
      <td>0.024939</td>
      <td>0.031442</td>
      <td>0.004444</td>
      <td>0</td>
      <td>T19</td>
    </tr>
    <tr>
      <th>T20</th>
      <td>0.000050</td>
      <td>0.058108</td>
      <td>0.051885</td>
      <td>0.186710</td>
      <td>0.012468</td>
      <td>0.000019</td>
      <td>0.553793</td>
      <td>0.035260</td>
      <td>0.099618</td>
      <td>0.002089</td>
      <td>2</td>
      <td>T20</td>
    </tr>
    <tr>
      <th>T21</th>
      <td>0.039673</td>
      <td>0.044623</td>
      <td>0.384067</td>
      <td>0.040921</td>
      <td>0.078062</td>
      <td>0.030982</td>
      <td>0.288701</td>
      <td>0.000003</td>
      <td>0.071898</td>
      <td>0.021070</td>
      <td>2</td>
      <td>T21</td>
    </tr>
    <tr>
      <th>T22</th>
      <td>0.012651</td>
      <td>0.035672</td>
      <td>0.121442</td>
      <td>0.013106</td>
      <td>0.246033</td>
      <td>0.105638</td>
      <td>0.464987</td>
      <td>0.000452</td>
      <td>0.000016</td>
      <td>0.000003</td>
      <td>2</td>
      <td>T22</td>
    </tr>
    <tr>
      <th>T23</th>
      <td>0.392302</td>
      <td>0.144848</td>
      <td>0.114833</td>
      <td>0.063183</td>
      <td>0.127046</td>
      <td>0.040140</td>
      <td>0.101926</td>
      <td>0.014310</td>
      <td>0.001060</td>
      <td>0.000352</td>
      <td>0</td>
      <td>T23</td>
    </tr>
    <tr>
      <th>T24</th>
      <td>0.175008</td>
      <td>0.162899</td>
      <td>0.047336</td>
      <td>0.274065</td>
      <td>0.059439</td>
      <td>0.000555</td>
      <td>0.123288</td>
      <td>0.018711</td>
      <td>0.104026</td>
      <td>0.034673</td>
      <td>0</td>
      <td>T24</td>
    </tr>
    <tr>
      <th>N1</th>
      <td>0.106625</td>
      <td>0.046408</td>
      <td>0.028343</td>
      <td>0.187387</td>
      <td>0.003192</td>
      <td>0.000003</td>
      <td>0.001076</td>
      <td>0.007439</td>
      <td>0.602166</td>
      <td>0.017359</td>
      <td>1</td>
      <td>N1</td>
    </tr>
    <tr>
      <th>N2</th>
      <td>0.086272</td>
      <td>0.027572</td>
      <td>0.054113</td>
      <td>0.161817</td>
      <td>0.003579</td>
      <td>0.000005</td>
      <td>0.000020</td>
      <td>0.016334</td>
      <td>0.384869</td>
      <td>0.265418</td>
      <td>1</td>
      <td>N2</td>
    </tr>
    <tr>
      <th>N3</th>
      <td>0.035202</td>
      <td>0.162609</td>
      <td>0.070341</td>
      <td>0.248312</td>
      <td>0.024191</td>
      <td>0.000020</td>
      <td>0.000086</td>
      <td>0.057123</td>
      <td>0.386722</td>
      <td>0.015394</td>
      <td>1</td>
      <td>N3</td>
    </tr>
    <tr>
      <th>N4</th>
      <td>0.026713</td>
      <td>0.017470</td>
      <td>0.000020</td>
      <td>0.207382</td>
      <td>0.000012</td>
      <td>0.000009</td>
      <td>0.001067</td>
      <td>0.007188</td>
      <td>0.740131</td>
      <td>0.000008</td>
      <td>1</td>
      <td>N4</td>
    </tr>
    <tr>
      <th>N5</th>
      <td>0.009044</td>
      <td>0.029328</td>
      <td>0.003403</td>
      <td>0.227718</td>
      <td>0.001141</td>
      <td>0.000010</td>
      <td>0.000044</td>
      <td>0.007892</td>
      <td>0.702250</td>
      <td>0.019170</td>
      <td>1</td>
      <td>N5</td>
    </tr>
    <tr>
      <th>N6</th>
      <td>0.111423</td>
      <td>0.071038</td>
      <td>0.009773</td>
      <td>0.459526</td>
      <td>0.005586</td>
      <td>0.000013</td>
      <td>0.001447</td>
      <td>0.008357</td>
      <td>0.325865</td>
      <td>0.006972</td>
      <td>1</td>
      <td>N6</td>
    </tr>
    <tr>
      <th>N7</th>
      <td>0.017955</td>
      <td>0.040369</td>
      <td>0.006294</td>
      <td>0.184749</td>
      <td>0.001804</td>
      <td>0.000008</td>
      <td>0.000932</td>
      <td>0.059184</td>
      <td>0.648346</td>
      <td>0.040358</td>
      <td>1</td>
      <td>N7</td>
    </tr>
    <tr>
      <th>N8</th>
      <td>0.117214</td>
      <td>0.088526</td>
      <td>0.082789</td>
      <td>0.294240</td>
      <td>0.006566</td>
      <td>0.000827</td>
      <td>0.000852</td>
      <td>0.035242</td>
      <td>0.308995</td>
      <td>0.064750</td>
      <td>1</td>
      <td>N8</td>
    </tr>
    <tr>
      <th>N9</th>
      <td>0.055111</td>
      <td>0.015404</td>
      <td>0.044575</td>
      <td>0.280784</td>
      <td>0.002841</td>
      <td>0.000004</td>
      <td>0.000826</td>
      <td>0.023095</td>
      <td>0.395444</td>
      <td>0.181917</td>
      <td>1</td>
      <td>N9</td>
    </tr>
    <tr>
      <th>N10</th>
      <td>0.003256</td>
      <td>0.011029</td>
      <td>0.012972</td>
      <td>0.168502</td>
      <td>0.000008</td>
      <td>0.000006</td>
      <td>0.000025</td>
      <td>0.003241</td>
      <td>0.643491</td>
      <td>0.157470</td>
      <td>1</td>
      <td>N10</td>
    </tr>
    <tr>
      <th>N11</th>
      <td>0.028233</td>
      <td>0.010143</td>
      <td>0.041976</td>
      <td>0.565064</td>
      <td>0.000009</td>
      <td>0.000007</td>
      <td>0.002922</td>
      <td>0.001449</td>
      <td>0.319083</td>
      <td>0.031115</td>
      <td>1</td>
      <td>N11</td>
    </tr>
  </tbody>
</table>
</div>



Here we can see that whole of the Normal samples are in cluster 1, so we can rename the name of clusters for future analysis. This step is optional!


```python
proportion_df.loc[proportion_df['Predicted_Labels']=='0', 'Predicted_Labels'] = 'Tumor 1'
proportion_df.loc[proportion_df['Predicted_Labels']=='1', 'Predicted_Labels'] = 'Normal'
proportion_df.loc[proportion_df['Predicted_Labels']=='2', 'Predicted_Labels'] = 'Tumor 2'
```

## Cell-type selection. 
<div class="alert alert-block alert-info"> 
Importantly, we determine which cell type derives the disease from the first stage to the second stage by detecting cell types having statistically significant changes between two sub-groups. This is done using Welch’s t-test, which is appropriate for unequal variances in two groups of samples.

Based on the adjusted p-value threshold you consider, you can choose how statistically significant the cell types you want to have.
  
</div>


```python
pl.pl.cell_type_diff_two_sub_patient_groups(proportion_df, proportion_df.columns[0:-2],
                                      group1 = 'Tumor 2', group2 = 'Tumor 1',
                                      pval_thr = 0.05, figsize = (15, 4))
```


    
![png](Patients_sub_group_detection_files/Patients_sub_group_detection_20_0.png)
    


<div class="alert alert-block alert-info"> 
Next, we can check the distribution of each patient sub-group in each cell type
</div>


```python
pl.pl.plot_cell_types_distributions(proportion_df, cell_types=['Stellate cell','Ductal cell type 2','Ductal cell type 1'],
                              figsize = (17,8),label_order=['Normal', 'Tumor 1', 'Tumor 2'],label_colors=['#FF7F00','#BCD2EE','#1874CD'])
```


    
![png](Patients_sub_group_detection_files/Patients_sub_group_detection_22_0.png)
    


<div class="alert alert-block alert-info"> 
In the statistical table generated in your results path, all cell types sorted based on the score and their adjusted p-value.
cell types with positive scores shows the most differences cell types in group 1 compared with group 2 and negative scores shows the other way. 
You can the results in 'Results_PILOT/Diff_Expressions_Results' folder.
</div>

<h2 style="font-weight: 300;font-size: 32px;"> Differential expression analysis</h2>

### Note:

<div class="alert alert-block alert-info"> 
This step needs the 'limma' package in R, You need to run this function to install it (if you already have not done)!
</div>


```python
#pl.tl.install_r_packages()
```

### Ductal cell type 1

<div class="alert alert-block alert-info"> 
In the next step, for specific cell types, we find the Differential genes between two interested patient sub-groups
Based on the fold change threshold, you can determine how much difference you want to see between two interested patients sub-groups. For the tutorial, we already saved the needed input for this function (2000 highly variable genes for each cell type). 
    
You can run this function for your own data, and it automatically extracts the gene expressions for each cell type. In case you want highly variable genes, please make the "highly_variable_genes_" parameter True.
</div>


```python
cell_type = "Ductal cell type 1" #look at the Cells folder
pl.tl.compute_diff_expressions(adata,cell_type, proportion_df,
                         fc_thr =  0.5, pval_thr = 0.05,
                         group1 = 'Tumor 1', group2 = 'Tumor 2',sample_col='sampleID',
                               col_cell='cell_types'
                            )
```

    run limma lmFit
    run limma eBayes



    
![png](Patients_sub_group_detection_files/Patients_sub_group_detection_29_1.png)
    


<h2 style="font-weight: 300;font-size: 32px;"> Gene Ontology analysis</h2>

<div class="alert alert-block alert-info">  Based on the adjusted p-value and fold change threshold, we select genes that are highly differentiated in each patient sub-groups and specify their annotation using the <a href="https://biit.cs.ut.ee/gprofiler/gost">gProfiler</a>
    
</div>


```python
pl.pl.gene_annotation_cell_type_subgroup(cell_type = cell_type, group = 'Tumor 1'
                                   ,num_gos=10)
```


    
![png](Patients_sub_group_detection_files/Patients_sub_group_detection_32_0.png)
    



```python
pl.pl.gene_annotation_cell_type_subgroup(cell_type = cell_type, group = 'Tumor 2',
                                   num_gos=10)
```


    
![png](Patients_sub_group_detection_files/Patients_sub_group_detection_33_0.png)
    


### Stellate cell


```python
cell_type = "Stellate cell" #look at the Cells folder
pl.tl.compute_diff_expressions(adata,cell_type, proportion_df,
                         fc_thr = 0.1, pval_thr = 0.05,
                         group1 = 'Tumor 1', group2 = 'Tumor 2',sample_col='sampleID',
                               col_cell='cell_types')
```

    run limma lmFit
    run limma eBayes



    
![png](Patients_sub_group_detection_files/Patients_sub_group_detection_35_1.png)
    



```python
pl.pl.gene_annotation_cell_type_subgroup(cell_type = cell_type, group = 'Tumor 1',
                                   num_gos=10)
```


    
![png](Patients_sub_group_detection_files/Patients_sub_group_detection_36_0.png)
    



```python
pl.pl.gene_annotation_cell_type_subgroup(cell_type = cell_type, group = 'Tumor 2',
                                  num_gos=10)
```


    
![png](Patients_sub_group_detection_files/Patients_sub_group_detection_37_0.png)
    


### Ductal cell type 2


```python
cell_type = "Ductal cell type 2" #look at the Cells folder
pl.tl.compute_diff_expressions(adata,cell_type, proportion_df,
                         fc_thr = 0.020, pval_thr = 0.05,
                         group1 = 'Tumor 1', group2 = 'Tumor 2',sample_col='sampleID',
                               col_cell='cell_types')
```

    run limma lmFit
    run limma eBayes



    
![png](Patients_sub_group_detection_files/Patients_sub_group_detection_39_1.png)
    



```python
pl.pl.gene_annotation_cell_type_subgroup(cell_type = cell_type, group = 'Tumor 1',
                                   num_gos=10)
```


    
![png](Patients_sub_group_detection_files/Patients_sub_group_detection_40_0.png)
    



```python
pl.pl.gene_annotation_cell_type_subgroup(cell_type = cell_type, group = 'Tumor 2',
                                   num_gos=10)
```


    
![png](Patients_sub_group_detection_files/Patients_sub_group_detection_41_0.png)
    



```python

```

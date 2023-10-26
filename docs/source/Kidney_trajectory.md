### Evaluate the presence of batch effects by PILOT

<div class="alert alert-block alert-info">
In this tutorial, we include statistical tests to evaluate association between the detected trajectories with any experimental or clinical variable provided in the Kidney dataset.
</div>


```python
import PILOT as pl
import scanpy as sc
from matplotlib.cm import get_cmap
import seaborn as sns
import pandas as pd
```

#### Reading the original Anndata (without filteration):
We consider the Kidney data without any filtering, we observe a high association with the tissue location: renal medulla, cortex of kidney, renal papilla or kidney).You can download the Anndata (h5ad) file from [here](https://costalab.ukaachen.de/open_data/PILOT/Kidney_ori.h5ad), and place it in the _Datasets_ folder.


```python
adata = sc.read_h5ad('/Datasets/Kidney_ori.h5ad')
```

##### Loading the required information and computing the Wasserstein distance:
<div class="alert alert-block alert-info"> In order to work with PILOT, ensure that your Anndata object is loaded and contains the required information.
    
Use the following parameters to configure PILOT for your analysis (Setting Parameters):
    
- adata: Pass your loaded Anndata object to PILOT.
    
- emb_matrix: Provide the name of the variable in the obsm level that holds the dimension reduction (PCA representation).
    
- clusters_col: Specify the name of the column in the observation level of your Anndata that corresponds to cell types or clusters.
    
- sample_col: Indicate the column name in the observation level of your Anndata that contains information about samples or patients.
    
- status: Provide the column name that represents the status or disease (e.g., "control" or "case").
       
</div>


```python
pl.tl.wasserstein_distance(
    adata,
    emb_matrix = 'X_pca',
    clusters_col = 'cell_type',
    sample_col = 'donor_id',
    status = 'disease'
    )
```

##### Trajectory:
<div class="alert alert-block alert-info"> 
 Here we show the Diffusion map of Wasserstein distance.  In the upcoming trajectory analysis, the labels "IZ" stands for the ischaemic zone tissue. 
</div>


```python
pl.pl.trajectory(adata, colors = ['red','orange','Blue'])
```


    
![png](Kidney_trajectory_files/Kidney_trajectory_8_0.png)
    


##### Fit a principal graph:
<div class="alert alert-block alert-info"> 
The difussion map creates an embedding that potentially reveals a trajectory in the data. Next, PILOT explores EIPLGraph to find the structure of the trajectory. An important parameter is the source_node, which indicates the start of the trajectory. Here, we selected a control sample (node with id = 7). This method returns ranked samples, which we define as a disease progression score (t = t1, ..., tn), where tl represents the ranking of the nth sample.
</div>



```python
pl.pl.fit_pricipla_graph(adata, source_node = 8)
```

##### Cell-type importance:
<div class="alert alert-block alert-info"> 
Next, we can use the robust regression model to find cells whose proportions change linearly or non-linearly with disease progression. As indicated in the paper, major hallmark of MI progression are detected, i.e., a decrease of cardiomyocyte cells (CM) and an increase of fibroblasts and myeloid cells.
</div>


```python
pl.tl.cell_importance(adata,heatmap_w = 20,height = 15,xlim = 30)
```

##### Statistical tests:
For categorical variables, this is based on ANOVA statistics on Trajectory analysis  while for numerical variables this is based on Spearman correlation. For these functions, provide the sample_col as the Sample/Patient column and your interested variables. Of note, these functions show just the significant variables (p-values) and ignore the insignificant ones.


```python
numeric = ['degen.score','aStr.score','aEpi.score','matrisome.score','collagen.score','glycoprotein.score','proteoglycan.score']
```


```python
categorical = ['BMI','hypertension','development_stage','sex','eGFR','diabetes_history','disease','tissue']
```

##### Categorical variables 


```python
pl.tl.correlation_categorical_with_trajectory(adata, sample_col = 'donor_id', features = categorical)
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
      <th>Feature</th>
      <th>ANOVA_FStatistic</th>
      <th>ANOVA_PValue</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>7</th>
      <td>tissue</td>
      <td>9.772883</td>
      <td>0.000022</td>
    </tr>
    <tr>
      <th>5</th>
      <td>diabetes_history</td>
      <td>3.429356</td>
      <td>0.038475</td>
    </tr>
  </tbody>
</table>
</div>



#####  Numerical variables : You can read the main reference of kidney data for more detailed info about these variables.


```python
pl.tl.correlation_numeric_with_trajectory(adata, sample_col = 'donor_id', features = numeric)
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
      <th>Feature</th>
      <th>Spearman_Correlation</th>
      <th>Spearman_PValue</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>degen.score</td>
      <td>0.262272</td>
      <td>0.032027</td>
    </tr>
  </tbody>
</table>
</div>



##### Visualizing Feature Distribution Within Patients sub-group 

Utilizing the 'trajectory' function, you can skillfully visualize the distribution of significant variables over the identified trajectory. Please configure the requisite columns as outlined at the outset of this tutorial.

###### Tissue


```python
palette = sns.color_palette("husl", len(adata.obs['tissue'].unique()))
colors = [palette[i] for i in range(len(adata.obs['tissue'].unique()))]
pl.tl.wasserstein_distance(
    adata,
    emb_matrix = 'X_pca',
    clusters_col = 'cell_type',
    sample_col = 'donor_id',
    status = 'tissue'
    )
```


```python
pl.pl.trajectory(adata, colors = colors)
```


    
![png](Kidney_trajectory_files/Kidney_trajectory_24_0.png)
    


#### Investigating after filtration

Here, we use the all previous steps for the filtered data set (only samples associated with the kidney (or whole kidney) location).You can download the Anndata (h5ad) file from [here](https://costalab.ukaachen.de/open_data/PILOT/Kidney_filtered.h5ad), and place it in the _Datasets_ folder.


```python
adata_filtered=sc.read_h5ad('/Datasets/Kidney_filtered.h5ad')
```


```python
pl.tl.wasserstein_distance(
    adata_filtered,
    emb_matrix  ='X_pca',
    clusters_col='cell_type',
    sample_col='donor_id',
    status='disease'
    )
```


```python
pl.pl.trajectory(adata_filtered, colors = ['red','orange','Blue'])
```


    
![png](Kidney_trajectory_files/Kidney_trajectory_29_0.png)
    



```python
pl.pl.fit_pricipla_graph(adata_filtered, source_node = 2)

```


```python
pl.tl.cell_importance(adata_filtered,heatmap_w=20,height=15,xlim=30)
```

##### Categorical variables 


```python
pl.tl.correlation_categorical_with_trajectory(adata_filtered, sample_col='donor_id', features=['BMI','hypertension','development_stage','sex','eGFR','diabetes_history','disease'])
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
      <th>Feature</th>
      <th>ANOVA_FStatistic</th>
      <th>ANOVA_PValue</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>6</th>
      <td>disease</td>
      <td>9.457675</td>
      <td>0.000566</td>
    </tr>
    <tr>
      <th>5</th>
      <td>diabetes_history</td>
      <td>12.517843</td>
      <td>0.001189</td>
    </tr>
    <tr>
      <th>0</th>
      <td>BMI</td>
      <td>3.969681</td>
      <td>0.016358</td>
    </tr>
  </tbody>
</table>
</div>



###### Diabetes_history


```python
palette = sns.color_palette("husl", len(adata_filtered.obs['diabetes_history'].unique()))
colors = [palette[i] for i in range(len(adata_filtered.obs['diabetes_history'].unique()))]
pl.tl.wasserstein_distance(
    adata_filtered,
    emb_matrix='X_pca',
    clusters_col='cell_type',
    sample_col='donor_id',
    status='diabetes_history'
    )
```


```python
pl.pl.trajectory(adata_filtered, colors = colors)
```


    
![png](Kidney_trajectory_files/Kidney_trajectory_36_0.png)
    


###### BMI


```python
palette = sns.color_palette("husl", len(adata_filtered.obs['BMI'].unique()))
colors = [palette[i] for i in range(len(adata_filtered.obs['BMI'].unique()))]
pl.tl.wasserstein_distance(
    adata_filtered,
    emb_matrix='X_pca',
    clusters_col='cell_type',
    sample_col='donor_id',
    status='BMI'
    )
```


```python
pl.pl.trajectory(adata_filtered, colors = colors)
```


    
![png](Kidney_trajectory_files/Kidney_trajectory_39_0.png)
    


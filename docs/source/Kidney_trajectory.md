### Evaluate the presence of sample level batch effects by PILOT (Trajectory)

<div class="alert alert-block alert-info">
In this tutorial, we will demonstrate how to use statistical tests to evaluate the potential association between the detected trajectories with any experimental or clinical variable present in a data set. This example will be based on the kidney single cell data.
</div>


```python
import PILOT as pl
import scanpy as sc
from matplotlib.cm import get_cmap
import seaborn as sns
import pandas as pd
```

#### Reading the original Anndata:
First, we consider the original kidney single cell data. You can download the Anndata (h5ad) file from [here](https://costalab.ukaachen.de/open_data/PILOT/Kidney_ori.h5ad), and place it in the _Datasets_ folder.


```python
adata = sc.read_h5ad('/Datasets/kidney_ori.h5ad')
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

##### Trajectory estimation:
<div class="alert alert-block alert-info"> 
As a next step, we estimate a diffusion map from the Wasserstein distance
</div>


```python
pl.pl.trajectory(adata, colors = ['red','orange','Blue'])
```


    
![png](Kidney_trajectory_files/Kidney_trajectory_8_0.png)
    


Then estimate the trajectory by fitting a principle graph and calculte the order of samples by 'cell_importance' functions. See main [tutorial](https://pilot.readthedocs.io/en/latest/Myocardial_infarction.html), for details. 


```python
#pl.pl.fit_pricipla_graph(adata, source_node = 8)
#pl.tl.cell_importance(adata,heatmap_w = 20,height = 15,xlim = 30)
```

##### Evaluation of the association of estimated disease progression with experimental factor:
A very important question is if PILOT analysis is affected by experimental artefacts (location of tissues, batch processing). To evaluate this, we use ANOVA statistics and Spearman correlation to check any association between any variables describing the disease cohort and the predicted variables. 

To run these functions, provide the sample_col as the Sample/Patient column and your interested variables. Of note, these functions show just the significant variables (p-value < 0.05).


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



#####  Numerical variables : Similarly, you can do the same analysis for numerical variables. 


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



We observe that there is a clear association of the trajectory with the origin of the biopsies. This is not surprising as PILOT uses cellular composition information, as samples at distinct regions do have distinct cells. Regarding continuous variable, we find an association with the degen.score, which is a gene signature score estimated from the comparisons of controls and disease samples (see data manuscript for more details https://doi.org/10.1038/s41586-023-05769-3). 

To double check these results, one can of course plot the trajectory by showing the location as labels. 

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


    
![png](Kidney_trajectory_files/Kidney_trajectory_21_0.png)
    


#### Filtering of samples

We focus therefore only on biopsies from sample locationn 'kidney', which were used in our benchmarking. You can download the Anndata (h5ad) file from [here](https://costalab.ukaachen.de/open_data/PILOT/Kidney_filtered.h5ad), and place it in the _Datasets_ folder.


```python
#adata_filtered=sc.read_h5ad('/Datasets/Kidney_filtered.h5ad')
adata_filtered=sc.read_h5ad('/data/scRNA/For_Mina/batch_PILOT/filter_data/Kidney.h5ad')
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


    
![png](Kidney_trajectory_files/Kidney_trajectory_26_0.png)
    


##### Categorical variables 

We next recheck the association of variables with the estimated trajectory. 


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



We observed data the true label (disease) has highest association followed by BMI or Diabetes History. While the latter are two clinical factors with a impact on kidney disease, as closer investigation of the variables indicates these association as spurious, as only controls have BMI values and most disease are diabetic. 

##### Diabetes_history


```python
contingency_table = pd.crosstab(adata_filtered.obs['disease'], adata_filtered.obs['diabetes_history'],values=adata_filtered.obs['donor_id'], aggfunc=pd.Series.nunique, margins=True, margins_name='Total Unique Patients')
contingency_table
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
      <th>diabetes_history</th>
      <th>No</th>
      <th>Yes</th>
      <th>Total Unique Patients</th>
    </tr>
    <tr>
      <th>disease</th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Normal</th>
      <td>18.0</td>
      <td>NaN</td>
      <td>18</td>
    </tr>
    <tr>
      <th>acute kidney failure</th>
      <td>NaN</td>
      <td>5.0</td>
      <td>5</td>
    </tr>
    <tr>
      <th>chronic kidney disease</th>
      <td>NaN</td>
      <td>13.0</td>
      <td>13</td>
    </tr>
    <tr>
      <th>Total Unique Patients</th>
      <td>18.0</td>
      <td>18.0</td>
      <td>36</td>
    </tr>
  </tbody>
</table>
</div>



##### BMI


```python
contingency_table = pd.crosstab(adata_filtered.obs['disease'], adata_filtered.obs['BMI'],values=adata_filtered.obs['donor_id'], aggfunc=pd.Series.nunique, margins=True, margins_name='Total Unique Patients')
contingency_table
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
      <th>BMI</th>
      <th>18.5—24.9</th>
      <th>25.0—29.9</th>
      <th>30.0 and Above</th>
      <th>unknown</th>
      <th>Total Unique Patients</th>
    </tr>
    <tr>
      <th>disease</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Normal</th>
      <td>2.0</td>
      <td>11.0</td>
      <td>4.0</td>
      <td>1.0</td>
      <td>18</td>
    </tr>
    <tr>
      <th>acute kidney failure</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>5.0</td>
      <td>5</td>
    </tr>
    <tr>
      <th>chronic kidney disease</th>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>13.0</td>
      <td>13</td>
    </tr>
    <tr>
      <th>Total Unique Patients</th>
      <td>2.0</td>
      <td>11.0</td>
      <td>4.0</td>
      <td>19.0</td>
      <td>36</td>
    </tr>
  </tbody>
</table>
</div>



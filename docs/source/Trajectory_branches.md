### Analysing branches of the trajectory

<div class="alert alert-block alert-info">
<b>PILOT</b>

A trajectory can have several parts or branches, e.g. one branch towards a disease sub-type and another one towards a different sub-type, or a healthy part and a disease part. The cell types that change can be different in each branch. In this tutorial we show how to split the trajectory by choosing nodes of the principal graph, and how to run the cell type importance (and the rest of the analysis) for each branch separately.

We use two examples: first the Myocardial Infarction data (a healthy part and a disease part), then, at the end, a lung cancer atlas with two disease branches that start from the normal samples (normal to COPD and normal to lung adenocarcinoma), together with the clinical variables of each branch. For the Myocardial Infarction data, you can download the Anndata (h5ad) file from [here](https://costalab.ukaachen.de/open_data/PILOT/myocardial_infarction.h5ad), and place it in the _Datasets_ folder. The first steps are the same as in the [Myocardial Infarction tutorial](https://pilot.readthedocs.io/en/latest/Myocardial_infarction.html).
</div>

```python
import pilotpy as pl
import scanpy as sc
import matplotlib.pyplot as plt
```

##### Reading Anndata

```python
adata = sc.read_h5ad('Datasets/myocardial_infarction.h5ad')
```

##### Wasserstein distance and trajectory:
<div class="alert alert-block alert-info">
We compute the Wasserstein distances between the samples and the trajectory, as in the Myocardial Infarction tutorial.
</div>

```python
pl.tl.wasserstein_distance(
    adata,
    emb_matrix = 'PCA',
    clusters_col = 'cell_subtype',
    sample_col = 'sampleID',
    status = 'Status',
    )
```

```python
pl.pl.trajectory(adata)
```

![png](Trajectory_branches_files/cell07_1.png)

##### Principal graph and its nodes:
<div class="alert alert-block alert-info">
PILOT fits a principal graph (a tree) through the samples and orders them along it, starting from the root (red star). The numbers in the plot are the node numbers. We use them to choose the branches.

Note: the node numbers can change when the graph is fitted again (e.g. in a new Python session or with other settings), so always look at your own plot before choosing the nodes.
</div>

```python
pl.pl.fit_pricipla_graph(adata, source_node = 'auto', root_status = 'control')
```

```
Root of the trajectory: node 6 (end node where the control samples come first).
```

![png](Trajectory_branches_files/cell09_1.png)

##### Which samples belong to which node:
<div class="alert alert-block alert-info">
Each sample belongs to its nearest node of the graph. `pl.tl.sample_nodes` shows the node, the status and the pseudotime of every sample, which helps to choose the nodes of a branch.
</div>

```python
pl.tl.sample_nodes(adata)
```

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>sampleID</th>
      <th>node</th>
      <th>status</th>
      <th>pseudotime</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>RZ_P9</td>
      <td>0</td>
      <td>control</td>
      <td>0.141687</td>
    </tr>
    <tr>
      <th>1</th>
      <td>control_P8</td>
      <td>1</td>
      <td>control</td>
      <td>0.759073</td>
    </tr>
    <tr>
      <th>2</th>
      <td>RZ_BZ_P12</td>
      <td>1</td>
      <td>control</td>
      <td>0.773940</td>
    </tr>
    <tr>
      <th>3</th>
      <td>IZ_P16</td>
      <td>2</td>
      <td>IZ</td>
      <td>1.617639</td>
    </tr>
    <tr>
      <th>4</th>
      <td>RZ_FZ_P5</td>
      <td>4</td>
      <td>control</td>
      <td>0.811143</td>
    </tr>
    <tr>
      <th>5</th>
      <td>control_P1</td>
      <td>4</td>
      <td>control</td>
      <td>0.828181</td>
    </tr>
    <tr>
      <th>6</th>
      <td>RZ_P6</td>
      <td>5</td>
      <td>control</td>
      <td>0.672575</td>
    </tr>
    <tr>
      <th>7</th>
      <td>control_P7</td>
      <td>6</td>
      <td>control</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>8</th>
      <td>GT_IZ_P13</td>
      <td>7</td>
      <td>IZ</td>
      <td>1.534083</td>
    </tr>
    <tr>
      <th>9</th>
      <td>RZ_P11</td>
      <td>8</td>
      <td>control</td>
      <td>0.410313</td>
    </tr>
    <tr>
      <th>10</th>
      <td>control_P17</td>
      <td>8</td>
      <td>control</td>
      <td>0.445190</td>
    </tr>
    <tr>
      <th>11</th>
      <td>GT_IZ_P9_rep2</td>
      <td>9</td>
      <td>IZ</td>
      <td>1.183608</td>
    </tr>
    <tr>
      <th>12</th>
      <td>IZ_P10</td>
      <td>10</td>
      <td>IZ</td>
      <td>1.081853</td>
    </tr>
    <tr>
      <th>13</th>
      <td>IZ_P15</td>
      <td>12</td>
      <td>IZ</td>
      <td>1.726663</td>
    </tr>
    <tr>
      <th>14</th>
      <td>GT_IZ_P9</td>
      <td>13</td>
      <td>IZ</td>
      <td>0.927182</td>
    </tr>
    <tr>
      <th>15</th>
      <td>RZ_P3</td>
      <td>14</td>
      <td>control</td>
      <td>0.990906</td>
    </tr>
    <tr>
      <th>16</th>
      <td>RZ_BZ_P3</td>
      <td>15</td>
      <td>control</td>
      <td>0.567245</td>
    </tr>
    <tr>
      <th>17</th>
      <td>RZ_BZ_P2</td>
      <td>15</td>
      <td>control</td>
      <td>0.584196</td>
    </tr>
    <tr>
      <th>18</th>
      <td>RZ_GT_P2</td>
      <td>16</td>
      <td>control</td>
      <td>0.506837</td>
    </tr>
    <tr>
      <th>19</th>
      <td>GT_IZ_P15</td>
      <td>18</td>
      <td>IZ</td>
      <td>1.309857</td>
    </tr>
  </tbody>
</table>

##### Branch 1, the healthy part:
<div class="alert alert-block alert-info">
`pl.tl.select_branch` keeps only the samples of the given nodes. All other samples are ignored in the next steps, so you can also use it to leave out samples you do not want to analyse. Parameters:

- `nodes`: the node numbers of the branch (any order), here [6, 0, 17, 11, 8, 16, 15, 5, 1, 4]: from the root to the end of the control samples.
- `name`: the name of the branch, used in the plot titles and in the names of the result folders.

The Wasserstein distances, the trajectory and the graph are not changed. `pl.pl.principal_graph` now shows the branch in orange; the other samples are faded.
</div>

```python
pl.tl.select_branch(adata, nodes = [6, 0, 17, 11, 8, 16, 15, 5, 1, 4], name = 'healthy part')
pl.pl.principal_graph(adata)
```

```
Branch 'healthy part': 12 of 20 samples (control: 12, IZ: 0).
```

![png](Trajectory_branches_files/cell13_1.png)

##### Cell type importance of the healthy part:
<div class="alert alert-block alert-info">
The next steps work as usual, but only use the samples of the branch, ordered along the trajectory. The results are saved in sub-folders named after the branch (e.g. 'Results_PILOT/plots/Cell_type_Report/branch_healthy_part').
</div>

```python
pl.tl.cell_importance_calculation(adata)
```

```python
pl.pl.heatmap_proportions(adata)
```

![png](Trajectory_branches_files/cell16_1.png)

```python
pl.pl.cell_importance_grid(adata)
```

![png](Trajectory_branches_files/cell17_1.png)

##### Branch 2, the disease part:
<div class="alert alert-block alert-info">
Now we choose the second branch, here the nodes [1, 4, 13, 14, 10, 9, 19, 18, 3, 7, 2, 12]: from the end of the healthy part to the end of the trajectory. We repeat the analysis. Nodes can be part of both branches (here the last two nodes of the healthy part).
</div>

```python
pl.tl.select_branch(adata, nodes = [1, 4, 13, 14, 10, 9, 19, 18, 3, 7, 2, 12], name = 'disease part')
pl.pl.principal_graph(adata)
```

```
Branch 'disease part': 12 of 20 samples (IZ: 7, control: 5).
```

![png](Trajectory_branches_files/cell19_1.png)

```python
pl.tl.cell_importance_calculation(adata)
```

```python
pl.pl.heatmap_proportions(adata)
```

![png](Trajectory_branches_files/cell21_1.png)

```python
pl.pl.cell_importance_grid(adata)
```

![png](Trajectory_branches_files/cell22_1.png)

##### Both branches side by side:
<div class="alert alert-block alert-info">
With `ax`, the principal graphs of the branches can be drawn next to each other.
</div>

```python
fig, axes = plt.subplots(1, 2, figsize = (16, 7))
pl.tl.select_branch(adata, nodes = [6, 0, 17, 11, 8, 16, 15, 5, 1, 4], name = 'healthy part')
pl.pl.principal_graph(adata, ax = axes[0])
pl.tl.select_branch(adata, nodes = [1, 4, 13, 14, 10, 9, 19, 18, 3, 7, 2, 12], name = 'disease part')
pl.pl.principal_graph(adata, ax = axes[1])
plt.tight_layout()
plt.show()
```

```
Branch 'healthy part': 12 of 20 samples (control: 12, IZ: 0).
Branch 'disease part': 12 of 20 samples (IZ: 7, control: 5).
```

![png](Trajectory_branches_files/cell24_1.png)

##### Back to all samples:
<div class="alert alert-block alert-info">
`nodes = None` removes the selection, so all samples are used again (for example before the gene analysis of the Myocardial Infarction tutorial).
</div>

```python
pl.tl.select_branch(adata, nodes = None)
pl.tl.cell_importance_calculation(adata)
```

```
No branch selected: all samples are used.
```

##### Summary:
<div class="alert alert-block alert-info">
Splitting the trajectory shows which cell types change in which part of the disease progression. In this data, only a few cell types change in the healthy part, while in the disease part the healthy cardiomyocytes (`healthy_CM`) decrease and myofibroblasts (`Myofib`), `SPP1_Macrophages` and proliferating cells (`prolif`) increase.

The same works for any trajectory: look at the principal graph, choose the nodes of each branch (or leave out the nodes of samples you do not want to analyse), and run the analysis for each branch.
</div>

### Second example: lung cancer atlas

<div class="alert alert-block alert-info">
<b>Lung cancer atlas</b>

The second example is the lung cancer single-cell atlas (Salcher et al., 2022; data described in [Zenodo 15575593](https://zenodo.org/records/15575593), atlas in [Zenodo 7227571](https://zenodo.org/records/7227571)). We use the lung tissue samples measured with 10x: 941,504 cells of 165 donors (51 normal, 18 chronic obstructive pulmonary disease (COPD), 71 lung adenocarcinoma, 12 squamous cell lung carcinoma, 13 non-small cell lung carcinoma), in 12 cell types (`ann_coarse`).

Here the trajectory does not go in one direction: one branch goes from the normal samples towards COPD and another one towards lung adenocarcinoma. We analyse both branches, and we use the clinical variables (2.1.2) to see which variables, and which batch, follow each branch.

Download the data from the links above and place it in `Datasets/` (here as `lung_atlas.h5ad`).
</div>

```python
import pandas as pd

adata_lung = sc.read_h5ad('Datasets/lung_atlas.h5ad')
```

##### Wasserstein distance:
<div class="alert alert-block alert-info">
The data comes from 11 studies, so we use the batch-corrected scVI embedding of the atlas (`X_scVI`) instead of the PCA. The samples are the donors (`donor_id`), the status is the disease.
</div>

```python
pl.tl.wasserstein_distance(
    adata_lung,
    emb_matrix = 'X_scVI',
    clusters_col = 'ann_coarse',
    sample_col = 'donor_id',
    status = 'disease',
    )
```

##### Trajectory:
<div class="alert alert-block alert-info">
With many samples, `n_neighbors` (the number of neighbor samples used for the trajectory, default 64) changes the shape. With 32 neighbors the three groups are separated: lung adenocarcinoma on the left, normal and COPD samples on the right, and a third arm on top.
</div>

```python
pl.pl.trajectory(adata_lung, n_neighbors = 32)
```

![png](Trajectory_branches_files/cell34_1.png)

##### Batch on the trajectory:
<div class="alert alert-block alert-info">
Before we interpret the branches, we check the batch. `pl.pl.trajectory_variable` colors the samples by any column, here the study. All COPD samples and 28 of the 51 normal samples come from one study (Adams_Kaminski_2020); they form the right arm. So the difference between this arm and the rest is partly a study (batch) effect. Many normal samples of the other studies are in the top arm, together with tumor samples of the same studies.
</div>

```python
pl.pl.trajectory_variable(adata_lung, 'study', show_graph = False)
```

![png](Trajectory_branches_files/cell36_1.png)

##### Principal graph:
<div class="alert alert-block alert-info">
The root is the end node with the most normal samples. The numbers are the node numbers; we use them for the branches (look at your own plot, the numbers can change between sessions).
</div>

```python
pl.pl.fit_pricipla_graph(adata_lung, source_node = 'auto', root_status = 'normal')
```

```
Root of the trajectory: node 1 (end node where the control samples come first).
```

![png](Trajectory_branches_files/cell38_1.png)

```python
nodes = pl.tl.sample_nodes(adata_lung)
pd.crosstab(nodes['node'], nodes['status'])
```

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>status</th>
      <th>chronic obstructive pulmonary disease</th>
      <th>lung adenocarcinoma</th>
      <th>squamous cell lung carcinoma</th>
      <th>non-small cell lung carcinoma</th>
      <th>normal</th>
    </tr>
    <tr>
      <th>node</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>20</td>
    </tr>
    <tr>
      <th>2</th>
      <td>0</td>
      <td>2</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0</td>
      <td>15</td>
      <td>1</td>
      <td>6</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0</td>
      <td>9</td>
      <td>4</td>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>5</th>
      <td>4</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>3</td>
    </tr>
    <tr>
      <th>6</th>
      <td>0</td>
      <td>6</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>7</th>
      <td>1</td>
      <td>8</td>
      <td>2</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>8</th>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>1</td>
      <td>2</td>
    </tr>
    <tr>
      <th>9</th>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>4</td>
    </tr>
    <tr>
      <th>10</th>
      <td>1</td>
      <td>5</td>
      <td>1</td>
      <td>0</td>
      <td>9</td>
    </tr>
    <tr>
      <th>11</th>
      <td>0</td>
      <td>2</td>
      <td>0</td>
      <td>2</td>
      <td>0</td>
    </tr>
    <tr>
      <th>12</th>
      <td>0</td>
      <td>2</td>
      <td>1</td>
      <td>0</td>
      <td>3</td>
    </tr>
    <tr>
      <th>13</th>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>14</th>
      <td>0</td>
      <td>9</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
    </tr>
    <tr>
      <th>15</th>
      <td>2</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>2</td>
    </tr>
    <tr>
      <th>16</th>
      <td>4</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>17</th>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>18</th>
      <td>0</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
    </tr>
    <tr>
      <th>19</th>
      <td>0</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>2</td>
    </tr>
  </tbody>
</table>

##### Clinical variables on the whole trajectory:
<div class="alert alert-block alert-info">
`pl.tl.correlation_with_trajectory` tests all sample variables at once (see the [Kidney tutorial](https://pilot.readthedocs.io/en/latest/Kidney_trajectory.html)). The study comes first: the trajectory is strongly linked to the batch, as the plot above shows.
</div>

```python
pl.tl.cell_importance_calculation(adata_lung)
res = pl.tl.correlation_with_trajectory(adata_lung)
res.head(10)
```

```
Tested 17 variables (12 categorical, 5 numeric); skipped 34 columns.
  skipped (differs between the cells of a sample (cell-level)): ann_fine, cell_type_predicted, leiden, ann_coarse, cell_type_tumor, cell_type_major, cell_type_neutro, cell_type_neutro_coarse, ...
  skipped (only one value): doublet_status, platform, suspension_type, organism, tissue, organism_ontology_term_id, tissue_ontology_term_id
  skipped (same groups as 'assay'): assay_ontology_term_id
  skipped (same groups as 'development_stage'): development_stage_ontology_term_id
  skipped (same groups as 'disease'): accLabel, disease_ontology_term_id
  skipped (same groups as 'sex'): sex_ontology_term_id
  skipped (too many different values (165 for 165 samples, e.g. an ID)): sample, accSample
  skipped (too many missing values (146 of 165)): EGFR_mutation, TP53_mutation, ALK_mutation, BRAF_mutation, ERBB2_mutation, KRAS_mutation, ROS_mutation
  skipped (too many missing values (165 of 165)): self_reported_ethnicity, self_reported_ethnicity_ontology_term_id
```

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Feature</th>
      <th>Type</th>
      <th>Level</th>
      <th>Test</th>
      <th>Statistic</th>
      <th>P_value</th>
      <th>Adjusted_P_value</th>
      <th>Kruskal_P_value</th>
      <th>N_samples</th>
      <th>N_groups</th>
      <th>Note</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>study</td>
      <td>categorical</td>
      <td>sample</td>
      <td>ANOVA</td>
      <td>21.118918</td>
      <td>2.578121e-24</td>
      <td>4.382805e-23</td>
      <td>5.841455e-16</td>
      <td>165</td>
      <td>11</td>
      <td></td>
    </tr>
    <tr>
      <th>1</th>
      <td>dataset</td>
      <td>categorical</td>
      <td>sample (most frequent)</td>
      <td>ANOVA</td>
      <td>14.563837</td>
      <td>5.533410e-23</td>
      <td>4.703398e-22</td>
      <td>3.542170e-14</td>
      <td>164</td>
      <td>17</td>
      <td>1 group(s) with &lt; 2 samples not tested</td>
    </tr>
    <tr>
      <th>2</th>
      <td>origin</td>
      <td>categorical</td>
      <td>sample (most frequent)</td>
      <td>ANOVA</td>
      <td>45.644463</td>
      <td>2.111309e-21</td>
      <td>1.196409e-20</td>
      <td>3.009417e-16</td>
      <td>165</td>
      <td>4</td>
      <td></td>
    </tr>
    <tr>
      <th>3</th>
      <td>tumor_stage</td>
      <td>categorical</td>
      <td>sample</td>
      <td>ANOVA</td>
      <td>63.220490</td>
      <td>5.448963e-21</td>
      <td>2.315809e-20</td>
      <td>2.512663e-16</td>
      <td>164</td>
      <td>3</td>
      <td></td>
    </tr>
    <tr>
      <th>4</th>
      <td>disease</td>
      <td>categorical</td>
      <td>sample</td>
      <td>ANOVA</td>
      <td>32.245659</td>
      <td>1.058288e-19</td>
      <td>3.598179e-19</td>
      <td>4.789446e-15</td>
      <td>165</td>
      <td>5</td>
      <td></td>
    </tr>
    <tr>
      <th>5</th>
      <td>uicc_stage</td>
      <td>categorical</td>
      <td>sample</td>
      <td>ANOVA</td>
      <td>31.564839</td>
      <td>2.388143e-19</td>
      <td>6.766405e-19</td>
      <td>7.487959e-15</td>
      <td>164</td>
      <td>5</td>
      <td></td>
    </tr>
    <tr>
      <th>6</th>
      <td>origin_fine</td>
      <td>categorical</td>
      <td>sample (most frequent)</td>
      <td>ANOVA</td>
      <td>17.112955</td>
      <td>7.580832e-13</td>
      <td>1.841059e-12</td>
      <td>6.453232e-10</td>
      <td>130</td>
      <td>6</td>
      <td></td>
    </tr>
    <tr>
      <th>7</th>
      <td>n_genes_by_counts</td>
      <td>numeric</td>
      <td>cell (mean)</td>
      <td>Spearman</td>
      <td>-0.498926</td>
      <td>9.070189e-12</td>
      <td>1.927415e-11</td>
      <td>NaN</td>
      <td>165</td>
      <td>165</td>
      <td></td>
    </tr>
    <tr>
      <th>8</th>
      <td>age</td>
      <td>numeric</td>
      <td>sample</td>
      <td>Spearman</td>
      <td>0.472102</td>
      <td>1.140865e-07</td>
      <td>2.154966e-07</td>
      <td>NaN</td>
      <td>114</td>
      <td>48</td>
      <td></td>
    </tr>
    <tr>
      <th>9</th>
      <td>pct_counts_mito</td>
      <td>numeric</td>
      <td>cell (mean)</td>
      <td>Spearman</td>
      <td>0.360428</td>
      <td>1.980500e-06</td>
      <td>3.366850e-06</td>
      <td>NaN</td>
      <td>165</td>
      <td>36</td>
      <td></td>
    </tr>
  </tbody>
</table>

```python
pl.pl.trajectory_associations(adata_lung, top = 6)
```

![png](Trajectory_branches_files/cell42_1.png)

##### Branch 1, normal to COPD:
<div class="alert alert-block alert-info">
The nodes [1, 9, 5, 8, 16, 15]: from the root (normal) along the right arm to the COPD samples. Along this branch the macrophages/monocytes decrease, while T cells, NK cells, plasma cells and mast cells increase.
</div>

```python
pl.tl.select_branch(adata_lung, nodes = [1, 9, 5, 8, 16, 15], name = 'normal to COPD')
pl.pl.principal_graph(adata_lung)
```

```
Branch 'normal to COPD': 53 of 165 samples (normal: 31, chronic obstructive pulmonary disease: 14, lung adenocarcinoma: 4, squamous cell lung carcinoma: 2, non-small cell lung carcinoma: 2).
```

![png](Trajectory_branches_files/cell44_1.png)

```python
pl.tl.cell_importance_calculation(adata_lung)
```

```python
pl.pl.cell_importance_grid(adata_lung)
```

![png](Trajectory_branches_files/cell46_1.png)

##### Clinical variables of branch 1:
<div class="alert alert-block alert-info">
The same tests, now only with the samples of the branch. Besides the disease, smoking (`ever_smoker`) follows this branch: COPD is mainly caused by smoking. The study is skipped here because it has the same groups as `dataset`.
</div>

```python
res_copd = pl.tl.correlation_with_trajectory(adata_lung)
res_copd.head(10)
```

```
Tested 13 variables (10 categorical, 3 numeric); skipped 38 columns.
  skipped (differs between the cells of a sample (cell-level)): ann_fine, cell_type_predicted, leiden, ann_coarse, cell_type_tumor, cell_type_major, cell_type_neutro, cell_type_neutro_coarse, ...
  skipped (fewer than 2 groups with >= 2 samples): total_counts_mito, pct_counts_mito
  skipped (only one value): doublet_status, platform, suspension_type, organism, tissue, organism_ontology_term_id, tissue_ontology_term_id
  skipped (same groups as 'assay'): assay_ontology_term_id
  skipped (same groups as 'dataset'): study
  skipped (same groups as 'disease'): accLabel, disease_ontology_term_id
  skipped (same groups as 'sex'): sex_ontology_term_id
  skipped (too many different values (32 for 48 samples, e.g. an ID)): development_stage, development_stage_ontology_term_id
  skipped (too many different values (53 for 53 samples, e.g. an ID)): sample, accSample
  skipped (too many missing values (51 of 53)): EGFR_mutation, TP53_mutation, ALK_mutation, BRAF_mutation, ERBB2_mutation, KRAS_mutation, ROS_mutation
  skipped (too many missing values (53 of 53)): self_reported_ethnicity, self_reported_ethnicity_ontology_term_id
```

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Feature</th>
      <th>Type</th>
      <th>Level</th>
      <th>Test</th>
      <th>Statistic</th>
      <th>P_value</th>
      <th>Adjusted_P_value</th>
      <th>Kruskal_P_value</th>
      <th>N_samples</th>
      <th>N_groups</th>
      <th>Note</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>disease</td>
      <td>categorical</td>
      <td>sample</td>
      <td>ANOVA</td>
      <td>5.122053</td>
      <td>0.001613</td>
      <td>0.020974</td>
      <td>0.003677</td>
      <td>53</td>
      <td>5</td>
      <td></td>
    </tr>
    <tr>
      <th>1</th>
      <td>ever_smoker</td>
      <td>categorical</td>
      <td>sample</td>
      <td>ANOVA</td>
      <td>8.099187</td>
      <td>0.006543</td>
      <td>0.042528</td>
      <td>0.008793</td>
      <td>49</td>
      <td>2</td>
      <td></td>
    </tr>
    <tr>
      <th>2</th>
      <td>origin</td>
      <td>categorical</td>
      <td>sample (most frequent)</td>
      <td>ANOVA</td>
      <td>4.396237</td>
      <td>0.017537</td>
      <td>0.046094</td>
      <td>0.021217</td>
      <td>52</td>
      <td>3</td>
      <td>1 group(s) with &lt; 2 samples not tested</td>
    </tr>
    <tr>
      <th>3</th>
      <td>dataset</td>
      <td>categorical</td>
      <td>sample (most frequent)</td>
      <td>ANOVA</td>
      <td>3.740497</td>
      <td>0.017690</td>
      <td>0.046094</td>
      <td>0.025798</td>
      <td>48</td>
      <td>4</td>
      <td>5 group(s) with &lt; 2 samples not tested</td>
    </tr>
    <tr>
      <th>4</th>
      <td>origin_fine</td>
      <td>categorical</td>
      <td>sample (most frequent)</td>
      <td>ANOVA</td>
      <td>5.811754</td>
      <td>0.019972</td>
      <td>0.046094</td>
      <td>0.025347</td>
      <td>48</td>
      <td>2</td>
      <td>2 group(s) with &lt; 2 samples not tested</td>
    </tr>
    <tr>
      <th>5</th>
      <td>tumor_stage</td>
      <td>categorical</td>
      <td>sample</td>
      <td>ANOVA</td>
      <td>4.162573</td>
      <td>0.021274</td>
      <td>0.046094</td>
      <td>0.024449</td>
      <td>53</td>
      <td>3</td>
      <td></td>
    </tr>
    <tr>
      <th>6</th>
      <td>uicc_stage</td>
      <td>categorical</td>
      <td>sample</td>
      <td>ANOVA</td>
      <td>3.135819</td>
      <td>0.033892</td>
      <td>0.062941</td>
      <td>0.039493</td>
      <td>52</td>
      <td>4</td>
      <td>1 group(s) with &lt; 2 samples not tested</td>
    </tr>
    <tr>
      <th>7</th>
      <td>age</td>
      <td>numeric</td>
      <td>sample</td>
      <td>Spearman</td>
      <td>0.283207</td>
      <td>0.051113</td>
      <td>0.083058</td>
      <td>NaN</td>
      <td>48</td>
      <td>32</td>
      <td></td>
    </tr>
    <tr>
      <th>8</th>
      <td>is_primary_data</td>
      <td>categorical</td>
      <td>sample</td>
      <td>ANOVA</td>
      <td>3.154673</td>
      <td>0.081676</td>
      <td>0.117976</td>
      <td>0.081781</td>
      <td>53</td>
      <td>2</td>
      <td></td>
    </tr>
    <tr>
      <th>9</th>
      <td>sex</td>
      <td>categorical</td>
      <td>sample</td>
      <td>ANOVA</td>
      <td>2.633029</td>
      <td>0.110951</td>
      <td>0.144236</td>
      <td>0.110280</td>
      <td>52</td>
      <td>2</td>
      <td></td>
    </tr>
  </tbody>
</table>

```python
pl.pl.trajectory_associations(adata_lung, top = 6)
```

![png](Trajectory_branches_files/cell49_1.png)

##### Branch 2, normal to lung adenocarcinoma:
<div class="alert alert-block alert-info">
The nodes [1, 0, 17, 11, 19, 18, 6, 4, 7, 14, 3]: the normal samples at the root, then the centre of the graph and the left arm to the lung adenocarcinoma samples. The COPD nodes are left out, so the branch goes directly from normal to cancer. Here the T cells increase strongly and the macrophages/monocytes decrease; the B cells also increase. In the clinical variables, the tumor stage comes first, but the study is still strongly linked.
</div>

```python
pl.tl.select_branch(adata_lung, nodes = [1, 0, 17, 11, 19, 18, 6, 4, 7, 14, 3], name = 'normal to adenocarcinoma')
pl.pl.principal_graph(adata_lung)
```

```
Branch 'normal to adenocarcinoma': 107 of 165 samples (lung adenocarcinoma: 57, normal: 28, non-small cell lung carcinoma: 9, squamous cell lung carcinoma: 7, chronic obstructive pulmonary disease: 6).
```

![png](Trajectory_branches_files/cell51_1.png)

```python
pl.tl.cell_importance_calculation(adata_lung)
```

```python
pl.pl.cell_importance_grid(adata_lung)
```

![png](Trajectory_branches_files/cell53_1.png)

```python
res_luad = pl.tl.correlation_with_trajectory(adata_lung)
res_luad.head(10)
```

```
Tested 16 variables (11 categorical, 5 numeric); skipped 35 columns.
  skipped (differs between the cells of a sample (cell-level)): ann_fine, cell_type_predicted, leiden, ann_coarse, cell_type_tumor, cell_type_major, cell_type_neutro, cell_type_neutro_coarse, ...
  skipped (only one value): doublet_status, platform, suspension_type, organism, tissue, organism_ontology_term_id, tissue_ontology_term_id
  skipped (same groups as 'assay'): assay_ontology_term_id
  skipped (same groups as 'disease'): accLabel, disease_ontology_term_id
  skipped (same groups as 'sex'): sex_ontology_term_id
  skipped (too many different values (107 for 107 samples, e.g. an ID)): sample, accSample
  skipped (too many different values (39 for 73 samples, e.g. an ID)): development_stage, development_stage_ontology_term_id
  skipped (too many missing values (107 of 107)): self_reported_ethnicity, self_reported_ethnicity_ontology_term_id
  skipped (too many missing values (96 of 107)): EGFR_mutation, TP53_mutation, ALK_mutation, BRAF_mutation, ERBB2_mutation, KRAS_mutation, ROS_mutation
```

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Feature</th>
      <th>Type</th>
      <th>Level</th>
      <th>Test</th>
      <th>Statistic</th>
      <th>P_value</th>
      <th>Adjusted_P_value</th>
      <th>Kruskal_P_value</th>
      <th>N_samples</th>
      <th>N_groups</th>
      <th>Note</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>tumor_stage</td>
      <td>categorical</td>
      <td>sample</td>
      <td>ANOVA</td>
      <td>48.781292</td>
      <td>1.245221e-15</td>
      <td>1.545800e-14</td>
      <td>6.876270e-12</td>
      <td>106</td>
      <td>3</td>
      <td></td>
    </tr>
    <tr>
      <th>1</th>
      <td>study</td>
      <td>categorical</td>
      <td>sample</td>
      <td>ANOVA</td>
      <td>16.783528</td>
      <td>1.932250e-15</td>
      <td>1.545800e-14</td>
      <td>2.617695e-10</td>
      <td>107</td>
      <td>9</td>
      <td></td>
    </tr>
    <tr>
      <th>2</th>
      <td>dataset</td>
      <td>categorical</td>
      <td>sample (most frequent)</td>
      <td>ANOVA</td>
      <td>11.871177</td>
      <td>4.416711e-15</td>
      <td>2.355579e-14</td>
      <td>4.719579e-09</td>
      <td>106</td>
      <td>15</td>
      <td>1 group(s) with &lt; 2 samples not tested</td>
    </tr>
    <tr>
      <th>3</th>
      <td>origin</td>
      <td>categorical</td>
      <td>sample (most frequent)</td>
      <td>ANOVA</td>
      <td>31.790541</td>
      <td>1.252541e-14</td>
      <td>5.010162e-14</td>
      <td>4.984414e-11</td>
      <td>107</td>
      <td>4</td>
      <td></td>
    </tr>
    <tr>
      <th>4</th>
      <td>uicc_stage</td>
      <td>categorical</td>
      <td>sample</td>
      <td>ANOVA</td>
      <td>24.949905</td>
      <td>2.214768e-14</td>
      <td>7.087258e-14</td>
      <td>1.070309e-10</td>
      <td>106</td>
      <td>5</td>
      <td></td>
    </tr>
    <tr>
      <th>5</th>
      <td>disease</td>
      <td>categorical</td>
      <td>sample</td>
      <td>ANOVA</td>
      <td>24.456651</td>
      <td>3.309169e-14</td>
      <td>8.824450e-14</td>
      <td>1.452318e-10</td>
      <td>107</td>
      <td>5</td>
      <td></td>
    </tr>
    <tr>
      <th>6</th>
      <td>n_genes_by_counts</td>
      <td>numeric</td>
      <td>cell (mean)</td>
      <td>Spearman</td>
      <td>-0.529713</td>
      <td>4.441277e-09</td>
      <td>1.015149e-08</td>
      <td>NaN</td>
      <td>107</td>
      <td>107</td>
      <td></td>
    </tr>
    <tr>
      <th>7</th>
      <td>origin_fine</td>
      <td>categorical</td>
      <td>sample (most frequent)</td>
      <td>ANOVA</td>
      <td>11.973986</td>
      <td>2.332729e-08</td>
      <td>4.665459e-08</td>
      <td>1.298342e-06</td>
      <td>75</td>
      <td>6</td>
      <td></td>
    </tr>
    <tr>
      <th>8</th>
      <td>age</td>
      <td>numeric</td>
      <td>sample</td>
      <td>Spearman</td>
      <td>0.507435</td>
      <td>4.608664e-06</td>
      <td>8.193180e-06</td>
      <td>NaN</td>
      <td>73</td>
      <td>39</td>
      <td></td>
    </tr>
    <tr>
      <th>9</th>
      <td>ever_smoker</td>
      <td>categorical</td>
      <td>sample</td>
      <td>ANOVA</td>
      <td>11.443066</td>
      <td>1.043643e-03</td>
      <td>1.669829e-03</td>
      <td>1.112707e-03</td>
      <td>97</td>
      <td>2</td>
      <td></td>
    </tr>
  </tbody>
</table>

```python
pl.pl.trajectory_associations(adata_lung, top = 6)
```

![png](Trajectory_branches_files/cell55_1.png)

```python
fig, axes = plt.subplots(1, 2, figsize = (16, 7))
pl.tl.select_branch(adata_lung, nodes = [1, 9, 5, 8, 16, 15], name = 'normal to COPD')
pl.pl.principal_graph(adata_lung, ax = axes[0])
pl.tl.select_branch(adata_lung, nodes = [1, 0, 17, 11, 19, 18, 6, 4, 7, 14, 3], name = 'normal to adenocarcinoma')
pl.pl.principal_graph(adata_lung, ax = axes[1])
plt.show()
```

```
Branch 'normal to COPD': 53 of 165 samples (normal: 31, chronic obstructive pulmonary disease: 14, lung adenocarcinoma: 4, squamous cell lung carcinoma: 2, non-small cell lung carcinoma: 2).
Branch 'normal to adenocarcinoma': 107 of 165 samples (lung adenocarcinoma: 57, normal: 28, non-small cell lung carcinoma: 9, squamous cell lung carcinoma: 7, chronic obstructive pulmonary disease: 6).
```

![png](Trajectory_branches_files/cell56_1.png)

##### Summary of the lung example:
<div class="alert alert-block alert-info">
The lung trajectory has two disease branches that start from the normal samples: one towards COPD and one towards lung adenocarcinoma, with different cell types changing in each branch. The clinical variables show what else follows each branch. Here the study is strongly linked to the trajectory: the COPD branch comes almost only from one study, so part of its changes can be a batch effect. Always check the batch (e.g. with `pl.pl.trajectory_variable` or `pl.tl.correlation_with_trajectory`) before interpreting a branch.
</div>

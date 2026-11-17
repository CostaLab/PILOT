### Analysing branches of the trajectory

<div class="alert alert-block alert-info">
<b>PILOT</b>

A trajectory can have several parts or branches, e.g. one branch towards a disease sub-type and another one towards a different sub-type, or a healthy part and a disease part. The cell types that change can be different in each branch. In this tutorial we show how to split the trajectory by choosing nodes of the principal graph, and how to run the cell type importance (and the rest of the analysis) for each branch separately.

We use the Myocardial Infarction data. You can download the Anndata (h5ad) file from [here](https://costalab.ukaachen.de/open_data/PILOT/myocardial_infarction.h5ad), and place it in the _Datasets_ folder. The first steps are the same as in the [Myocardial Infarction tutorial](https://pilot.readthedocs.io/en/latest/Myocardial_infarction.html).
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


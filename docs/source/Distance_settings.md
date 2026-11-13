### Settings of the Wasserstein distance

<div class="alert alert-block alert-info">
<b>PILOT</b>

In this tutorial we show the settings of `pl.tl.wasserstein_distance` and their effect on the trajectory, the UMAP of the samples and the computation time, using the Myocardial Infarction data. You can download the Anndata (h5ad) file from [here](https://costalab.ukaachen.de/open_data/PILOT/myocardial_infarction.h5ad), and place it in the _Datasets_ folder.

PILOT describes every sample by its cell type proportions. The Wasserstein distance between two samples is the cheapest way to turn the proportions of one sample into those of the other. Two choices define it:

1. <b>The cost matrix</b>: how expensive it is to move cells from one cell type to another, i.e. how different two cell types are (`cost_method`, `metric`, `n_pcs`, `centroid`).
2. <b>The transport method</b>: how the cells are moved between two samples (`regularized`).

The default settings are the ones of the PILOT paper and work well for this data. The other settings are useful for other data sets, e.g. with many cell types, rare cell types or strong batch effects.
</div>

```python
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import pilotpy as pl
from sklearn.metrics import silhouette_score
```

##### Reading Anndata

```python
adata = sc.read_h5ad('Datasets/myocardial_infarction.h5ad')
```

##### A small helper to compare the settings:
<div class="alert alert-block alert-info">
For each setting, the helper below computes the Wasserstein distance, measures the time, and shows the trajectory and the UMAP of the samples colored by status side by side.

It also computes the Silhouette score of the status on the distances: 1 means that control and IZ samples are perfectly separated, 0 that they are mixed. It is a rough check only. A setting should not be chosen only because it separates the known groups best; the goal of PILOT is to find the disease progression, also within the groups.
</div>

```python
results = []

def run_setting(name, **settings):
    start = time.perf_counter()
    pl.tl.wasserstein_distance(
        adata,
        emb_matrix = 'PCA',
        clusters_col = 'cell_subtype',
        sample_col = 'sampleID',
        status = 'Status',
        **settings)
    seconds = time.perf_counter() - start
    distances = (adata.uns['EMD'] + adata.uns['EMD'].T) / 2
    np.fill_diagonal(distances, 0)
    silhouette = silhouette_score(distances, adata.uns['real_labels'], metric = 'precomputed')
    results.append({'setting': name, 'seconds': round(seconds, 1), 'silhouette': round(silhouette, 3)})

    fig, axes = plt.subplots(1, 2, figsize = (13, 5))
    pl.pl.trajectory(adata, ax = axes[0], plot_titel = f'Trajectory: {name}')
    pl.pl.umap_samples(adata, ax = axes[1], title = f'UMAP: {name}')
    plt.tight_layout()
    plt.show()
```

##### Default settings:
<div class="alert alert-block alert-info">
The default cost is the distance between the centers (median) of the cell types in the PCA space, with the cosine metric, and exact optimal transport.
</div>

```python
run_setting('default')
```

![png](Distance_settings_files/cell08_1.png)

##### Number of principal components (`n_pcs`):
<div class="alert alert-block alert-info">
`n_pcs` uses only the first principal components to compute the cost. Few components keep only the main differences between the cell types; more components also keep the finer ones. By default all components of `emb_matrix` are used (30 in this data).
</div>

```python
for n_pcs in [2, 5, 10, 25, 30]:
    run_setting(f'n_pcs = {n_pcs}', n_pcs = n_pcs)
```

![png](Distance_settings_files/cell10_1.png)

![png](Distance_settings_files/cell10_2.png)

![png](Distance_settings_files/cell10_3.png)

![png](Distance_settings_files/cell10_4.png)

![png](Distance_settings_files/cell10_5.png)

##### Center of a cell type (`centroid`):
<div class="alert alert-block alert-info">
With the default `cost_method = 'centroid'`, each cell type is represented by one point: the `'median'` (default, robust to outlier cells) or the `'mean'` of its cells.
</div>

```python
run_setting("centroid = 'mean'", centroid = 'mean')
```

![png](Distance_settings_files/cell12_1.png)

##### Metric (`metric`):
<div class="alert alert-block alert-info">
The metric measures the distance between two cell types. All possible values:

- `'cosine'` (default): angle between the profiles, ignores their size.
- `'correlation'`: 1 - Pearson correlation (cosine after centering).
- `'euclidean'`: straight-line distance; also compares the size.
- `'sqeuclidean'`: squared euclidean; far cell types cost much more.
- `'cityblock'` (Manhattan): sum of the absolute differences.
- `'chebyshev'`: largest difference in one dimension.
- `'minkowski'`: generalisation of euclidean and cityblock.
- `'canberra'`: sum of relative differences; sensitive to small values.
- `'braycurtis'`: relative sum of differences; meant for non-negative values.
- `'seuclidean'`: euclidean with every dimension scaled by its variance.
- `'mahalanobis'`: euclidean that takes correlated dimensions into account.
- `'jensenshannon'`: for non-negative values that are like proportions, e.g. pathomics features or `cost_method = 'expression'`; not for PCA, so it is not shown here.
</div>

```python
for metric in ['correlation', 'euclidean', 'sqeuclidean', 'cityblock', 'chebyshev', 'minkowski',
               'canberra', 'braycurtis', 'seuclidean', 'mahalanobis']:
    run_setting(f"metric = '{metric}'", metric = metric)
```

![png](Distance_settings_files/cell14_1.png)

![png](Distance_settings_files/cell14_2.png)

![png](Distance_settings_files/cell14_3.png)

![png](Distance_settings_files/cell14_4.png)

![png](Distance_settings_files/cell14_5.png)

![png](Distance_settings_files/cell14_6.png)

![png](Distance_settings_files/cell14_7.png)

![png](Distance_settings_files/cell14_8.png)

![png](Distance_settings_files/cell14_9.png)

![png](Distance_settings_files/cell14_10.png)

##### Other ways to compute the cost (`cost_method`):
<div class="alert alert-block alert-info">
Besides the distance between the centers of the cell types (`'centroid'`, default), PILOT can compute the cost in other ways. All of them give a small cost for similar cell types and a large cost for different ones, so they can be used for the Wasserstein distance:

- `'distribution'`: compares all cells of two cell types, not only their centers (energy distance), so the spread and the sub-groups of a cell type count too. `max_cells` (default 500) random cells per cell type are used.
- `'graph'`: two cell types are close when many of their cells are neighbors in the k-nearest-neighbor graph of the cells (like PAGA), also through intermediate cell types (e.g. healthy, then damaged cardiomyocytes). `n_neighbors_cost` (default 15) sets the number of neighbors.
- `'expression'`: distance between the average gene expression (pseudobulk) of the cell types on the `n_genes` (default 2000) most variable genes. It needs counts in `adata.X` (or in `layer`).
- `'jaccard'`: Jaccard distance between the marker genes of the cell types: 1 - (shared markers / all markers). `n_markers` (default 100) markers per cell type are used. It needs counts in `adata.X` (or in `layer`).
- `'precomputed'`: your own cost matrix, given with `cost`.
</div>

```python
run_setting("cost_method = 'distribution'", cost_method = 'distribution', metric = 'euclidean')
```

![png](Distance_settings_files/cell16_1.png)

```python
run_setting("cost_method = 'graph'", cost_method = 'graph')
```

![png](Distance_settings_files/cell17_1.png)

```python
run_setting("cost_method = 'expression'", cost_method = 'expression', metric = 'correlation')
```

![png](Distance_settings_files/cell18_1.png)

```python
run_setting("cost_method = 'jaccard'", cost_method = 'jaccard')
```

![png](Distance_settings_files/cell19_1.png)

As an example of your own cost matrix, we use the cell type hierarchy of the data: two subtypes of the same cell type (e.g. two fibroblast subtypes) cost 0.5, two different cell types cost 1. Any matrix with the cell types as rows and columns, 0 on the diagonal, can be used.

```python
major = adata.obs.groupby('cell_subtype', observed = True)['cell_type'].agg(lambda x: x.value_counts().index[0])
cost = pd.DataFrame([[0 if a == b else (0.5 if major[a] == major[b] else 1) for b in major.index]
                     for a in major.index], index = major.index, columns = major.index)
run_setting('own cost: cell type hierarchy', cost_method = 'precomputed', cost = cost)
```

![png](Distance_settings_files/cell21_1.png)

##### Transport method (`regularized`):
<div class="alert alert-block alert-info">
- `'unreg'` (default): exact optimal transport.
- `'sinkhorn'`: entropic (Sinkhorn) optimal transport with strength `reg` (default 0.1). It is smoother, but the distance of a sample to itself is not 0.
- `'sinkhorn_divergence'`: Sinkhorn, corrected so that the distance of a sample to itself is 0.
- `'unbalanced'`: the proportions do not have to be matched exactly; mass that is not moved is penalised with strength `reg_m` (default 1). Smaller values make the distance less sensitive to rare cell types and outliers; large values give almost the same result as `'unreg'`.
</div>

```python
run_setting("regularized = 'sinkhorn_divergence'", regularized = 'sinkhorn_divergence')
```

![png](Distance_settings_files/cell23_1.png)

```python
for reg_m in [0.1, 1]:
    run_setting(f"regularized = 'unbalanced', reg_m = {reg_m}", regularized = 'unbalanced', reg_m = reg_m)
```

![png](Distance_settings_files/cell24_1.png)

![png](Distance_settings_files/cell24_2.png)

##### Summary:
<div class="alert alert-block alert-info">
The table shows the computation time of the Wasserstein distance and the Silhouette score of the status for every setting.
</div>

```python
summary = pd.DataFrame(results)
summary
```

<div>
<table border="1" class="dataframe dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>setting</th>
      <th>seconds</th>
      <th>silhouette</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>default</td>
      <td>0.3</td>
      <td>0.492</td>
    </tr>
    <tr>
      <th>1</th>
      <td>n_pcs = 2</td>
      <td>0.2</td>
      <td>0.572</td>
    </tr>
    <tr>
      <th>2</th>
      <td>n_pcs = 5</td>
      <td>0.2</td>
      <td>0.513</td>
    </tr>
    <tr>
      <th>3</th>
      <td>n_pcs = 10</td>
      <td>0.2</td>
      <td>0.512</td>
    </tr>
    <tr>
      <th>4</th>
      <td>n_pcs = 25</td>
      <td>0.3</td>
      <td>0.493</td>
    </tr>
    <tr>
      <th>5</th>
      <td>n_pcs = 30</td>
      <td>0.3</td>
      <td>0.492</td>
    </tr>
    <tr>
      <th>6</th>
      <td>centroid = 'mean'</td>
      <td>0.2</td>
      <td>0.484</td>
    </tr>
    <tr>
      <th>7</th>
      <td>metric = 'correlation'</td>
      <td>0.3</td>
      <td>0.492</td>
    </tr>
    <tr>
      <th>8</th>
      <td>metric = 'euclidean'</td>
      <td>0.3</td>
      <td>0.462</td>
    </tr>
    <tr>
      <th>9</th>
      <td>metric = 'sqeuclidean'</td>
      <td>0.3</td>
      <td>0.525</td>
    </tr>
    <tr>
      <th>10</th>
      <td>metric = 'cityblock'</td>
      <td>0.3</td>
      <td>0.448</td>
    </tr>
    <tr>
      <th>11</th>
      <td>metric = 'chebyshev'</td>
      <td>0.3</td>
      <td>0.492</td>
    </tr>
    <tr>
      <th>12</th>
      <td>metric = 'minkowski'</td>
      <td>0.3</td>
      <td>0.462</td>
    </tr>
    <tr>
      <th>13</th>
      <td>metric = 'canberra'</td>
      <td>0.3</td>
      <td>0.391</td>
    </tr>
    <tr>
      <th>14</th>
      <td>metric = 'braycurtis'</td>
      <td>0.3</td>
      <td>0.453</td>
    </tr>
    <tr>
      <th>15</th>
      <td>metric = 'seuclidean'</td>
      <td>0.3</td>
      <td>0.440</td>
    </tr>
    <tr>
      <th>16</th>
      <td>metric = 'mahalanobis'</td>
      <td>0.3</td>
      <td>0.376</td>
    </tr>
    <tr>
      <th>17</th>
      <td>cost_method = 'distribution'</td>
      <td>1.6</td>
      <td>0.402</td>
    </tr>
    <tr>
      <th>18</th>
      <td>cost_method = 'graph'</td>
      <td>7.9</td>
      <td>0.345</td>
    </tr>
    <tr>
      <th>19</th>
      <td>cost_method = 'expression'</td>
      <td>3.5</td>
      <td>0.394</td>
    </tr>
    <tr>
      <th>20</th>
      <td>cost_method = 'jaccard'</td>
      <td>3.5</td>
      <td>0.385</td>
    </tr>
    <tr>
      <th>21</th>
      <td>own cost: cell type hierarchy</td>
      <td>0.2</td>
      <td>0.405</td>
    </tr>
    <tr>
      <th>22</th>
      <td>regularized = 'sinkhorn_divergence'</td>
      <td>1.1</td>
      <td>0.506</td>
    </tr>
    <tr>
      <th>23</th>
      <td>regularized = 'unbalanced', reg_m = 0.1</td>
      <td>5.4</td>
      <td>0.599</td>
    </tr>
    <tr>
      <th>24</th>
      <td>regularized = 'unbalanced', reg_m = 1</td>
      <td>5.3</td>
      <td>0.582</td>
    </tr>
  </tbody>
</table>
</div>

##### Conclusion:
<div class="alert alert-block alert-info">
For this data, all settings place the control samples at one end and the IZ samples at the other end of the trajectory. The default settings are a good choice here, and they are also among the fastest.

- Fewer principal components (`n_pcs`) or unbalanced transport separate control and IZ a bit more.
- `'distribution'`, `'graph'`, `'expression'` and `'jaccard'` take a few seconds and give a lower score here. They are useful when cell types have sub-groups, when the PCA is affected by batch effects, or to use gene expression directly.
- The computation times are for 115,517 cells, 33 cell types and 20 samples.
</div>


### Export and import the PILOT results with MuData

<div class="alert alert-block alert-info">
<b>MuData</b>

[MuData](https://mudata.readthedocs.io) keeps several AnnData objects ("modalities") in one object and one file (`.h5mu`). PILOT can

- **export** its results as a MuData with two modalities: `cells` (your single-cell data) and `samples` (one row per patient with the PILOT results), with `pl.tl.to_mudata`;
- **import** such a file again and continue the analysis without computing the distances again, with `pl.tl.from_mudata`;
- **use a MuData as input**, e.g. a multi-modal object with an `rna` modality, in `pl.tl.wasserstein_distance`.

The `samples` modality is a normal AnnData at the patient level, so it can be used directly by other tools (scanpy, scikit-learn, ...).

Here we use the Myocardial Infarction data of the [Myocardial Infarction tutorial](https://pilot.readthedocs.io/en/latest/Myocardial_infarction.html) as an example; the same code works for any data. MuData is an optional package: `pip install mudata`.
</div>

```python
import pilotpy as pl
import scanpy as sc
import mudata as md
```

##### PILOT analysis (as in the Myocardial Infarction tutorial)

```python
adata = sc.read_h5ad('Datasets/myocardial_infarction.h5ad')
pl.tl.wasserstein_distance(
    adata,
    emb_matrix = 'PCA',
    clusters_col = 'cell_subtype',
    sample_col = 'sampleID',
    status = 'Status',
    )
pl.pl.trajectory(adata)
```

![png](MuData_files/cell04_1.png)

```python
pl.pl.fit_pricipla_graph(adata, source_node = 'auto', root_status = 'control')
```

```
Root of the trajectory: node 7 (end node where the control samples come first).
```

![png](MuData_files/cell05_1.png)

```python
pl.tl.cell_importance_calculation(adata)
```

##### Export:
<div class="alert alert-block alert-info">
`pl.tl.to_mudata` puts the results in a MuData object and, with `path`, also writes it to a `.h5mu` file. The `samples` modality has one row per sample:

- `X`: the cell type proportions (samples x cell types);
- `obs`: the status and the clinical variables, the pseudotime, the position on the trajectory and the principal graph node (and the clusters, if computed);
- `obsm`: the trajectory (`X_trajectory`) and the sample UMAP (`X_umap`, if computed);
- `obsp['wasserstein']`: the Wasserstein distances between the samples;
- `varp['cost']`: the cost matrix between the cell types;
- `uns['pilot']`: the settings and the principal graph.
</div>

```python
mdata = pl.tl.to_mudata(adata, path = 'Results_PILOT/MI_PILOT.h5mu')
mdata
```

```
MuData object with n_obs × n_vars = 115537 × 28441
  uns:	'pilot_modalities'
  2 modalities
    cells:	115517 × 28408
      obs:	'orig_ident', 'nCount_RNA', 'nFeature_RNA', 'percent_mt', 'doublet_score', 'doublet', 'dissociation_s1', 'opt_clust', 'patient', 'batch', 'opt_clust_integrated', 'cell_type', 'ident', 'n_genes', 'cell_subtype', 'cell_subtype_available', 'cell_subtype2', 'sampleID', 'patient_group', 'Status'
      obsm:	'HARMONY', 'PCA', 'UMAP_HARMONY'
      layers:	'logcounts'
    samples:	20 × 33
      obs:	'status', 'doublet', 'patient', 'batch', 'cell_subtype_available', 'patient_group', 'Status', 'pseudotime', 'graph_node', 'trajectory_order'
      uns:	'pilot'
      obsm:	'X_trajectory'
      obsp:	'wasserstein'
      varp:	'cost'
```

```python
mdata['samples'].obs.head()
```

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>status</th>
      <th>doublet</th>
      <th>patient</th>
      <th>batch</th>
      <th>cell_subtype_available</th>
      <th>patient_group</th>
      <th>Status</th>
      <th>pseudotime</th>
      <th>graph_node</th>
      <th>trajectory_order</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>control_P1</th>
      <td>control</td>
      <td>singlet</td>
      <td>P1</td>
      <td>A</td>
      <td>yes</td>
      <td>myogenic</td>
      <td>control</td>
      <td>0.823557</td>
      <td>13</td>
      <td>12</td>
    </tr>
    <tr>
      <th>RZ_FZ_P5</th>
      <td>control</td>
      <td>singlet</td>
      <td>P5</td>
      <td>A</td>
      <td>yes</td>
      <td>myogenic</td>
      <td>control</td>
      <td>0.810914</td>
      <td>13</td>
      <td>11</td>
    </tr>
    <tr>
      <th>RZ_BZ_P3</th>
      <td>control</td>
      <td>singlet</td>
      <td>P3</td>
      <td>A</td>
      <td>yes</td>
      <td>myogenic</td>
      <td>control</td>
      <td>0.567274</td>
      <td>15</td>
      <td>6</td>
    </tr>
    <tr>
      <th>RZ_BZ_P2</th>
      <td>control</td>
      <td>singlet</td>
      <td>P2</td>
      <td>A</td>
      <td>yes</td>
      <td>myogenic</td>
      <td>control</td>
      <td>0.584191</td>
      <td>15</td>
      <td>7</td>
    </tr>
    <tr>
      <th>RZ_P3</th>
      <td>control</td>
      <td>singlet</td>
      <td>P3</td>
      <td>A</td>
      <td>yes</td>
      <td>myogenic</td>
      <td>control</td>
      <td>0.988766</td>
      <td>14</td>
      <td>14</td>
    </tr>
  </tbody>
</table>

```python
mdata['samples'].obsp['wasserstein'].shape
```

```
(20, 20)
```

##### Only the patient level:
<div class="alert alert-block alert-info">
With `include_cells = False` only the `samples` modality is kept. The file is then very small and can be shared easily (no single-cell data in it). You can also get the patient-level AnnData directly with `pl.tl.sample_adata(adata)`.
</div>

```python
samples_only = pl.tl.to_mudata(adata, include_cells = False, path = 'Results_PILOT/MI_PILOT_samples.h5mu')
samples_only
```

```
MuData object with n_obs × n_vars = 20 × 33
  uns:	'pilot_modalities'
  1 modality
    samples:	20 × 33
      obs:	'status', 'doublet', 'patient', 'batch', 'cell_subtype_available', 'patient_group', 'Status', 'pseudotime', 'graph_node', 'trajectory_order'
      uns:	'pilot'
      obsm:	'X_trajectory'
      obsp:	'wasserstein'
      varp:	'cost'
```

##### Import:
<div class="alert alert-block alert-info">
`pl.tl.from_mudata` reads the file and returns the cells with the PILOT results restored (Wasserstein distances, trajectory, principal graph, pseudotime). The next steps work directly; for example the trajectory is not computed again.
</div>

```python
adata_loaded = pl.tl.from_mudata('Results_PILOT/MI_PILOT.h5mu')
pl.pl.trajectory(adata_loaded)
```

![png](MuData_files/cell14_1.png)

```python
pl.pl.principal_graph(adata_loaded)
```

![png](MuData_files/cell15_1.png)

##### A MuData as input:
<div class="alert alert-block alert-info">
If your data is already a MuData (e.g. RNA and protein of CITE-seq, or RNA and ATAC), give it to `pl.tl.wasserstein_distance` and choose the modality with the cells with `modality`. Clinical variables in the global `mdata.obs` are copied to the modality. The results are stored in that modality, so continue with `mdata['rna']`.
</div>

```python
rna = sc.read_h5ad('Datasets/myocardial_infarction.h5ad')
mdata_input = md.MuData({'rna': rna})

pl.tl.wasserstein_distance(
    mdata_input,
    modality = 'rna',
    emb_matrix = 'PCA',
    clusters_col = 'cell_subtype',
    sample_col = 'sampleID',
    status = 'Status',
    )
adata_rna = mdata_input['rna']
pl.pl.trajectory(adata_rna)
```

![png](MuData_files/cell17_1.png)


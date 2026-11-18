```{eval-rst}
.. module:: pilotpy
```

```{eval-rst}
.. automodule:: pilotpy
   :noindex:
```

# API

Import pilotpy as:

```
import pilotpy
```

## tools: `tl`

```{eval-rst}
.. module:: pilotpy.tl
```

```{eval-rst}
.. currentmodule:: pilotpy
```


### Trajectory & Clusters Analysis

```{eval-rst}
.. autosummary::
   :toctree: generated/

   tl.wasserstein_distance
   tl.cost_centroid
   tl.cost_distribution
   tl.cost_graph
   tl.cost_expression
   tl.cost_jaccard
   tl.cost_precomputed
   tl.trajectory_calculation
   tl.umap_calculation
   tl.fit_principal_graph
   tl.compute_pseudotime
   tl.sample_nodes
   tl.select_branch
   tl.cell_importance
   tl.cell_importance_calculation
   tl.genes_importance
   tl.gene_cluster_differentiation
   tl.results_gene_cluster_differentiation
   tl.compute_diff_expressions
```

### Clinical and technical variables

```{eval-rst}
.. autosummary::
   :toctree: generated/

   tl.sample_variables
   tl.correlation_with_trajectory
   tl.correlation_with_clustering
   tl.name_clusters
   tl.correlation_categorical_with_trajectory
   tl.correlation_numeric_with_trajectory
   tl.correlation_categorical_with_clustering
   tl.correlation_numeric_with_clustering
```

### MuData

```{eval-rst}
.. autosummary::
   :toctree: generated/

   tl.sample_adata
   tl.to_mudata
   tl.from_mudata
```
## plot: `pl`

```{eval-rst}
.. module:: pilotpy.pl
```

```{eval-rst}
.. currentmodule:: pilotpy
```
### plots

```{eval-rst}
.. autosummary::
   :toctree: generated/

   pl.heatmaps
   pl.heatmap_cost
   pl.heatmap_wasserstein
   pl.trajectory
   pl.trajectory_plotting
   pl.umap_samples
   pl.trajectory_variable
   pl.umap_variable
   pl.trajectory_associations
   pl.select_best_sil
   pl.clustering_emd
   pl.cluster_associations
   pl.plot_trajectory_continuous
   pl.fit_pricipla_graph
   pl.principal_graph
   pl.heatmap_proportions
   pl.cell_importance_grid
   pl.cell_importance_plot
   pl.go_enrichment
   pl.exploring_specific_genes
   pl.genes_selection_heatmap
   pl.plot_hallmark_genes_clusters
   pl.plot_gene_list_pattern
  
```


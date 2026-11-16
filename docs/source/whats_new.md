# What's new

## Version 2.1.0

**Distance settings.** You can now choose how PILOT measures the distance between samples (`pl.tl.wasserstein_distance`). The default settings are unchanged. See the [Distance settings tutorial](Distance_settings.md).

- `cost_method`: how the cost between two cell types is computed: `'centroid'` (default), `'distribution'`, `'graph'`, `'expression'`, `'jaccard'` or your own matrix (`'precomputed'`).
- `n_pcs` and `centroid` (`'median'` or `'mean'`) for the default cost.
- All 12 values of `metric` are documented.
- `regularized`: new transport methods `'sinkhorn_divergence'` and `'unbalanced'` (with `reg_m`).
- The run time of each step is stored.

**Trajectory and plots.**

- `pl.pl.trajectory_variable` and `pl.pl.umap_variable`: color the patients by any clinical or technical variable (e.g. age, sex, number of cells, pseudotime).
- New trajectory parameters: `n_neighbors` (new name for `knn`), `bandwidth_type` and `components`.
- `ax=` for the trajectory and UMAP plots, to draw several settings side by side.

## Version 2.0.19

- `pl.pl.umap_samples`: UMAP of the patients from the Wasserstein distances.
- Stable trajectory plots, automatic colors, and the trajectory is computed once and reused.
- `source_node='auto'`: the principal graph starts at the control samples.
- Separate steps for the principal graph (`tl.fit_principal_graph`, `pl.principal_graph`, `tl.compute_pseudotime`) and the cell type importance (`tl.cell_importance_calculation`, `pl.heatmap_proportions`, `pl.cell_importance_grid`, `pl.cell_importance_plot`). The old one-call functions still work.
- Separate heatmaps: `pl.heatmap_cost`, `pl.heatmap_wasserstein`, `pl.heatmap_proportions`.
- Clear error messages, full docstrings and bug fixes.

# What's new

## Version 2.1.2

**Clinical and technical variables, automatically.**

- `pl.tl.correlation_with_trajectory(adata)`: tests **all** columns of `adata.obs` against the trajectory in one table. The usable variables are found automatically, their type (categorical or numeric) is detected, and columns that cannot be used (one value only, IDs, too many missing values, cell-level labels, duplicates such as 'sex' and 'sex_ontology_term_id') are skipped with a reason. Categorical: ANOVA (and Kruskal-Wallis); numeric: Spearman; p-values corrected (Benjamini-Hochberg). See the [Batch effects: trajectory](Kidney_trajectory.md) tutorial.
- `pl.pl.trajectory_associations(adata)`: the significant variables side by side on the trajectory (or the sample UMAP with `basis = 'umap'`).
- `pl.tl.correlation_with_clustering(adata, proportion_df)` and `pl.pl.cluster_associations(adata)`: the same for the sub-groups of patients (Chi-squared / ANOVA).
- `pl.pl.clustering_emd(..., rename_clusters = True)` and `pl.tl.name_clusters`: the clusters are named after the majority status of their samples (e.g. 'control', 'disease_1', 'disease_2') instead of 0, 1, 2. See the [Batch effects: sub-groups](Kidney_clusters.md) tutorial.
- `pl.tl.sample_variables`: the table of usable sample variables and their types.

**Faster gene analysis.** `pl.tl.genes_importance` saves the expression of each cell type as a small sparse `.h5ad` file in `Results_PILOT/cells` instead of a large CSV file (option `cells_format = 'csv'` for the old format). All functions that read this folder accept both formats, so results of older versions keep working without running anything again.

**MuData.** `pl.tl.to_mudata` saves the results as a MuData with a patient-level modality (one row per sample: proportions, clinical variables, pseudotime, clusters, trajectory, Wasserstein distances); `pl.tl.from_mudata` loads it again and continues the analysis; `pl.tl.wasserstein_distance` accepts a MuData (`modality = 'rna'`). See the [MuData tutorial](MuData.md).

**New tutorial.** [COVID-19 from raw data](COVID_from_scratch.md): quality control, normalization, PCA, Harmony and Leiden, then PILOT and the clinical variables.

**Fixes.**

- `pl.tl.gene_cluster_differentiation`, `pl.pl.exploring_specific_genes` and `pl.pl.genes_selection_heatmap` work again with current scikit-learn and NumPy (2.x).
- `pl.tl.extract_cells_from_gene_expression` now extracts all given cell types (it stopped after the first one).

## Version 2.1.1

**Trajectory branches.** Split the trajectory and analyse each branch separately. See the [Trajectory branches tutorial](Trajectory_branches.md).

- `pl.tl.select_branch(adata, nodes=[...], name=...)`: keep only the samples of the chosen nodes of the principal graph. All other samples are ignored by the next steps (cell type importance, heatmap of the proportions, cell type plots), whose results are saved in a folder per branch. `nodes=None` uses all samples again.
- `pl.tl.sample_nodes(adata)`: the nearest graph node, status and pseudotime of every sample, to help choose the nodes.
- `pl.pl.principal_graph` highlights the selected branch and accepts `ax=` to show branches side by side.
- Without a selected branch, all results are exactly the same as in 2.1.0.

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

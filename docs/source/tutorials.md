# Tutorials

Each tutorial is a Jupyter notebook in the [Tutorial](https://github.com/CostaLab/PILOT/tree/main/Tutorial) folder on GitHub. New users should start with the Myocardial infarction tutorial.

## Disease trajectories

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} 🫀 Myocardial infarction (scRNA-seq)
:link: Myocardial_infarction
:link-type: doc

The full PILOT workflow: Wasserstein distances, trajectory, patient UMAP, clinical variables, principal graph and cell type importance.

**👉 Start here**
:::

:::{grid-item-card} ⚙️ Distance settings
:link: Distance_settings
:link-type: doc

How the metric, the cost between cell types and the transport method change the trajectory, the UMAP and the run time.

**✨ New in 2.1.0**
:::

:::{grid-item-card} 🔬 Kidney IgAN (pathomics)
:link: Combination_Kidney_IgAN
:link-type: doc

Trajectories from tissue structures (glomeruli and tubuli), and how to combine two modalities.
:::
::::

## Patient sub-groups

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} 👥 Pancreas cancer sub-groups
:link: Patients_sub_group_detection
:link-type: doc

Find sub-groups of patients by clustering their distances, then rank the cell types and genes of each group.
:::
::::

## Batch effects

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} 📈 Batch effects in a trajectory
:link: Kidney_trajectory
:link-type: doc

Test whether the trajectory is related to a technical or clinical variable (kidney scRNA-seq).
:::

:::{grid-item-card} 🧩 Batch effects in sub-groups
:link: Kidney_clusters
:link-type: doc

Test whether the patient sub-groups are related to a technical or clinical variable (kidney scRNA-seq).
:::
::::

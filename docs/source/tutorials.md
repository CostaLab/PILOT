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

:::{grid-item-card} 🦠 COVID-19 from raw data
:link: COVID_from_scratch
:link-type: doc

Start from the raw counts: quality control, normalization, PCA, Harmony and Leiden, then PILOT and the clinical variables of the patients.

**✨ New in 2.1.2**
:::

:::{grid-item-card} ⚙️ Distance settings
:link: Distance_settings
:link-type: doc

How the metric, the cost between cell types and the transport method change the trajectory, the UMAP and the run time.

**✨ New in 2.1.0**
:::

:::{grid-item-card} 🌿 Trajectory branches
:link: Trajectory_branches
:link-type: doc

Split the trajectory by choosing nodes of the principal graph, and find the cell types that change in each branch. Examples: Myocardial Infarction, and a lung cancer atlas with two branches (normal to COPD, normal to cancer) and the clinical variables of each branch.

**✨ New in 2.1.1**
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

Test all technical and clinical variables against the trajectory at once; the variables and their types are found automatically (kidney scRNA-seq).

**✨ Updated in 2.1.2**
:::

:::{grid-item-card} 🧩 Batch effects in sub-groups
:link: Kidney_clusters
:link-type: doc

Name the sub-groups after their samples (e.g. control_1, control_2) and test all variables against them at once (kidney scRNA-seq).

**✨ Updated in 2.1.2**
:::
::::

## Saving and sharing results

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} 📦 Export and import with MuData
:link: MuData
:link-type: doc

Save the PILOT results with one row per patient (MuData), load them again, or start from a MuData object.

**✨ New in 2.1.2**
:::
::::

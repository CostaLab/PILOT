# PILOT

```{div} pilot-tagline
**P**at**I**ent-**L**evel analysis with **O**ptimal **T**ransport for single-cell genomics and pathomics data.
```

PILOT compares patients (samples) by their cell type proportions. It computes the Wasserstein distance between every pair of samples with optimal transport, and uses these distances to find disease trajectories and patient sub-groups. It then shows which cell types, genes and tissue structures change along the trajectory or between the groups.

```{image} _static/pilot_overview.png
:alt: Overview of PILOT
:width: 100%
```

## Key features

::::{grid} 1 2 2 2
:gutter: 2

:::{grid-item}
🧭 **Disease trajectories**: order patients from healthy to disease (pseudotime) with diffusion maps and a principal graph.
:::

:::{grid-item}
👥 **Patient sub-groups**: cluster patients by their Wasserstein distances.
:::

:::{grid-item}
🔬 **Interpretation**: find the cell types, genes and tissue structures that change along the trajectory or between groups.
:::

:::{grid-item}
⚙️ **Flexible distances**: choose the metric, the cost between cell types and the transport method.
:::
::::

## Get started

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} 📦 Installation
:link: installation
:link-type: doc

Install PILOT with pip and get the tutorial data.
:::

:::{grid-item-card} 📘 Tutorials
:link: tutorials
:link-type: doc

Step-by-step examples on single-cell and pathomics data.
:::

:::{grid-item-card} 📖 API reference
:link: api
:link-type: doc

All functions of `pl.tl` (tools) and `pl.pl` (plots) with their parameters.
:::

:::{grid-item-card} ✨ What's new
:link: whats_new
:link-type: doc

New in version 2.1.0: distance settings and clinical variables on the trajectory.
:::
::::

Questions, bugs or ideas? Please open an issue on [GitHub](https://github.com/CostaLab/PILOT/issues).

## Citation

If you use PILOT, please cite:

> Joodaki, M., Shaigan, M., Parra, V. *et al.* Detection of PatIent-Level distances from single cell genomics and pathomics data with Optimal Transport (PILOT). *Molecular Systems Biology* 20, 57–74 (2024). [doi:10.1038/s44320-023-00003-8](https://link.springer.com/article/10.1038/s44320-023-00003-8)

```
@article{joodaki2024detection,
  title={Detection of PatIent-Level distances from single cell genomics and pathomics data with Optimal Transport (PILOT)},
  author={Joodaki, Mehdi and Shaigan, Mina and Parra, Victor and B{\"u}low, Roman D and Kuppe, Christoph and H{\"o}lscher, David L and Cheng, Mingbo and Nagai, James S and Goedertier, Micha{\"e}l and Bouteldja, Nassim and others},
  journal={Molecular systems biology},
  volume={20},
  number={2},
  pages={57--74},
  year={2024},
  publisher={Nature Publishing Group UK London}
}
```

```{toctree}
:hidden:
:caption: Getting started

installation
Tutorial overview <tutorials>
whats_new
```

```{toctree}
:hidden:
:caption: Tutorials

Myocardial infarction <Myocardial_infarction>
Distance settings <Distance_settings>
Trajectory branches <Trajectory_branches>
Kidney IgAN (pathomics) <Combination_Kidney_IgAN>
Patient sub-groups (pancreas) <Patients_sub_group_detection>
Batch effects: trajectory <Kidney_trajectory>
Batch effects: sub-groups <Kidney_clusters>
```

```{toctree}
:hidden:
:caption: Reference

api
```


# Welcome to the exciting world of data analysis with PILOT!

```{eval-rst}
..  figure:: logo.png
    :scale: 20%
    :align: right
```

**PILOT** is a Python library for Detection of **P**at**I**ent-**L**evel distances from single cell genomics and pathomics data with **O**ptimal **T**ransport.

🚀 In these three comprehensive tutorials, we'll guide you through a fascinating journey across diverse datasets. 📊✨ Uncover the intricate details of cellular behavior in Myocardial Infarction single cell data, explore the intricate landscape of Kidney IgAN(Glomeruli) & Kidney IgAN(Tubule) with pathomics data, and unravel the complexities of Patients sub-group detection while ranking cells and genes in Pancreas data. Whether you're a beginner or an experienced data enthusiast, our step-by-step guides will empower you to harness the power of PILOT to derive insights and make impactful discoveries from these intricate datasets. Let's dive in and unlock the hidden insights within the data together! 🧬🔍💡


## Installation Guide

Follow these steps to install and set up PILOT:

```bash
conda create --name PILOT python=3.11.5 r-base
conda activate PILOT
pip install pilotpy
```
Once you've completed these steps, you can proceed to run the tutorials and explore the features of PILOT. 
When doing so, remember to move to the tutorial folder, as all the work will be performed there:
```bash

cd Tutorial


```


## Citation
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
---
caption: scRNA-seq Analysis
maxdepth: 2
---
Myocardial_infarction
```

```{toctree}
---
maxdepth: 2
caption: Pathomics data Analysis
---
Combination_Kidney_IgAN
```

```{toctree}
---
maxdepth: 2
caption: Patients sub-group detection
---
Patients_sub_group_detection
```

```{toctree}
---
maxdepth: 2
caption: Evaluation of the presence of
batch effects(Trajectory)
---
Kidney_trajectory
```
```{toctree}
---
maxdepth: 2
caption: Evaluation of the presence of
batch effects(Clusters)
---
Kidney_clusters
```

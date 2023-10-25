
# Welcome to the exciting world of data analysis with PILOT!

**PILOT** is a Python library for Detection of **P**at**I**ent-**L**evel distances from single cell genomics and pathomics data with **O**ptimal **T**ransport.

🚀 In these three comprehensive tutorials, we'll guide you through a fascinating journey across diverse datasets. 📊✨ Uncover the intricate details of cellular behavior in Myocardial Infarction single cell data, explore the intricate landscape of Kidney IgAN(Glomeruli) & Kidney IgAN(Tubule) with pathomics data, and unravel the complexities of Patients sub-group detection while ranking cells and genes in Pancreas data. Whether you're a beginner or an experienced data enthusiast, our step-by-step guides will empower you to harness the power of PILOT to derive insights and make impactful discoveries from these intricate datasets. Let's dive in and unlock the hidden insights within the data together! 🧬🔍💡


## Installation Guide

Follow these steps to install and set up PILOT:

```bash
git clone https://github.com/CostaLab/PILOT
cd PILOT
conda create --name PILOT r-base
conda activate PILOT
conda install -c conda-forge rpy2
conda install jupyter
pip install .
```
Once you've completed these steps, you can proceed to run the tutorials and explore the features of PILOT. 
When doing so, remember to move to the tutorial folder, as all the work will be performed there:
```bash

cd Tutorial


```



```{toctree}
---
caption: Tutorial for scRNA-seq Analysis
maxdepth: 2
---
Myocardial_infarction
```

```{toctree}
---
maxdepth: 2
caption: Tutorial for pathomics data Analysis
---
Combination_Kidney_IgAN
```

```{toctree}
---
maxdepth: 2
caption: Tutorial for Patients sub-group detection
---
Patients_sub_group_detection
```

```{toctree}
---
maxdepth: 2
caption: Tutorial for evaluation of the presence of batch effects(Trajectory)
---
Kidney_trajectory
```
```{toctree}
---
maxdepth: 2
caption: Tutorial for evaluation of the presence of batch effects(Clusters)
---
Kidney_clusters
```
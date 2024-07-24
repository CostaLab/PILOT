# PILOT  <img src="img/logo.png" align="right" width="300" />

[![GitHub license](https://img.shields.io/github/license/CostaLab/PILOT.svg)](https://github.com/CostaLab/PILOT?tab=MIT-1-ov-file#MIT-1-ov-file)

Although clinical applications represent the next challenge in single-cell genomics and digital pathology, we still lack computational methods to analyze single-cell or pathomics data to find sample-level trajectories or clusters associated with diseases. This remains challenging as single-cell/pathomics data are multi-scale, i.e., a sample is represented by clusters of cells/structures, and samples cannot be easily compared with each other. Here we propose PatIent Level analysis with Optimal Transport (PILOT). PILOT uses optimal transport to compute the Wasserstein distance between two individual single-cell samples. This allows us to perform unsupervised analysis at the sample level and uncover trajectories or cellular clusters associated with disease progression. We evaluate PILOT and competing approaches in single-cell genomics or pathomics studies involving various human diseases with up to 600 samples/patients and millions of cells or tissue structures. Our results demonstrate that PILOT detects disease-associated samples from large and complex single-cell or pathomics data. Moreover, PILOT provides a statistical approach to find changes in cell populations, gene expression, and tissue structures related to the trajectories or clusters supporting interpretation of predictions.

![plot](./img/plot.png)


Current version for PILOT is 2.0.6

## Installation
The easiest way to install PILOT and the required packages is using the following way:

```terminal

conda create --name PILOT python=3.11.5 r-base
conda activate PILOT
pip install pilotpy
```
Once you've completed these steps, you can proceed to run the tutorials and explore the features of PILOT. 
When doing so, remember to move to the tutorial folder, as all the work will be performed there:

```terminal
git clone https://github.com/CostaLab/PILOT.git
cd PILOT/Tutorial
```

## [Tutorial&Data sets](https://pilot.readthedocs.io/en/latest/index.html)
There are five tutorials, one for [Myocardial Infarction (single cell data)](https://pilot.readthedocs.io/en/latest/Myocardial_infarction.html), and the second tutorial for [pathomics data, the combination of Kidney IgAN(G) & Kidney IgAN(T)](https://pilot.readthedocs.io/en/latest/Combination_Kidney_IgAN.html), and the third one  for [Patients sub-group detection and then ranking cells/genes (Pancreas data)](https://pilot.readthedocs.io/en/latest/Patients_sub_group_detection.html) and the forth one for 
[evaluation of the presence of batch effects in Trajectory](https://pilot.readthedocs.io/en/latest/Kidney_trajectory.html) and the last one for [evaluation of the presence of batch effects in detected sub-groups](https://pilot.readthedocs.io/en/latest/Kidney_clusters.html).


You can access the used data sets by PILOT in Part 1 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4740646.svg)](https://zenodo.org/records/8370081) and  Part 2 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4740646.svg)](https://zenodo.org/records/7957118)


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





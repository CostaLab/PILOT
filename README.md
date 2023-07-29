# PILOT

PILOT uses optimal transport to compute the Wasserstein distance between two single single-cell experiments. This allows us to perform unsupervised analysis at the sample level and to uncover trajectories associated with disease progression. Moreover, PILOT provides a statistical approach to delineate non-linear changes in cell populations, gene expression and tissues structures related to the disease trajectories.  We evaluate PILOT and competing approaches in  disease single-cell genomics and pathomics studies with up to 1.000 patients/donors and millions of cells or structures. Results demonstrate that PILOT detects disease-associated samples, cells, and genes from large and complex single-cell and pathomics data.


![plot](./img/plot.png)


Current version for PILOT is 1.1.0

## User installation
The easiest way to install PILOT and the required packages is using pip,

```terminal
git clone https://github.com/CostaLab/PILOT
cd PILOT
pip install .
```
## [Tutorial](https://github.com/CostaLab/PILOT/tree/main/Tutorial)
There are two tutorials, one for [Myocardial Infarction (single cell data)](https://github.com/CostaLab/PILOT/blob/main/Tutorial/%20Myocardial%20infarction.ipynb) and the second tutorial for [pathomics data, the combination of Kidney IgAN(G) & Kidney IgAN(T)](https://github.com/CostaLab/PILOT/blob/main/Tutorial/Combination_Kidney_IgAN.ipynb).




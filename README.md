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
For applying PILOT to your data, we provide tutorials here, you just need to after installation.

There are two tutorials, one for myocardial infarction data (single cell data) and the second tutorial is for pathomics data (Kidney IgAN(G) & Kidney IgAN(T)).
For [Tutorial of single-cell processing](https://github.com/CostaLab/PILOT/blob/main/Tutorial/%20Myocardial%20infarction.ipynb) and for [Tutorial of pathomics data processing]](https://github.com/CostaLab/PILOT/blob/main/Tutorial/%20Myocardial%20infarction.ipynb)

You can see the required data and its structure in [Datasets](https://github.com/CostaLab/PILOT/tree/main/Tutorial/Datasets).



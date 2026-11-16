# Installation

## Install PILOT

PILOT needs Python 3.9 or newer. We recommend a new conda environment:

```bash
conda create -n pilot python=3.11
conda activate pilot
pip install pilotpy
```

Check the installation:

```python
import pilotpy as pl
print(pl.__version__)
```

## Get the tutorials

The tutorials are Jupyter notebooks in the `Tutorial` folder of the GitHub repository. Run them from that folder, because PILOT writes its results there (in `Results_PILOT`):

```bash
git clone https://github.com/CostaLab/PILOT.git
cd PILOT/Tutorial
```

## Data

Each tutorial says which data it uses and where to download it. Place the files in the `Tutorial/Datasets` folder. All data sets used in the PILOT paper are on Zenodo: [Part 1](https://zenodo.org/records/8370081) and [Part 2](https://zenodo.org/records/7957118).

## Import

```python
import pilotpy as pl
```

- `pl.tl`: tools (calculations), e.g. `pl.tl.wasserstein_distance`.
- `pl.pl`: plots, e.g. `pl.pl.trajectory`.

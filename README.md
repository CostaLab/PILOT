# PILOT  <img src="img/logo.png" align="right" width="300" />

[![GitHub license](https://img.shields.io/github/license/CostaLab/PILOT.svg)](https://github.com/CostaLab/PILOT?tab=MIT-1-ov-file#MIT-1-ov-file)

**Authors:**
 Mehdi Joodaki<sup>[*1]</sup>
 ,Mina Shaigan<sup>[*1]</sup>
 ,Victor Parra<sup>[1]</sup>
 ,Roman D. Bülow<sup>[2]</sup>
 ,Christoph Kuppe<sup>[3]</sup>
 ,David L. Hölscher<sup>[2]</sup>
 ,Mingbo Cheng<sup>[1]</sup>
 ,James S. Nagai<sup>[1]</sup>
 ,Michaël Goedertier<sup>[1,2]</sup>
 ,Nassim Bouteldja<sup>[2]</sup>
 ,Vladimir Tesar<sup>[5]</sup> 
,Jonathan Barratt<sup>[6,7]</sup>
 ,Ian S.D. Roberts<sup>[8]</sup>
 ,Rosanna Coppo<sup>[9]</sup>
 ,Rafael Kramann<sup>[3,4]</sup>
 ,Peter Boor<sup>[2,@]</sup>
 ,Ivan G. Costa<sup>[1,@]</sup>

**Affiliations:**
- [1] Institute for Computational Genomics, Joint Research Center for Computational Biomedicine, RWTH Aachen University Medical School
- [2] Institute of Pathology, Laboratory of Nephropathology, RWTH Aachen University Medical School
- [3] Institute of Experimental Medicine and Systems Biology, RWTH Aachen University
- [4] Department of Internal Medicine, Nephrology and Transplantation, Erasmus Medical Center
- [5] Department of Nephrology, $1^{st}$ Faculty of Medicine and General University Hospital, Charles University, Prague, Czech Republic
- [6] John Walls Renal Unit, University Hospital of Leicester National Health Service Trust, Leicester, United Kingdom
- [7] Department of Cardiovascular Sciences, University of Leicester, Leicester, United Kingdom
- [8] Department of Cellular Pathology, Oxford University Hospitals National Health Services Foundation Trust, Oxford, United Kingdom
- [9] Fondazione Ricerca Molinette. Regina Margherita Children's University Hospital, Torino, Italy
  
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

# References
 PILOT: Detection of PatIent-Level distances from single cell genomics and pathomics data with Optimal Transport (PILOT) [link](https://www.embopress.org/doi/full/10.1038/s44320-023-00003-8)








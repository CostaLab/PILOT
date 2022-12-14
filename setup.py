from setuptools import setup

setup(
name='PILOT',
version='0.1.0',
author='Mehdi Joodaki',
author_email='judakimehdi@gmail.com',
packages=['PILOT'],
scripts=['PILOT/Cell_gene_selection.py','PILOT/Trajectory.py'],
install_requires=[
"cycler",
"gprofiler-official==0.3.5",
"joypy",
"leidenalg",
"numpy",
"matplotlib",
"pandas",
"plotly",
"plotnine",
"POT",
"pydiffmap",
"scanpy",
"scikit_learn",
"scikit_network",
"scipy",
"seaborn",
"shap",
"statsmodels"
],
)

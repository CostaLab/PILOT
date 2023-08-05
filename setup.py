from setuptools import setup

setup(
    name='PILOT',
    version='1.1.0',
    author='Mehdi Joodaki',
    author_email='judakimehdi@gmail.com',
    packages=['PILOT'],
    scripts=['PILOT/Cell_gene_selection.py','PILOT/Trajectory.py','PILOT/patients_sub_clustering.py'],
    install_requires=[
        "cycler",
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
        "statsmodels",
        "elpigraph-python",
        "adjustText",
        "gprofiler-official",
    ],
)


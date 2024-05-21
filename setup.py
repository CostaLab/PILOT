from setuptools import find_packages, setup

setup(
    name='PILOT',
    version='2.0.4',
    author='Mehdi Joodaki',
    author_email='judakimehdi@gmail.com',
    url='https://github.com/CostaLab/PILOT',
    
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
    packages=find_packages()
)


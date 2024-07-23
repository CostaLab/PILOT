from setuptools import find_packages, setup

setup(
    name='pilotpy',
    version='2.0.6',
    author='Mehdi Joodaki',
    author_email='judakimehdi@gmail.com',
    url='https://github.com/CostaLab/PILOT',
    python_requires='>=3.11.5,<3.12',
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
            "statsmodels",
            "elpigraph-python",
            "adjusttext",
            "gprofiler-official",
            "rpy2",
  
        ],
        packages=find_packages()
)



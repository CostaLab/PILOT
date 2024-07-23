from setuptools import find_packages, setup

setup(
    name='pilotpy',
    version='2.0.6',
    author='Mehdi Joodaki',
    author_email='judakimehdi@gmail.com',
    url='https://github.com/CostaLab/PILOT',
    install_requires=[
          "rpy2>=3.5.11",
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
           
  
        ],
        packages=find_packages()
)



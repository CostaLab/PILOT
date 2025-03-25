from setuptools import find_packages, setup

setup(
    name='pilotpy',
    version='2.0.6',
    author='Mehdi Joodaki',
    author_email='judakimehdi@gmail.com',
    url='https://github.com/CostaLab/PILOT',
    python_requires='>=3.11.5,<3.12',
    install_requires=[
            "cycler>=0.11.0,<0.12.0",
            "joypy>=0.2.6,<0.3.0",
            "leidenalg>=0.10.1,<0.11.0",
            "numpy>=1.24.4,<1.25.0",
            "matplotlib==3.8.0",
            "pandas>=2.0.3,<2.1.0",
            "plotly>=5.22.0,<5.23.0",
            "plotnine>=0.12.3,<0.13.0",
            "pot>=0.9.1,<0.10.0",
            "pydiffmap>=0.2.0.1,<0.3.0",
            "scanpy>=1.9.5,<1.10.0",
            "scikit-learn>=1.3.0,<1.4.0",
            "scikit-network>=0.31.0,<0.32.0",
            "scipy>=1.11.2,<1.12.0",
            "seaborn>=0.12.2,<0.13.0",
            "shap>=0.42.1,<0.43.0",
            "statsmodels>=0.14.0,<0.15.0",
            "elpigraph-python>=0.3.1,<0.4.0",
            "adjusttext>=0.8,<0.9",
            "gprofiler-official>=1.0.0,<1.1.0",
            "rpy2>=3.5.11",
            "gseapy>=1.1.7"
        ],
        packages=find_packages()
)


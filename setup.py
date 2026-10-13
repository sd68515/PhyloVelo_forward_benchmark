from setuptools import find_packages, setup

setup(
    name="phylovelo",
    version="0.1.1",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=[
        "biopython>=1.78",
        "matplotlib>=3.4",
        "numpy>=1.21",
        "pandas>=1.2",
        "scikit-learn>=0.24",
        "scipy>=1.6",
        "statsmodels>=0.12",
        "tqdm>=4.60",
    ],
    author="Kun Wang and contributors",
    license="BSD-3-Clause",
)

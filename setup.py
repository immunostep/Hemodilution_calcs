from setuptools import setup, find_packages

setup(
    name="hemodilution_calcs",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "fcsparser==0.2.8",
        "hdbscan==0.8.39",
        "matplotlib==3.10.7",
        "numpy==1.26.4",
        "pandas==2.3.3",
        "scikit-learn==1.7.2",
        "scipy==1.16.3",
        "seaborn==0.13.2"
    ],
    python_requires=">=3.12",
)

from setuptools import setup, find_packages

setup(
    name="metacontam",
    version="0.0.1",
    description="Contamination detection tool for metagenomic sequencing data",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Lifemining-lab/Metacontam",
    packages=find_packages(),
    package_data={
        "Metacontam": [
            "species_under_blacklist.pkl",
            "Netcomi.R",
        ]
    },
    scripts=["metacontam"],
    python_requires=">=3.10",
    install_requires=[
        "pandas>=2.0.3",
        "numpy>=1.25.2",
        "networkx==3.4.2",
        "biopython>=1.74",
        "matplotlib>=3.8.2",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)

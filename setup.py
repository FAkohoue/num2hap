from setuptools import setup, find_packages

setup(
    name="num2hap",
    version="0.1",
    author="Félicien Akohoue",
    author_email="akohoue.f@gmail.com",
    description="A package for converting numeric genotype data to hapmap diploid format",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/FAkohoue/num2hap",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "tqdm"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
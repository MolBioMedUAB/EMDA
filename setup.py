from setuptools import setup, find_packages
import EMDA

with open("README.md", "r") as f:
    long_description = f.read()

#version='0.2.0a0'
#version='{{VERSION_PLACEHOLDER}}'

setup(
    name="EasyMDA",
    version=EMDA.__version__,
    description="Easy MD Analysis (EMDA) is a Python package created with the aim of providing an easy, yet powerful way to perform analysis of MD simulations.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Miquel Canyelles Niño",
    author_mail="mcanyellesnino@gmail.com",
    packages=find_packages(where='.'),
    include_package_data=True,
    install_requires=["MDAnalysis", "tqdm", "matplotlib"],
    keywords="biochemistry, simulations, MDAnalysis, molecular dynamics",
    url="https://github.com/MolBioMedUAB/EMDA",
    download_url=f"https://github.com/MolBioMedUAB/EMDA/archive/refs/tags/{EMDA.__version__}.tar.gz",
)

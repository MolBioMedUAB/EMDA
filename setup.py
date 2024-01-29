from setuptools import setup, find_packages

with open("README.md", "r") as f:
    long_description = f.read()

version='0.1.0'
#version='{{VERSION_PLACEHOLDER}}'

setup(
    name="EasyMDA",
    vesion=version,
    description="Easy MD Analysis (EMDA) is a Python package created with the aim of providing an easy, yet powerful way to perform analysis of MD simulations.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Miquel Canyelles Niño",
    author_mail="mcanyellesnino@gmail.com",
    packages=find_packages(),
    include_package_data=True,
    install_requires=["numpy", "pyyaml", "MDAnalysis", "tqdm"],
    keywords="biochemistry, simulations, MDAnalysis, molecular dynamics",
    url="https://github.com/MolBioMedUAB/EMDA",
    download_url=f"https://github.com/MolBioMedUAB/EMDA/archive/refs/tags/{version}.tar.gz",
)

# Easy MD Analysis (EMDA)

## Description

Easy MD Analysis (EMDA) is a Python package based on [MDAnalysis](http://mdanalysis.org) created with the aim of providing an easy, yet powerful way to perform analysis of MD simulations. 

## Features

### Measurers

Several kinds of measures have been implemented in EMDA. Each of them has its own adder function and calculator function. Adder functions are the functions through which the user is able to request a certain kind of measure, while calculators are the engines that perform the measure when the run function is executed. 

The available measures are listed below:
- Distance
- Angle
- Dihedral
- Planar angle
- Distance of bridging waters between two sets of atoms
- RMSD
- Contacts, both of a group of atoms and of a whole protein

## Installation

### Using pip

```bash
pip install EasyMDA
```

### From source code

1. Clone the GitHub repository in your local machine:
```bash
git clone https://github.com/MolBioMedUAB/EMDA
```

2. Move to the repository direcctory:
```bash
cd EMDA
```

3. Install using pip:
```bash
pip install .
```


## How to use EMDA

EMDA has been design to perform analysis of MD trajectories in three steps: measure, analyse, and plot. Nonetheless, measuring is the most time-consuming task in an analysis. Thus, the code has been structured in such a way that all the requested measures are firstly added using the adders and the ran with the run method. On the other hand, analysis of the measures and plotting tasks are executed at the moment. An example Jupyter notebook showing how to perform an analysis of a sample trajectory can be found (here)[https://github.com/MolBioMedUAB/EMDA/blob/main/example/example.ipynb].



# Easy MD Analysis (EMDA)

## Description

Easy MD Analysis (EMDA) is a Python package based on [MDAnalysis](http://mdanalysis.org) created with the aim of providing an easy, yet powerful way to perform analysis of MD simulations. 

## Table of contents

- [Features](#features)
- [Installation](#installation)
- [How to use EMDA](#how-to-use-emda)
- [Future implementations](#future-implementations)
- [Changelog](#changelog)

## Features

### Measurers

Several kinds of measures have been implemented in EMDA. Each of them has its own adder function and calculator function. Adder functions are the functions through which the user is able to request a certain kind of measure, while calculators are the engines that perform the measure when the run function is executed. 

The available measures are listed below:
- __Distance__: measures the distance between two sets of atoms
- __Angle__: measures the angle between three atoms
- __Dihedral__: measures the dihedral angle between four atoms
- __Planar angle__: measures the angle between the closest planes to two sets of at least three atoms
- __Distance of bridging waters between two sets of atoms__: identifies the closest water that is bridging between two sets of atoms and measures the distances to each
- __RMSD__: measures the RMSD of a set of atoms (or the whole system) in reference of a frame of the structure
- __Contacts__, both of a group of atoms and of a whole protein: identifies the contacts stablished by a selection in a given radius or the contacts of each residue.

### Analysers

Some analysis can be performed from previous measures and they are stored as Analysis classes. In this case, each analysis has its own analyser. In opposition with measures, analysis are executed when requested. All the analysers functions' names start with the 'analyse_' string.

The available analyis are listed below:
- __value__: analyses the value of a frame-wise measure (like distance, for instance) and returns frame-wise list containing True if the value is between the given values or False if it is not.
- __contacts_frequency__: analyses the contacts and returns a dictionary containing the contacts that take place and how many times it takes place (in an absolute or relative number).
- __contacts_amounts__: analyses the contacts and returns a frame-wise list containing how many contacts a selection (or a residue) stablishes in each frame.
- __NACs__ (near-attack conformations): analyses two or more analysed values (so a frame-wise boolean list) and returns the combination of all the values as a boolean frame-wise list.


### Plotters

Some analysis or measures can be plotted. The plotters functions (named with the plot_ prefix) take the analysis or measures' result and returns a plot depending on the type of data.

The available plotters are listed below:
- __values__: plots a float-containing frame-wise list. A similar method has been implemented inside the Measure class.
- __contacts_frequencies_diff__: external plotter (so it is not a method of the EMDA class). It takes two contacts_frequency (analyser) results (or two lists of), compares them so a bar plot is returned containing the contacts that are the most different between the two sets.

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

EMDA has been design to perform analysis of MD trajectories in three steps: measure, analyse, and plot. Nonetheless, measuring is the most time-consuming task in an analysis. Thus, the code has been structured in such a way that all the requested measures are firstly added using the adders and the ran with the run method. On the other hand, analysis of the measures and plotting tasks are executed at the moment. 

An example Jupyter notebook showing how to perform an analysis of a sample trajectory can be found [here](https://github.com/MolBioMedUAB/EMDA/blob/main/example/example.ipynb).



## Future implementations

Future features that will be added in the future can be found in the [TO-DO](https://github.com/MolBioMedUAB/EMDA/blob/main/TO-DO.md) file


## Changelog

Additions or fixed bugs in each version can be found in the [CHANGELOG](https://github.com/MolBioMedUAB/EMDA/blob/main/CHANGELOG.md) file.
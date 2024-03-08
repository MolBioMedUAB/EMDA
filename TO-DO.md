# TO-DO

## GENERAL
- [ ] Add async features to run method
- [ ] Add print statistics' method to Measure to print avg of value (for dist and similar) and sd.
- [ ] Add H-bond as measure
- [ ] Add selection exporter
- [ ] Add pymol script creator for exporting loaded variants/replicas and selections.
- [ ] Add RMSF, HELANAL and other analysers from MDAnalysis
- [ ] Add trajectory exporter
- [ ] Add saver for specific results, not the whole EMDA status.


## PER SECTION
### CONTACTS/PROTEIN_CONTACTS
- [ ] Add the possibility to give a selection instead of the whole protein. It has to have more than one residue.
- [ ] Add other options besides distances to be calculated between conctacting groups such as contacting area, etc.

### PLOTTERS
- [ ] Move ext_plot_contacts_frequencies_differences to internal.
- [ ] Add multiple values in plot_value so several RMSDs, distances, etc can be plotted in the same plot (per variant and replica)
- [ ] :exclamation: Fix global colorbar in plot_probability_densities
- [ ] :exclamation: Add axis label builder for specifying distance, angle, etc labels of PDF plots.

### RMSD ADDER AND CALCULATOR
- [ ] Add the option to set the average structure as reference

### ANALYSE_CONTACTS_FREQUENCY
- [ ] Allow normalisation by variant and not only by variant & replica.
- [ ] Allow global normalisation by variant & replica, or variant.

## PER VERSION
 
### 1.1.0
- [ ] Fix distWATbridge
- [ ] Fix pKa

---------------------------------------------------------------------

# DONE

## GENERAL
- [X] Add some code previous to the running cycle in the run method to extend selections to these variants/replicas that have been added after the selection
- [X] Add option to load universes in memory.
- [X] Add multi-variant and multi-replica function to EMDA class, so more than one variant (parameters and trajectory) and/or more than one replica (trajectory) can be loaded into the EMDA class. Thus, all the measures, analysis and plots are added once and executed for all.
- [X] Add a method to EMDA class to check if selections are correct within the different vairants/replicas
- [X] :exclamation: Add unwrap option when a variant/replica is loaded
- [X] Add saving and loading options.

## PER-SECTION

### ANALYSIS
- [X] :exclamation: Add analysis to check if a contact is present or not
- [X] Contacts amounts to analyse_value and NACs

### ADD_CONTACTS
- [X] Split selection and whole protein contacts into to different adders
- [X] Remove out_format-related features

### EMDA.RUN
- [X] Update adders before running
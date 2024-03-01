# TO-DO

## GENERAL
- [ ] Add a method to EMDA class to check if selections are correct within the different vairants/replicas
- [ ] Add async features to run method
- [ ] Add print statistics' method to Measure to print avg of value (for dist and similar) and sd.
- [ ] :exclamation: Add unwrap option when a variant/replica is loaded
- [ ] Add saving and loading options.
- [ ] Add H-bond as measure
- [ ] :exclamation: Add anaylsis to check if a contact is present or not
- [ ] Contacts amounts to analyse_value and NACs

## PER SECTION
### CONTACTS/PROTEIN_CONTACTS
- [ ] Add the possibility to give a selection instead of the whole protein. It has to have more than one residue.

### PLOTTERS
- [ ] Move ext_plot_contacts_frequencies_differences to internal.
- [ ] Add multiple values in plot_value

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

## PER-SECTION

### ADD_CONTACTS
- [X] Split selection and whole protein contacts into to different adders
- [X] Remove out_format-related features

### EMDA.RUN
- [X] Update adders before running
# TO-DO

## GENERAL

- [ ] Add multi-variant and multi-replica function to EMDA class, so more than one variant (parameters and trajectory) and/or more than one replica (trajectory) can be loaded into the EMDA class. Thus, all the measures, analysis and plots are added once and executed for all.
- [ ] Add a method to EMDA class to check if selections are correct within the different vairants/replicas
- [ ] Add async features to run method
- [ ] Add print statistic method to Measure to print avg of value (for dist and similar) and sd.

## PER SECTION
### ADD_CONTACTS


### RMSD ADDER AND CALCULATOR
- [ ] Add the option to set the average structure as reference

### ANALYSE_CONTACTS_FREQUENCY
- [ ] Allow normalisation by variant and not only by variant & replica.
- [ ] Allow global normalisation by variant & replica, or variant.

### EMDA.RUN
- [ ] Update adders before running so 

###

## PER VERSION
 
### 1.1.0
- [ ] Fix distWATbridge
- [ ] Fix pKa

---------------------------------------------------------------------

# DONE
- [X] Add some code previous to the running cycle in the run method to extend selections to these variants/replicas that have been added after the selection


## ADD_CONTACTS
- [X] Split selection and whole protein contacts into to different adders
- [X] Remove out_format-related features


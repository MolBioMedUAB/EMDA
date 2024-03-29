# CHANGELOG

## 0.3.0
- Now only parameters are mandatory, so multiframe PDBs are accepted.
- Plotter for plotting a contacts_frequency has been added. 
- Normalisation to the most frequent contact in analyse_contacts_frequency has been added as an option.


## 0.2.0
- analyse_NACs analyser has been created to combine 2 or more analyse_value results. get_most_frequent function has been added to tools.py. NotEqualLenghtsError and NotEnoughDataError exceptions have been created for the analyser.
- analyse_contacts_amount has been created to calculate the number of contacts in each frame.
- save_result and read_result methods of EMDA class have been added to store data from measures or analysis' result attribute. 
- plot method of Measure class has been added to plot frame-wise results such as distance, angle, dihedral, RMSD and planar_angle.
- ext_plot_contacts_frequencies_differences plotter has been created. It takes two (or two sets of) frequency of contacts' analysis result and returns a bar graph showing the contacts that are the most different between the two input lists depending on the given threshold.
- exclude option of EMDA's run method has been fixed. It was ignored if a string was given as input instead of a list.
- Classes, methods and functions' descriptions have been added or improved.


## 0.1.0

- First analysers have been created for evaluating values within ranges and to calculate frequencies of contacts.
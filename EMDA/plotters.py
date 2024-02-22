import matplotlib.pyplot as plt
from numpy import absolute as abs

"""
TO BUILD:
    - [] Distance vs contacts frequency/amounts as scatter
    - [] 
"""

def plot_measure(self, measure_name, same_y : bool = True, same_x : bool = True, axis_label_everywhere=False, merge_replica_plots=False, out_name=False):
    """
    DESCRIPTION:
        Function for plotting all measures in a Measure class. Only suitable for frame-wise lists 
        like the obtained with distance, angle, dihedral, planar angle and RMSD.

    OPTIONS:
        - measure_name:             [str]        name of the measure listed as EMDA.measures key
        - same_y, same_x:           [bool]       shares the y and/or x among all plots in the same row/column, so plots have the same axis dimensions
        - axis_label_everywhere:    [bool]       Adds the x and y axis labels to all the subplots instead of only to the ones at the left and bottom
        - out_name:                 [pseudobool] False by default. If a string is given, the plot will be saved.
        - 
    """
    print('merge', merge_replica_plots)

    # Plotting from a plotter
    if measure_name != None:
        variants = len(self.measures[measure_name].result)

        if merge_replica_plots:
            max_replicas = 1
        else :
            max_replicas = max([ len(self.measures[measure_name].result[variant]) for variant in list(self.measures[measure_name].result) ])

        fig, axs = plt.subplots(ncols=variants, nrows=max_replicas, sharey=same_y, sharex=same_x)

        for v_num, variant in enumerate(list(self.measures[measure_name].result.keys())):
            for r_num, replica in enumerate(list(self.measures[measure_name].result[variant].keys())):

                if merge_replica_plots:
                    axs[v_num].plot(
                    range(1, len(self.measures[measure_name].result[variant][replica])+1),
                    self.measures[measure_name].result[variant][replica],
                    c = f"C{v_num}",
                    label=replica
                    )

                    if v_num == 0:
                        if self.measures[measure_name].type in ("distance"):
                            axs[v_num].set_ylabel("Distance (Å)")

                        elif self.measures[measure_name].type in ("RMSD"):
                            axs[v_num].set_ylabel("RMSD (Å)")

                        elif self.measures[measure_name].type in ("angle"):
                            axs[v_num].set_ylabel("Angle (°)")

                        elif self.measures[measure_name].type in ("planar_angle"):
                            axs[v_num].set_ylabel("Planar angle (°)")
                        
                        elif self.measures[measure_name].type in ("dihedral"):
                            axs[v_num].set_ylabel("Dihedral angle (°)")

                    if r_num == max_replicas-1 or axis_label_everywhere:
                        axs[v_num].set_xlabel("Frame")
                        
                    axs[v_num].set_title(f"{variant}, replicas {', '.join(list(self.measures[measure_name].result[variant].keys()))}")  

                else :
                    axs[r_num, v_num].plot(
                        range(1, len(self.measures[measure_name].result[variant][replica])+1),
                        self.measures[measure_name].result[variant][replica],
                        c = f"C{v_num}"
                        )
                    
                    if v_num == 0 or axis_label_everywhere:
                        if self.measures[measure_name].type in ("distance"):
                            axs[r_num, v_num].set_ylabel("Distance (Å)")

                        elif self.measures[measure_name].type in ("RMSD"):
                            axs[r_num, v_num].set_ylabel("RMSD (Å)")

                        elif self.measures[measure_name].type in ("angle"):
                            axs[r_num, v_num].set_ylabel("Angle (°)")

                        elif self.measures[measure_name].type in ("planar_angle"):
                            axs[r_num, v_num].set_ylabel("Planar angle (°)")
                        
                        elif self.measures[measure_name].type in ("dihedral"):
                            axs[r_num, v_num].set_ylabel("Dihedral angle (°)")

                    if r_num == max_replicas-1 or axis_label_everywhere:
                        axs[v_num, r_num].set_xlabel("Frame")
                    
                    axs[v_num, r_num].set_title(f"{variant}, {replica}")    

        fig.suptitle("Plots for " + r"$\bf{%s}$" % measure_name.replace('_', '\_') +  " Measure")


    # code for Measure's method
    elif measure_name == None:
        variants = len(self.result)

        if merge_replica_plots:
            max_replicas = 1
        else :
            max_replicas = max([ len(self.result[variant]) for variant in list(self.result) ])

        fig, axs = plt.subplots(ncols=variants, nrows=max_replicas, sharey=same_y, sharex=same_x)

        for v_num, variant in enumerate(list(self.result.keys())):
            for r_num, replica in enumerate(list(self.result[variant].keys())):

                if merge_replica_plots:
                    axs[v_num].plot(
                        range(1, len(self.result[variant][replica])+1),
                        self.result[variant][replica],
                        c = f"C{v_num}",
                        label=replica
                    )

                    if v_num == 0:
                        if self.type in ("distance"):
                            axs[v_num].set_ylabel("Distance (Å)")

                        elif self.type in ("RMSD"):
                            axs[v_num].set_ylabel("RMSD (Å)")

                        elif self.type in ("angle"):
                            axs[v_num].set_ylabel("Angle (°)")

                        elif self.type in ("planar_angle"):
                            axs[v_num].set_ylabel("Planar angle (°)")
                        
                        elif self.type in ("dihedral"):
                            axs[v_num].set_ylabel("Dihedral angle (°)")

                    if r_num == max_replicas-1 or axis_label_everywhere:
                        axs[v_num].set_xlabel("Frame")
                    
                    axs[v_num].set_title(f"{variant}, replicas {', '.join(list(self.result[variant].keys()))}")  

                else :
                    axs[r_num, v_num].plot(
                        range(1, len(self.result[variant][replica])+1),
                        self.result[variant][replica],
                        c = f"C{v_num}"
                        )
                    
                    if v_num == 0 or axis_label_everywhere:
                        if self.type in ("distance"):
                            axs[r_num, v_num].set_ylabel("Distance (Å)")

                        elif self.type in ("RMSD"):
                            axs[r_num, v_num].set_ylabel("RMSD (Å)")

                        elif self.type in ("angle"):
                            axs[r_num, v_num].set_ylabel("Angle (°)")

                        elif self.type in ("planar_angle"):
                            axs[r_num, v_num].set_ylabel("Planar angle (°)")
                        
                        elif self.type in ("dihedral"):
                            axs[r_num, v_num].set_ylabel("Dihedral angle (°)")

                    if r_num == variants-1 or axis_label_everywhere:
                        axs[r_num, v_num].set_xlabel("Frame")
                    
                    axs[r_num, v_num].set_title(f"{variant}, {replica}")    
                    
        fig.suptitle("Plots for " + r"$\bf{%s}$" % self.name.replace('_', '\_') +  " Measure")
        

    fig.tight_layout()

    #if merge_replica_plots:
    #    plt.legend()

    if out_name != False and isinstance(out_name, str):
        if not out_name.endswith((".png", ".jpg", ".jpeg", ".tiff")):
            out_name = ".".join(out_name.split('.')[:-1]) + ".png"

        plt.savefig(out_name, dpi=300, bbox_inches='tight')

    plt.show()
    plt.close()


def plot_contacts_frequency(self, analysis_name, out_name=None, fill_empty=False, width_plot=0.5):
    """
    DESCRIPTION:
        Plotter that takes the result of a contacts_frequency analysis and plots each interaction as a bar plot
    """

    if self.analyses[analysis_name].type != "contacts_frequency":
        #raise 
        pass

    if fill_empty:
        all_protein = self.universe.select_atoms('protein')
        for resname, resid in zip(list(all_protein.residues.resnames), list(all_protein.residues.resids)):
            if f"{str(resname)}{str(resid)}" not in list(self.analyses[analysis_name].result.keys()):
                self.analyses[analysis_name].result[f"{str(resname)}{str(resid)}"] = 0

        plt.figure(figsize=(len(self.analyses[analysis_name].result) * width_plot, 5))

        plt.bar(
            [int(res[3:]) for res in list(self.analyses[analysis_name].result.keys())],
            #list(self.analyses[analysis_name].result.keys()),
            [self.analyses[analysis_name].result[res] for res in list(self.analyses[analysis_name].result.keys())]
            #list(self.analyses[analysis_name].result.values())
        )
    
        plt.xticks(
        [int(res[3:]) for res in list(self.analyses[analysis_name].result.keys()) if self.analyses[analysis_name].result[res] != 0],
        [res for res in list(self.analyses[analysis_name].result.keys()) if self.analyses[analysis_name].result[res] != 0],
        rotation=90
        )

    elif not fill_empty:
        plt.figure(figsize=(len(self.analyses[analysis_name].result) * width_plot, 5))

        plt.bar(
            list(self.analyses[analysis_name].result.keys()),
            list(self.analyses[analysis_name].result.values())
        )


        plt.xticks(rotation=45, ha="right")

    

    plt.show()
    plt.close()




def ext_plot_contacts_frequencies_differences(
    contacts_ref,
    contacts_tgt,
    threshold=1,
    testing=False,
    remove_consequent=False,
    return_results=False,
    return_labels=False,
    width_plot=0.5,
    save_plot=False,
):
    """
    DESCRIPTION:
        Plotter that compares two different lists of protein contacts' frequencies and plots only those that are not similar based on a threshold.


    OPTIONS:
        - threshold:            sets the minimum difference that has to be between the reference and the target.
        - testing:              returns only the number of obtained contacts instead of the plot.
        - remove_consequent:    removes contacts between a residue and its following and previous.
        - return_results:       returns the dictionary with the compared results.
        - return_labels:        returns the list of labels of the obtained results.
        - width_plot:           sets the width per bar of the output plot.
        - save_plot:            pseudoboolean value to save the generated plot. If True, it will be stored as 'contacts_frequencies_diffs.png'. 
                                If a string is given, it will be used as the file name's with the .png extension.

    BUILDING NOTES:
        - < 0 --> more in tgt
        - > 0 --> more in ref
    """

    def calculate_average(contacts_dicts):
        # creates contacts_avg list with all keys and value 0
        contacts_avg = {}

        for contacts_dict in contacts_dicts:
            for k, value in contacts_dict.items():
                for k_, value_ in value.items():
                    if k in list(contacts_avg.keys()):
                        if k_ in list(contacts_avg[k].keys()):
                            contacts_avg[k][k_] += value_ / len(contacts_dicts)
                        else:
                            contacts_avg[k][k_] = value_ / len(contacts_dicts)
                    else:
                        contacts_avg[k] = {}
                        contacts_avg[k][k_] = value_ / len(contacts_dicts)

        return contacts_avg

    # calculate averages if more than one result is given as input
    if isinstance(contacts_ref, list):
        contacts_ref = calculate_average(contacts_ref)

    if isinstance(contacts_tgt, list):
        contacts_tgt = calculate_average(contacts_tgt)

    max_contact = 0
    labels_translator = {}
    contacts_ref_total = {}
    for k in list(contacts_ref.keys()):
        labels_translator[k[3:]] = k
        for k_, contacts in contacts_ref[k].items():
            contacts_ref_total[f"{k[3:]}-{k_[3:]}"] = contacts
            if max_contact < contacts:
                max_contact = contacts

    contacts_tgt_total = {}
    for k in list(contacts_tgt.keys()):
        for k_, contacts in contacts_tgt[k].items():
            contacts_tgt_total[f"{k[3:]}-{k_[3:]}"] = contacts
            if max_contact < contacts:
                max_contact = contacts

    important_contacts = {}
    for contact in list(
        set(contacts_ref_total.keys()) | set(contacts_tgt_total.keys())
    ):
        if contact not in list(contacts_ref_total.keys()) and contact in list(
            contacts_tgt_total.keys()
        ):  # and contacts_tgt_total[contact] > threshold:
            if contacts_tgt_total[contact] > threshold:
                important_contacts[contact] = -contacts_tgt_total[contact]

        elif contact in list(contacts_ref_total.keys()) and contact not in list(
            contacts_tgt_total.keys()
        ):  #
            if contacts_ref_total[contact] > threshold:
                important_contacts[contact] = contacts_ref_total[contact]

        else:
            if (
                abs(contacts_tgt_total[contact] - contacts_ref_total[contact])
                > threshold
            ):
                #                important_contacts[contact] = + contacts_tgt_total[contact] - contacts_ref_total[contact] # if more in tgt, positive
                important_contacts[contact] = (
                    -contacts_tgt_total[contact] + contacts_ref_total[contact]
                )  # if more in ref, positive

    if remove_consequent:
        for k in list(important_contacts.keys()):
            k_ = k.split("-")
            if int(k_[0]) + 1 == int(k_[1]):
                del important_contacts[k]
            elif int(k_[0]) == int(k_[1]) - 1:
                del important_contacts[k]


    # Remove redundant
    for k in list(important_contacts.keys()):
        for k_ in list(important_contacts.keys()):
            if k.split('-')[0] == k_.split('-')[1] and k.split('-')[1] == k_.split('-')[0]:
                del important_contacts[k]

    if len(important_contacts) == 0:
        print(
            "No contact has been found to be different between the two input sets. Change the threshold for plotting or redo the contacts' measure changing the radius."
        )
        return

    if testing:
        print(
            f"{len(important_contacts)} contacts have been detected with a threshold of {threshold}."
        )
        return

    col = [
        "red" if value < 0 else "blue" for value in list(important_contacts.values())
    ]

    plt.figure(figsize=(len(important_contacts) * width_plot, 5))

    plt.bar(
        [
            f"{labels_translator[k.split('-')[0]]}-{labels_translator[k.split('-')[1]]}"
            for k in list(important_contacts.keys())
        ],
        list(important_contacts.values()),
        color=col,
    )

    plt.xticks(rotation=45, ha="right")

    if max_contact == 100:
        plt.ylabel("Frequency (%)")

    else:
        plt.ylabel("Frequency (number of)")

    plt.show()

    if save_plot:
        plt.savefig('contacts_frequencies_diffs.png', dpi=300, bbox_inches='tight')

    elif isinstance(save_plot, str):
        plt.savefig(save_plot + '.png', dpi=300, bbox_inches='tight')

        

    plt.close()

    if return_results and return_labels:
        return important_contacts, return_labels
    elif return_results:
        return important_contacts
    elif return_labels:
        return return_labels


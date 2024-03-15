import matplotlib.pyplot as plt
from numpy import absolute as abs
from numpy import average, std, arange, array
import numpy as np

from .exceptions import NotCompatibleAnalysisForPlotterError, NotCompatibleMeasureForPlotterError
from .tools import get_dictionary_structure

"""
TO BUILD:
    - [] Distance vs contacts frequency/amounts as scatter
    - [] 
"""

def plot_measure(self, measure_name, same_y : bool = True, same_x : bool = True, axis_label_everywhere : bool =False, combine_replicas : bool =False, width_per_replica : float = 4, height_per_variant : float = 4, out_name=False):
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

    y_labels = {
        "distance" : "Distance (Å)",
        "RMSD" : "RMSD (Å)",
        "angle" : "Angle (°)",
        "planar_angle" : "Planar angle (°)",
        "dihedral" : "Dihedral angle (°)",
        "contacts_amount": "Number of contacts"
    }

    # Check if plotting as plotter or as class' method
    if measure_name == None:
        measure_obj = self
    else :
        measure_obj = self.measures[measure_name]

    
    if measure_obj.type not in ("distance", "angle", "dihedral", "RMSD", "planar_angle", "contacts_amount"):
        raise NotCompatibleMeasureForPlotterError
    if measure_obj.type == "contacts_amounts" and measure_obj.options["mode"] not in ("contacts"):
        raise NotCompatibleMeasureForPlotterError
    
    variants = len(measure_obj.result)

    if combine_replicas:
        max_replicas = 1
    else :
        max_replicas = max([ len(measure_obj.result[variant]) for variant in list(measure_obj.result) ])

    #fig, axs = plt.subplots(ncols=variants, nrows=max_replicas, sharey=same_y, sharex=same_x) --> axs[r_num, v_num]
    # plotting replicas in X axis and variant in Y axis
    fig, axs = plt.subplots(
        ncols=max_replicas, nrows=variants, 
        sharey=same_y, sharex=same_x, 
        figsize=(max_replicas*width_per_replica, variants*height_per_variant)
    )

    # Check if only one variant

    if variants == 1 and max_replicas == 1:
        variant = list(measure_obj.result.keys())[0]
        replica = list(measure_obj.result[variant].keys())[0]

        axs.plot(
            range(1, len(measure_obj.result[variant][replica])+1),
            measure_obj.result[variant][replica],
            c = f"C0",
            label=replica
            )
        
        axs.set_ylabel(y_labels[measure_obj.type])
        axs.set_xlabel("Frame")
    
    else :
        # If there's only one replica, treat it as merged
        if max_replicas == 1:
            combine_replicas = True

        for v_num, variant in enumerate(list(measure_obj.result.keys())):
            for r_num, replica in enumerate(list(measure_obj.result[variant].keys())):

                if combine_replicas:
                    axs[r_num].plot(
                    range(1, len(measure_obj.result[variant][replica])+1),
                    measure_obj.result[variant][replica],
                    c = f"C{v_num}",
                    label=replica
                    )

                    if r_num == 0 or axis_label_everywhere:
                        axs[r_num].set_ylabel(y_labels[measure_obj.type])

                    if v_num == variants-1 or axis_label_everywhere:
                        axs[r_num].set_xlabel("Frame")
                        
                    axs[r_num].set_title(f"{variant}, replicas {', '.join(list(measure_obj.result[variant].keys()))}")

                elif variants == 1:
                    axs[r_num].plot(
                    range(1, len(measure_obj.result[variant][replica])+1),
                    measure_obj.result[variant][replica],
                    c = f"C{v_num}",
                    label=replica
                    )

                    if r_num == 0 or axis_label_everywhere:
                        axs[r_num].set_ylabel(y_labels[measure_obj.type])

                    if v_num == variants-1 or axis_label_everywhere:
                        axs[r_num].set_xlabel("Frame")
                        
                    axs[r_num].set_title(f"{variant}, {replica}")  

                else :
                    axs[v_num, r_num].plot(
                        range(1, len(measure_obj.result[variant][replica])+1),
                        measure_obj.result[variant][replica],
                        c = f"C{v_num}"
                        )
                    
                    if r_num == 0 or axis_label_everywhere:
                        axs[v_num, r_num].set_ylabel(y_labels[measure_obj.type])
                    
                    if v_num == variants-1 or axis_label_everywhere:
                        axs[v_num, r_num].set_xlabel("Frame")
                    
                    axs[v_num, r_num].set_title(f"{variant}, {replica}")

            if r_num != max_replicas-1:
                for r_num_ in range(r_num, max_replicas):
                    if v_num == variants-1 or axis_label_everywhere:
                        axs[v_num, r_num_].plot()
                        axs[v_num, r_num_].set_xlabel("Frame")

            


    if measure_name == None:
        fig.suptitle("Plots for " + r"$\bf{%s}$" % self.name.replace('_', '\_') +  " Measure")

    else :
        fig.suptitle("Plots for " + r"$\bf{%s}$" % measure_name.replace('_', '\_') +  " Measure")

    
    fig.tight_layout()

    #if merge_replica_plots:
    #    plt.legend()

    if out_name != False and isinstance(out_name, str):
        if not out_name.endswith((".png", ".jpg", ".jpeg", ".tiff")):
            out_name = ".".join(out_name.split('.')[:-1]) + ".png"

        plt.savefig(out_name, dpi=300, bbox_inches='tight')

    plt.show()
    plt.close()



def plot_NACs(self, analysis_name, merge_replicas=False, percentage=False, error_bar=True, bar_width=0.1, width=None, title=None, out_name=False):
    """
    DESCRIPTION:
        Function for plotting NACs (or value-type Analysis) as a bar plot where all the variants are compared
    """

    # Check if plotting as plotter or as class' method
    if analysis_name == None:
        analysis_obj = self
    else :
        analysis_obj = self.analyses[analysis_name]

    if analysis_obj.type not in ('value', 'NACs'):
        raise NotCompatibleMeasureForPlotterError

    # create pyplot subplot
    fig, ax = plt.subplots(ncols=1,nrows=1)
    
    # Check width
    if width == None:
        fig.set_figwidth(10*bar_width*len(analysis_obj.result))
    elif isinstance(width, (float, int)):
        fig.set_figwidth(width)
    elif width.lower() in ('auto', 'automatic'):
        pass
    else :
        print("width value ({width}) cannot be understood, using default value.")

    # run code if merging replicas
    if merge_replicas:
        max_replicas = max([ len(analysis_obj.result[variant]) for variant in list(analysis_obj.result) ])
        for v_num, variant in enumerate(list(analysis_obj.result.keys())):
            # calculate the avg value for the bar (it height)
            if percentage:
                avgs = average([analysis_obj.result[variant][replica].count(True)*100/len(analysis_obj.result[variant][replica]) for replica in list(analysis_obj.result[variant].keys())])
            elif not percentage:
                avgs = average([analysis_obj.result[variant][replica].count(True) for replica in list(analysis_obj.result[variant].keys())])

            # plot a bar for each replica
            ax.bar(variant, avgs, bar_width*max_replicas, color = f"C{v_num}")
            if error_bar and max_replicas != 1:
                ax.errorbar(variant, avgs, 
                        yerr=std([analysis_obj.result[variant][replica].count(True) for replica in list(analysis_obj.result[variant].keys())]), 
                        color = f"k",
                        solid_capstyle='butt',
                        capsize=bar_width*72 # in to pt is 72, 5 is to make it wider
                        )

    elif not merge_replicas: 
        # Step 1: Process data
        # Counting True values for each replica in each variant
        true_counts = {variant: [sum(replica) for replica in replicas.values()] for variant, replicas in analysis_obj.result.items()}
        all_counts  = {variant: [len(replica) for replica in replicas.values()] for variant, replicas in analysis_obj.result.items()}

        # Find the maximum number of replicas to standardize the data structure
        max_replicas = max(len(replicas) for replicas in true_counts.values())

        # Initialize a list to hold the count of Trues for each replica in each variant
        counts_by_replica = {f'R{i+1}': [] for i in range(max_replicas)}

        # Populate the counts by iterating through each variant and replica
        v_num = 0
        for variant in true_counts.values():
            for i, count in enumerate(variant):
                if percentage:
                    counts_by_replica[f'R{i+1}'].append(count*100/all_counts[list(true_counts.keys())[v_num]][i])
                else :
                    counts_by_replica[f'R{i+1}'].append(count)
            v_num += 1
            # If a variant has fewer replicas, append 0 to ensure equal length lists
            for i in range(len(variant), max_replicas):
                counts_by_replica[f'R{i+1}'].append(0)

        # Step 2: Organize data for plotting
        labels = list(true_counts.keys())  # Variant names
        num_variants = len(true_counts)
        x = arange(num_variants)  # the label locations
    

        for i, (replica, counts) in enumerate(counts_by_replica.items()):
            ax.bar(x + i*bar_width, counts, bar_width, label=replica)

        # Add some text for labels, title, and custom x-axis tick labels, etc.
            
        ax.set_xticks(x + bar_width * (max_replicas / 2 - 0.5))
        ax.set_xticklabels(labels)
       
        ax.legend()


    if percentage:
        ax.set_ylim([0,100])

    
    if title == None:
        ax.set_title('Plot for ' + r"$\bf{%s}$" % analysis_obj.name.replace('_', '\_') + " Analysis")
    elif title == '':
        pass
    else :
        ax.set_title(title)

    if percentage:
        ax.set_ylabel('Percentage of NACs (%)')
        ax.set_yticks(range(0,101, 10))

    elif not percentage:
        ax.set_ylabel('Number of NACs')
    
#    if not merge_replicas:
#        ax.set_xticks(x + bar_width * (max_replicas / 2 - 0.5))
#        ax.set_xticklabels(labels)
        
    ax.grid(axis='x')
    ax.set_xlabel('Variant')
        
    if out_name != False and isinstance(out_name, str):
        if not out_name.endswith((".png", ".jpg", ".jpeg", ".tiff")):
            out_name = ".".join(out_name.split('.')[:-1]) + ".png"

        plt.savefig(out_name, dpi=300, bbox_inches='tight')

    plt.show()
    plt.close()
            


def plot_contacts_frequency(
        #self, analysis_name, fill_empty=False, width_plot=0.5, out_name=None
        self, analysis_name, same_y : bool = True, same_x : bool = True, axis_label_everywhere : bool =False, merge_replicas : bool =False, error_bar=True, bar_width=0.8, errorbar_width=5 , width_per_replica : float = 4, height_per_variant : float = 4, residue_label_rotation=45, out_name=False
        ):
    """
    DESCRIPTION:
        Plotter that takes the result of a contacts_frequency analysis and plots each interaction as a bar plot

    TODO:
        - [ ] Add normalisation to max per replica
    """


    #y_labels = {
        #"distance" : "Distance (Å)",
        #"RMSD" : "RMSD (Å)",
        #"angle" : "Angle (°)",
        #"planar_angle" : "Planar angle (°)",
        #"dihedral" : "Dihedral angle (°)"
    #}

    # Check if plotting as plotter or as class' method
    if analysis_name == None:
        analysis_obj = self
    else :
        analysis_obj = self.analyses[analysis_name]

    if analysis_obj.type not in ("contacts_frequency") and analysis_obj.options['mode'] not in ('contacts'):
        raise NotCompatibleAnalysisForPlotterError

    variants = len(analysis_obj.result)
    max_replicas = max([ len(analysis_obj.result[variant]) for variant in list(analysis_obj.result) ])

    if merge_replicas:
        # if only one replica, no merge is possible.
        if max_replicas == 1:
            merge_replicas = False

        else:
            max_replicas = 1
            merged_replicas = {}
            for variant in list(analysis_obj.result.keys()):
                merged_replicas[variant] = {}
                for replica in list(analysis_obj.result[variant].keys()):
                    for residue in list(analysis_obj.result[variant][replica].keys()):
                        merged_replicas[variant][residue] = []

            for variant in list(merged_replicas.keys()):
                for residue in list(merged_replicas[variant].keys()):
                    for replica in list(analysis_obj.result[variant].keys()):
                        if residue not in list(analysis_obj.result[variant][replica].keys()):
                            merged_replicas[variant][residue].append(0)
                        else :
                            merged_replicas[variant][residue].append(analysis_obj.result[variant][replica][residue])

            avgs, errs = {}, {}
            for variant in list(merged_replicas.keys()):
                avgs[variant], errs[variant] = {}, {}

                for residue in merged_replicas[variant]:
                    avgs[variant][residue] = average(merged_replicas[variant][residue])
                    errs[variant][residue] =     std(merged_replicas[variant][residue])

    #fig, axs = plt.subplots(ncols=variants, nrows=max_replicas, sharey=same_y, sharex=same_x) --> axs[r_num, v_num]
    # plotting replicas in X axis and variant in Y axis
    fig, axs = plt.subplots(
        ncols=max_replicas, nrows=variants, 
        sharey=same_y, sharex=same_x, 
        figsize=(max_replicas*width_per_replica, variants*height_per_variant)
    )

    # Check if only one variant
    if variants == 1 and max_replicas == 1:
        variant = list(analysis_obj.result.keys())[0]
        replica = list(analysis_obj.result[variant].keys())[0]

        axs.bar(
            list(analysis_obj.result[variant][replica].keys()),
            list(analysis_obj.result[variant][replica].values()),
            color = f"C0",
            #label=replica
            )
        
        #axs.set_ylabel(y_labels[analysis_obj.type])
        axs.set_xlabel("Residue")
        #axs.set_xticks(rotation=45)
    
    else :

        for v_num, variant in enumerate(list(analysis_obj.result.keys())):
            for r_num, replica in enumerate(list(analysis_obj.result[variant].keys())):

                if merge_replicas:
                    axs[v_num].bar(
                        list(avgs[variant].keys()),
                        list(avgs[variant].values()),
                        bar_width,
                        color = f"C{v_num}",
                    )

                    if error_bar:
                        #if bar_width < 1:
                        #    capsize = bar_width*15
                        #elif bar_width >= 1:
                        #    capsize = bar_width*2
                        
                        axs[v_num].errorbar(
                            list(avgs[variant].keys()),
                            list(avgs[variant].values()),
                            yerr=list(errs[variant].values()),
                            xerr=None,
                            ls='none',
                            color = f"k",
                            solid_capstyle='butt',
                            capsize=errorbar_width #capsize #bar_width*5 #72/2 # in to pt is 72, 5 is to make it wider
                        )

                    #if r_num == 0 or axis_label_everywhere:
                    #    axs[v_num].set_ylabel(y_labels[analysis_obj.type])

                    if v_num == variants-1 or axis_label_everywhere:
                        axs[v_num].set_xlabel("Residue")
                        
                    #axs[v_num].set_xticks(rotation=45)
                    axs[v_num].set_title(f"{variant}, average of replicas {', '.join(list(analysis_obj.result[variant].keys()))}")  

                elif max_replicas == 1:
                    axs[v_num].bar(
                    list(analysis_obj.result[variant][replica].keys()),
                    list(analysis_obj.result[variant][replica].values()),
                    bar_width,
                    color = f"C{v_num}",
                    label=replica
                    )

                    #if r_num == 0 or axis_label_everywhere:
                    #    axs[v_num].set_ylabel(y_labels[analysis_obj.type])

                    if v_num == variants-1 or axis_label_everywhere:
                        axs[v_num].set_xlabel("Residue")
                        
                    #axs[v_num].set_xticks(rotation=45)
                    axs[v_num].set_title(f"{variant}, replicas {', '.join(list(analysis_obj.result[variant].keys()))}")  

                elif variants == 1:
                    axs[r_num].bar(
                        list(analysis_obj.result[variant][replica].keys()),
                        list(analysis_obj.result[variant][replica].values()),
                        bar_width,
                        color = f"C{v_num}",
                        #label=replica
                    )

                    #if r_num == 0 or axis_label_everywhere:
                    #    axs[r_num].set_ylabel(y_labels[analysis_obj.type])

                    if v_num == variants-1 or axis_label_everywhere:
                        axs[r_num].set_xlabel("Residue")
                        
                    axs[r_num].set_title(f"{variant}, {replica}")  

                else :
                    axs[v_num, r_num].bar(
                        list(analysis_obj.result[variant][replica].keys()),
                        list(analysis_obj.result[variant][replica].values()),
                        bar_width,
                        color = f"C{v_num}"
                        )
                    
                    #if r_num == 0 or axis_label_everywhere:
                    #    axs[v_num, r_num].set_ylabel(y_labels[analysis_obj.type])
                    
                    if v_num == variants-1 or axis_label_everywhere:
                        axs[v_num, r_num].set_xlabel("Residue")
                    
                    axs[v_num, r_num].set_title(f"{variant}, {replica}")

            if r_num != max_replicas-1:
                for r_num_ in range(r_num, max_replicas):
                    if v_num == variants-1 or axis_label_everywhere:
                        axs[v_num, r_num_].plot()
                        axs[v_num, r_num_].set_xlabel("Residue")

            
    fig.autofmt_xdate(rotation=residue_label_rotation)


    if analysis_name == None:
        fig.suptitle("Plots for " + r"$\bf{%s}$" % self.name.replace('_', '\_') +  " Measure")

    else :
        fig.suptitle("Plots for " + r"$\bf{%s}$" % analysis_name.replace('_', '\_') +  " Measure")

    if analysis_obj.options['percentage']:
        plt.setp(axs, ylim=(0,100))
    
    fig.tight_layout()

    #if merge_replica_plots:
    #    plt.legend()

    if out_name != False and isinstance(out_name, str):
        if not out_name.endswith((".png", ".jpg", ".jpeg", ".tiff")):
            out_name = ".".join(out_name.split('.')[:-1]) + ".png"

        plt.savefig(out_name, dpi=300, bbox_inches='tight')

    plt.show()
    plt.close()



def plot_probability_densities(self, analysis_name, plot_minima : bool = True, plot_measures : bool = True, same_y : bool = True, same_x : bool = True, axis_label_everywhere : bool = False, colorbar_everywhere : bool = True, width_per_replica : float = 4, height_per_variant : float = 4, color_map='RdBu_r', show_contour_lines : bool = True, set_names_in_axis = False,  levels_lines : int = 25, levels_fill : int = 25, out_name=False):
    """
    DESCRIPTION:
        Plotter for probability density maps.

    ARGUMENTS:
        - plot_minima:  adds minima to contour plots.

    TODO:
        - [ ] fix global colorbar
    """
    def lighten_color(color, amount=0.5):
        """
        Lightens the given color by multiplying (1-luminosity) by the given amount.
        Input can be matplotlib color string, hex string, or RGB tuple.

        Examples:
        >> lighten_color('g', 0.3)
        >> lighten_color('#F034A3', 0.6)
        >> lighten_color((.3,.55,.1), 0.5)
        """
        import matplotlib.colors as mc
        import colorsys
        try:
            c = mc.cnames[color]
        except:
            c = color
        c = colorsys.rgb_to_hls(*mc.to_rgb(c))
        return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

    axis_labels = {
        "distance" : "Distance (Å)",
        "RMSD" : "RMSD (Å)",
        "angle" : "Angle (°)",
        "planar_angle" : "Planar angle (°)",
        "dihedral" : "Dihedral angle (°)",
        "contacts_amount": "Number of contacts"
    }

    if plot_minima and plot_measures:
        scatter_colors = {
            'minima'   : 'w',
            'measures' : 'k'
        }

    elif plot_minima and not plot_measures:
        scatter_colors = {
            'minima'   : 'k',
        }
    
    elif not plot_minima and plot_measures:
        scatter_colors = {
            'measures' : 'k'
        }

    # Check if plotting as plotter or as class' method
    if analysis_name == None:
        analysis_obj = self
    else :
        analysis_obj = self.analyses[analysis_name]
    
    if analysis_obj.type not in ("pdf"):
        raise NotCompatibleAnalysisForPlotterError
    

    if isinstance(set_names_in_axis, bool):
        if set_names_in_axis:
            x_axis = analysis_obj.options["selection_names"][0] + ' ' + axis_labels[analysis_obj.options["measure_types"][0]]
            y_axis = analysis_obj.options["selection_names"][1] + ' ' + axis_labels[analysis_obj.options["measure_types"][1]]

        else :
            x_axis = axis_labels[analysis_obj.options["measure_types"][0]]
            y_axis = axis_labels[analysis_obj.options["measure_types"][1]]
    
    elif isinstance(set_names_in_axis, list):
        x_axis = set_names_in_axis[0] + ' ' + axis_labels[analysis_obj.options["measure_types"][0]]
        y_axis = set_names_in_axis[1] + ' ' + axis_labels[analysis_obj.options["measure_types"][1]]
        

    variants = len(analysis_obj.result)

    max_replicas = max([ len(analysis_obj.result[variant]) for variant in list(analysis_obj.result) ])

    if self.analyses[analysis_name].options['merge_replicas']:
        gray_scale = 1/max_replicas
        
        max_replicas = 1



    #fig, axs = plt.subplots(ncols=variants, nrows=max_replicas, sharey=same_y, sharex=same_x) --> axs[r_num, v_num]
    # plotting replicas in X axis and variant in Y axis
    fig, axs = plt.subplots(
        ncols=max_replicas, nrows=variants, 
        sharey=same_y, sharex=same_x, 
        figsize=(max_replicas*width_per_replica, variants*height_per_variant)
    )

    # Check if only one variant

    if variants == 1 and max_replicas == 1:
        variant = list(analysis_obj.result.keys())[0]
        replica = list(analysis_obj.result[variant].keys())[0]

        if show_contour_lines:
            axs.tricontour(
                np.array(analysis_obj.result[variant][replica]['lscape'])[:,0],
                np.array(analysis_obj.result[variant][replica]['lscape'])[:,1],
                np.array(analysis_obj.result[variant][replica]['lscape'])[:,4],
                levels = levels_lines, linewidths=0.5, colors='k')
        cntr = axs.tricontourf(
            np.array(analysis_obj.result[variant][replica]['lscape'])[:,0],
            np.array(analysis_obj.result[variant][replica]['lscape'])[:,1],
            np.array(analysis_obj.result[variant][replica]['lscape'])[:,4],
            levels = levels_fill, cmap=color_map)
        
        if plot_minima:
            axs.scatter(
                np.array(analysis_obj.result[variant][replica]['mins'])[:,1],
                np.array(analysis_obj.result[variant][replica]['mins'])[:,2],
                color=scatter_colors['minima'])

        if plot_measures:
            if self.analyses[analysis_name].options['merge_replicas']:
                for r_num, replica in enumerate(list(analysis_obj.result[variant].keys())):
                    axs.scatter(
                        self.measures[analysis_obj.measure_name[0]].result[variant][replica],
                        self.measures[analysis_obj.measure_name[1]].result[variant][replica],
                        color=lighten_color('black', amount=(1-r_num*gray_scale)))

            else :
                axs.scatter(
                    self.measures[analysis_obj.measure_name[0]].result[variant][replica],
                    self.measures[analysis_obj.measure_name[1]].result[variant][replica],
                    color=scatter_colors['measures'])
                        
        axs.set_xlabel(x_axis)
        axs.set_ylabel(y_axis)

        cbar = fig.colorbar(cntr, ax=axs, extend='max')
        cbar.set_label('E/RT')


    elif variants == 1:

        for r_num, replica in enumerate(list(analysis_obj.result[variant].keys())):
            if show_contour_lines:
                axs[r_num].tricontour(
                    np.array(analysis_obj.result[variant][replica]['lscape'])[:,0],
                    np.array(analysis_obj.result[variant][replica]['lscape'])[:,1],
                    np.array(analysis_obj.result[variant][replica]['lscape'])[:,4],
                    levels = levels_lines, linewidths=0.5, colors='k')
            cntr = axs[r_num].tricontourf(
                np.array(analysis_obj.result[variant][replica]['lscape'])[:,0],
                np.array(analysis_obj.result[variant][replica]['lscape'])[:,1],
                np.array(analysis_obj.result[variant][replica]['lscape'])[:,4],
                levels = levels_fill, cmap=color_map)
            
            if plot_measures:
                axs[r_num].scatter(
                    self.measures[analysis_obj.measure_name[0]].result[variant][replica],
                    self.measures[analysis_obj.measure_name[1]].result[variant][replica],
                    color=scatter_colors['measures'])

            if plot_minima:
                axs[r_num].scatter(
                    np.array(analysis_obj.result[variant][replica]['mins'])[:,1],
                    np.array(analysis_obj.result[variant][replica]['mins'])[:,2],
                    color=scatter_colors['minima'])

            if r_num == 0 or axis_label_everywhere:
                axs[r_num].set_ylabel(y_axis)

            if v_num == variants-1 or axis_label_everywhere:
                axs[r_num].set_xlabel(x_axis)
                
            axs[r_num].set_title(f"{variant}, {replica}")  

            if colorbar_everywhere:
                cbar = fig.colorbar(cntr, ax=axs[r_num], extend='max')
                cbar.set_label('E/RT')


    elif max_replicas == 1:
        r_num, replica = 0, 'R1'

        for v_num, variant in enumerate(list(analysis_obj.result.keys())):
            if show_contour_lines:
                axs[v_num].tricontour(
                    np.array(analysis_obj.result[variant][replica]['lscape'])[:,0],
                    np.array(analysis_obj.result[variant][replica]['lscape'])[:,1],
                    np.array(analysis_obj.result[variant][replica]['lscape'])[:,4],
                    levels = levels_lines, linewidths=0.5, colors='k')
            cntr = axs[v_num].tricontourf(
                np.array(analysis_obj.result[variant][replica]['lscape'])[:,0],
                np.array(analysis_obj.result[variant][replica]['lscape'])[:,1],
                np.array(analysis_obj.result[variant][replica]['lscape'])[:,4],
                levels = levels_fill, cmap=color_map)
            
            if plot_measures:
                axs[v_num].scatter(
                    self.measures[analysis_obj.measure_name[0]].result[variant][replica],
                    self.measures[analysis_obj.measure_name[1]].result[variant][replica],
                    color=scatter_colors['measures'])

            if plot_minima:
                axs[v_num].scatter(
                    np.array(analysis_obj.result[variant][replica]['mins'])[:,1],
                    np.array(analysis_obj.result[variant][replica]['mins'])[:,2],
                    color=scatter_colors['minima'])

            if r_num == 0 or axis_label_everywhere:
                axs[v_num].set_ylabel(y_axis)

            if v_num == variants-1 or axis_label_everywhere:
                axs[v_num].set_xlabel(x_axis)
                
            axs[v_num].set_title(f"{variant}, {replica}")  

            if colorbar_everywhere:
                cbar = fig.colorbar(cntr, ax=axs[v_num], extend='max')
                cbar.set_label('E/RT')




    else :
        for v_num, variant in enumerate(list(analysis_obj.result.keys())):
            for r_num, replica in enumerate(list(analysis_obj.result[variant].keys())):
                if show_contour_lines:
                    axs[v_num, r_num].tricontour(
                        np.array(analysis_obj.result[variant][replica]['lscape'])[:,0],
                        np.array(analysis_obj.result[variant][replica]['lscape'])[:,1],
                        np.array(analysis_obj.result[variant][replica]['lscape'])[:,4],
                        levels = levels_lines, linewidths=0.5, colors='k')
                    
                cntr = axs[v_num, r_num].tricontourf(
                    np.array(analysis_obj.result[variant][replica]['lscape'])[:,0],
                    np.array(analysis_obj.result[variant][replica]['lscape'])[:,1],
                    np.array(analysis_obj.result[variant][replica]['lscape'])[:,4],
                    levels = levels_fill, cmap=color_map)
                    
                if plot_measures:
                    axs[v_num, r_num].scatter(
                        self.measures[analysis_obj.measure_name[0]].result[variant][replica],
                        self.measures[analysis_obj.measure_name[1]].result[variant][replica],
                        color=scatter_colors['measures'])

                if plot_minima:
                    axs[v_num, r_num].scatter(
                        np.array(analysis_obj.result[variant][replica]['mins'])[:,1],
                        np.array(analysis_obj.result[variant][replica]['mins'])[:,2],
                        color=scatter_colors['minima'])

                if r_num == 0 or axis_label_everywhere:
                    axs[v_num, r_num].set_ylabel(y_axis)
                
                if v_num == variants-1 or axis_label_everywhere:
                    axs[v_num, r_num].set_xlabel(x_axis)
                
                if colorbar_everywhere:
                    cbar = fig.colorbar(cntr, ax=axs[v_num, r_num], extend='max')
                    cbar.set_label('E/RT')
                
                axs[v_num, r_num].set_title(f"{variant}, {replica}")


    if not colorbar_everywhere:
        cbar = fig.colorbar(cntr, ax=axs[v_num, r_num], extend='max')
        cbar.set_label('E/RT')

    if analysis_name == None:
        fig.suptitle("Plots for " + r"$\bf{%s}$" % self.name.replace('_', '\_') +  " Measure")

    else :
        fig.suptitle("Plots for " + r"$\bf{%s}$" % analysis_name.replace('_', '\_') +  " Measure")

    
    fig.tight_layout()

    if out_name != False and isinstance(out_name, str):
        if not out_name.endswith((".png", ".jpg", ".jpeg", ".tiff")):
            out_name = ".".join(out_name.split('.')[:-1]) + ".png"

        plt.savefig(out_name, dpi=300, bbox_inches='tight')

    plt.show()
    plt.close()



################# EXTERNAL PLOTTERS #################

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
            if (
                k.split("-")[0] == k_.split("-")[1]
                and k.split("-")[1] == k_.split("-")[0]
            ):
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

    if save_plot:
        plt.savefig("contacts_frequencies_diffs.png", dpi=300, bbox_inches="tight")

    elif isinstance(save_plot, str):
        plt.savefig(save_plot + ".png", dpi=300, bbox_inches="tight")

    plt.show()
    plt.close()

    if return_results and return_labels:
        return important_contacts, return_labels
    elif return_results:
        return important_contacts
    elif return_labels:
        return return_labels


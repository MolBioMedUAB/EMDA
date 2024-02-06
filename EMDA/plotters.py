import matplotlib.pyplot as plt
from numpy import absolute as abs

"""
TO BUILD:
    - [] Distance vs contacts frequency/amounts as scatter
    - [] 
"""

def plot_values(self, measure_name, out_name=None):

    plt.plot(self.measure[measure_name].result)

    plt.show()

    if out_name != None:
        pass


def ext_plot_contacts_frequencies_differences(contacts_ref, contacts_tgt, threshold=1, testing=False, remove_consequent=False, return_results=False, return_labels=False, width_plot=0.5):
    """
    DESCRIPTION:
        Plotter that compares two different lists of protein contacts' frequencies and plots only those that are not similar based on a threshold.


    OPTIONS:
        - threshold: sets the minimum difference that has to be between the reference and the target.

        
    TO-DO:
        - [x] add option to only print z

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
                            contacts_avg[k][k_] += value_/len(contacts_dicts)
                        else :
                            contacts_avg[k][k_] =  value_/len(contacts_dicts)
                    else :
                        contacts_avg[k] = {}
                        contacts_avg[k][k_] =  value_/len(contacts_dicts)


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
    for contact in list(set(contacts_ref_total.keys()) | set(contacts_tgt_total.keys())):
        if (contact not in list(contacts_ref_total.keys()) and contact in list(contacts_tgt_total.keys())):# and contacts_tgt_total[contact] > threshold:
            if contacts_tgt_total[contact] > threshold:
                important_contacts[contact] = -contacts_tgt_total[contact]

        elif (contact in list(contacts_ref_total.keys()) and contact not in list(contacts_tgt_total.keys())):# 
            if contacts_ref_total[contact] > threshold:
                important_contacts[contact] = contacts_ref_total[contact]
        
        else :
            if abs(contacts_tgt_total[contact] - contacts_ref_total[contact]) > threshold:
#                important_contacts[contact] = + contacts_tgt_total[contact] - contacts_ref_total[contact] # if more in tgt, positive
                important_contacts[contact] = - contacts_tgt_total[contact] + contacts_ref_total[contact] #if more in ref, positive

    if remove_consequent:
        for k in list(important_contacts.keys()):
            k_ = k.split('-')
            if int(k_[0])+1 == int(k_[1]):
                del important_contacts[k]
            elif int(k_[0]) == int(k_[1])-1:
                del important_contacts[k]

    if len(important_contacts) == 0:
        print("No contact has been found to be different between the two input sets. Change the threshold for plotting or redo the contacts' measure changing the radius.")
        return


    if testing:
        print(f"{len(important_contacts)} contacts have been detected with a threshold of {threshold}.")
        return

    col = ['red' if value < 0 else 'blue' for value in list(important_contacts.values())]


    plt.figure(figsize=(len(important_contacts)*width_plot, 5))

    plt.bar(
        [f"{labels_translator[k.split('-')[0]]}-{labels_translator[k.split('-')[1]]}" for k in list(important_contacts.keys())],
        list(important_contacts.values()),
        color = col,
    )
    
    plt.xticks(rotation=45, ha='right')

    if max_contact == 100:
        plt.ylabel('Frequency (%)')

    else :
        plt.ylabel('Frequency (number of)')

    plt.show()
    plt.close()

    if return_results and return_labels:
        return important_contacts, return_labels
    elif return_results:
        return important_contacts
    elif return_labels:
        return return_labels

    

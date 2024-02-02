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


def ext_plot_contacts_frequencies_differences(contacts_ref, contacts_tgt, threshold=1, testing=False, remove_consequent=False):
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

    labels_translator = {}
    contacts_ref_total = {}
    for k in list(contacts_ref.keys()):
        for k_, contacts in contacts_ref[k].items():
            contacts_ref_total[f"{k[3:]}-{k_[3:]}"] = contacts
            labels_translator[f"{k[3:]}-{k_[3:]}"] = f"{k}-{k_}"

    contacts_tgt_total = {}
    for k in list(contacts_tgt.keys()):
        for k_, contacts in contacts_tgt[k].items():
            contacts_tgt_total[f"{k[3:]}-{k_[3:]}"] = contacts_tgt[k][k_]

    important_contacts = {}
    for contact in list(set(contacts_ref_total.keys()) | set(contacts_tgt_total.keys())):
        if (contact not in list(contacts_ref_total.keys()) and contact in list(contacts_tgt_total.keys())) and contacts_tgt_total[contact] > threshold:
            important_contacts[contact] = -contacts_tgt_total[contact]

        elif (contact in list(contacts_ref_total.keys()) and contact not in list(contacts_tgt_total.keys())) and contacts_ref_total[contact] > threshold:
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


    if testing:
        print(f"{len(important_contacts)} contacts have been detected with a threshold of {threshold}.")
        return


    col = ['red' if value < 0 else 'blue' for value in list(important_contacts.values())]

    plt.bar(
        [labels_translator[k] for k in list(important_contacts.keys())],
        list(important_contacts.values()),
        color = col,
    )
    
    plt.xticks(rotation=45, ha='right')
    
    plt.show()
    plt.close()


    

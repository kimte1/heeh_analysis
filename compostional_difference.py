import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import collections



p1 = 'SLEELLKLAEEALKRGKTIRILGFEISSEEALRRFEEWLRRFI' # heeh420
p2 = 'DIEEIEKKARKILEKGDSIEIAGFEVRDEEDLKKILEWLRRHG' # heeh958




def residue_counts(string):
    ''' return dict of residue and counts '''
    unique = dict(zip(list(string),[list(string).count(i) for i in list(string)]))
    return(unique)


res_counts_p1=residue_counts(p1)
res_counts_p2=residue_counts(p2)



def fill_in_missing_res(p1, p2):
    ''' '''
    all_res = list(p1.keys())+list(p2.keys())
    unique_res = np.unique(all_res)

    missing_p1 = set(unique_res)^set(p1)
    missing_p2 = set(unique_res)^set(p2)

    for missing in missing_p1:
        p1[missing]=0
    for missing in missing_p2:
        p2[missing]=0

    p1 = collections.OrderedDict(sorted(p1.items()))
    p2 = collections.OrderedDict(sorted(p2.items()))
    return p1, p2


dict_res_p1, dict_res_p2=fill_in_missing_res(p1=res_counts_p1, p2=res_counts_p2)

dict_res_p1, dict_res_p2





def plot_counts():
    p1, p2 = fill_in_missing_res(p1=res_counts_p1, p2=res_counts_p2)


    plt.subplot(2, 1, 1)
    plt.bar(p1.keys(), p1.values())
    plt.subplot(2, 1, 2)
    plt.bar(p2.keys(), p2.values())
    plt.ylim(0, 10,2)
    #plt.savefig('../figures/figure_components/total_counts.svg')
    #return p1
plot_counts()



def abs_diff():
    p1, p2 = fill_in_missing_res(p1=res_counts_p1, p2=res_counts_p2)

    abs_diff_dict = {}

    for key in p1.keys():
        abs_diff_dict[key] = abs(p1[key]-p2[key])
    return abs_diff_dict

abs_diff = abs_diff()
abs_diff

def plot_abs_diff():

    plt.bar(abs_diff.keys(), abs_diff.values())
    plt.ylim(0,10)
    #plt.savefig('../figures/figure_components/abs_diff.svg')
plot_abs_diff()

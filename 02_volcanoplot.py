''' This script performs a binomial test for which residues are more likely to be found in rd5 stable designs '''



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import math


''' load files '''
def import_data(metrics_data, scores_data):
    ''' returns df with designs' metrics and stabilityscore data '''
    metrics = pd.read_csv('datasets/' + metrics_data)
    scores = pd.read_csv('datasets/' + scores_data)
    df = pd.merge(left = metrics, right = scores, left_on = 'name', right_on = 'name', how = 'inner')
    df = df.rename(columns = {'sequence_x': 'sequence'} )
    return(df)

df = import_data('rd5_metrics.csv', 'rd5_stability_scores.csv')



def volcano_plot():
    # make dict of aa seq
    aa_list = ['A','F','I','L','M','V','W','Y','G','P','Q','N','C','S','T','H','R','K','D','E']
    # make sub-df separating ss >= 1.0 and ss < 1.0
    s = df.query('stabilityscore > 1.0')
    no = df.query('stabilityscore <= 1.0')
    # make lists of sequences
    seq_s = list(s['sequence'])
    seq_no = list(no['sequence'])
    # create empty list with each element being a position's dictionary { aa : frequency}
    pos_s = []
    pos_no = []
    for i in range(43):
        pos_s.append(dict.fromkeys(aa_list, 0)) # create keys
        pos_no.append(dict.fromkeys(aa_list, 0))
    # populate empty list with aa frequencies
    for seq in seq_s: # for each design
        for i in range(43): # go through each position of the design
            aa = seq[i] # id the aa at each position
            pos_s[i][aa] +=1 # add the aa to the position list
    for seq in seq_no:
        for i in range(43):
            aa = seq[i]
            pos_no[i][aa] +=1
    to_add = []
    for position in range(len(pos_s)): # go through each position down the heatmap
        for residue in range(len(aa_list)): # go through each aa on the heatmap for a specific position (ie look at each square)
            name = aa_list[residue] + str(position+1) # create a name for each 'square' on the heatmap
            k = pos_s[position].get(aa_list[residue]) # get the number of success for each pos_aa
            k_no = pos_no[position].get(aa_list[residue]) # get the number of no success for each pos_aa
            k_total = k + k_no
            array = (k, k_no)
            if k_total == 0:
                pass
            else:
                obs_ss = k/(k_total) # probability of observed success rate = ( # of ss designs with aa X / total # designs with aa X )
                theor_ss = len(s)/len(df) # probability of theoretical success 2157/5799
                ''' fc is based on overall theoretical success '''
                fc = obs_ss/theor_ss
                try:
                    log2fc = math.log2(fc)
                    pv = stats.binom_test(x = array, p = theor_ss) # binomial test
                    log10pv = math.log(pv) * -1
                    row = [name, k, k_total, obs_ss, theor_ss, pv, log10pv, fc, log2fc]
                    to_add.append(row) # append list
                except ValueError:
                    pass
    # create df
    binom_df = pd.DataFrame(to_add, columns =[ 'name', 'k', 'total', 'obs_ss', 'theor_ss', 'pv', 'log10pv', 'fc', 'log2fc'])

    # plot
    plt.scatter(x = binom_df['log2fc'], y = binom_df['log10pv'], alpha = 0.3, edgecolors = 'none' ), plt.grid(axis = 'y', linewidth = 1)
    plt.xlabel('log2 fold change')
    plt.ylabel('-log10 pvalue')
    plt.savefig('volcanoplot.svg')


volcano_plot()

''' this script analyzes loops constriction to stability '''


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from scipy.stats import ranksums
from scipy.stats import ks_2samp
from matplotlib import rcParams
rcParams['font.family'] = 'arial'

pwd

''' load files '''
def load_data():
    comp = pd.read_csv('../../datasets/HEEH_KT_rd6_lowenergy_comp_all_filters_addtlfeatures_210402.csv')
    exp = pd.read_csv('../../datasets/HEEH_KT_rd6_lowenergy_exp_all_filters_addtlfeatures_210402.csv')
    df = pd.merge(left=comp, right=exp, on='name', how='inner')
    return df
df=load_data()




def middle_loop_information():
    ''' add len_L2 and residue_L2 info to df '''
    sec_structures = [  sec_str for sec_str in df['dssp'] ]
    sequences = [ seq for seq in df['sequence']]
    #stabilityscores = [ score for score in df['stabilityscore']]
    all_residues = []
    len_residues = []

    for structure, sequence in zip(sec_structures, sequences):
        loop2_start = structure.find('EL')+1 # find start position for the middle loop
        loop2_end = structure.find('LE', loop2_start)+1 # find end position for middle loop
        residues = sequence[loop2_start:loop2_end] # find residues for middle loop
        all_residues.append(residues)
        len_residues.append(len(residues))
    df['len_L2'] = len_residues
    df['res_L2'] = all_residues

    sub_df = df.query('len_L2 == 2')
    sub_df['hphob_1'] = [ 1 if res in 'AFILMWVY' else 0 for res in sub_df['res_L2'].str[0]]
    sub_df['hphob_2'] = [ 2 if res in 'AFILMWVY' else 0 for res in sub_df['res_L2'].str[1]]
    sub_df['hphob_1_2'] = [ hphob_1+hphob_2 for hphob_1, hphob_2 in zip(sub_df['hphob_1'], sub_df['hphob_2']) ]

    # there's only 23 designs whose hphob value is 3 which is both loop positions are hphob, so filter those out
    sub_df = sub_df.query('hphob_1_2 !=3')
    return sub_df

sub_df = middle_loop_information()

len(sub_df)
3772/3994

sns.distplot(sub_df.query('hphob_1==1'))
sub_df

sns.countplot(sub_df['hphob_1'])
len(sub_df.query('hphob_1_2==3'))
sns.countplot(sub_df['hphob_1_2'])



def boxplot_middle_loop():
    plt.rcParams['axes.axisbelow'] = True
    plt.grid(axis = 'y', zorder=-1)
    flierprops = dict(marker='o', markersize=2)
    #colors=sns.color_palette("seagreen")
    colors=sns.cubehelix_palette(start=.5, rot=-.75,)# as_cmap=True)

    sns.boxplot(x='hphob_1_2', y='stabilityscore', data=sub_df.query('len_L2 ==2'),
                palette=colors, flierprops=flierprops, zorder=10,
                showmeans=True,
                meanprops={'marker':'o', 'markeredgecolor':'black', 'markersize':5})

    # wilcoxon rank sums test
    zero = sub_df.query('hphob_1_2 == 0')
    one = sub_df.query('hphob_1_2 == 1')
    two = sub_df.query('hphob_1_2 == 2')

    #three = sub_df.query('hphob_1_2 == 3')
    rs_01 = ks_2samp(zero['stabilityscore'], one['stabilityscore'], )
    rs_02 = ks_2samp(zero['stabilityscore'], two['stabilityscore'], )
    #rs_03 = ranksums(zero['stabilityscore'], three['stabilityscore'], )
    rs_12 = ks_2samp(one['stabilityscore'], two['stabilityscore'], )
    #rs_13 = ranksums(one['stabilityscore'], three['stabilityscore'], )
    #rs_23 = ranksums(two['stabilityscore'], three['stabilityscore'], )
    #plt.savefig('../figures/rd6_loop_constriction_boxplot.svg')
    return rs_01, rs_02, rs_12
boxplot_middle_loop()


np.median(sub_df.query('hphob_1_2==0')['stabilityscore'])
np.median(sub_df.query('hphob_1_2==1')['stabilityscore'])
np.median(sub_df.query('hphob_1_2==2')['stabilityscore'])
#np.median(sub_df.query('hphob_1_2==3')['stabilityscore'])

sns.distplot(sub_df.query('hphob_1_2 == 0')['stabilityscore']), sns.distplot(sub_df.query('hphob_1_2 == 1')['stabilityscore']), sns.distplot(sub_df.query('hphob_1_2 == 2')['stabilityscore'])
np.mean(sub_df.query('hphob_1_2==0')['stabilityscore'])
np.mean(sub_df.query('hphob_1_2==1')['stabilityscore'])
np.mean(sub_df.query('hphob_1_2==2')['stabilityscore'])
#np.mean(sub_df.query('hphob_1_2==3')['stabilityscore'])
2/9

def loop_contriction(n_iterations):
    ''' bootstrap 1000 times and output success ratio and stabilityscores for design has 0, 1, or 2 hphob resiudes for designs with len_L2 == 2 , '''
    n_designs = len(sub_df)
    # bootstrapping
    iterations = []
    hphob_category = []
    hphob_counts = []
    ss_ratios = []
    mean_stabilities = []
    for iteration in range(n_iterations):
        sample = sub_df.sample(n=n_designs, replace=True)

        for i in range(3):
            sub_sample = sample.query('hphob_1_2 == @i ')
            mean_stabilities.append(np.mean(sample['stabilityscore']))
            stable_sub_sample = sub_sample.query('stabilityscore > 1')
            iterations.append(iteration+1)
            hphob_category.append(i)
            ss_ratios.append(len(stable_sub_sample)/len(sub_sample))
    bootstrap_df = pd.DataFrame()
    bootstrap_df['iteration'] = iterations
    bootstrap_df['hphob_category'] = hphob_category
    bootstrap_df['ss_ratio'] = ss_ratios
    bootstrap_df['mean_stability'] = mean_stabilities
    return bootstrap_df

bootstrap_df = loop_contriction(n_iterations=1000)

bootstrap_df
# save data
bootstrap_df.to_csv('../datasets/rd6_bootstrap_loop_constriction_data_no_bothhphob.csv', index=False, header=True )

#conf_int_df = pd.read_csv('../datasets/rd6_bootstrap_loop_constriction_data.csv')

bootstrap_df.head()

def plot_average_stability(n_iterations):

    # keep 95% conf int
    conf_int = 0.95
    margin = (1-conf_int)/2*n_iterations
    upper = (n_iterations-margin)
    lower = margin

    conf_int_df = pd.DataFrame()
    for hphob_cat in range(3):
        sub = bootstrap_df.query('hphob_category == @hphob_cat').sort_values(by = 'ss_ratio').reset_index().loc[lower:upper+1]
        conf_int_df = pd.concat([conf_int_df, sub])

    for i in range(3):
        ratio_min = min(conf_int_df.query('hphob_category==@i')['ss_ratio'])
        ratio_max = max(conf_int_df.query('hphob_category==@i')['ss_ratio'])
        plt.plot( [i, i], [ratio_min, ratio_max], color='black', zorder=10)
    plt.grid(axis='y')
    sns.barplot(x='hphob_category', y='ss_ratio', data=bootstrap_df, palette=['#D5DBC4', '#A0C29F','#2CA02C', '#4F7881'])
    plt.savefig('../figures/rd6_loopconstriction_barplot.svg')


    # ks 2 sample test
    #zero = conf_int_df.query('hphob_category == 0')
    #one = conf_int_df.query('hphob_category == 1')
    #two = conf_int_df.query('hphob_category == 2')
    #three = conf_int_df.query('hphob_category == 3')
    #ks_01 = ranksums(zero['mean_stability'], one['mean_stability'], )
    #ks_02 = ranksums(zero['mean_stability'], two['mean_stability'], )
    #ks_03 = ranksums(zero['mean_stability'], three['mean_stability'], )
    #ks_12 = ranksums(one['mean_stability'], two['mean_stability'], )
    #ks_13 = ranksums(one['mean_stability'], three['mean_stability'], )
    #ks_23 = ranksums(two['mean_stability'], three['mean_stability'], )
    #print(ks_01, ks_02, ks_03, ks_12, ks_13, ks_23)
    #plt.savefig('../figures/FigS16_rd6_loop_constriction_average_bootstrap.svg')
plot_average_stability(n_iterations=1000)



len(sub_df.query('hphob_1_2==0 & stabilityscore >1'))/len(sub_df.query('hphob_1_2==0'))
len(sub_df.query('hphob_1_2==1 & stabilityscore >1'))/len(sub_df.query('hphob_1_2==1'))
len(sub_df.query('hphob_1_2==2 & stabilityscore >1'))/len(sub_df.query('hphob_1_2==2'))
len(sub_df.query('hphob_1_2==3 & stabilityscore >1'))/len(sub_df.query('hphob_1_2==3'))
len(sub_df.query('hphob_1_2==3'))


def plot_success_ratio():
    #plt.grid(axis='y', zorder=1)

    for i in range(4):
        ratio_min = min(conf_int_df.query('hphob_category==@i')['ss_ratio'])
        ratio_max = max(conf_int_df.query('hphob_category==@i')['ss_ratio'])
        plt.plot( [i, i], [ratio_min, ratio_max], color='black', zorder=10)

    sns.barplot(x='hphob_category', y='ss_ratio', data=conf_int_df, palette=['#D5DBC4', '#A0C29F','#2CA02C', '#4F7881'])
    plt.savefig('../figures/FigS16_rd6_loopconstriction_barplot.svg')
plot_success_ratio()


np.mean(conf_int_df.query('hphob_category==0')['ss_ratio'])
np.mean(conf_int_df.query('hphob_category==1')['ss_ratio'])
np.mean(conf_int_df.query('hphob_category==2')['ss_ratio'])
np.mean(conf_int_df.query('hphob_category==3')['ss_ratio'])

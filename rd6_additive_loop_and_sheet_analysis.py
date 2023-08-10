''' analyze whether having a specific ABEGO pattern at each loop impacts stability in an additive way
 - develop a regression model analzing loop abegos and strand legnths '''


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import Lasso
from sklearn import metrics
from sklearn.utils import resample

from matplotlib import cm
from matplotlib.colors import DivergingNorm
import itertools
from itertools import product
from itertools import permutations
from itertools import combinations
import statistics

from matplotlib import rcParams
rcParams['font.family'] = 'arial'

pwd
''' load files '''
def load_data():
    comp = pd.read_csv('datasets/rd6_metrics_fig4b-g.csv')
    exp = pd.read_csv('datasets/rd6_stability_scores_fig4b-g.csv')
    df = pd.merge(left=comp, right=exp, on='name', how='inner')
    return df
df=load_data()




def filter_abegos():
    ''' identify and keep ABEGO patterns that are represented > 300 times '''
    n_loops = 3
    # identify all unique ABEGO sequences at each loop
    unique_abegos_L1 = np.unique(df['L1'])
    unique_abegos_L2 = np.unique(df['L2'])
    unique_abegos_L3 = np.unique(df['L3'])
    unique_abegos_all_loops = [unique_abegos_L1, unique_abegos_L2, unique_abegos_L3]
    loop_numbers = ['L1', 'L2', 'L3']

    # create list with all unique abego names and the how many of them there are
    all_abego_names = [ [] for i in range(n_loops)]
    n_all_abegos = [ [] for i in range(n_loops)]
    for i, abegos_at_particular_loop in zip(range(n_loops), unique_abegos_all_loops):
        loop_number = 'L'+str(i+1)
        for abego in (abegos_at_particular_loop):
            n_abegos = len(df [df[loop_number]==abego])
            n_all_abegos[i].append(n_abegos)
            all_abego_names[i].append(abego)

    # identify unique strands
    unique_strands = np.unique(df['len_strand1'])
    # keep strands if there are > 300 of them
    len_all_strands = [ [] for i in range(len(unique_strands))]
    filter_strands = []
    for len_strand in unique_strands:
        n_strand = len(df.query('len_strand1 == @len_strand'))
        if n_strand > 300:
            filter_strands.append(len_strand)

    # keep ABEGOs if there are > 300 of them
    filter_abego_names = [ [] for i in range(n_loops)]
    for i in range(n_loops):
        for abego, n in zip(all_abego_names[i], n_all_abegos[i]):
            if n > 300:
                filter_abego_names[i].append(abego)

    return filter_abego_names, filter_strands, df
filter_abego_names, filter_strands, df = filter_abegos()


def get_dummies():
    dummy_df = pd.get_dummies(df, columns = ['L1', 'L2', 'L3', 'len_strand1'])
    return filter_abego_names, pd.concat([df[['L1', 'L2', 'L3', 'len_strand1']], dummy_df], axis = 1)
filter_abego_names, df = get_dummies()







def regression_model_loops(n_iterations):
    ''' reg model with filtered abegos as variables '''
    n_loops = 3
    # get abego columns to analyze
    abegos_to_analyze = []
    for i in range(n_loops):
        for abego in filter_abego_names[i]:
            abegos_to_analyze.append('L'+str(i+1)+'_'+abego)

    strands_to_analyze = [ 'len_strand1_'+str(strand_num) for strand_num in filter_strands]
    feat_to_analyze = abegos_to_analyze + strands_to_analyze

    # bootstrap linear regresesion
    n_designs = len(df)
    summary_df = pd.DataFrame()

    for iteration in range(n_iterations):
        print('iteration ', iteration)
        sample = df.sample(n=n_designs, replace=True)

        X = sample[feat_to_analyze]
        y = sample['stabilityscore']
        lr = Lasso(alpha = 0.00000001)
        lr.fit(X, y)
        pred = lr.predict(X)

        sub_df = pd.DataFrame()
        sub_df['iteration'] = [iteration]*len(feat_to_analyze)
        sub_df['feature'] = feat_to_analyze
        sub_df['coef'] = lr.coef_
        sub_df['rmse'] = np.sqrt(metrics.mean_squared_error(pred, y))

        # adjust coefficient magnitude
        l1_coef_values = [ coef for coef, feat in zip(sub_df['coef'], sub_df['feature']) if feat[:2] == 'L1' ]
        l2_coef_values = [ coef for coef, feat in zip(sub_df['coef'], sub_df['feature']) if feat[:2] == 'L2' ]
        l3_coef_values = [ coef for coef, feat in zip(sub_df['coef'], sub_df['feature']) if feat[:2] == 'L3' ]
        strand_coef_values = [ coef for coef, feat in zip(sub_df['coef'], sub_df['feature']) if feat[:3] == 'len' ]

        avg_l1 = np.average(l1_coef_values)
        avg_l2 = np.average(l2_coef_values)
        avg_l3 = np.average(l3_coef_values)
        avg_strand = np.average(strand_coef_values)

        new_l1_coef_val = [ val-avg_l1 for val in l1_coef_values ]
        new_l2_coef_val = [ val-avg_l2 for val in l2_coef_values ]
        new_l3_coef_val = [ val-avg_l3 for val in l3_coef_values ]
        new_strand_coef_val = [ val-avg_strand for val in strand_coef_values ]

        all_new_val = new_l1_coef_val + new_l2_coef_val + new_l3_coef_val + new_strand_coef_val
        sub_df['adj_coef'] = all_new_val

        summary_df = pd.concat([summary_df, sub_df])
    return summary_df, feat_to_analyze,

summary_df, feat_to_analyze, = regression_model_loops(n_iterations=1000)


def conf_int(n_iterations):
    # keep 95% conf int
    conf_int = 0.95
    margin = (1-conf_int)/2*n_iterations
    upper = (n_iterations-margin)
    lower = margin

    conf_int_df = pd.DataFrame()
    for feat in feat_to_analyze:
        sub = summary_df.query('feature == @feat').sort_values(by = 'adj_coef').reset_index().loc[lower:upper+1]
        conf_int_df = pd.concat([conf_int_df, sub])
    return conf_int_df
conf_int_df = conf_int(n_iterations=1000)



def barplot_magnitudes_coeff(n_iterations):
    ''' barplot showing the top topological features from the regression model '''
    plt.grid(axis = 'x', zorder=-10)
    plt.xticks(rotation = 45)
    orders = ['L1_GBB', 'L1_AGBB', 'L2_GG', 'L2_EA', 'L3_AB', 'L3_BAAB', 'L3_BBGB', 'len_strand1_3', 'len_strand1_5']
    for feat, i in zip(orders, range(len(orders))):
        feat_min = min(conf_int_df.query('feature==@feat')['adj_coef'])
        feat_max = max(conf_int_df.query('feature==@feat')['adj_coef'])
        plt.plot([feat_min, feat_max], [i, i], color='black',  zorder=10)

    sns.barplot(x='adj_coef', y='feature', order=orders, color='#5AA095', data=summary_df, zorder=10, ci=None) # flierprops=flierprops,

barplot_magnitudes_coeff(1000)









def coeff_difference(n_iterations):

    summary_df.sort_values(by=['iteration', 'feature'], inplace=True)
    L1_diff = []
    L2_diff = []
    L3_diff1 = []
    L3_diff2 = []
    L3_diff3 = []
    strand_diff1 = []
    strand_diff2 = []
    strand_diff3 = []
    iteration = []
    feature=[]
    var=[]

    for i in range(n_iterations):
        sub = summary_df.query('iteration==@i')
        print(i)

        L1_diff = float(sub.query('feature=="L1_GBB"')['adj_coef']) - float(sub.query('feature=="L1_AGBB"')['adj_coef'])
        L2_diff = float(sub.query('feature=="L2_GG"')['adj_coef']) - float(sub.query('feature=="L2_EA"')['adj_coef'])
        L3_diff1 = float(sub.query('feature=="L3_AB"')['adj_coef']) - float(sub.query('feature=="L3_BAAB"')['adj_coef']) # AB - BAAB
        L3_diff2 = float(sub.query('feature=="L3_AB"')['adj_coef']) - float(sub.query('feature=="L3_BBGB"')['adj_coef']) # AB - BBGB
        L3_diff3 = float(sub.query('feature=="L3_BAAB"')['adj_coef']) - float(sub.query('feature=="L3_BBGB"')['adj_coef']) # BAAB - BBGB
        strand_diff1 = float(sub.query('feature=="len_strand1_3"')['adj_coef']) - float(sub.query('feature=="len_strand1_4"')['adj_coef']) # strand 3 - 4
        strand_diff2 = float(sub.query('feature=="len_strand1_3"')['adj_coef']) - float(sub.query('feature=="len_strand1_5"')['adj_coef']) # strand 3 - 5
        strand_diff3 = float(sub.query('feature=="len_strand1_4"')['adj_coef']) - float(sub.query('feature=="len_strand1_5"')['adj_coef']) # strand 4 - 5

        iteration.extend([i, i, i, i, i, i, i, i])
        feature.extend(['L1_GBB-L1_AGBB', 'L2_GG-L2_EA', 'L3_AB-L3_BAAB', 'L3_AB-L3_BBGB', 'L3_BAAB-L3_BBGB',
                        'strand3-4', 'strand3-5', 'strand4-5'])
        var.extend([L1_diff, L2_diff, L3_diff1, L3_diff2, L3_diff3, strand_diff1, strand_diff2, strand_diff3])

    newdf = pd.DataFrame()
    newdf['iteration'] = iteration
    newdf['feature'] = feature
    newdf['coef_diff'] = var
    return newdf

newdf = coeff_difference(n_iterations=1000)



def conf_int(n_iterations):
    # keep 95% conf int
    conf_int = 0.95
    margin = (1-conf_int)/2*n_iterations
    upper = (n_iterations-margin)
    lower = margin

    unique_feats = np.unique(newdf['feature'])
    conf_int_df = pd.DataFrame()
    for feat in unique_feats:
        sub = newdf.query('feature == @feat').sort_values(by = 'coef_diff').reset_index().loc[lower:upper+1]
        conf_int_df = pd.concat([conf_int_df, sub])
    return conf_int_df
conf_int_df = conf_int(n_iterations=1000)




def plot_coeff_magnitudes():

    orders = ['L1_GBB-L1_AGBB', 'L2_GG-L2_EA', 'L3_AB-L3_BBGB', 'L3_AB-L3_BAAB', 'L3_BAAB-L3_BBGB', 'strand3-5', 'strand4-5', 'strand3-4']
    #95ci
    for feat, i in zip(orders, range(len(orders))):
        feat_min = min(conf_int_df.query('feature==@feat')['coef_diff'])
        feat_max = max(conf_int_df.query('feature==@feat')['coef_diff'])
        plt.plot([feat_min, feat_max], [i, i], color='black', linewidth=1,  zorder=5)

    # 1 std
    for feat, i in zip(orders, range(len(orders))):
        sd = statistics.stdev(newdf.query('feature==@feat')['coef_diff'])
        feat_min = np.mean(newdf.query('feature==@feat')['coef_diff'])-sd
        feat_max = np.mean(newdf.query('feature==@feat')['coef_diff'])+sd
        plt.plot([feat_min, feat_max], [i, i], color='black', linewidth=3,  zorder=10)


    plt.rcParams['axes.axisbelow'] = True
    sns.barplot(x='coef_diff', y='feature', data=newdf, order=orders, color='#5AA095', ci=None)#, zorder=10)
    #sns.boxplot(x='feature', y='coef_diff', data=conf_int_df, color='#5AA095', zorder=10)
    plt.grid(axis='x', zorder=-1)

    plt.xticks(rotation=90,)
    #plt.savefig('../figures/rd6_abego_strand_bootstrap_barplot_diff.svg')
plot_coeff_magnitudes()



def keep_designs_abego_combinations():
    ''' keep designs whose abego pattern combination and strand length > 100 examples '''
    all_abego_combinations = list(product(*filter_abego_names))
    df2 = pd.DataFrame()
    for abego_combo in all_abego_combinations:
        for strand in filter_strands:
            abego = abego_combo[0] + '-' + abego_combo[1] + '-' + abego_combo[2]

            sub_df = df.query('abego_loop ==@abego & len_strand1==@strand')
            if len(sub_df) > 100:
                df2 = pd.concat([df2, sub_df])
    df2['abego_strand'] = [ abego+'-'+str(len_strand) for abego, len_strand in zip(df2['abego_loop'], df2['len_strand1'])]
    return df2
sub_df = keep_designs_abego_combinations()



def boxplot():
    ''' outputs boxplot of abego patterns vs exp score overlayed with a scatterplot of abego patterns vs pred score '''
    # create Lasso regression and add pred_ss data
    X = sub_df[feat_to_analyze]
    y = sub_df['stabilityscore']
    lr = Lasso(alpha = 0.00000001)
    lr.fit(X, y)
    pred = lr.predict(X)
    sub_df['pred_ss'] = pred

    pred_ss_df = pd.DataFrame()
    pred_ss_df['abego_strand'] = np.unique(sub_df['abego_strand'])
    pred_ss_df['pred_ss'] = np.unique(sub_df['pred_ss'])

    # boxplot
    cm = 1/2.54  # centimeters in inches
    plt.figure(figsize=(18*cm, 19.68*cm))
    plt.rcParams['axes.axisbelow'] = True
    plt.grid(axis = 'x', zorder=-1)
    flierprops = dict(marker='o', markersize=1)
    sns.boxplot(x='stabilityscore', y='abego_strand', data=sub_df.sort_values('pred_ss', ascending=True), color='#63978F',  flierprops=flierprops, zorder=10)
    sns.scatterplot(x='pred_ss', y='abego_strand', data=pred_ss_df.sort_values('pred_ss', ascending=True), color='#FFC107', edgecolor=None, zorder=10)
    plt.xlim([0,4])
    #plt.savefig('../figures/rd6_abego_strand_bootstrap_barplot.svg')
boxplot()




def success_rate_btsp(n_iterations):
    ''' determine success ratio by bootstrapping and keep 95% conf int '''
    #lr_results_df, feat_to_analyze, filter_abego_names, unique_strand_lengths, df = regression_model_loops()

    # identify unique abego-strands
    unique_abego_strands = np.unique(sub_df['abego_strand'])
    n_designs = len(sub_df)
    iterations = []
    abegos = []
    strands = []
    abego_strands = []
    n_subsets = []
    n_subsets_stable = []
    success_ratios = []

    for iteration in range(n_iterations):
        #sample = resample(df, n_samples=n_designs, replace=True)
        sample = sub_df.sample(n=n_designs, replace=True)
        for abego_strand in unique_abego_strands:
            subset = sample.query('abego_strand==@abego_strand')

            if len(subset) >= 1:
                subset_stable = subset.query('stabilityscore > 1')
                ss_ratio = len(subset_stable)/len(subset)

                iterations.append(iteration+1)
                abego_strands.append(abego_strand)
                success_ratios.append(ss_ratio)
                n_subsets.append(len(subset))
                n_subsets_stable.append(len(subset_stable))
                print(iteration+1, abego_strand, ss_ratio)
    bootstrap_df = pd.DataFrame()
    bootstrap_df['iteration'] = iterations
    bootstrap_df['abego_strand'] = abego_strands
    bootstrap_df['success_ratio'] = success_ratios
    bootstrap_df['n_subset'] = n_subsets
    bootstrap_df['n_ss_subset'] = n_subsets_stable
    return bootstrap_df
btsp = success_rate_btsp(n_iterations = 1000)



def plot_barplot_successratio(n_iterations):

    conf_int = 0.95
    margin = (1-conf_int)/2*n_iterations
    upper = (n_iterations-margin)
    lower = margin

    conf_int_df = pd.DataFrame()
    for abego_strand in np.unique(btsp['abego_strand']):
        sub_df = btsp.query('abego_strand == @abego_strand').sort_values(by = 'success_ratio').reset_index().loc[lower:upper+1]
        conf_int_df = pd.concat([conf_int_df, sub_df])
    #return conf_int_df

    # order that corresponds to the boxplot
    orders = ['GBB-EA-BBGB-5', 'AGBB-EA-AB-5', 'GBB-GG-BAAB-5',
             'GBB-GG-BBGB-5', 'GBB-EA-AB-5',  'GBB-GG-AB-5',
            'AGBB-EA-AB-3', 'GBB-GG-BAAB-3', 'GBB-GG-BBGB-3',
             'GBB-EA-AB-3', 'AGBB-GG-AB-3','GBB-GG-AB-3' ]
    plt.figure(figsize=(10, 9))
    plt.grid(axis = 'x', zorder=-1)

    for order, i in zip(orders, range(len(orders))):
        ratio_min = min(conf_int_df.query('abego_strand==@order')['success_ratio'])
        ratio_max = max(conf_int_df.query('abego_strand==@order')['success_ratio'])
        plt.plot( [ratio_min, ratio_max],[i, i], color='black', zorder=10)

    sns.barplot(y='abego_strand',x='success_ratio', data=btsp, order=orders, color='#5AA095', zorder=5, ci=None)

    #plt.xticks(size=18)
    #plt.yticks(size=18)
    #plt.savefig('../figures/rd6_abego_strand_bootstrap_ssrate_barplot.svg')
plot_barplot_successratio(n_iterations=1000)




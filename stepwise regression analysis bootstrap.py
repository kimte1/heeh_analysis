
''' This script generates a 10-feature linear regression model '''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression

from sklearn import metrics
from sklearn.metrics import r2_score
from matplotlib import cm
from matplotlib.colors import DivergingNorm


def load_data():
    ''' returns dataframe with combined metrics and stability scores, and a list of quantitative features for analysis '''
    metrics = pd.read_csv('datasets/rd5_metrics.csv')
    scores = pd.read_csv('.datasets/rd5_stability_scores.csv')
    df = pd.merge(left = metrics, right = scores, on = 'name', how = 'inner')
    features = [ x for x in metrics.columns if ((metrics[x].dtype in [np.dtype('int64'), np.dtype('float64')]) and (False not in list(metrics[x].values == metrics[x].values))) ]
    # remove features that are difficult to generalize or is an artifcat
    features_subset = [feat for feat in features if 'site' not in feat if 'mer' not in feat if 'nearest' not in feat if 'tryp' not in feat if 'score' not in feat if 'intra' not in feat if 'frags' not in feat if 'helices_freq_R' not in feat if 'helices_freq_K' not in feat ]
    return df, features_subset

df, quant_feat = load_data()


''' +++++++++ identify top features by stepwise regression model +++++++++'''
def lr_find_feat(df, *args):
    ''' forward selection linear regression to find the next best feature '''
    feat_list = [ feat for feat in args[0] ]
    y = df['stabilityscore']
    data = []
    for feature in quant_feat:
        feat_list.append(feature)
        X = tx_df[feat_list]
        lr = LinearRegression()
        lr.fit(X, y)
        pred = lr.predict(X)
        coef = lr.coef_
        r = np.corrcoef(pred, y)[0][1]
        abs_coef = np.abs(coef)
        rmse = np.sqrt(metrics.mean_squared_error(pred, y))
        to_add = feature, r, coef[-1], abs_coef[-1], rmse
        data.append(to_add)
        del feat_list[-1]
    top_feat = pd.DataFrame(data, columns = ['feature', 'r', 'coef', 'abs_coef', 'rmse'])
    return top_feat.sort_values('rmse', ascending = True).head(40)



''' +++++++++ select 10 features +++++++++'''
def modify_features(df):
    ''' returns a datframe with re-calculated features '''
    df['freq_FILMWY'] = df['freq_F'] + df['freq_I'] + df['freq_L'] + df['freq_M'] + df['freq_W'] + df['freq_Y']
    df['n_FILMWY'] = df['freq_FILMWY'] * 43
    df['ER_contacts'] = df['n_ER_3d_contacts_3A'] / 2
    df['EE_contacts'] = df['n_EE_3d_contacts_3A'] / 2
    df['netcharge_squared_2'] = df['netcharge_squared'] / 12
    df['abego_res_profile_std'] = df['abego_res_profile'] / np.std(df['abego_res_profile'])

    feat_list = ['n_FILMWY',
                'abego_res_profile_std',
               'ST_helixcaps',
               'hphob_pos1_pos2_pos42_pos43',
               'surface_freq_hphob',
                'hphob_sc_contacts',
                'ER_contacts',
                'Tend_netq',
                'netcharge_squared_2',
                'EE_contacts',
                'buns_sc_heavy', ]
    return df, feat_list
df, feat_list = modify_features(df)


''' +++++++++ regression model by boostrapping +++++++++'''
def fwd_select_model(n_iterations, feat_list, transform = True):

    n_designs=len(df)
    summary_df = pd.DataFrame()
    pred_df = pd.DataFrame()

    for iteration in range(n_iterations):
        print(iteration)
        sample = df.sample(n=n_designs, replace=True)
        if transform:
            ss = StandardScaler()
            X = ss.fit_transform(sample[feat_list])
        else:
            X = sample[feat_list]

        y = sample['stabilityscore']
        lr = LinearRegression()
        lr.fit(X, y)
        # add pred_ss to a sub df
        pred_col = 'pred_ss_'+str(iteration)
        pred_df[pred_col] = lr.predict(X)
        # add regression results in sub df
        sub_df = pd.DataFrame()
        sub_df['iteration'] = [iteration]*len(feat_list)
        sub_df['feature'] = feat_list
        sub_df['coef'] = lr.coef_
        sub_df['rmse'] = np.sqrt(metrics.mean_squared_error(lr.predict(X), y))
        summary_df = pd.concat([summary_df, sub_df])

    pred_df['avg_pred_ss'] = pred_df.mean(axis=1)

    return summary_df, pred_df

summary_tx_true, pred_tx_true  = fwd_select_model(n_iterations=1000, feat_list=feat_list, transform=True)
summary_tx_false, pred_tx_false  = fwd_select_model(n_iterations=1000, feat_list=feat_list, transform=False)




''' test how much the model's correlation coeff improves when adding Rosetta score terms'''

score_terms = ['fa_atr', 'fa_rep', 'fa_elec', 'fa_sol',
            'lk_ball', 'lk_ball_iso', 'lk_ball_bridge', 'lk_ball_bridge_uncpl',
            'rama_prepro', 'p_aa_pp', 'omega', 'pro_close',
            'fa_dun_rot', 'fa_dun_dev', 'fa_dun_semi',
            'hbond_sr_bb', 'hbond_lr_bb', 'hbond_bb_sc', 'hbond_sc',
            'fa_intra_atr_xover4', 'fa_intra_rep_xover4', 'fa_intra_elec', 'fa_intra_sol_xover4',
            'hxl_tors', 'ref']


feat_list_2 = feat_list + score_terms
feat_list_2_no_proclose = feat_list2.remove('pro_close') # remove 1 score term


summary_tx_false_scoreterms, pred_tx_false_scoreterms  = fwd_select_model(n_iterations=1000, feat_list=feat_list_2, transform=False)
summary_tx_true_scoreterms, pred_tx_true_scoreterms  = fwd_select_model(n_iterations=1000, feat_list=feat_list_2_no_proclose, transform=True)



def regplot(df, pred):
    ''' regplot exp_ss vs. pred_ss '''
    sns.regplot(x=pred['pred_ss_0'], y=df['stabilityscore'], scatter_kws = {'alpha': 0.1}),
    plt.ylabel('Experimental Stability Score', fontsize = 15),
    plt.xlabel('Predicted Stability Score', fontsize = 15),
    print(np.corrcoef(pred['avg_pred_ss'], df['stabilityscore'])[0][1])
regplot(df=df, pred=pred_tx_true_scoreterms)


def barplot_coefficients(n_iterations, summary_df, orders, feat_list):
    ''' create barplot showing the top topological features '''
    # keep 95% conf int
    conf_int = 0.95
    margin = (1-conf_int)/2*n_iterations
    upper = (n_iterations-margin)
    lower = margin
    conf_int_df = pd.DataFrame()
    for feat in feat_list:
        sub = summary_df.query('feature == @feat').sort_values(by = 'coef').reset_index().loc[lower:upper+1]
        conf_int_df = pd.concat([conf_int_df, sub])
    # identify error bars for each filter
    plt.grid(axis = 'x', zorder=-1)
    for feat, i in zip(orders, range(len(orders))):
        feat_min = min(conf_int_df.query('feature==@feat')['coef'])
        feat_max = max(conf_int_df.query('feature==@feat')['coef'])
        plt.plot([feat_min, feat_max], [i, i], color='black',  zorder=10)
    sns.barplot(x='coef', y='feature', data=summary_df, order=orders, saturation = 0.75, edgecolor='black', ci=None, zorder = 2)


barplot_coefficients(n_iterations=1000, summary_df=summary_tx_false_scoreterms, orders=feat_list_2_no_proclose, feat_list=feat_list_2_no_proclose)
barplot_coefficients(n_iterations=1000, summary_df=summary_tx_false_scoreterms, orders=feat_list_2, feat_list=feat_list_2)
barplot_coefficients(n_iterations=1000, summary_df=summary_tx_true_scoreterms, orders=feat_list_2, feat_list=feat_list_2)

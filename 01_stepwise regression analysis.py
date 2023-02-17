
'''
perform stepwise linear regression 

'''


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression, LassoCV
from sklearn import metrics


# load data
metrics = pd.read_csv('datasets/rd5_metrics.csv')
scores = pd.read_csv('datasets/rd5_stability_scores.csv')

df = pd.merge(left = metrics, right = scores, on = 'name', how = 'inner')




def quant_feat(comp_data):
    ''' returns a list of quantitative features for analysis '''

    # list of all features
    features = [ x for x in metrics.columns if ((metrics[x].dtype in [np.dtype('int64'), np.dtype('float64')]) and (False not in list(metrics[x].values == metrics[x].values))) ]
    # remove features that are difficult to generalize
    features_subset = [feat for feat in features if 'site' not in feat if '3A' not in feat if '5A' not in feat if 'energies' not in feat if 'mer' not in feat if 'nearest' not in feat if 'tryp' not in feat if 'score' not in feat if 'intra' not in feat if 'frags' not in feat ]
    #remove individual residues since we're interested in general features
    indiv_res = [ 'G', 'A', 'V', 'L', 'I', 'P', 'M', 'W', 'F', 'Q', 'N', 'C', 'Y', 'S', 'T', 'H', 'R', 'K', 'D', 'E' ]
    indiv_freq_res = [ 'freq_%s' %res for res in indiv_res]
    all_indiv_freq_res = [feat for feat in features_subset for res in indiv_freq_res if res in feat]
    features_subset2 = list(set(all_indiv_freq_res) ^ set(features_subset))
    return features_subset2


quant_feat = quant_feat('metrics.csv')


''' Forward Selection Linear Regression '''
 ''' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ '''
def transform_data(df):
    # standardize data
    ss = StandardScaler()
    X = ss.fit_transform(df[quant_feat])
    tx_df = pd.DataFrame(X, columns = quant_feat)
    return tx_df


def lr_find_feat(df, *args):
    ''' stepwise linear regression to find the next best feature '''
    feat_list = [ feat for feat in args[0] ]
    y = df['stabilityscore']
    data = []
    for feature in quant_feat:
        lr = LinearRegression()
        feat_list.append(feature)
        X = tx_df[feat_list]
        lr.fit(X, y)
        pred = lr.predict(X)
        r = np.corrcoef(pred, y)[0][1]
        coef = lr.coef_
        abs_coef = np.abs(coef)
        rmse = np.sqrt(metrics.mean_squared_error(pred, y))
        to_add = feature, r, coef[-1], abs_coef[-1], rmse
        data.append(to_add)
        del feat_list[-1]
    top_feat = pd.DataFrame(data, columns = ['feature', 'r', 'coef', 'abs_coef', 'rmse'])
    return top_feat.sort_values('rmse', ascending = True).head(10)



def lr_stepwise(df, *args):
    ''' stepwise linear regression using top features'''
    feat_list = [ feat for feat in args[0] ]
    lr = LinearRegression()
    X = tx_df[feat_list]
    y = df['stabilityscore']
    lr.fit(X, y)
    pred = lr.predict(X)
    return feat_list, pred, y, lr.coef_


def barplot_stepwise(df, *args):
    ''' create barplot showing the top topological features '''
    feat_list, pred, y, coef = lr_stepwise(df, *args)
    feat_df = pd.DataFrame()
    feat_df = feat_df.assign(feature = feat_list, coef = coef, abs_coef = np.abs(coef))
    plt.figure(figsize = (6, 5))
    sns.barplot(x = 'coef', y = 'feature', data = feat_df.head(30)), plt.xlabel('Coefficient', fontsize = 15), plt.ylabel('Features', fontsize = 15)



def regplot(df, *args):
    ''' regplot exp_ss vs. pred_ss '''
    feat_list, pred, y, coef = lr_stepwise(df, *args)
    reg = sns.regplot(x = pred, y = y, scatter_kws = {'alpha': 0.1}), plt.ylabel('Observed Stability Score', fontsize = 30), plt.xlabel('Predicted Stability Score', fontsize = 30),
    plt.xticks(fontsize = 20), plt.yticks(fontsize = 20)
    plt.savefig('stepwise_regplot.svg')
    np.corrcoef(pred, y)[0][1]

    return reg


#########






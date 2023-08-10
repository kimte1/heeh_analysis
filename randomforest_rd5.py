import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

def import_merge_data(comp_data, scores_data):
    comp = pd.read_csv('../datasets/' + comp_data)
    scores = pd.read_csv('../datasets/' + scores_data)
    # merge all rd5 data
    df = pd.merge(left = comp, right = scores, left_on = 'name', right_on = 'name', how = 'inner')
    df = df.rename(columns = {'sequence_x': 'sequence'} )
    return(df)


def quant_feat_adj(metrics_data, quant_feat_list):
    comp = pd.read_csv('../datasets/' + metrics_data)
    with open('../datasets/' + quant_feat_list, 'r') as file:
        quant_feat = [line.replace('\n', '') for line in file]
    quant_feat = [ x for x in quant_feat if ((comp[x].dtype in [np.dtype('int64'), np.dtype('float64')]) and (False not in list(comp[x].values == comp[x].values))) ]
    return quant_feat

rd5 = import_merge_data('TK_HEEH_rd5_metrics_updated.csv', 'TK_HEEH_rd5_stabilityscore2.csv' )
quant_feat3 = quant_feat_adj('TK_HEEH_rd5_metrics_updated.csv', 'quant_feat_list3' )


quant_feat4 = [ feature for feature in quant_feat3 if 'basic' not in feature and 'charged' not in feature]


X = rd5[quant_feat4]
y = rd5['stabilityscore']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=42)


# Scale the data to be between -1 and 1
scaler = StandardScaler()
scaler.fit(X_train)
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)
# Establish model
model = RandomForestRegressor(n_jobs=-1)
# Try different numbers of n_estimators - this will take a minute or so
estimators = np.arange(10, 200, 10)
scores = []
for n in estimators:
    model.set_params(n_estimators=n)
    model.fit(X_train, y_train)
    scores.append(model.score(X_test, y_test))
plt.title("Effect of n_estimators"), plt.xlabel("n_estimator"), plt.ylabel("score"), plt.plot(estimators, scores)

# get importance
importance = model.feature_importances_
importance
# summarize feature importance
for i,v in enumerate(importance):
	print('Feature: %0d, Score: %.5f' % (i,v))
# plot feature importance
importance_df = pd.DataFrame()
importance_df['score'] = importance
importance_df['feature'] = quant_feat4
importance_df = importance_df.sort_values(by = 'score', ascending = False)

sns.barplot(x = 'feature', y = 'score', data = importance_df.head(15)), plt.xticks(rotation = 90)





rd5.values[2]

''' ECLAT '''
#Importing Libraries

transactions = []
for i in range(0, len(rd5)):
    transactions.append([str(rd5.values[i,j]) for j in range(0, 20)])
transactions[1]

from apyori import apriori

rules = apriori(transactions, min_support = 0.003, min_confidence = 0.2, min_lift = 3, min_length = 2)
result = list(rules)
result










import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
''' rd6 '''
rd6 = pd.read_csv('../datasets/rd6_kt_metrics_score_fixed_200427_updated.csv')

sns.distplot(rd6.query('stabilityscore < 1.0')['stabilityscore'], label = 'Unstable'), sns.distplot(rd6.query('stabilityscore >= 1.0')['stabilityscore'], label = 'Stable'), plt.xlabel('Stability Score', fontsize = 20), plt.legend(fontsize = 15)
len(rd6.query('stabilityscore >= 1.0')['stabilityscore'])
rd5 = pd.read_csv('../datasets/TK_HEEH_rd5_metrics_updated.csv')
rd5.head(20)
rd5.query('name == "HEEH_TK_rd5_3711.pdb" ')['hydrophobicity']

exp = pd.read_csv('../datasets/TK_HEEH_rd5_stabilityscore2.csv')
exp.query('name == "HEEH_TK_rd5_3711.pdb" ')['stabilityscore']


hydro = [ 1174, 968, 1258, 1018, 647, 842 ]
score = [2.277102, 2.361214, 2.664394, 2.451481, 1.852318, 1.600967]

sns.scatterplot(rd5['hydrophobicity'], rd5['stabilityscore'], alpha= 0.2), plt.xlabel('Hydrophobicity', fontsize = 20), plt.ylabel('Stability Score', fontsize = 20), plt.scatter(hydro, score, )

sns.boxplot(rd5['n_hydrophobic'], rd5['stabilityscore']), plt.xlabel('Number of Hydrophobic Residues', fontsize = 20), plt.ylabel('Stability Score', fontsize = 20)

rd5['stabilityscore_bool'] = [ 1 if score >= 1.0 else 0 for score in rd5['stabilityscore'] ]
sns.scatterplot(rd5['percent_hydrophobic'], rd5['netcharge_squared'], hue = rd5['stabilityscore_bool'], palette = 'colorblind', sizes = (5, 100))
sns.set_palette("colorblind")
sns.scatterplot(rd5['percent_hydrophobic'], rd5['netcharge_squared'], hue = rd5['stabilityscore'],palette = 'rainbow', alpha=0.2,), plt.xlabel('Percent Hydrophobic', fontsize = 15), plt.ylabel('Netcharge Squared', fontsize = 15)

sns.scatterplot(rd5['percent_hydrophobic'], rd5['score_per_res'], hue = rd5['stabilityscore'], palette = 'rainbow',)



sns.scatterplot(rd5['hydrophobicity'],rd5['netcharge_squared'], hue = rd5['stabilityscore'], palette = 'rainbow' ), plt.xlabel('Hydrophobicity', fontsize = 15), plt.ylabel('Netcharge Squared', fontsize = 15)






sns.regplot(rd5['pred'])


reg = pd.read_csv('../datasets/TK_HEEH_rd5_reg.csv')
rd5 = pd.merge(left = rd5, right = reg,  left_on = 'name', right_on = 'name', how = 'inner')

len(rd5)

sns.regplot(rd5['pred_ss_4'], rd5['stabilityscore'], scatter_kws = {'alpha':0.1}), plt.xlabel('Predicted Stability', fontsize = 15), plt.ylabel('Observed Stability', fontsize = 15)


np.corrcoef(rd5['pred_ss_5'], rd5['stabilityscore'])


rd5['hphob_n22'] = [1 if residue in 'AFILMVWY' else 0 for sequence in rd5['sequence'] for residue in sequence ]

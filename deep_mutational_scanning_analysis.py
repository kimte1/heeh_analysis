''' this script analyzes how different residue type mutants impact the average stability score '''


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#metrics = pd.read_csv('../datasets/TK_HEEH_rd6_metrics.csv')


df = pd.read_csv('../datasets/HEEH_KT_rd6_lowenergy_exp_210331_2.csv')


# add category
def cat(x):
    for i in wts:
        if 'pdb_' in x:
            return 'mutant'



def add_columns():
    wts=['HEEH_TK_rd5_0018','HEEH_TK_rd5_0341','HEEH_TK_rd5_0420','HEEH_TK_rd5_0614','HEEH_TK_rd5_0958','HEEH_TK_rd5_3711']

    df['cat'] = [cat(x) for x in df['name'].values ]

    mutants_df = df.query('cat == "mutant"')

    wt_residues = []
    positions = []
    for name in mutants_df['name'].values:
        wt_residue = name[-4]
        if wt_residue == '_':
            wt_residue = name[-3]
        try:
            int_name = int(name[22:24])
        except:
            int_name = int(name[22:23])
        positions.append(int_name)
        wt_residues.append(wt_residue)
    mutants_df['pos'] = positions
    mutants_df['wt'] = wt_residues
    mutants_df['mutant'] = [ name[-1] for name in mutants_df['name'].values]


    helix_core = [ 5, 9, 12, 13, 32, 34, 35, 36, 39]
    helix_solv = [3, 4, 6, 7, 8, 10, 11, 14, 15, 30, 31, 33, 37, 38, 40, 41]
    loop1_3 = [ 16, 17, 18, 27, 28, 29] #, 22, 23, 27, 28, 29]
    loop2 = [22, 23]

    strand_core = [19, 21, 24, 26]
    strand_solv = [20, 25 ]
    ends = [1, 2, 42, 43]
    topology = []
    for pos in mutants_df['pos']:

        if pos in helix_core:
            topology.append('helix_core')
        if pos in helix_solv:
            topology.append('helix_solv')
        if pos in loop1_3:
            topology.append('loop1_3')
        if pos in loop2:
            topology.append('loop2')
        if pos in strand_core:
            topology.append('strand_core')
        if pos in strand_solv:
            topology.append('strand_solv')
        if pos in ends:
            topology.append('ends')
    mutants_df['topology'] = topology
    #print(len(topology))

    branch_hphob = 'AILMV'
    aromatic_hphob = 'FYW'
    polar = 'QNCST'
    charged = 'HRKDE'
    other = 'GP'
    wt_type = []
    mut_type = []
    for wt, mutant in zip(mutants_df['wt'].values, mutants_df['mutant'].values):
        if wt in branch_hphob:
            wt_type.append('br_hphob')
        if wt in aromatic_hphob:
            wt_type.append('ar_hphob')
        if wt in polar:
            wt_type.append('polar')
        if wt in charged:
            wt_type.append('charged')
        if wt in other:
            wt_type.append('other')

        if mutant in branch_hphob:
            mut_type.append('br_hphob')
        if mutant in aromatic_hphob:
            mut_type.append('ar_hphob')
        if mutant in polar:
            mut_type.append('polar')
        if mutant in charged:
            mut_type.append('charged')
        if mutant in other:
            mut_type.append('other')

    mutants_df['wt_type'] = wt_type
    mutants_df['mut_type'] = mut_type

    return mutants_df
mutants_df = add_columns()

len(mutants_df)
mutants_df.head()


def plot():

    all_topology = []
    all_wt_types = []
    all_mut = []
    all_scores = []
    counts=[]
    for topology in np.unique(mutants_df['topology']):
        #print(topology)
        wt_types = np.unique(mutants_df['wt_type'])
        mut_types = np.unique(mutants_df['mut_type'])
        count = 0
        for wt_type in wt_types:

            #wt = mutants_df.query('wt_type == @wt_type')
            for mut_type in mut_types:
                subdf = mutants_df.query('topology == @topology & wt_type == @wt_type & mut_type == @mut_type')
                avg_score = np.mean(subdf['stabilityscore'])
                all_topology.append(topology)
                all_wt_types.append(wt_type)
                all_mut.append(mut_type)
                all_scores.append(avg_score)
                count+=1
                counts.append(count)
                #print(wt_type, mut_type, avg_score)
    subdf2 = pd.DataFrame()
    subdf2['count'] = counts
    subdf2['topology'] = all_topology
    subdf2['wt_type'] = all_wt_types
    subdf2['mut_type'] = all_mut
    subdf2['avg_score'] = all_scores

    return subdf2
subdf = plot()

subdf



def plot2():
    n_topologies = len(np.unique(subdf['topology']))

    fig, axs = plt.subplots(n_topologies, 1, figsize=(13,50))
    for topo, i in zip(np.unique(subdf['topology']), range(n_topologies)):
        sub = subdf.query('topology == @topo')

        sns.barplot(x='count', y='avg_score', data=sub, ax=axs[i], color='#95C8E1')
        axs[i].set_title(topo, size = 20)
plot2()










==== ===== ====

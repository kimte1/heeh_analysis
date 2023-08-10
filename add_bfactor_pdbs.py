
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


df=pd.read_csv('../datasets/HEEH_KT_rd6_lowenergy_exp_210331.csv')

df['sequence'] = df['protein_sequence_c']
protein_names=['HEEH_TK_rd5_0018','HEEH_TK_rd5_0341','HEEH_TK_rd5_0420','HEEH_TK_rd5_0614','HEEH_TK_rd5_0958','HEEH_TK_rd5_3711']

def cat(x):
    for i in protein_names:
        if i in x:
            return i
    if 'PG_hp' in x:
        return 'PG_hp'
    else:
        return 'design'
df['cat'] = [cat(x) for x in df['name'].values]



aas='QENHDRKTSAGMLVIWYFP'
def makessm(c="HEEH_TK_rd5_0018.pdb",value='stabilityscore',center=1,annot=True):
    wtseq = df.query('name == "%s.pdb"' % c)['sequence'].values[0]
    wt=df.query('name == "%s.pdb"' % c)[value].values[0]

    out=np.zeros((43,19))
    out[:] = wt

    subdf=df.query('cat == "%s"' % c)
    for name, v in zip(subdf['name'], subdf[value]):
        if name.count('_') == 4:
            pos = int(name.split('_')[-1][1:-1])
            aa = name[-1]
            out[pos-1, aas.index(aa)] = v

    plt.figure(figsize=(20,7))

    if annot:
        annotation = np.where(out > 1, out, '').T
    else:
        annotation=False

    if center=='wt': center=wt

    hm = sns.heatmap(out.T, xticklabels=['%s%s' % (aa, pos) for aa, pos in zip(wtseq, range(1,44))], yticklabels=[x for x in aas],
               cmap='bwr',center=center,annot=annotation,fmt='.3s', )
    plt.plot([0,43],[5,5],color='black')
    plt.plot([0,43],[7,7],color='black')
    plt.plot([0,43],[11,11],color='black')
    plt.plot([0,43],[15,15],color='black')
    plt.plot([0,43],[18,18],color='black')



    for i in range(43):
        plt.scatter(i+0.5,aas.index(wtseq[i])+0.5, s=30,color='black')

    cax = plt.gcf().axes[-1]
    cax.plot([0,2],[wt,wt],color='black',linewidth=3)
    return hm, out


sensitivity={}
for wt in protein_names:
    hm, out = makessm(wt,annot=False,center='wt')
    sensitivity[wt] = np.average(out, axis=1)
    plt.title('%s.pdb Stability Score' % wt,fontsize=15)
    plt.tight_layout()
    plt.yticks(rotation=0)
    plt.savefig('../figures/%s.pdb_new.svg' % wt)
    plt.figure()

plt.figure(figsize=(24,7))
for name in protein_names:


    plt.plot(sensitivity[name]  - np.max(sensitivity[name]),linewidth=2,label=name)
plt.legend()
plt.grid()



for name in protein_names:
    avg_muts = sensitivity[name]  - np.max(sensitivity[name])

    plt.figure(figsize=(24,2))
    plt.title(name)
    plt.plot(avg_muts,linewidth=2)
    plt.grid()

    newlines = []
    with open('../rd5_pdbs/'+name + '.pdb') as file:
        lines = file.readlines()


    for line in lines:
        if 'ATOM' not in line:
            newlines.append(line)
        else:
            resid = int(line[24:27])
            newlines.append(line[0:61] + ('%5s' % ('%.2f' % avg_muts[resid - 1])) + line[66:])

    with open('../rd5_pdbs/'+name+'.pdb' + '_avgmut.pdb','w') as file:
        for line in newlines:
            file.write(line)

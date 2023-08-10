

### compare the

#import Bio.PDB
#from Bio.PDB import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
# Import SciPi Library
from scipy.spatial import distance
from math import dist


''' find average coordinate positions for 40 nmr pdbs '''

def load_and_clean():
    chains = ['chainA', 'chainB']
    allchains_df = pd.DataFrame()
    # load data
    for chain in chains:
        df=pd.read_csv('../data/NMR/heeh_341_'+chain+'.pdb', delim_whitespace=True)
        # get rid of NaN and re-index
        df_atoms=df.dropna().reset_index()
        # convert atom num and residue position to integers
        # convert coordinates to float
        cols_int = ['level_1', 'level_5']
        cols_float = ['level_6', 'level_7', 'level_8', 'level_9', 'MODEL']
        for col in cols_int:
            df_atoms[col] = df_atoms[col].astype(int)
        for col in cols_float:
            df_atoms[col] = df_atoms[col].astype(float)
        allchains_df = pd.concat([allchains_df, df_atoms])
    return allchains_df
    allchains_df


def average_nmr_df():
    ''' create one df that takes the average coordinates for all 20 nmr structures '''
    allchains_df = load_and_clean()

    all_new_rows = []
    highest_atom_num = allchains_df['level_1'].max()
    for atom_num in range(1, highest_atom_num+1 ): #739
        sub_df = allchains_df.query('level_1==@atom_num').reset_index()
        level_6_avg = '%.3f'%sub_df['level_6'].mean()
        level_7_avg = '%.3f'%sub_df['level_7'].mean()
        level_8_avg = '%.3f'%sub_df['level_8'].mean()
        level_9 = '%.2f'%sub_df['level_9'].mean()
        MODEL_avg = '%.2f'%sub_df['MODEL'].mean()
        last = sub_df['1'].iloc[0]

        new_row1=[]
        new_row2=[]
        combine_new_row=[]

        new_row1 = list(sub_df.iloc[0][['level_0', 'level_1', 'level_2', 'level_3', 'level_4', 'level_5']])
        new_row2= [level_6_avg, level_7_avg, level_8_avg, level_9, MODEL_avg]

        combine_new_row = new_row1 + new_row2 + [last]
        all_new_rows.append(combine_new_row)

    df=pd.DataFrame(all_new_rows, columns=['atom', 'atom_num', 'atom_type', 'residue', 'res_abbrev', 'res_pos', 'x', 'y', 'z', 'occupancy', 'b-factor', 'atom2'])

    return df



def load_design():
    heeh_341=pd.read_csv('../data/pdbs/HEEH_TK_rd5_0341.pdb',  delim_whitespace=True)
    heeh_341=heeh_341.dropna()
    heeh_341.columns =['atom', 'atom_num', 'atom_type', 'residue', 'res_abbrev', 'res_pos', 'x', 'y', 'z', 'occupancy', 'b-factor', 'atom2']
    heeh_341['atom_num'] = heeh_341['atom_num'].astype(int)
    heeh_341['res_pos'] = heeh_341['res_pos'].astype(int)
    return heeh_341

def carbon_distance(carbon):
    ''' make df with just a_carbon or b_carbon info and calculate distance '''
    nmr_df = average_nmr_df()
    design_df = load_design()

    # subset b-carbon data
    nmr_bc=nmr_df.query('atom_type==@carbon')
    design_bc=design_df.query('atom_type==@carbon')
    # renumber residue position for nmr df
    nmr_bc['res_pos']=nmr_bc['res_pos']-21


    merge = pd.merge(left=nmr_bc[['residue', 'res_pos', 'x', 'y', 'z']], right=design_bc[['residue', 'res_pos', 'x', 'y', 'z']], on='res_pos', how='inner')

    all_distances = []
    for i in range(len(merge)):
        nmr_x = float(merge.iloc[i, 2])
        nmr_y = float(merge.iloc[i, 3])
        nmr_z = float(merge.iloc[i, 4])
        nmr_coord = (nmr_x, nmr_y, nmr_z)

        design_x = merge.iloc[i, 6]
        design_y = merge.iloc[i, 7]
        design_z = merge.iloc[i, 8]
        design_coord = (design_x, design_y, design_z)

        distance = dist(nmr_coord, design_coord)
        all_distances.append(distance)

    merge['dist']=all_distances
    return merge


carbon_a_df=carbon_distance(carbon='CA')
carbon_b_df=carbon_distance(carbon='CB')

carbon_b_df.sort_values('dist')


carbon_a_df.sort_values('dist', ascending=False).head(10)

carbon_b_df.sort_values('dist', ascending=False).head(10)

def plot_distance(df, title):
    plt.figure(figsize=(10,4))
    sns.lineplot(x='res_pos', y='dist', data=df)
    plt.title('distance bw ' + title + ' design model vs nmr structure')
    #plt.savefig('../figures/figure_components/distance_'+title+'.svg')
plot_distance(carbon_a_df, title='a-carbons')

plot_distance(carbon_b_df, title='b-carbons')

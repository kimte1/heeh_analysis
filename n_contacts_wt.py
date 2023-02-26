### goal: compare change in stability (due to mutation) vs. number of contacts for wt hphob residues



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import distance
from math import dist
import itertools
import glob as glob
import os

### chdir to where the dms pdbs are 
os.chdir('datasets/pdbs/rd5_dms_lowest_models')



def get_wt_filename(pdbs):
    ''' get wt pdb file '''
    wt = []
    for file in pdbs:
        if file.split('_')[4].split('.')[0] == 'fasta': # find wt file
            wt.append(file)
    return wt[0]

def wt_bcarbon(pdb):
    ''' return df with atomic coordinates of b-carbon for wt '''
    with open(pdb, 'r') as file:
        data = file.read().splitlines()
    new_data = []
    for line in data[1:]:
        new_data.append(line.split())
    df = pd.DataFrame(new_data, columns=['atom', 'atom_num', 'atom_type', 'residue', 'chain', 'res_pos', 'x', 'y', 'z', 'occupancy', 'b-factor', 'atom2'])
    df=df.dropna()
    df['atom_num'] = df['atom_num'].astype(int)
    df['res_pos'] = df['res_pos'].astype(int)
    df['x'] = df['x'].astype(float)
    df['y'] = df['y'].astype(float)
    df['z'] = df['z'].astype(float)
    df.drop(columns=['b-factor', 'atom2', 'occupancy', 'chain'], inplace=True)
    df=df.query('atom_type=="CB"')
    return df




def res_type(df):
    ''' add column with 1-letter residue abbreviation'''
    aa_dict = {'GLY':'G',
               'ALA':'A',
              'VAL':'V',
              'LEU':'L',
              'ILE':'I',
              'PRO':'P',
              'MET':'M',
              'TRP':'W',
              'PHE':'F',
              'GLN':'Q',
              'ASN':'N',
              'TYR':'Y',
              'SER':'S',
              'THR':'T',
              'HIS':'H',
              'ARG':'R',
              'LYS':'K',
              'ASP':'D',
              'GLU':'E',}

    all_abbrev = []
    for res in list(df['residue']):
        abbrev = aa_dict[res]
        all_abbrev.append(abbrev)
    df['res'] = all_abbrev
    return df


def pairwise_dist(df, design):
    ''' returns data on the euclidean distance for all possible beta carbon positions in rosetta model) '''
    pairwise = list(itertools.permutations(df['res_pos'], 2)) # get all pairwise position combinations
    
    pairs = []
    pairwise_dist = []
    res1_type = []
    res2_type = []
    
    for pair in pairwise:
        first_pair = pair[0]
        second_pair = pair[1]
        
        res1 = df.query('res_pos==@first_pair')['res'].tolist()[0]
        res2 = df.query('res_pos==@second_pair')['res'].tolist()[0]
        # get coordinates
        first_x = float(df.query('res_pos==@first_pair')['x'])
        first_y = float(df.query('res_pos==@first_pair')['y'])
        first_z = float(df.query('res_pos==@first_pair')['z'])
        first_coord = np.array((first_x, first_y, first_z), dtype=float)
        second_x = float(df.query('res_pos==@second_pair')['x'])
        second_y = float(df.query('res_pos==@second_pair')['y'])
        second_z = float(df.query('res_pos==@second_pair')['z'])
        second_coord = np.array((second_x, second_y, second_z), dtype=float)
        # calculate distance
        distance = np.linalg.norm(first_coord-second_coord)
        # add results to lists
        pairs.append(pair)
        pairwise_dist.append(distance)
        res1_type.append(res1)
        res2_type.append(res2)
        
    newdf = pd.DataFrame(pairs, columns=['pos1', 'pos2'])
    newdf['pos_pair'] = pairs
    newdf['dist'] = pairwise_dist
    newdf['res1'] = res1_type
    newdf['res2'] = res2_type
    newdf['design'] = 'HEEH_TK_rd5_'+ design
    return newdf

def keep_dist_close_contact(df):
    ''' distance < 7A'''
    return df.query('dist < 8')


def n_contacts(df):
    ''' quantify number of contacts and add it to df '''
    newdf=pd.DataFrame()
    for design in ['HEEH_TK_rd5_0018', 'HEEH_TK_rd5_0341', 'HEEH_TK_rd5_0420', 'HEEH_TK_rd5_0614',
                   'HEEH_TK_rd5_0958', 'HEEH_TK_rd5_3711']: 
        
        subdf=df.query('design==@design')
        contact_counts = subdf['pos1']
        df_contacts = contact_counts.value_counts().rename_axis('pos1').reset_index(name='n_contacts') # turn contacts freq into df
        
        subdf=pd.merge(left=subdf, right=df_contacts, on='pos1', how='inner')
        
        newdf = pd.concat([newdf, subdf])
    return newdf 


def final(design):
    list_pdb_filenames = subset_pdbs(design)
    wt_file = get_wt_filename(pdbs=list_pdb_filenames)
    bcarbon_coord = wt_bcarbon(pdb=wt_file)
    bcarbon_coord = res_type(df=bcarbon_coord)
    pairwise = pairwise_dist(df=bcarbon_coord, design=design)
    close_contacts_df = keep_dist_close_contact(df=pairwise)
    n_contacts_df = n_contacts(df=close_contacts_df)
    return n_contacts_df
    

pairwise_df=pd.DataFrame()

for design in ['0018', '0341', '0420', '0614', '0958', '3711']:
    subdf = final(design=design)
    pairwise_df = pd.concat([pairwise_df, subdf])





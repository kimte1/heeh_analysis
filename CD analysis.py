''' this script analyzes CD data '''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import savgol_filter
from scipy import interpolate



def import_data(metrics_data, scores_data):
    ''' returns df with designs' metrics and stabilityscore data '''

    metrics = pd.read_csv('../datasets/' + metrics_data)
    scores = pd.read_csv('../datasets/' + scores_data)
    df = pd.merge(left = metrics, right = scores, left_on = 'name', right_on = 'name', how = 'inner')
    df = df.rename(columns = {'sequence_x': 'sequence'} )
    return(df)
df_5 = import_data('TK_HEEH_rd5_addtl_metrics.csv', 'TK_HEEH_rd5_stabilityscore.csv')



def scattterplot_hphob():
    ''' returns a scatterplot of score vs hydrophobicity, highlighting 6 designs for CD analysis '''

    plt.figure(dpi = 300)
    sns.scatterplot(x = 'hydrophobicity', y = 'stabilityscore', data = df_5, color = '#88CCEE', alpha = 0.2, linewidth = 0)
    plt.xlabel('Hydrophobicity', fontsize = 30), plt.xticks(fontsize = 20)
    plt.ylabel('Stability Score', fontsize = 30), plt.yticks( fontsize = 20),

    select_hphob_designs = ['0018.pdb', '0341.pdb', '0420.pdb', '0614.pdb', '0958.pdb', '3711.pdb']
    for file in select_hphob_designs:
        design = 'HEEH_TK_rd5_' + file
        hphob = float(df_5.query('name == @design')['hydrophobicity'])
        score = float(df_5.query('name == @design')['stabilityscore'])
        plt.scatter(x = hphob, y = score, color = '#882255')
#    plt.savefig('Fig3.tiff')
scattterplot_hphob()


#def load_cd_melt_data():
    # load df with CD data
    cd_df = pd.read_csv('../data_CD/CD_data_analysis.csv')
    melt_df = pd.read_csv('../data_CD/melt_data_analysis.csv')
#    return cd_df, melt_df
#cd_df, melt_df = load_cd_melt_data()

def add_molar_ellipticy(cd_df):
    ''' returns df with CD data that substracts background noise and converts to residue ellipticity '''

    # pick columns with data
    #col = df.columns[1:19]
    # remove pbs background for each wavescan
    #for c in col:
    #    name = c + '-pbs'
    #    df[name] = df[c] - df['pbs']
    # design names, molecular weights, concentrations, and the 3 temperature ranges in which scans were done
    design = ['18', '341', '420', '614', '958', '3711' ]
    mw = [ 7440, 7464, 7472, 7447, 7430, 7554]
    conc = [ 0.32, 0.326, 0.136, 0.112, 0.0975, 0.174 ]
    wavescan = ['CD_25', 'CD_95', 'CD_25return']
    # create new colums that converts wavescan data to molar ellipticity (me)
    for i in range(len(design)):
        for scan in wavescan:
            col_name = design[i] + '_' + scan
            new_col_name = design[i] + '_' + scan + '_me'
            # conversion = scan * mw / n_residues / (10 * len_path_cell * concentration) / factor of 3
            # units = mdeg * g/mol / 64 / (10 * cm * g/L) / 1000
            me_conversion = mw[i] / 64 / (10 * 0.1 * conc[i]) / 1000
            cd_df[new_col_name] = (cd_df[col_name]-cd_df['pbs']) * me_conversion

add_molar_ellipticy(cd_df)


def plot_CD(df, design):
    ''' return CD plot for a specific design '''

    line_25 = design + '_CD_25_me'
    line_95 = design + '_CD_95_me'
    line_25_return = design + '_CD_25return_me'

    sns.lineplot(df['wavelength'], df[line_25], color = 'black', label = '25C'),
    sns.lineplot(df['wavelength'], df[line_95], color = '#D81B60', label = '95C'),
    sns.lineplot(df['wavelength'], df[line_25_return], color = '#1E88E5', label = '25C return'),
    plt.xlim(200)
    plt.ylim(-13,5)
    #plt.savefig('cdplot_' + design + '.svg')

plot_CD(cd_df, '18')
plot_CD(cd_df, '341')
plot_CD(cd_df, '420')
plot_CD(cd_df, '614')
plot_CD(cd_df, '958')
plot_CD(cd_df, '3711')





def plot_thermal_melt(melt_df):
    ''' returns df with thermal melt data that converts to residue ellipticity '''

    # design names, molecular weights, and concentrations
    design = ['18', '341', '420', '614', '958', '3711']
    mw = [ 7440, 7464, 7472, 7447, 7430, 7554]
    conc = [ 0.32, 0.326, 0.136, 0.112, 0.0975, 0.133 ]
    # create new colums that converts wavescan data to molar ellipticity (me)
    for i in range(len(design)):
        col_name = design[i] + '_melt'
        new_col_name = design[i] + '_melt_me'
        # conversion = scan * mw / n_residues / (10 * len_path_cell * concentration)
        # units = mdeg * g/mol / 64 / (10 * cm * g/L)
        me_conversion = mw[i] / 64 / (10 * 0.1 * conc[i]) / 1000
        melt_df[new_col_name] = melt_df[col_name] * me_conversion
plot_thermal_melt(melt_df)



def plot_melt(df, design):
    ''' return melt plot for a specific design '''

    col_name = design + '_melt_me'
    df['smooth_'+col_name] = savgol_filter(df[col_name], window_length=5,  polyorder=3)
    sns.lineplot(df['temp'], df['smooth_'+col_name], color = '#D81B60'),
    #sns.lineplot(df['temp'], df[col_name], color = '#D81B60'),
    plt.ylim(df[col_name].min()-1, 0)

    #plt.savefig('../figures/meltplot_' + design + '.svg')

plot_melt(melt_df, '18')
plot_melt(melt_df, '341')
plot_melt(melt_df, '420')
plot_melt(melt_df, '614')
plot_melt(melt_df, '958')
plot_melt(melt_df, '3711')





#def plot_rate_change(df, design):
#    designs = melt_df.columns[13:]

#    col_name = design + '_melt_me'
#    sns.regplot(df['temp'], df[col_name].diff(), order = 2, color = '#D81B60'),
    #plt.ylim(df[col_name].diff().min()-1, 0)
    #plt.savefig('meltplot_' + design + '.svg')
#plot_rate_change(melt_df, '18')



plot_rate_change(melt_df, '341')
plot_rate_change(melt_df, '420')
plot_rate_change(melt_df, '614')
plot_rate_change(melt_df, '958')
plot_rate_change(melt_df, '3711')

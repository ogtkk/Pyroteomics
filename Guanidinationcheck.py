import Functions as func
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotter as plot
import os
import re
from datetime import datetime as dt

filenames = ['160306ko-Guan.txt']
data = [func.ReadFile(filename) for filename in filenames]
data = func.MergeFilesNoRemove(data)

time = dt.now()
f_time = time.strftime('%Y%m%d%H%M%S')
f_name = f_time + re.compile(u".txt").sub("", 'Guanidination') + "_result"
if not os.path.isdir(f_name):
    os.mkdir(f_name)
os.chdir(f_name)
homedir = os.getcwd()
data = func.getSeqProperties(data)
data = func.GuanidinationEfficiency(data)
data = func.extractPhospho(data)
dfs_temp = func.DivideMergeFile(data)
plot.chargeStateID_barplot(dfs_temp)
dfs_temp = {k: func.RemoveDupli(v) for k, v in dfs_temp.items()}
dfs = dfs_temp

plot.Guanidination_efficiency_barplot(dfs)
dic = {'Dataname': [],
       'Fullmeanarea': [],
       'Partialmeanarea': [],
       'Nonemeanarea': [],
       'N-termmeanarea': [],
       'ID': [],
       'FullID': [],
       'PartialID': [],
       'NoneID': [],
       'N-termID': [],
       'Guanidyl_1': [],
       'Guanidyl_2': [],
       'Guanidyl_3': [],
       'Guanidyl_4': [],
       'Guanidyl_1_area': [],
       'Guanidyl_2_area': [],
       'Guanidyl_3_area': [],
       'Guanidyl_4_area': [],
       'phospho_1': [],
       'phospho_2': [],
       'phospho_3': [],
       'phospho_4': [],
       'phospho_>4': [],
       'phospho_area_1': [],
       'phospho_area_2': [],
       'phospho_area_3': [],
       'phospho_area_4': [],
       'phospho_area_>4': [],
       'terminal_K': [],
       'terminal_R': [],
       'average_number_of_K': [],
       'missclvg_0': [],
       'missclvg_1': [],
       'missclvg_2': [],
       'missclvg_0_area': [],
       'missclvg_1_area': [],
       'missclvg_2_area': [],
       '..K..K': [],
       '..K..R': [],
       '..R..K': [],
       '..R..R': [],
       '..K..Cterm': [],
       '..R..Cterm': [],
       '..K..': [],
       '..R..': [],
       '..K..K_area': [],
       '..K..R_area': [],
       '..R..K_area': [],
       '..R..R_area': [],
       '..K..Cterm_area': [],
       '..R..Cterm_area': [],
       '..K.._score': [],
       '..R.._score': [],
       'K_0': [],
       'K_1': [],
       'K_2': [],
       'K_3': [],
       'K=1-terminalK': [],
       'K=1-terminalR': []
       }
i = 0
for k, v in dfs.items():
    area = np.sum(v[u'PepArea'])
    farea_list = v[v[u'Guanidination_efficiency']
        .str.contains('Full')][u'PepArea'].values
    parea_list = v[v[u'Guanidination_efficiency']
        .str.contains('Partial')][u'PepArea'].values
    narea_list = v[v[u'Guanidination_efficiency']
        .str.contains('None')][u'PepArea'].values
    oarea_list = v[v[u'Guanidination_efficiency']
        .str.contains('Overlabel')][u'PepArea'].values
    dic['Dataname'].append(k)
    dic['Fullmeanarea'].append(np.nanmean(farea_list))
    dic['Partialmeanarea'].append(np.nanmean(parea_list))
    dic['Nonemeanarea'].append(np.nanmean(narea_list))
    dic['N-termmeanarea'].append(np.nanmean(oarea_list))
    dic['ID'].append(len(v))
    dic['FullID'].append(len(v[v[u'Guanidination_efficiency'].str.contains('Full')]))
    dic['PartialID'].append(len(v[v[u'Guanidination_efficiency'].str.contains('Partial')]))
    dic['NoneID'].append(len(v[v[u'Guanidination_efficiency'].str.contains('None')]))
    dic['N-termID'].append(len(v[v[u'Guanidination_efficiency'].str.contains('Overlabel')]))
    dic['Guanidyl_1'].append(len(v[v[u'number_of_Guanidination'] == 1]))
    dic['Guanidyl_2'].append(len(v[v[u'number_of_Guanidination'] == 2]))
    dic['Guanidyl_3'].append(len(v[v[u'number_of_Guanidination'] == 3]))
    dic['Guanidyl_4'].append(len(v[v[u'number_of_Guanidination'] == 4]))
    dic['Guanidyl_1_area'].append(np.sum(v[v[u'number_of_Guanidination'] == 1][u'PepArea'].values))
    dic['Guanidyl_2_area'].append(np.sum(v[v[u'number_of_Guanidination'] == 2][u'PepArea'].values))
    dic['Guanidyl_3_area'].append(np.sum(v[v[u'number_of_Guanidination'] == 3][u'PepArea'].values))
    dic['Guanidyl_4_area'].append(np.sum(v[v[u'number_of_Guanidination'] == 4][u'PepArea'].values))
    dic['phospho_1'].append(len(v[v[u'number_of_phospho'] == 1]))
    dic['phospho_2'].append(len(v[v[u'number_of_phospho'] == 2]))
    dic['phospho_3'].append(len(v[v[u'number_of_phospho'] == 3]))
    dic['phospho_4'].append(len(v[v[u'number_of_phospho'] == 4]))
    dic['phospho_>4'].append(len(v[v[u'number_of_phospho'] > 4]))
    dic['phospho_area_1'].append(np.sum(v[v[u'number_of_phospho'] == 1][u'PepArea'].values))
    dic['phospho_area_2'].append(np.sum(v[v[u'number_of_phospho'] == 2][u'PepArea'].values))
    dic['phospho_area_3'].append(np.sum(v[v[u'number_of_phospho'] == 3][u'PepArea'].values))
    dic['phospho_area_4'].append(np.sum(v[v[u'number_of_phospho'] == 4][u'PepArea'].values))
    dic['phospho_area_>4'].append(np.sum(v[v[u'number_of_phospho'] > 4][u'PepArea'].values))
    dic['terminal_K'].append(len(v[v[u'Seq'].str.endswith('K')]))
    dic['terminal_R'].append(len(v[v[u'Seq'].str.endswith('R')]))
    dic['average_number_of_K'].append(np.nanmean(v[u'number_of_k'].values))
    dic['missclvg_0'].append(len(v[v[u'MissClvg'] == 0]))
    dic['missclvg_1'].append(len(v[v[u'MissClvg'] == 1]))
    dic['missclvg_2'].append(len(v[v[u'MissClvg'] == 2]))
    dic['missclvg_0_area'].append(v[v[u'MissClvg'] == 0][u'PepArea'].sum())
    dic['missclvg_1_area'].append(v[v[u'MissClvg'] == 1][u'PepArea'].sum())
    dic['missclvg_2_area'].append(v[v[u'MissClvg'] == 2][u'PepArea'].sum())
    dic['K_0'].append(len(v[v[u'number_of_k'] == 0]))
    dic['K_1'].append(len(v[v[u'number_of_k'] == 1]))
    dic['K_2'].append(len(v[v[u'number_of_k'] == 2]))
    dic['K_3'].append(len(v[v[u'number_of_k'] == 3]))
    vk = v[v['number_of_k'] == 1]
    dic['K=1-terminalK'].append(len(vk[vk[u'Seq'].str.endswith('K')]))
    dic['K=1-terminalR'].append(len(vk[vk[u'Seq'].str.endswith('R')]))
    dic['..K..'].append(len(v[v[u'Seq'].str.contains('.*K.+')]))
    dic['..R..'].append(len(v[v[u'Seq'].str.contains('.*R.+')]))
    v = v[v[u'MissClvg'] == 1]
    dic['..K..K'].append(len(v[v[u'Seq'].str.contains('.*K.*K$')]))
    dic['..K..R'].append(len(v[v[u'Seq'].str.contains('.*K.*R$')]))
    dic['..R..K'].append(len(v[v[u'Seq'].str.contains('.*R.*K$')]))
    dic['..R..R'].append(len(v[v[u'Seq'].str.contains('.*R.*R$')]))
    dic['..K..Cterm'].append(len(v[v[u'Seq'].str.contains('K.*[^KR]$')]))
    dic['..R..Cterm'].append(len(v[v[u'Seq'].str.contains('R.*[^KR]$')]))
    dic['..K..K_area'].append(np.median(v[v[u'Seq'].str.contains('.*K.*K$')][u'PepArea']))
    dic['..K..R_area'].append(np.median(v[v[u'Seq'].str.contains('.*K.*R$')][u'PepArea']))
    dic['..R..K_area'].append(np.median(v[v[u'Seq'].str.contains('.*R.*K$')][u'PepArea']))
    dic['..R..R_area'].append(np.median(v[v[u'Seq'].str.contains('.*R.*R$')][u'PepArea']))
    dic['..K..Cterm_area'].append(np.median(v[v[u'Seq'].str.contains('K.*[^KR]$')][u'PepArea']))
    dic['..R..Cterm_area'].append(np.median(v[v[u'Seq'].str.contains('R.*[^KR]$')][u'PepArea']))
    dic['..K.._score'].append(np.median(v[v[u'Seq'].str.contains('.*K.+')][u'PeptScore'].values))
    dic['..R.._score'].append(np.median(v[v[u'Seq'].str.contains('.*R.+')][u'PeptScore'].values))
result = pd.DataFrame(dic, columns=['Dataname', 'Fullmeanarea', 'Partialmeanarea',
                                    'Nonemeanarea', 'N-termmeanarea',
                                    'ID', 'FullID', 'PartialID', 'NoneID',
                                    'N-termID', 'Guanidyl_1', 'Guanidyl_2',
                                    'Guanidyl_3', 'Guanidyl_4', 'Guanidyl_1_area',
                                    'Guanidyl_2_area', 'Guanidyl_3_area', 'Guanidyl_4_area',
                                    'phospho_1', 'phospho_2', 'phospho_3', 'phospho_4',
                                    'phospho_>4', 'phospho_area_1', 'phospho_area_2',
                                    'phospho_area_3', 'phospho_area_4', 'phospho_area_>4',
                                    'terminal_K', 'terminal_R', 'average_number_of_K',
                                    'missclvg_0', 'missclvg_1', 'missclvg_2',
                                    'missclvg_0_area', 'missclvg_1_area', 'missclvg_2_area',
                                    '..K..K', '..K..R', '..R..K', '..R..R',
                                    '..K..Cterm', '..R..Cterm', '..K..', '..R..',
                                    '..K..K_area', '..K..R_area', '..R..K_area',
                                    '..R..R_area', '..K..Cterm_area', '..R..Cterm_area',
                                    '..K.._score', '..R.._score',
                                    'K_0', 'K_1', 'K_2', 'K_3',
                                    'K=1-terminalK', 'K=1-terminalR'])
result = result.sort(columns='Dataname')
result = result.T
result.to_csv('Summary_table.csv')

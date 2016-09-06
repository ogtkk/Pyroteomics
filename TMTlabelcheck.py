import Functions as func
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotter as plot
import os
import re
from datetime import datetime as dt

filename = '160319ko-TMT.txt'
data = func.ReadFile(filename)

time = dt.now()
f_time = time.strftime('%Y%m%d%H%M%S')
f_name = f_time + re.compile(u".txt").sub("", filename) + "_result"
if not os.path.isdir(f_name):
    os.mkdir(f_name)
os.chdir(f_name)
homedir = os.getcwd()

data = func.getSeqProperties(data)
data = func.TMTlabelEfficiency(data)
data = func.GRAVYScoreCalculator(data)
data = func.pICalculator(data)
data = func.extractPhospho(data)
dfs_temp = func.DivideMergeFile(data)
dfs_temp = {k: func.RemoveDupli(v) for k, v in dfs_temp.items()}
dfs = dfs_temp

#dfs['Sep-Pak_pH_4.5'] = func.MergeFiles([dfs_temp['160312ko02-SepPak_pH45.raw'], dfs_temp['160312ko08-SepPak_pH45.raw'], dfs_temp['160312ko14-SepPak_pH45.raw']])
#dfs['Sep-Pak_pH_6.5'] = func.MergeFiles([dfs_temp['160312ko03-SepPak_pH65.raw'], dfs_temp['160312ko09-SepPak_pH65.raw'], dfs_temp['160312ko15-SepPak_pH65.raw']])
#dfs['RP-C18_TFA'] = func.MergeFiles([dfs_temp['160312ko04-RP-C18_TFA.raw'], dfs_temp['160312ko10-RP-C18_TFA.raw'], dfs_temp['160312ko16-RP-C18_TFA.raw']])
#dfs['RP-C18_HFBA'] = func.MergeFiles([dfs_temp['160312ko05-RP-C18_HFBA.raw'], dfs_temp['160312ko11-RP-C18_HFBA.raw'], dfs_temp['160312ko17-RP-C18_HFBA.raw']])
#dfs['Solution'] = func.MergeFiles([dfs_temp['160312ko06-ref_solution.raw'], dfs_temp['160312ko12-ref_solution.raw'], dfs_temp['160312ko18-ref_solution.raw']])

plot.chargeStateID_barplot(dfs)
plot.plot_all(dfs)
plot.meanPeptideIDwithSD(dfs)
plot.mzViolinPlot(data)

dic = {'Dataname': [],
       'Fullarea': [],
       'Partialarea': [],
       'Nonearea': [],
       'Fullmeanarea': [],
       'Partialmeanarea': [],
       'Nonemeanarea': [],
       'N-termarea': [],
       'ID': [],
       'FullID': [],
       'PartialID': [],
       'NoneID': [],
       'N-termID': [],
       'S_area': [],
       'T_area': [],
       'Y_area': [],
       'H_area': [],
       'S_ID': [],
       'T_ID': [],
       'Y_ID': [],
       'H_ID': [],
       'TMT_1': [],
       'TMT_2': [],
       'TMT_3': [],
       'TMT_4': [],
       'TMT_1_area': [],
       'TMT_2_area': [],
       'TMT_3_area': [],
       'TMT_4_area': [],
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
    farea_list = v[v[u'efficiency'].str.contains('Full')][u'PepArea'].values
    parea_list = v[v[u'efficiency'].str.contains('Partial')][u'PepArea'].values
    narea_list = v[v[u'efficiency'].str.contains('None')][u'PepArea'].values
    ntarea_list = v[v[u'Mod'].str.contains('N-term')][u'PepArea'].values
    dic['Dataname'].append(k)
    dic['Fullarea'].append(np.sum(farea_list))
    dic['Partialarea'].append(np.sum(parea_list))
    dic['Nonearea'].append(np.sum(narea_list))
    dic['Fullmeanarea'].append(np.nanmean(farea_list))
    dic['Partialmeanarea'].append(np.nanmean(parea_list))
    dic['Nonemeanarea'].append(np.nanmean(narea_list))
    dic['N-termarea'].append(np.sum(ntarea_list))
    dic['ID'].append(len(v))
    dic['FullID'].append(len(v[v[u'efficiency'].str.contains('Full')]))
    dic['PartialID'].append(len(v[v[u'efficiency'].str.contains('Partial')]))
    dic['NoneID'].append(len(v[v[u'efficiency'].str.contains('None')]))
    dic['N-termID'].append(len(v[v[u'Mod'].str.contains('N-term')]))
    dic['S_area'].append(np.sum(v[v[u'number_of_S_TMT'] == 1][u'PepArea'].values))
    dic['T_area'].append(np.sum(v[v[u'number_of_T_TMT'] == 1][u'PepArea'].values))
    dic['Y_area'].append(np.sum(v[v[u'number_of_Y_TMT'] == 1][u'PepArea'].values))
    dic['H_area'].append(np.sum(v[v[u'number_of_H_TMT'] == 1][u'PepArea'].values))
    dic['S_ID'].append(len(v[v[u'number_of_S_TMT'] == 1]))
    dic['T_ID'].append(len(v[v[u'number_of_T_TMT'] == 1]))
    dic['Y_ID'].append(len(v[v[u'number_of_Y_TMT'] == 1]))
    dic['H_ID'].append(len(v[v[u'number_of_H_TMT'] == 1]))
    dic['TMT_1'].append(len(v[v[u'number_of_TMT'] == 1]))
    dic['TMT_2'].append(len(v[v[u'number_of_TMT'] == 2]))
    dic['TMT_3'].append(len(v[v[u'number_of_TMT'] == 3]))
    dic['TMT_4'].append(len(v[v[u'number_of_TMT'] == 4]))
    dic['TMT_1_area'].append(np.sum(v[v[u'number_of_TMT'] == 1][u'PepArea'].values))
    dic['TMT_2_area'].append(np.sum(v[v[u'number_of_TMT'] == 2][u'PepArea'].values))
    dic['TMT_3_area'].append(np.sum(v[v[u'number_of_TMT'] == 3][u'PepArea'].values))
    dic['TMT_4_area'].append(np.sum(v[v[u'number_of_TMT'] == 4][u'PepArea'].values))
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
result = pd.DataFrame(dic, columns=['Dataname', 'Fullarea', 'Partialarea',
                                    'Nonearea', 'N-termarea', 'Fullmeanarea',
                                    'Partialmeanarea', 'Nonemeanarea',
                                    'ID', 'FullID', 'PartialID', 'NoneID',
                                    'N-termID', 'S_area', 'T_area', 'Y_area',
                                    'H_area', 'S_ID', 'T_ID', 'Y_ID', 'H_ID',
                                    'TMT_1', 'TMT_2', 'TMT_3',
                                    'TMT_4', 'TMT_1_area', 'TMT_2_area',
                                    'TMT_3_area', 'TMT_4_area', 'phospho_1',
                                    'phospho_2', 'phospho_3', 'phospho_4', 'phospho_>4',
                                    'phospho_area_1', 'phospho_area_2', 'phospho_area_3',
                                    'phospho_area_4', 'phospho_area_>4',
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
result = result.sort_values(by='Dataname')
result = result.T
result.to_csv('Summary_table.csv')

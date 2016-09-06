import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st
import os
from collections import Counter
import re
import matplotlib as mpl
import Functions as func

data_TFA = func.ReadFile('160207ko14-RPsol_6plex_itraq.txt')
data_HFBA = func.ReadFile('160227ko27-9plex_full.txt')

reporter_values_TFA = (
    u'iTRAQ:126.127725', u'iTRAQ:127.131079',
    u'iTRAQ:128.128114', u'iTRAQ:129.131468',
    u'iTRAQ:130.141141', u'iTRAQ:131.138176')

reporter_values_HFBA = (u'iTRAQ:129.137787', u'iTRAQ:130.134822', u'iTRAQ:130.141141', 
                        u'iTRAQ:128.128114', u'iTRAQ:128.134433', u'iTRAQ:129.131468')

reporters_HFBA = {'RP-C18': [u'iTRAQ:129.137787', u'iTRAQ:130.134822', u'iTRAQ:130.141141'],
             'Solution': [u'iTRAQ:128.128114', u'iTRAQ:128.134433', u'iTRAQ:129.131468']}
reporters_TFA = {'RP-C18': [u'iTRAQ:126.127725', u'iTRAQ:127.131079', u'iTRAQ:128.128114'],
             'Solution': [u'iTRAQ:129.131468', u'iTRAQ:130.141141', u'iTRAQ:131.138176']}

df_TFA_temp = pd.DataFrame()
if bool(reporter_values_TFA) is True:
    for cn in reporter_values_TFA:
        series = data_TFA.ix[:, data_TFA.columns.map(
            lambda x: x.startswith(cn))]
        df_TFA_temp = pd.concat([df_TFA_temp, series], axis=1)
df_TFA_temp = df_TFA_temp.join(data_TFA[u'SeqModModDetail'])
df_TFA_temp = df_TFA_temp.groupby([u'SeqModModDetail']).mean()
df_TFA_temp = df_TFA_temp.replace(0, np.nan)
df_TFA_temp = df_TFA_temp.dropna()
for title, reporter in reporters_TFA.items():
    df_TFA_Int = pd.DataFrame()
    index_TFA = []
    average_TFA = []

    for cn in reporter:
        df_TFA_Int[cn] = df_TFA_temp[cn]
    for i in df_TFA_Int.index:
        index_TFA.append(i)
        average_TFA.append(np.nanmean(df_TFA_Int.loc[i]))
        df_TFA_temp[u'AVERAGE: ' + title] = pd.Series(average_TFA, index=index_TFA)


df_TFA_temp['log2'] = df_TFA_temp[u'AVERAGE: ' + 'RP-C18'] / df_TFA_temp[u'AVERAGE: ' + 'Solution']
df_TFA_temp['SeqModModDetail'] = df_TFA_temp.index
df_TFA = pd.merge(func.RemoveDupli(data_TFA), df_TFA_temp,
              on='SeqModModDetail', how='left')

df_HFBA_temp = pd.DataFrame()
if bool(reporter_values_HFBA) is True:
    for cn in reporter_values_HFBA:
        series = data_HFBA.ix[:, data_HFBA.columns.map(
            lambda x: x.startswith(cn))]
        df_HFBA_temp = pd.concat([df_HFBA_temp, series], axis=1)
df_HFBA_temp = df_HFBA_temp.join(data_HFBA[u'SeqModModDetail'])
df_HFBA_temp = df_HFBA_temp.groupby([u'SeqModModDetail']).mean()
df_HFBA_temp = df_HFBA_temp.replace(0, np.nan)
df_HFBA_temp = df_HFBA_temp.dropna()
for title, reporter in reporters_HFBA.items():
    df_HFBA_Int = pd.DataFrame()
    index_HFBA = []
    average_HFBA = []

    for cn in reporter:
        df_HFBA_Int[cn] = df_HFBA_temp[cn]
    for i in df_HFBA_Int.index:
        index_HFBA.append(i)
        average_HFBA.append(np.nanmean(df_HFBA_Int.loc[i]))
        df_HFBA_temp[u'AVERAGE: ' + title] = pd.Series(average_HFBA, index=index_HFBA)


df_HFBA_temp['log2'] = df_HFBA_temp[u'AVERAGE: ' + 'RP-C18'] / df_HFBA_temp[u'AVERAGE: ' + 'Solution']
median = np.median(df_HFBA_temp['log2'])
df_HFBA_temp['SeqModModDetail'] = df_HFBA_temp.index
df_HFBA = pd.merge(func.RemoveDupli(data_HFBA), df_HFBA_temp,
              on='SeqModModDetail', how='left')

df_merge = pd.merge(df_HFBA, df_TFA, on='SeqModModDetail', how='left', suffixes=('_HFBA', '_TFA'))
df_merge = df_merge.dropna(axis=0, subset=['log2_TFA', 'log2_HFBA'])
df_merge['HFBA/TFA'] = df_merge['log2_HFBA'] / df_merge['log2_TFA']
bins = map(lambda x: float(x) * 10, range(11))
df_merge['RetTime_bin'] = pd.cut(df_merge.RetTime_HFBA, bins)
fig = plt.figure()
fig.patch.set_color('white')
ax1 = fig.add_subplot(111)
whisprop = dict(linestyle='--', linewidth=1.25)
capprop = dict(linestyle='-', linewidth=1.5)
boxprop = dict(linestyle='-', linewidth=1.25)
medprop = dict(linestyle='-', linewidth=2.0)
flprop = dict(markersize=8, color='black')
bp = df_merge.boxplot(column='HFBA/TFA', by='RetTime_bin', grid=False, whis=[5,95],
           ax=ax1, positions=map(lambda x: float(x) + 0.5, range(10)),
           whiskerprops=whisprop, capprops=capprop, boxprops=boxprop,
           medianprops=medprop, flierprops=flprop, return_type='dict')
[[item.set_color('black') for item in bp[key]['medians']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['boxes']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['whiskers']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['caps']] for key in bp.keys()]

ax1.set_xticks(range(len(bins)))
ax1.set_xticklabels([0, 10,20,30,40,50,60,70,80,90,100])
ax1.set_xlim(2,9)
fig.suptitle('')
ax1.set_title('')
ax1.set_xlabel('Retention time (min)', labelpad=20)
ax1.set_ylabel('Relative reporter ion intensity\n(HFBA / TFA)')
ax1.axhline(1, linestyle='--', color='black')
plt.grid(False)
plt.grid(axis='y')
plt.subplots_adjust(bottom=0.2, left=0.2, right=0.8)
plt.savefig('boxplot_xRetTyInt.png')
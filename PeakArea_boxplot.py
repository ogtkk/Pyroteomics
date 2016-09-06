import Functions as func
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
from datetime import datetime as dt
import scipy.stats as st
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator
from beeswarm import *

plt.rcParams['font.family'] = 'sans-serif' 
plt.rcParams['font.size'] = 24
plt.rcParams['xtick.major.pad'] = 12
plt.rcParams['xtick.major.size'] = 0 
plt.rcParams['xtick.major.width'] = 0
plt.rcParams['ytick.major.pad'] = 12
plt.rcParams['ytick.major.size'] = 6 
plt.rcParams['ytick.major.width'] = 1.0
plt.rcParams['xtick.minor.size'] = 0 
plt.rcParams['xtick.minor.width'] = 0
plt.rcParams['xtick.minor.pad'] = 60
plt.rcParams['figure.figsize'] = [30/2.54, 30/2.54]
plt.rcParams['mathtext.default'] = 'regular'


filename = '160312ko-tipTMT_full.txt'
data = func.ReadFile(filename)
time = dt.now()
f_time = time.strftime('%Y%m%d%H%M%S')
f_name = 'For_paper'
if not os.path.isdir(f_name):
    os.mkdir(f_name)
os.chdir(f_name)
homedir = os.getcwd()

data = func.extractPhospho(data)
data = func.getSeqProperties(data)
dfs_temp = func.DivideMergeFile(data)
dfs_temp = {k: func.RemoveDupli(v) for k, v in dfs_temp.items()}
dfs = {}

dfs['Sep-Pak_pH_6.5'] = func.MergeFiles([dfs_temp['160312ko03-SepPak_pH65.raw'], dfs_temp['160312ko09-SepPak_pH65.raw'], dfs_temp['160312ko15-SepPak_pH65.raw']])
dfs['RP-C18_TFA'] = func.MergeFiles([dfs_temp['160312ko04-RP-C18_TFA.raw'], dfs_temp['160312ko10-RP-C18_TFA.raw'], dfs_temp['160312ko16-RP-C18_TFA.raw']])
dfs['RP-C18_HFBA'] = func.MergeFiles([dfs_temp['160312ko05-RP-C18_HFBA.raw'], dfs_temp['160312ko11-RP-C18_HFBA.raw'], dfs_temp['160312ko17-RP-C18_HFBA.raw']])
dfs['Solution'] = func.MergeFiles([dfs_temp['160312ko06-ref_solution.raw'], dfs_temp['160312ko12-ref_solution.raw'], dfs_temp['160312ko18-ref_solution.raw']])
df_merge = pd.merge(dfs['RP-C18_TFA'], dfs['RP-C18_HFBA'], on='SeqModModDetail', how='inner', suffixes=['_TFA', '_HFBA'])
df_merge_2 = pd.merge(dfs['RP-C18_HFBA'], dfs['Sep-Pak_pH_6.5'], on='SeqModModDetail', how='inner', suffixes=['_HFBA', '_SepPak'])
df_merge_3 = pd.merge(dfs['RP-C18_TFA'], dfs['Sep-Pak_pH_6.5'], on='SeqModModDetail', how='inner', suffixes=['_TFA', '_SepPak'])
df_merge_4 = pd.merge(dfs['RP-C18_TFA'], dfs['Solution'], on='SeqModModDetail', how='inner', suffixes=['_TFA', '_Solution'])
df_merge_3['Arearatio'] = np.log2(df_merge_3['PepArea_TFA'] / df_merge_3['PepArea_SepPak'])
df_merge_4['Arearatio'] = np.log2(df_merge_4['PepArea_TFA'] / df_merge_4['PepArea_Solution'])

fig = plt.figure()
fig.patch.set_color('white')
ax1 = fig.add_subplot(111)
ax1.hist(df_merge_4['Arearatio'], bins=80, range=[-4,4], color='0.4', edgecolor='0.4')
plt.show()
plt.clf()

plt.rcParams['xtick.major.size'] = 6 
plt.rcParams['xtick.major.width'] = 1.0

bins = map(lambda x: float(x) * 5, range(11))
df_merge_4['RetTime_bin'] = pd.cut(df_merge_4.RetTime_TFA, bins)
fig = plt.figure()
fig.patch.set_color('white')
ax1 = fig.add_subplot(111)
whisprop = dict(linestyle='--', linewidth=1.25)
capprop = dict(linestyle='-', linewidth=1.5)
boxprop = dict(linestyle='-', linewidth=1.25)
medprop = dict(linestyle='-', linewidth=2.0)
flprop = dict(markersize=8, color='black')
bp = df_merge_4.boxplot(column='Arearatio', by='RetTime_bin', grid=False, whis=[5,95],
           ax=ax1, positions=map(lambda x: float(x) + 0.5, range(10)),
           whiskerprops=whisprop, capprops=capprop, boxprops=boxprop,
           medianprops=medprop, flierprops=flprop, return_type='dict')
[[item.set_color('black') for item in bp[key]['medians']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['boxes']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['whiskers']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['caps']] for key in bp.keys()]

ax1.set_xticks(range(len(bins)))
ax1.set_xticklabels([0,5,10,15,20,25,30,35,40,45,50])
ax1.set_xlim(4,9)
fig.suptitle('')
ax1.set_title('')
ax1.set_xlabel('Retention time (min)', labelpad=20)
ax1.set_ylabel('Relative Peak Area\n($Log_{2}$ TFA / Solution)')
plt.grid(False)
plt.grid(axis='y')
plt.subplots_adjust(bottom=0.2, left=0.2, right=0.8)
plt.show()
plt.clf()


df_merge_u = pd.merge(dfs['RP-C18_TFA'], dfs['Sep-Pak_pH_6.5'], on='SeqModModDetail', how='outer', suffixes=['_TFA', '_SepPak'])
df_merge_u.to_csv('testarea.csv')
df_merge_u = df_merge_u.dropna(subset=['Seq_SepPak'])
plt.hist(df_merge_u.RetTime_TFA, bins=10, range=[0,50], color='0.4', edgecolor='0.4')
plt.title('n=%d' % len(df_merge_u))
plt.show()
plt.clf()
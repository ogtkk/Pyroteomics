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
import plotter as plot
plt.rcParams['font.family'] = 'sans-serif' 
plt.rcParams['font.size'] = 24
plt.rcParams['xtick.major.pad'] = 12
plt.rcParams['xtick.major.size'] = 6 
plt.rcParams['xtick.major.width'] = 1.0
plt.rcParams['ytick.major.pad'] = 12
plt.rcParams['ytick.major.size'] = 6 
plt.rcParams['ytick.major.width'] = 1.0
plt.rcParams['xtick.minor.size'] = 0 
plt.rcParams['xtick.minor.width'] = 0
plt.rcParams['xtick.minor.pad'] = 60
plt.rcParams['figure.figsize'] = [30/2.54, 30/2.54]
plt.rcParams['mathtext.default'] = 'regular'


filename = '160227ko27-9plex_full.txt'
data = func.ReadFile(filename)
time = dt.now()
f_time = time.strftime('%Y%m%d%H%M%S')
f_name = 'For_paper'
if not os.path.isdir(f_name):
    os.mkdir(f_name)
os.chdir(f_name)
homedir = os.getcwd()

data = func.extractPhospho(data)
#plot.reporterionIntensity_boxplot(data)
data = func.getSeqProperties(data)
df = data
#data = func.GRAVYScoreCalculator(data)
#data = func.pICalculator(data)
plot.reproducibilityCheck(data, tagname='TMT10plex')
#plot.log2histID(data, tagname='160227')
#data.to_csv('160308ko-testA.csv')
#data[['iTRAQ:126.127725', 'iTRAQ:127.124760', 'iTRAQ:127.131079']] =\
#    data[['iTRAQ:126.127725', 'iTRAQ:127.124760', 'iTRAQ:127.131079']].replace(0, np.nan)
#data = data.dropna(axis=0, subset=['iTRAQ:126.127725', 'iTRAQ:127.124760', 'iTRAQ:127.131079'])
#data.to_csv('160308ko-test.csv')

reporter_values = (
    u'iTRAQ:126.127725', u'iTRAQ:127.124760',
    u'iTRAQ:127.131079', u'iTRAQ:128.128114',
    u'iTRAQ:128.134433', u'iTRAQ:129.131468',
    u'iTRAQ:129.137787', u'iTRAQ:130.134822',
    u'iTRAQ:130.141141', u'iTRAQ:131.138176')

reporters = {'RP-C18': [u'iTRAQ:129.137787', u'iTRAQ:130.134822', u'iTRAQ:130.141141'],
             'Solution': [u'iTRAQ:128.128114', u'iTRAQ:128.134433', u'iTRAQ:129.131468']}

df_temp = pd.DataFrame()
df_temp = df.iloc[:, df.columns.isin(reporter_values)]
df_temp = df_temp.join(df[u'SeqModModDetail'])
df_temp = df_temp.groupby([u'SeqModModDetail']).mean()
df_temp = df_temp.replace(0, np.nan)
df_temp = df_temp.dropna(subset=[u'iTRAQ:128.128114',
            u'iTRAQ:128.134433', u'iTRAQ:129.131468',
            u'iTRAQ:129.137787', u'iTRAQ:130.134822',
            u'iTRAQ:130.141141'])
for title, reporter in reporters.items():
    df_Int = pd.DataFrame()
    index = []
    average = []
    for cn in reporter:
        df_Int[cn] = df_temp[cn]
    for i in df_Int.index:
        index.append(i)
        average.append(np.nanmean(df_Int.loc[i]))
        df_temp[u'AVERAGE: ' + title] = pd.Series(average, index=index)
df_temp['log2'] = np.log2(df_temp[u'AVERAGE: ' + 'RP-C18'] / df_temp[u'AVERAGE: ' + 'Solution'])
median = np.median(df_temp['log2'])
df_temp['SeqModModDetail'] = df_temp.index
df = pd.merge(func.RemoveDupli(df), df_temp,
              on='SeqModModDetail', how='left')
bins = map(lambda x: float(x) * 10, range(11))
df['RetTime_bin'] = pd.cut(df.RetTime, bins)
fig = plt.figure()
fig.patch.set_color('white')
ax1 = fig.add_subplot(111)
whisprop = dict(linestyle='--', linewidth=1.25)
capprop = dict(linestyle='-', linewidth=1.5)
boxprop = dict(linestyle='-', linewidth=1.25)
medprop = dict(linestyle='-', linewidth=2.0)
flprop = dict(markersize=8, color='black')
bp = df.boxplot(column='log2', by='RetTime_bin', grid=False,
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
ax1.set_ylabel('Relative reporter ion intensity\n($Log_{2}$ RP-C18 / Solution)')
plt.grid(False)
plt.grid(axis='y')
plt.subplots_adjust(bottom=0.2, left=0.2, right=0.8)
plt.savefig('boxplot_xRetTyInt.png')
plt.clf()

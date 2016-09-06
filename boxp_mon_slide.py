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
#from beeswarm import *

plt.rcParams['font.family'] = 'sans-serif' 
plt.rcParams['font.size'] = 40
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

fnamepath = re.compile(u"D.+ManualInputFile\\\\|D.+ManualInputFile/")
fnameext = re.compile(u".raw")
filename_1 = '160403ko-tipTMT-full.txt'
filename_2 = '160403ko-tipTMT-full_2.txt'
data_1 = func.ReadFile(filename_1)
data_2 = func.ReadFile(filename_2)
data = data_2
#pd.concat([data_1, data_2], ignore_index=True)
data[u'RawdataFile'] = data[u'RawdataFile'].apply(
    lambda x: fnamepath.sub("", x))
data[u'RawdataFile'] = data[u'RawdataFile'].apply(
    lambda x: fnameext.sub("", x))
data['RetTime'] = data['RetTime'] + 135 * (data['RawdataFile'].str[-1:].astype('int') - 1)
data.loc[data['RawdataFile'].str[-1:].astype('int') != 1, 'RetTime']\
    = data.loc[data['RawdataFile'].str[-1:].astype('int') != 1, 'RetTime'] - 15
time = dt.now()
f_time = time.strftime('%Y%m%d%H%M%S')
f_name = 'For_paper'
if not os.path.isdir(f_name):
    os.mkdir(f_name)
os.chdir(f_name)
os.chdir('160403')
homedir = os.getcwd()

data = func.extractPhospho(data)
data = func.getSeqProperties(data)
df = data

reporter_values = (
    u'iTRAQ:126.127725', u'iTRAQ:127.124760',
    u'iTRAQ:127.131079', u'iTRAQ:128.128114',
    u'iTRAQ:128.134433', u'iTRAQ:129.131468',
    u'iTRAQ:129.137787', u'iTRAQ:130.134822',
    u'iTRAQ:130.141141', u'iTRAQ:131.138176')

reporters = {'Solution': [u'iTRAQ:126.127725', u'iTRAQ:127.124760',
                          u'iTRAQ:127.131079', u'iTRAQ:128.128114', u'iTRAQ:128.134433'],
             'Tip': [u'iTRAQ:129.131468', u'iTRAQ:129.137787', u'iTRAQ:130.134822',
                     u'iTRAQ:130.141141', u'iTRAQ:131.138176']}

df_temp = pd.DataFrame()
df_temp = df.iloc[:, df.columns.isin(reporter_values)]
df_temp = df_temp.join(df[u'SeqModModDetail'])
df_temp = df_temp.groupby([u'SeqModModDetail']).mean()
df_temp = df_temp.replace(0, np.nan)
dtA = pd.DataFrame(df_temp[reporters['Solution']])
dtB = pd.DataFrame(df_temp[reporters['Tip']])
df_temp['Solution_N'] = 5 - dtA.isnull().sum(axis=1)
df_temp['Tip_N'] = 5 - dtB.isnull().sum(axis=1)
df_temp = df_temp[df_temp['Solution_N'] >= 3]
df_temp = df_temp[df_temp['Tip_N'] >= 3]
#df_temp = df_temp.dropna(subset=[u'iTRAQ:126.127725', u'iTRAQ:127.124760',
#                                 u'iTRAQ:127.131079', u'iTRAQ:128.128114',
#                                 u'iTRAQ:128.134433', u'iTRAQ:129.131468',
#                                 u'iTRAQ:129.137787', u'iTRAQ:130.134822',
#                                 u'iTRAQ:130.141141', u'iTRAQ:131.138176'])
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
df_temp['log2'] = np.log2(df_temp[u'AVERAGE: ' + 'Tip'] / df_temp[u'AVERAGE: ' + 'Solution'])
median = np.median(df_temp['log2'])
df_temp['SeqModModDetail'] = df_temp.index
df = pd.merge(func.RemoveDupli(df), df_temp,
              on='SeqModModDetail', how='left')

fig = plt.figure()
fig.patch.set_color('white')
ax1 = fig.add_subplot(111)
ax1.hist(df['log2'], bins=80, range=[-6,6], color='0.4', edgecolor='0.4')
ax1.axvline(median, color='black', linestyle='--', linewidth=1.5)
ax1.text(median + 0.5, 400, 'median = %f' % median)
ax1.set_xlabel('Relative reporter ion intensity\nLog$_2$ (RP-C18 / Solution)')
ax1.set_ylabel('Number of phosphopeptides')
plt.subplots_adjust(bottom=0.2, left=0.2, right=0.95)
plt.savefig('hist_HFBA160403_2.png')
plt.clf()

plt.rcParams['xtick.major.size'] = 6 
plt.rcParams['xtick.major.width'] = 1.0

bins = list(map(lambda x: float(x) * 40, range(14)))
df['RetTime_bin'] = pd.cut(df.RetTime, bins)
fig = plt.figure()
fig.patch.set_color('white')
ax1 = fig.add_subplot(111)
whisprop = dict(linestyle='--', linewidth=1.25)
capprop = dict(linestyle='-', linewidth=1.5)
boxprop = dict(linestyle='-', linewidth=1.25)
medprop = dict(linestyle='-', linewidth=2.0)
flprop = dict(markersize=8, color='black')
bp = df.boxplot(column='log2', by='RetTime_bin', grid=False, whis=[5,95],
           ax=ax1, positions=list(map(lambda x: float(x) + 0.5, range(13))),
           whiskerprops=whisprop, capprops=capprop, boxprops=boxprop,
           medianprops=medprop, flierprops=flprop, return_type='dict')
[[item.set_color('black') for item in bp[key]['medians']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['boxes']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['whiskers']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['caps']] for key in bp.keys()]

ax1.set_xticks(range(len(bins)))
ax1.set_xticklabels([0,40,80,120,160,200,240,280,320,360,400,440,480,520,560], size=20)
ax1.set_xlim(1,13)
ax1.set_ylim(-6, 6)
fig.suptitle('')
ax1.set_title('(n=%d)' % len(df_temp))
ax1.set_xlabel('Retention time (min)', labelpad=20)
ax1.set_ylabel('Relative reporter ion intensity\n($Log_{2}$ RP-C18 / Solution)')
plt.grid(False)
plt.grid(axis='y')
plt.subplots_adjust(bottom=0.2, left=0.2, right=0.8)
plt.savefig('boxplot_xRetTyInt-monHFBA160403_2.png')
plt.clf()

os.chdir('../')
os.chdir('../')
data_TFA = func.ReadFile('160319ko-solTFA.txt')
data_HFBA = func.ReadFile('160319ko-solHFBA.txt')
os.chdir('For_paper')
os.chdir('test')
data_TFA[u'RawdataFile'] = data_TFA[u'RawdataFile'].apply(
    lambda x: fnamepath.sub("", x))
data_TFA[u'RawdataFile'] = data_TFA[u'RawdataFile'].apply(
    lambda x: fnameext.sub("", x))
data_TFA['RetTime'] = data_TFA['RetTime'] + 120 * (data_TFA['RawdataFile'].str[-1:].astype('int') - 1)
data_HFBA[u'RawdataFile'] = data_HFBA[u'RawdataFile'].apply(
    lambda x: fnamepath.sub("", x))
data_HFBA[u'RawdataFile'] = data_HFBA[u'RawdataFile'].apply(
    lambda x: fnameext.sub("", x))
data_HFBA['RetTime'] = data_HFBA['RetTime'] + 120 * (data_HFBA['RawdataFile'].str[-1:].astype('int') - 1)


reporter_values_TFA = (
    u'iTRAQ:126.127725', u'iTRAQ:127.124760',
    u'iTRAQ:127.131079', u'iTRAQ:128.128114',
    u'iTRAQ:128.134433', u'iTRAQ:129.131468',
    u'iTRAQ:129.137787', u'iTRAQ:130.134822',
    u'iTRAQ:130.141141', u'iTRAQ:131.138176')

reporters_TFA = {'Solution': [u'iTRAQ:126.127725', u'iTRAQ:127.124760',
                              u'iTRAQ:127.131079', u'iTRAQ:128.128114', u'iTRAQ:128.134433'],
                 'Tip': [u'iTRAQ:129.131468', u'iTRAQ:129.137787', u'iTRAQ:130.134822',
                         u'iTRAQ:130.141141', u'iTRAQ:131.138176']}

reporter_values_HFBA = (
    u'iTRAQ:126.127725', u'iTRAQ:127.124760',
    u'iTRAQ:127.131079', u'iTRAQ:128.128114',
    u'iTRAQ:128.134433', u'iTRAQ:129.131468',
    u'iTRAQ:129.137787', u'iTRAQ:130.134822',
    u'iTRAQ:130.141141', u'iTRAQ:131.138176')

reporters_HFBA = {'Solution': [u'iTRAQ:126.127725', u'iTRAQ:127.124760',
                               u'iTRAQ:127.131079', u'iTRAQ:128.128114', u'iTRAQ:128.134433'],
                  'Tip': [u'iTRAQ:129.131468', u'iTRAQ:129.137787', u'iTRAQ:130.134822',
                          u'iTRAQ:130.141141', u'iTRAQ:131.138176']}


df_TFA_temp = pd.DataFrame()
if bool(reporter_values_TFA) is True:
    for cn in reporter_values_TFA:
        series = data_TFA.ix[:, data_TFA.columns.map(
            lambda x: x.startswith(cn)))]
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


df_TFA_temp['log2'] = np.log2(df_TFA_temp[u'AVERAGE: ' + 'Tip'] / df_TFA_temp[u'AVERAGE: ' + 'Solution'])
df_TFA_temp['ratio'] = df_TFA_temp[u'AVERAGE: ' + 'Tip'] / df_TFA_temp[u'AVERAGE: ' + 'Solution']
df_TFA_temp['SeqModModDetail'] = df_TFA_temp.index
df_TFA = pd.merge(func.RemoveDupli(data_TFA), df_TFA_temp,
              on='SeqModModDetail', how='left')

df_HFBA_temp = pd.DataFrame()
if bool(reporter_values_HFBA) is True:
    for cn in reporter_values_HFBA:
        series = data_HFBA.ix[:, data_HFBA.columns.map(
            lambda x: x.startswith(cn)))]
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


df_HFBA_temp['log2'] = np.log2(df_HFBA_temp[u'AVERAGE: ' + 'Tip'] / df_HFBA_temp[u'AVERAGE: ' + 'Solution'])
df_HFBA_temp['ratio'] = df_HFBA_temp[u'AVERAGE: ' + 'Tip'] / df_HFBA_temp[u'AVERAGE: ' + 'Solution']
median = np.median(df_HFBA_temp['log2'])
df_HFBA_temp['SeqModModDetail'] = df_HFBA_temp.index
df_HFBA = pd.merge(func.RemoveDupli(data_HFBA), df_HFBA_temp,
              on='SeqModModDetail', how='left')

df_merge = pd.merge(df_HFBA, df_TFA, on='SeqModModDetail', how='inner', suffixes=('_HFBA', '_TFA'))
df_merge = df_merge.dropna(axis=0, subset=['log2_TFA', 'log2_HFBA'])
df_merge['HFBA/TFA'] = np.log2(df_merge['ratio_HFBA'] / df_merge['ratio_TFA'])
plt.rcParams['xtick.major.size'] = 6 
plt.rcParams['xtick.major.width'] = 1.0

bins = list(map(lambda x: float(x) * 40, range(14))
df_merge['RetTime_bin'] = pd.cut(df_merge.RetTime_HFBA, bins)
fig = plt.figure()
fig.patch.set_color('white')
ax1 = fig.add_subplot(111)
whisprop = dict(linestyle='--', linewidth=2.25)
capprop = dict(linestyle='-', linewidth=2.5)
boxprop = dict(linestyle='-', linewidth=2.25)
medprop = dict(linestyle='-', linewidth=3.0)
flprop = dict(markersize=8, color='black')
bp = df_merge.boxplot(column='HFBA/TFA', by='RetTime_bin', grid=False, whis=[5,95],
           ax=ax1, positions=list(map(lambda x: float(x) + 0.5, range(13))),
           whiskerprops=whisprop, capprops=capprop, boxprops=boxprop,
           medianprops=medprop, flierprops=flprop, return_type='dict')
[[item.set_color('black') for item in bp[key]['medians']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['boxes']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['whiskers']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['caps']] for key in bp.keys()]

ax1.set_xticks(range(len(bins)))
ax1.set_xticklabels([0,'',80,'',160,'',240,'',320,'',400,'',480,'',560], size=32)
ax1.set_xlim(3,11)
ax1.set_ylim(-2,3)
fig.suptitle('')
ax1.set_title('')
ax1.set_xlabel('Retention time (min)', labelpad=20)
ax1.set_ylabel('Relative reporter ion intensity\n$Log_2$(Tip_HFBA / Tip_TFA)')
ax1.set_title('')
plt.grid(False)
plt.grid(axis='y')
plt.subplots_adjust(bottom=0.2, left=0.2, right=0.8)
plt.savefig('boxplot_xRetTyInt-monHFBA-TFA.png')
plt.clf()

bins = list(map(lambda x: float(x) * 40, range(14)))
df_TFA['RetTime_bin'] = pd.cut(df_TFA.RetTime, bins)
fig = plt.figure()
fig.patch.set_color('white')
ax1 = fig.add_subplot(111)
whisprop = dict(linestyle='--', linewidth=2.25)
capprop = dict(linestyle='-', linewidth=2.5)
boxprop = dict(linestyle='-', linewidth=2.25)
medprop = dict(linestyle='-', linewidth=3.0)
flprop = dict(markersize=8, color='black')
bp = df_TFA.boxplot(column='log2', by='RetTime_bin', grid=False, whis=[5,95],
           ax=ax1, positions=list(map(lambda x: float(x) + 0.5, range(13))),
           whiskerprops=whisprop, capprops=capprop, boxprops=boxprop,
           medianprops=medprop, flierprops=flprop, return_type='dict')
[[item.set_color('black') for item in bp[key]['medians']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['boxes']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['whiskers']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['caps']] for key in bp.keys()]

ax1.set_xticks(range(len(bins)))
ax1.set_xticklabels([0,'',80,'',160,'',240,'',320,'',400,'',480,'',560], size=32)
ax1.set_xlim(3,11)
ax1.set_ylim(-4, 4)
fig.suptitle('')
ax1.set_title('')
ax1.set_xlabel('Retention time (min)', labelpad=20)
ax1.set_ylabel('Relative reporter ion intensity\n($Log_{2}$ Tip / Solution)')
plt.grid(False)
plt.grid(axis='y')
plt.subplots_adjust(bottom=0.2, left=0.2, right=0.8)
plt.savefig('boxplot_xRetTyInt-monTFA.png')
plt.clf()

bins = list(map(lambda x: float(x) * 40, range(14)))
df_HFBA['RetTime_bin'] = pd.cut(df_HFBA.RetTime, bins)
fig = plt.figure()
fig.patch.set_color('white')
ax1 = fig.add_subplot(111)
whisprop = dict(linestyle='--', linewidth=2.25)
capprop = dict(linestyle='-', linewidth=2.5)
boxprop = dict(linestyle='-', linewidth=2.25)
medprop = dict(linestyle='-', linewidth=3.0)
flprop = dict(markersize=8, color='black')
bp = df_HFBA.boxplot(column='log2', by='RetTime_bin', grid=False, whis=[5,95],
           ax=ax1, positions=list(map(lambda x: float(x) + 0.5, range(13))),
           whiskerprops=whisprop, capprops=capprop, boxprops=boxprop,
           medianprops=medprop, flierprops=flprop, return_type='dict')
[[item.set_color('black') for item in bp[key]['medians']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['boxes']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['whiskers']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['caps']] for key in bp.keys()]

ax1.set_xticks(range(len(bins)))
ax1.set_xticklabels([0,'',80,'',160,'',240,'',320,'',400,'',480,'',560], size=32)
ax1.set_xlim(3,11)
ax1.set_ylim(-4, 4)
fig.suptitle('')
ax1.set_title('')
ax1.set_xlabel('Retention time (min)', labelpad=20)
ax1.set_ylabel('Relative reporter ion intensity\n($Log_{2}$ Method 3 / Method 2)')
plt.grid(False)
plt.grid(axis='y')
plt.subplots_adjust(bottom=0.2, left=0.2, right=0.8)
plt.savefig('boxplot_xRetTyInt-monHFBA_ps.png')
plt.clf()

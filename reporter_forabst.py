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

plt.rcParams['font.family'] = 'sans-serif' 
plt.rcParams['font.size'] = 24
plt.rcParams['xtick.major.pad'] = 16
plt.rcParams['xtick.major.size'] = 0 
plt.rcParams['xtick.major.width'] = 0
plt.rcParams['xtick.minor.size'] = 0 
plt.rcParams['xtick.minor.width'] = 0
plt.rcParams['xtick.minor.pad'] = 60
plt.rcParams['figure.figsize'] = [30/2.54, 30/2.54]


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
data = func.getSeqProperties(data)
data = func.GRAVYScoreCalculator(data)
#data = func.pICalculator(data)

reporters = {'RP-C18': ['iTRAQ:129.137787', 'iTRAQ:130.134822', 'iTRAQ:130.141141'],
             'Solution': ['iTRAQ:128.128114', 'iTRAQ:128.134433', 'iTRAQ:129.131468'],
             'TMT-HM': ['iTRAQ:126.127725', 'iTRAQ:127.124760', 'iTRAQ:127.131079']}

data = data.replace(0, np.nan)
df_Average = pd.DataFrame()
df = func.RemoveDupli(data)
for title, reporter in reporters.items():
    cn = 'AVERAGE'
    df_temp = data.dropna(subset=reporter)
    df_temp[cn] = df_temp.loc[:, reporter].mean(axis=1)
    mean_s = df_temp[[cn, 'SeqModModDetail']].groupby(['SeqModModDetail']).mean()
    mean_s['SeqModModDetail'] = mean_s.index
    df = pd.merge(df, mean_s, on='SeqModModDetail', how='left', suffixes=['','_'+title])
df.dropna(subset=['AVERAGE', 'AVERAGE_Solution'])
df['log2'] = np.log2(df['AVERAGE'] / df['AVERAGE_Solution'])

fig = plt.figure()
fig.patch.set_color('white')
ax1 = fig.add_subplot(111)
ax1.hist(df['log2'], bins=40, range=[-6,2], color='0.4', edgecolor='0.4')
plt.show()
plt.clf()


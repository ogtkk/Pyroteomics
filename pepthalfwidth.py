
import Functions as func
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import lines
import os
import re
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

fnamepath = re.compile(u"D.+ManualInputFile\\\\|D.+ManualInputFile/")
fnameext = re.compile(u".raw")
filename = '160403ko-batch2_merge.txt'
data = func.ReadFile(filename)
data = func.RemoveDupli(data)
data[u'RawdataFile'] = data[u'RawdataFile'].apply(
    lambda x: fnamepath.sub("", x))
data[u'RawdataFile'] = data[u'RawdataFile'].apply(
    lambda x: fnameext.sub("", x))
data['RetTime'] = data['RetTime'] + 135 * (data['RawdataFile'].str[-1:].astype('int') - 1)
data.loc[data['RawdataFile'].str[-1:].astype('int') != 1, 'RetTime']\
    = data.loc[data['RawdataFile'].str[-1:].astype('int') != 1, 'RetTime'] - 15
f_name = 'For_paper'
if not os.path.isdir(f_name):
    os.mkdir(f_name)
os.chdir(f_name)
os.chdir('160403')
homedir = os.getcwd()

data = func.extractPhospho(data)
data = func.getSeqProperties(data)
df = data
df = df.dropna(subset=['PepHalfWidth', 'RetTime'])
df['PepHalfWidth'] = df['PepHalfWidth'].convert_objects(convert_numeric=True)


fig = plt.figure()
fig.patch.set_color('white')
ax1 = fig.add_subplot(111)
ax1.hist(df['PepHalfWidth'], bins=20, range=[0,4], color='0.4', edgecolor='0.4')
ax1.set_xlabel('HalfWidth')
ax1.set_ylabel('Number of phosphopeptides')
plt.subplots_adjust(bottom=0.2, left=0.2, right=0.95)
plt.savefig('hist_Halfwidth_2.png')
plt.clf()

plt.rcParams['xtick.major.size'] = 6 
plt.rcParams['xtick.major.width'] = 1.0

bins = map(lambda x: float(x) * 40, range(14))
df['RetTime_bin'] = pd.cut(df.RetTime, bins)
fig = plt.figure()
fig.patch.set_color('white')
ax1 = fig.add_subplot(111)
whisprop = dict(linestyle='--', linewidth=1.25)
capprop = dict(linestyle='-', linewidth=1.5)
boxprop = dict(linestyle='-', linewidth=1.25)
medprop = dict(linestyle='-', linewidth=2.0)
flprop = dict(markersize=8, color='black')
bp = df.boxplot(column='PepHalfWidth', by='RetTime_bin', grid=False, whis=[5,95],
           ax=ax1, positions=map(lambda x: float(x) + 0.5, range(13)),
           whiskerprops=whisprop, capprops=capprop, boxprops=boxprop,
           medianprops=medprop, flierprops=flprop, return_type='dict')
[[item.set_color('black') for item in bp[key]['medians']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['boxes']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['whiskers']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['caps']] for key in bp.keys()]

ax1.set_xticks(range(len(bins)))
ax1.set_xticklabels([0,40,80,120,160,200,240,280,320,360,400,440,480,520,560], size=20)
ax1.set_xlim(1,13)
ax1.set_ylim(0, 4)
fig.suptitle('')
ax1.set_title('(n=%d)' % len(data))
ax1.set_xlabel('Retention time (min)', labelpad=20)
ax1.set_ylabel('Halfwidth')
plt.grid(False)
plt.grid(axis='y')
plt.subplots_adjust(bottom=0.2, left=0.2, right=0.8)
plt.savefig('boxplot_xRetTyHalfw_2.png')
plt.clf()

#df = df.dropna(subset=['PepHalfWidth', 'RetTime'])
#x = df['RetTime'].convert_objects(convert_numeric=True).values
#y = df['PepHalfWidth'].values
#slope, intercept, r_value, p_value, std_err = st.linregress(x,y)
#print slope, intercept, r_value, p_value, std_err, x, y
#func = lambda x: x * slope + intercept
#fig.patch.set_color('white')
#ax1 = fig.add_subplot(111)
#ax1.scatter(x, y, c='0.4', s=4)
#line = lines.Line2D([0, 570], [func(0), func(570)])
#ax1.add_line(line)
#ax1.set_xlim(0,570)
#ax1.set_ylim(0, 4)
#fig.suptitle('')
#ax1.set_title('(n=%d)' % len(data))
#ax1.set_xlabel('Retention time (min)', labelpad=20)
#ax1.set_ylabel('HalfWidth')
#plt.grid(False)
#plt.grid(axis='y')
#plt.subplots_adjust(bottom=0.2, left=0.2, right=0.8)
#plt.savefig('scatter_xRetTyHalfW_2.png')
#plt.clf()


fig.patch.set_color('white')
ax1 = fig.add_subplot(111)
whisprop = dict(linestyle='--', linewidth=1.25)
capprop = dict(linestyle='-', linewidth=1.5)
boxprop = dict(linestyle='-', linewidth=1.25)
medprop = dict(linestyle='-', linewidth=2.0)
flprop = dict(markersize=8, color='black')
df = df.dropna(subset=['PepSymmetryFacter'])
df['PepSymmetryFacter'] = df['PepSymmetryFacter'].convert_objects(convert_numeric=True) * 100
bp = df.boxplot(column='PepSymmetryFacter', by='RetTime_bin', grid=False, whis=[5,95],
           ax=ax1, positions=map(lambda x: float(x) + 0.5, range(13)),
           whiskerprops=whisprop, capprops=capprop, boxprops=boxprop,
           medianprops=medprop, flierprops=flprop, return_type='dict')
[[item.set_color('black') for item in bp[key]['medians']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['boxes']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['whiskers']] for key in bp.keys()]
[[item.set_color('black') for item in bp[key]['caps']] for key in bp.keys()]

ax1.set_xticks(range(len(bins)))
ax1.set_xticklabels([0,40,80,120,160,200,240,280,320,360,400,440,480,520,560], size=20)
ax1.set_xlim(1,13)
ax1.set_ylim(0, 400)
fig.suptitle('')
ax1.set_title('(n=%d)' % len(data))
ax1.set_xlabel('Retention time (min)', labelpad=20)
ax1.set_ylabel('SymmetryFacter')
plt.grid(False)
plt.grid(axis='y')
plt.subplots_adjust(bottom=0.2, left=0.2, right=0.8)
plt.savefig('boxplot_xRetTySymm_2.png')
plt.clf()

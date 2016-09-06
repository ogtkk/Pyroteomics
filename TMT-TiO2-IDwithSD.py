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
plt.rcParams['xtick.major.pad'] = 16
plt.rcParams['xtick.major.size'] = 0 
plt.rcParams['xtick.major.width'] = 0
plt.rcParams['xtick.minor.size'] = 0 
plt.rcParams['xtick.minor.width'] = 0
plt.rcParams['xtick.minor.pad'] = 60
plt.rcParams['figure.figsize'] = [30/2.54, 30/2.54]
plt.rcParams['mathtext.default'] = 'regular'


filename = '../ResultFiles/2016-08/160804/160319ko-TMT-TiO2-fix.txt'
data = func.ReadFile(filename)

time = dt.now()
f_time = time.strftime('%Y%m%d%H%M%S')
f_name = 'For_paper'
if not os.path.isdir(f_name):
    os.mkdir(f_name)
os.chdir(f_name)
homedir = os.getcwd()

data = func.getSeqProperties(data)
data = func.TMTlabelEfficiency(data)
data = func.GRAVYScoreCalculator(data)
#data = func.pICalculator(data)
data = func.extractPhospho(data)
dfs_temp = func.DivideMergeFile(data)
dfs_temp = {k: func.RemoveDupli(v) for k, v in dfs_temp.items()}
dfs = {}

dfs['TMT-TiO$_2$'] = func.MergeFilesNoRemove([dfs_temp['160319ko02-TMT-TiO2.raw'], dfs_temp['160319ko03-TMT-TiO2.raw'], dfs_temp['160319ko04-TMT-TiO2.raw']])
dfs['TiO$_2$-TMT'] = func.MergeFilesNoRemove([dfs_temp['160319ko05-TiO2-TMT.raw'], dfs_temp['160319ko06-TiO2-TMT.raw'], dfs_temp['160319ko07-TiO2-TMT.raw']])

fig = plt.figure()
ax1 = fig.add_subplot(111)
fig.patch.set_color('white')

plt.grid(False)
plt.grid(axis="y")
i = 0
keys = ['TMT-TiO$_2$', 'TiO$_2$-TMT']
for k in keys:
    v = dfs[k]
    li = v[u'RawdataFile'].drop_duplicates().values
    idnumbers = [len(v[v[u'RawdataFile'] == fn]) for fn in li]
    meanID = np.nanmean(idnumbers)
    sd = np.std(idnumbers)
    ax1.bar(float(i) + 0.2, meanID, yerr=sd, ecolor="black",
            width=0.6, color='black', edgecolor='black')
    i += 1

plt.subplots_adjust(bottom=0.15, left=0.2, right=0.8)
plt.ylabel('Number of phosphopeptides')
ax1.yaxis.set_label_coords(-0.25, 0.5)

plt.xlim(0,i)
plt.xticks(list(map(lambda x: float(x) + 0.5, range(len(keys)))),
           ['TMT-TiO$_2$', 'TiO$_2$-TMT'], ha='center')
#minor_locator = AutoMinorLocator(2)
#ax1.xaxis.set_minor_locator(minor_locator)
plt.savefig('p-peptidesID-TMT-TiO2.png')
plt.clf()

plt.rcParams['xtick.major.size'] = 4
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['xtick.direction'] = 'inout'

fig = plt.figure()
ax1 = fig.add_subplot(111)
fig.patch.set_color('white')

plt.grid(False)

i = 0
cols = ['0.4', 'white']
ecols = ['0.4', 'black']
style = ['stepfilled', 'step']
for k in keys:
    v = dfs[k]
    v = func.RemoveDupli(v)
    ax1.hist(v[u'GRAVY_score'].values, color=cols[i], edgecolor=ecols[i],
             histtype=style[i], range=[-4, 4], bins=40, label=k)
    i += 1
plt.legend(bbox_to_anchor=[1, 0.9], loc='upper right', frameon=False, fontsize=20)
ax1.set_xticks([-4, -2, 0, 2, 4])
plt.xlabel('GRAVY score')
plt.ylabel('Number of phosphopeptides')
plt.subplots_adjust(bottom=0.2, right=0.8, left=0.2)
plt.savefig("GRAVY_hist.png")
plt.clf()

radar_data = pd.DataFrame()

for k in keys:
    v = func.RemoveDupli(dfs[k])
    # Get average values and arrange these valeus for radar chart(0-4)
    meanAcid = v['number_of_acidic_residue'].mean() * 4 / 8
    meanBase = v['number_of_basic_residue'].mean()
    meanGRAVY = (v['GRAVY_score'].mean() + 1.8) * 10 * 4 / 8 
    meanLength = (v['peptide_length'].mean() - 6) * 4 / 24
    meanPhospho = v['number_of_phospho'].mean() * 4 / 2
    radar_data.loc['No. of D, E\n(0, 8)', k] = meanAcid
    radar_data.loc['No. of H, K, R\n(0, 4)', k] = meanBase
    radar_data.loc['GRAVY\n(-1.8, -1.0)', k] = meanGRAVY
    radar_data.loc['Length\n(6, 30)', k] = meanLength
    radar_data.loc['No. of phosphate\ngroups\n(0, 2)', k] = meanPhospho

import pylab as pl


pl.rcParams['xtick.major.size'] = 6 
pl.rcParams['xtick.major.width'] = 1
pl.rcParams['ytick.major.size'] = 6 
pl.rcParams['ytick.major.width'] = 1

class Radar(object):

    def __init__(self, fig, titles, labels, rect=None):
        if rect is None:
            rect = [0.05, 0.05, 0.95, 0.95]

        self.n = len(titles)
        self.angles = np.arange(90, 90+360, 360.0/self.n)
        self.axes = [fig.add_subplot(111, projection="polar", label="axes%d" % i) 
                         for i in range(self.n)]

        for i in range(self.n):
            self.ax = self.axes[i]
            self.ax.set_thetagrids(self.angles, labels=titles, fontsize=28, va='center', ha='center', frac=1.2)
            self.ax.xaxis.grid(ls='-', lw=1.5, color='#4e454a')
            self.ax.set_yticks([1,2,3,4])
            self.ax.yaxis.grid(lw=1.5, color='#4e454a', alpha=0.5)

        for ax in self.axes[1:]:
            ax.patch.set_visible(False)
            ax.grid("off")
            ax.xaxis.set_visible(False)

        for ax, angle, label in zip(self.axes, self.angles, labels):
            ax.set_rgrids(range(1, 6), angle=angle, labels=label)
            ax.spines["polar"].set_visible(False)
            ax.set_ylim(0, 4)

    def plot(self, values, *args, **kw):
        angle = np.deg2rad(np.r_[self.angles, self.angles[0]])
        values = np.r_[values, values[0]]
        self.ax.plot(angle, values, *args, **kw)

titles = radar_data.index

labels = [
    ['','','','',''],
    ['','','','',''],
    ['','','','',''],
    ['','','','',''],
    ['','','','','']
]

fig = pl.figure(figsize=(40/2.54, 40/2.54))
radar_data = radar_data * 5 / 4

radar = Radar(fig, titles, labels)
radar.plot(radar_data.loc[:,'TMT-TiO$_2$'].values, "-", lw=3, color="#4878CF", label='TMT-TiO$_2$')
radar.plot(radar_data.loc[:,'TiO$_2$-TMT'].values, "-", lw=3, color="#6ACC65", label='TiO$_2$-TMT')

radar.ax.legend(bbox_to_anchor=(1.25, 0.25), loc='center', frameon=False, fontsize=32)
radar.ax.xaxis.grid(ls='-', lw=1.0)
radar.ax.set_yticks([0,5])
radar.ax.yaxis.grid(lw=2.5)
pl.subplots_adjust(top=0.7, bottom=0.3, left=0.3, right=0.7)
pl.savefig('radar_test.png', transparent='true')


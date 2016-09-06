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
plt.rcParams['font.size'] = 24
plt.rcParams['xtick.major.pad'] = 16
plt.rcParams['xtick.major.size'] = 0 
plt.rcParams['xtick.major.width'] = 0
plt.rcParams['xtick.minor.size'] = 0 
plt.rcParams['xtick.minor.width'] = 0
plt.rcParams['xtick.minor.pad'] = 60
plt.rcParams['figure.figsize'] = [30/2.54, 30/2.54]


filename = '160312ko-tipTMT_full.txt'
data = func.ReadFile(filename)

time = dt.now()
f_time = time.strftime('%Y%m%d%H%M%S')
f_name = 'For_paper'
if not os.path.isdir(f_name):
    os.mkdir(f_name)
os.chdir(f_name)
os.chdir('test')
homedir = os.getcwd()

data = func.getSeqProperties(data)
data = func.TMTlabelEfficiency(data)
data = func.GRAVYScoreCalculator(data)
#data = func.pICalculator(data)
data = func.extractPhospho(data)
dfs_temp = func.DivideMergeFile(data)
dfs_temp = {k: func.RemoveDupli(v) for k, v in dfs_temp.items()}
dfs = {}

dfs['pH 4.5'] = func.MergeFilesNoRemove([dfs_temp['160312ko02-SepPak_pH45.raw'], dfs_temp['160312ko08-SepPak_pH45.raw'], dfs_temp['160312ko14-SepPak_pH45.raw']])
dfs['pH 6.5'] = func.MergeFilesNoRemove([dfs_temp['160312ko03-SepPak_pH65.raw'], dfs_temp['160312ko09-SepPak_pH65.raw'], dfs_temp['160312ko15-SepPak_pH65.raw']])
dfs['TFA'] = func.MergeFilesNoRemove([dfs_temp['160312ko04-RP-C18_TFA.raw'], dfs_temp['160312ko10-RP-C18_TFA.raw'], dfs_temp['160312ko16-RP-C18_TFA.raw']])
dfs['HFBA'] = func.MergeFilesNoRemove([dfs_temp['160312ko05-RP-C18_HFBA.raw'], dfs_temp['160312ko11-RP-C18_HFBA.raw'], dfs_temp['160312ko17-RP-C18_HFBA.raw']])
dfs['Solution'] = func.MergeFilesNoRemove([dfs_temp['160312ko06-ref_solution.raw'], dfs_temp['160312ko12-ref_solution.raw'], dfs_temp['160312ko18-ref_solution.raw']])

fig = plt.figure()
ax1 = fig.add_subplot(111)
fig.patch.set_color('white')

plt.grid(False)
plt.grid(axis="y")
i = 0
#keys = ['pH 4.5', 'pH 6.5', 'TFA', 'HFBA', 'Solution']
keys = ['pH 4.5', 'pH 6.5']
for k in keys:
    v = dfs[k]
    li = v[u'RawdataFile'].drop_duplicates().values
    idnumbers = [len(v[v[u'RawdataFile'] == fn]) for fn in li]
    meanID = np.nanmean(idnumbers)
    sd = np.std(idnumbers)
    ax1.bar(float(i) + 0.2, meanID, yerr=sd, error_kw={'ecolor':'black', 'linewidth':2},
            width=0.6, color='black', edgecolor='black')
    i += 1

plt.subplots_adjust(bottom=0.15, left=0.15, right=0.85)
plt.ylabel('Number of phosphopeptides')
plt.ylim(0, 1000)
plt.xlim(0,i)
plt.xticks(map(lambda x: float(x) + 0.5, range(len(keys))),
           ['pH 4.5', 'pH 6.5'], ha='center')
#minor_locator = AutoMinorLocator(2)
#ax1.xaxis.set_minor_locator(minor_locator)
#ax1.set_xticks([1.0, 2.0, 3.0, 4.0], minor=True)
#ax1.set_xticklabels(['SepPak tC18', '', 'RP-C18, pH 6.5', ''], minor=True)
plt.savefig('p-peptidesID-2.png')
plt.clf()

ratioData = pd.read_csv('160312ko-ratio_full.csv', header=0, index_col=0)
ratioData['Average'] = [np.nanmean(ratioData.loc[i, :]) for i in ratioData.index]
ratioData['SD'] = [np.std(ratioData.loc[i, :]) for i in ratioData.index]

fig = plt.figure()
ax1 = fig.add_subplot(111)
fig.patch.set_color('white')
ratio = ratioData.loc[:, ['n=1', 'n=2', 'n=3']].values
ax1.boxplot(ratio.T, meanline=True, positions=[0.5, 1.5],
            showbox=False, showfliers=False, showcaps=False, showmeans=True,
            whis=0, medianprops={'alpha': 0},
            meanprops={'linestyle': '-', 'linewidth': 3.0, 'color': 'black'})
b0 = beeswarm(ratio, positions=[0.5, 1.5], s=300,
         ylim=(65,100), ax=ax1, labels=ratioData.index)
#ax1.set_xticks([1.0, 2.0, 3.0, 4.0], minor=True)
#ax1.set_xticklabels(['SepPak tC18', '', 'RP-C18, pH 6.5', ''], minor=True)
plt.grid(False)
plt.grid(axis="y")
plt.yticks(map(lambda x: x * 5, range(22)))
plt.subplots_adjust(bottom=0.15, left=0.15, right=0.85)
plt.xlim(0,2)
plt.ylim(65,100)
plt.ylabel('Labeling rate (%)')
plt.savefig('Labeling_ratio-2.png')
plt.clf()

df = func.ReadFile('160313ko-tipTMT-labeling_rate.txt')
df = func.getSeqProperties(df)
df = func.TMTlabelEfficiency(df)
df['TotPepArea'] = df['PepArea'] * (df['number_of_k'] + 1)
df['TMTPepArea'] = df['PepArea'] * (df['number_of_TMT'])
df['number_of_RP'] = df['number_of_k'] + 1
dfs_temp = func.DivideMergeFile(df)
dfs_temp = {k: func.RemoveDupli(v) for k, v in dfs_temp.items()}
ratio = {}
for k, v in dfs_temp.items():
    #ratio[k] = v['TMTPepArea'].sum() / v['TotPepArea'].sum() * 100
    ratio[k] = float(v['number_of_TMT'].sum()) / float(v['number_of_RP'].sum()) * 100

ratioData = pd.DataFrame(index=['pH 4.5', 'pH 6.5', 'TFA', 'HFBA', 'Solution'],
                         columns=['n=1', 'n=2', 'n=3'])



ratioData.loc['pH 4.5',:] = [ratio['160312ko02-SepPak_pH45.raw'], ratio['160312ko08-SepPak_pH45.raw'], ratio['160312ko14-SepPak_pH45.raw']]
ratioData.loc['pH 6.5',:] = [ratio['160312ko03-SepPak_pH65.raw'], ratio['160312ko09-SepPak_pH65.raw'], ratio['160312ko15-SepPak_pH65.raw']]
ratioData.loc['TFA',:] = [ratio['160312ko04-RP-C18_TFA.raw'], ratio['160312ko10-RP-C18_TFA.raw'], ratio['160312ko16-RP-C18_TFA.raw']]
ratioData.loc['HFBA',:] = [ratio['160312ko05-RP-C18_HFBA.raw'], ratio['160312ko11-RP-C18_HFBA.raw'], ratio['160312ko17-RP-C18_HFBA.raw']]
ratioData.loc['Solution',:] = [ratio['160312ko06-ref_solution.raw'], ratio['160312ko12-ref_solution.raw'], ratio['160312ko18-ref_solution.raw']]




ratioData['Average'] = [np.nanmean(ratioData.loc[i, :]) for i in ratioData.index]
ratioData['SD'] = [np.std(ratioData.loc[i, :]) for i in ratioData.index]

fig = plt.figure()
ax1 = fig.add_subplot(111)
fig.patch.set_color('white')
ratio = ratioData.loc[:, ['n=1', 'n=2', 'n=3']].values
ax1.boxplot(ratio.T, meanline=True, positions=[0.5, 1.5, 2.5, 3.5, 4.5],
            showbox=False, showfliers=False, showcaps=False, showmeans=True,
            whis=0, medianprops={'alpha': 0},
            meanprops={'linestyle': '-', 'linewidth': 3.0, 'color': 'black'})
b0 = beeswarm(ratio, positions=[0.5, 1.5, 2.5, 3.5, 4.5], s=300,
         ylim=(85,100), ax=ax1, labels=ratioData.index)
ax1.set_xticks([1.0, 2.0, 3.0, 4.0], minor=True)
ax1.set_xticklabels(['SepPak tC18', '', 'RP-C18, pH 6.5', ''], minor=True)
plt.grid(False)
plt.grid(axis="y")
plt.yticks(map(lambda x: float(x) * 2.5, range(44)))
plt.subplots_adjust(bottom=0.15, left=0.15, right=0.85)
plt.xlim(0,5)
plt.ylim(85,100)
plt.ylabel('Labeling rate (%)')
plt.savefig('Labeling_ratio_SP.png')
plt.clf()

multiphosData = pd.read_csv('160312ko-multiphospho.csv', header=0, index_col=0)

fig = plt.figure()
ax1 = fig.add_subplot(111)
fig.patch.set_color('white')
multiphosData['multiID'] = multiphosData[['2','3','4']].sum(axis=1)
multiphosData['monoratio'] =  multiphosData['1'].astype(float) / multiphosData[['1','2','3','4']].astype(float).sum(axis=1) * 100
multiphosData['multiratio'] = multiphosData['multiID'].astype(float) / multiphosData[['1','2','3','4']].astype(float).sum(axis=1) * 100
ax1.bar([0.2,1.2,2.2,3.2,4.2], multiphosData['monoratio'], width=0.6, color='white')
ax1.bar([0.2,1.2,2.2,3.2,4.2], multiphosData['multiratio'], bottom=multiphosData['monoratio'], width= 0.6,
        color='black')
ax1.set_xticks([0.5, 1.5, 2.5, 3.5, 4.5])
ax1.set_xticklabels(multiphosData.index)
plt.grid(False)
plt.grid(axis="y")
plt.subplots_adjust(bottom=0.15, left=0.15, right=0.85)
plt.xlim(0,5)
plt.ylim(60,100)
plt.ylabel('rate (%)')
plt.savefig('multiphos_ratio.png')
plt.clf()

radar_data = pd.DataFrame()
keys = ['HFBA', 'TFA']

for k in keys:
    v = func.RemoveDupli(dfs[k])
    # Get average values and arrange these valeus for radar chart(0-8)
    meanAcid = v['number_of_acidic_residue'].mean()
    meanBase = v['number_of_basic_residue'].mean() * 2
    meanGRAVY = (v['GRAVY_score'].mean() + 1.5) * 10
    meanLength = (v['peptide_length'].mean() - 6) / 3 
    meanPhospho = v['number_of_phospho'].mean() * 4
    radar_data.loc['No. of D, E\n(0, 8)', k] = meanAcid
    radar_data.loc['No. of H, K, R\n(0, 4)', k] = meanBase
    radar_data.loc['GRAVY\n(-1.5, -0.7)', k] = meanGRAVY
    radar_data.loc['length\n(6, 24)', k] = meanLength
    radar_data.loc['No. of phosphate groups\n(0, 2)', k] = meanPhospho

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
plt.savefig("GRAVY_hist_HFBA.png")
plt.clf()


import pylab as pl


pl.rcParams['xtick.major.size'] = 6 
pl.rcParams['xtick.major.width'] = 1
pl.rcParams['ytick.major.size'] = 6 
pl.rcParams['ytick.major.width'] = 1

class Radar(object):

    def __init__(self, fig, titles, labels, rect=None):
        if rect is None:
            rect = [0.05, 0.05, 0.95, 0.95]
        fig.patch.set_color('white')

        self.n = len(titles)
        self.angles = np.arange(90, 90+360, 360.0/self.n)
        self.axes = [fig.add_subplot(111, projection="polar", label="axes%d" % i) 
                         for i in range(self.n)]

        for i in range(self.n):
            self.ax = self.axes[i]
            self.ax.set_thetagrids(self.angles, labels=titles, fontsize=20, va='center', ha='center')
            self.ax.xaxis.grid(ls='-', lw=1.5, color='0.6')
            self.ax.set_yticks([0,1,2,3,4])
            self.ax.yaxis.grid(lw=1.0, alpha=0.5)

        for ax in self.axes[1:]:
            ax.patch.set_visible(False)
            ax.grid("off")
            ax.xaxis.set_visible(False)

        for ax, angle, label in zip(self.axes, self.angles, labels):
            ax.set_rgrids(range(1, 6), angle=angle, labels=label)
            ax.spines["polar"].set_visible(False)
            ax.set_ylim(0, 5)

    def plot(self, values, *args, **kw):
        angle = np.deg2rad(np.r_[self.angles, self.angles[0]])
        values = np.r_[values, values[0]]
        self.ax.plot(angle, values, *args, **kw)

radar_data = radar_data / 8 * 5
titles = radar_data.index

labels = [
    ['','','','',''],
    ['','','','',''],
    ['','','','',''],
    ['','','','',''],
    ['','','','','']
]

fig = pl.figure(figsize=(30/2.54, 30/2.54))

radar = Radar(fig, titles, labels)
radar.plot(radar_data.loc[:,'TFA'].values, "-", lw=2, color="black", label=radar_data.columns[0])
radar.plot(radar_data.loc[:,'HFBA'].values, "--", lw=2, color="black", label=radar_data.columns[1])

radar.ax.legend(bbox_to_anchor=(1.05, 0.2), loc='center', frameon=False, fontsize=24)
radar.ax.xaxis.grid(ls='-', lw=1.0)
radar.ax.set_yticks([0,5])
radar.ax.yaxis.grid(lw=1.5)
pl.subplots_adjust(top=0.92, bottom=0.02, left=0.12, right=0.82)
pl.show()



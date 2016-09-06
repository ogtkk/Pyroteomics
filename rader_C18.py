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
plt.rcParams['xtick.major.pad'] = 16
plt.rcParams['xtick.major.size'] = 0 
plt.rcParams['xtick.major.width'] = 0
plt.rcParams['xtick.minor.size'] = 0 
plt.rcParams['xtick.minor.width'] = 0
plt.rcParams['xtick.minor.pad'] = 60
plt.rcParams['figure.figsize'] = [30/2.54, 30/2.54]
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20




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

keys = ['pH 6.5', 'Solution']
radar_data = pd.DataFrame()

for k in keys:
    v = func.RemoveDupli(dfs[k])
    # Get average values and arrange these valeus for radar chart(0-4)
    meanAcid = v['number_of_acidic_residue'].mean() * 4 / 10
    meanBase = v['number_of_basic_residue'].mean() * 4 / 3
    meanGRAVY = (v['GRAVY_score'].mean() + 1.6) * 10 * 4 / 10 
    meanLength = (v['peptide_length'].mean() - 6) * 4 / 24
    meanPhospho = v['number_of_phospho'].mean() * 4 / 2
    radar_data.loc['No. of D, E\n(0, 10)', k] = meanAcid
    radar_data.loc['No. of H, K, R\n(0, 3)', k] = meanBase
    radar_data.loc['GRAVY\n(-1.6, -0.6)', k] = meanGRAVY
    radar_data.loc['length\n(6, 30)', k] = meanLength
    radar_data.loc['No. of phosphate groups\n(0, 2)', k] = meanPhospho

print radar_data
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
            self.ax.set_yticks([1,2,3,4])
            self.ax.yaxis.grid(lw=1.0, alpha=0.5)

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

fig = pl.figure(figsize=(30/2.54, 30/2.54))
radar_data = radar_data * 5 / 4

radar = Radar(fig, titles, labels)
radar.plot(radar_data.loc[:,'pH 6.5'].values, "-", lw=2, color="black", label=radar_data.columns[0])
radar.plot(radar_data.loc[:,'Solution'].values, "--", lw=2, color="black", label=radar_data.columns[1])

radar.ax.legend(bbox_to_anchor=(1.05, 0.2), loc='center', frameon=False, fontsize=24)
radar.ax.xaxis.grid(ls='-', lw=1.0)
radar.ax.set_yticks([0,5])
radar.ax.yaxis.grid(lw=1.5)
pl.subplots_adjust(top=0.92, bottom=0.02, left=0.12, right=0.82)
pl.show()


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

filename = '160409ko-TiO2-TMT.txt'
data = func.ReadFile(filename)

time = dt.now()
f_time = time.strftime('%Y%m%d%H%M%S')
f_name = 'For_paper'
if not os.path.isdir(f_name):
    os.mkdir(f_name)
os.chdir(f_name)
os.chdir('160409')
#os.chdir('MIX')

homedir = os.getcwd()

data = func.getSeqProperties(data)
data = func.TMTlabelEfficiency(data)
data = func.GRAVYScoreCalculator(data)
#data = func.pICalculator(data)
data = func.extractPhospho(data)
dfs_temp = func.DivideMergeFile(data)
dfs_temp = {k: func.RemoveDupli(v) for k, v in dfs_temp.items()}
dfs = {}

dfs['TMT-TiO$_2$'] = dfs_temp['160409ko02-TMT-TiO2-3plex.raw']
dfs['TiO$_2$-TMT'] = dfs_temp['160409ko03-TiO2-TMT-3plex.raw']
dfs['MIX'] = dfs_temp['160409ko04-MIX-6plex.raw']
data = dfs['TiO$_2$-TMT']
#data = pd.concat([dfs['TMT-TiO$_2$'], dfs['TiO$_2$-TMT']])
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#fig.patch.set_color('white')
#
#plt.grid(False)
#plt.grid(axis="y")
#i = 0
#keys = ['TMT-TiO$_2$', 'TiO$_2$-TMT']
#for k in keys:
#    v = dfs[k]
#    li = v[u'RawdataFile'].drop_duplicates().values
#    idnumbers = [len(v[v[u'RawdataFile'] == fn]) for fn in li]
#    meanID = np.nanmean(idnumbers)
#    sd = np.std(idnumbers)
#    ax1.bar(float(i) + 0.2, meanID, yerr=sd, ecolor="black",
#            width=0.6, color='black', edgecolor='black')
#    i += 1
#
#plt.subplots_adjust(bottom=0.15, left=0.15, right=0.85)
#plt.ylabel('Number of phosphopeptides')
#ax1.yaxis.set_label_coords(-0.125, 0.5)
#plt.xlim(0,i)
#plt.xticks(map(lambda x: float(x) + 0.5, range(len(keys))),
#           keys, ha='center')
##minor_locator = AutoMinorLocator(2)
##ax1.xaxis.set_minor_locator(minor_locator)
#plt.savefig('p-peptidesID-TMT-TiO2.png')
#plt.clf()
#
#plt.rcParams['xtick.major.size'] = 4
#plt.rcParams['xtick.major.width'] = 1
#plt.rcParams['xtick.direction'] = 'inout'
#
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#fig.patch.set_color('white')
#
#plt.grid(False)
#
#i = 0
#cols = ['0.4', 'white']
#ecols = ['0.4', 'black']
#style = ['stepfilled', 'step']
#for k in keys:
#    v = dfs[k]
#    v = func.RemoveDupli(v)
#    ax1.hist(v[u'GRAVY_score'].values, color=cols[i], edgecolor=ecols[i],
#             histtype=style[i], range=[-4, 4], bins=40, label=k)
#    i += 1
#plt.legend(bbox_to_anchor=[1, 0.9], loc='upper right', frameon=False, fontsize=20)
#ax1.set_xticks([-4, -2, 0, 2, 4])
#plt.xlabel('GRAVY score')
#plt.ylabel('Number of phosphopeptides')
#plt.subplots_adjust(bottom=0.2, right=0.8, left=0.2)
#plt.savefig("GRAVY_hist.png")
#plt.clf()

reporter_values = (
    u'iTRAQ:126.127725', u'iTRAQ:127.124760',
    u'iTRAQ:127.131079', u'iTRAQ:128.128114',
    u'iTRAQ:128.134433', u'iTRAQ:129.131468',
    u'iTRAQ:129.137787', u'iTRAQ:130.134822',
    u'iTRAQ:130.141141', u'iTRAQ:131.138176')

reporters = {#'TMT-TiO$_2$': [u'iTRAQ:126.127725', u'iTRAQ:127.131079', u'iTRAQ:128.128114'],
             'TiO$_2$-TMT': [u'iTRAQ:129.137787', u'iTRAQ:130.134822', u'iTRAQ:131.138176']}

df_temp = pd.DataFrame()
df_temp = data.iloc[:, data.columns.isin(reporter_values)]
df_temp = df_temp.join(data[u'SeqModModDetail'])
df_temp = df_temp.groupby([u'SeqModModDetail']).mean()
df_temp = df_temp.replace(0, np.nan)
df_temp = df_temp.dropna(subset=[#u'iTRAQ:126.127725', u'iTRAQ:127.131079', u'iTRAQ:128.128114',
                                 u'iTRAQ:129.137787', u'iTRAQ:130.134822', u'iTRAQ:131.138176'])

#for title, reporter in reporters.items():
#    df_Int = pd.DataFrame()
#    index = []
#    average = []
#    for cn in reporter:
#        df_Int[cn] = df_temp[cn]
#    for i in df_Int.index:
#        index.append(i)
#        average.append(np.nanmean(df_Int.loc[i]))
#        df_temp[u'AVERAGE: ' + title] = pd.Series(average, index=index)
#
#df_temp['log2'] = np.log2(df_temp[u'AVERAGE: ' + 'TMT-TiO$_2$'] / df_temp[u'AVERAGE: ' + 'TiO$_2$-TMT'])
#median = np.median(df_temp['log2'])
#df_temp['SeqModModDetail'] = df_temp.index
#df = pd.merge(func.RemoveDupli(data), df_temp,
#              on='SeqModModDetail', how='left')
##df.to_csv('text.csv')
#fig = plt.figure()
#fig.patch.set_color('white')
#ax1 = fig.add_subplot(111)
#ax1.hist(df['log2'], bins=40, range=[-8,2], color='0.4', edgecolor='0.4')
#ax1.axvline(median, color='black', linestyle='--', linewidth=1.5)
#plt.show()
#plt.clf()
#
#fig = plt.figure()
#fig.patch.set_color('white')
#ax1 = fig.add_subplot(111)
#ax1.hist(df['log2'], bins=40, range=[-8,2], color='0.4', edgecolor='0.4')
#ax1.axvline(median, color='black', linestyle='--', linewidth=1.5)
#plt.show()
#plt.clf()
#
#plt.rcParams['xtick.major.size'] = 6 
#plt.rcParams['xtick.major.width'] = 1.0
#
#bins = map(lambda x: float(x) * 10, range(11))
#df['RetTime_bin'] = pd.cut(df.RetTime, bins)
#fig = plt.figure()
#fig.patch.set_color('white')
#ax1 = fig.add_subplot(111)
#whisprop = dict(linestyle='--', linewidth=1.25)
#capprop = dict(linestyle='-', linewidth=1.5)
#boxprop = dict(linestyle='-', linewidth=1.25)
#medprop = dict(linestyle='-', linewidth=2.0)
#flprop = dict(markersize=8, color='black')
#bp = df.boxplot(column='log2', by='number_of_acidic_residue', grid=False, whis=[5,95],
#           ax=ax1,
#           whiskerprops=whisprop, capprops=capprop, boxprops=boxprop,
#           medianprops=medprop, flierprops=flprop, return_type='dict')
#[[item.set_color('black') for item in bp[key]['medians']] for key in bp.keys()]
#[[item.set_color('black') for item in bp[key]['boxes']] for key in bp.keys()]
#[[item.set_color('black') for item in bp[key]['whiskers']] for key in bp.keys()]
#[[item.set_color('black') for item in bp[key]['caps']] for key in bp.keys()]
#
#ax1.set_xticks(range(len(bins)))
#ax1.set_xticklabels(range(26))
#ax1.set_xlim(0,26)
#fig.suptitle('')
#ax1.set_title('')
#ax1.set_xlabel('number of D, E', labelpad=20)
#ax1.set_ylabel('Relative reporter ion intensity\n($Log_{2}$ TMT-TiO2 / TiO2-TMT)')
#plt.grid(False)
#plt.grid(axis='y')
#plt.subplots_adjust(bottom=0.2, left=0.2, right=0.8)
#plt.savefig('boxplot_xDEyratio.png')
#plt.clf()
#
#plt.rcParams['xtick.major.size'] = 6 
#plt.rcParams['xtick.major.width'] = 1.0
#
#bins = map(lambda x: float(x) * 10, range(11))
#df['RetTime_bin'] = pd.cut(df.RetTime, bins)
#fig = plt.figure()
#fig.patch.set_color('white')
#ax1 = fig.add_subplot(111)
#whisprop = dict(linestyle='--', linewidth=1.25)
#capprop = dict(linestyle='-', linewidth=1.5)
#boxprop = dict(linestyle='-', linewidth=1.25)
#medprop = dict(linestyle='-', linewidth=2.0)
#flprop = dict(markersize=8, color='black')
#bp = df.boxplot(column='log2', by='RetTime_bin', grid=False, whis=[5,95],
#           ax=ax1, positions=map(lambda x: float(x) + 0.5, range(10)),
#           whiskerprops=whisprop, capprops=capprop, boxprops=boxprop,
#           medianprops=medprop, flierprops=flprop, return_type='dict')
#[[item.set_color('black') for item in bp[key]['medians']] for key in bp.keys()]
#[[item.set_color('black') for item in bp[key]['boxes']] for key in bp.keys()]
#[[item.set_color('black') for item in bp[key]['whiskers']] for key in bp.keys()]
#[[item.set_color('black') for item in bp[key]['caps']] for key in bp.keys()]
#
#ax1.set_xticks(range(len(bins)))
#ax1.set_xticklabels([0, 10,20,30,40,50,60,70,80,90,100])
#ax1.set_xlim(2,9)
#fig.suptitle('')
#ax1.set_title('')
#ax1.set_xlabel('Retention time (min)', labelpad=20)
#ax1.set_ylabel('Relative reporter ion intensity\n($Log_{2}$ RP-C18 / Solution)')
#plt.grid(False)
#plt.grid(axis='y')
#plt.subplots_adjust(bottom=0.2, left=0.2, right=0.8)
#plt.savefig('boxplot_xRetTyInt.png')
#plt.clf()


df_t = pd.DataFrame()
if bool(reporter_values) is True:
    for cn in reporter_values:
        series = data.ix[:, data.columns.map(
            lambda x: x.startswith(cn))]
        df_t = pd.concat([df_t, series], axis=1)
df_t = df_t.join(data[u'SeqModModDetail'])
df_t = df_t.replace(0, np.nan)
for title, reporter in reporters.items():
    df_temp = df_t.dropna(subset=reporter)
    df_temp = df_temp.groupby([u'SeqModModDetail']).mean()
    df_lnInt = pd.DataFrame()
    df_Int = df_temp.iloc[:, df_temp.columns.isin(reporter)]
    for cn in reporter:
        df_lnInt[cn] = np.log10(df_Int[cn])
    Rep1 = df_lnInt[reporter[0]]
    Rep2 = df_lnInt[reporter[1]]
    Rep3 = df_lnInt[reporter[2]]
    fig_title = title + ' (n=%d)' % len(df_lnInt)
    fig = plt.figure()
    fig.text(0.5, 0.85, fig_title, size=24, ha='center')
    ax1 = fig.add_subplot(334)
    ax2 = fig.add_subplot(338)
    ax3 = fig.add_subplot(337)
    ax1.scatter(Rep1, Rep2, color="gray", s=8, alpha=0.6)
    ax1.set_xlim(0, 10)
    ax1.set_ylim(0, 10)
    ax1.set_xticklabels('')
    ax1.set_ylabel('Rep. 2\nlog10 (Intensity)', size=12)
    ax2.scatter(Rep2, Rep3, color="gray", s=8, alpha=0.6)
    ax2.set_xlim(0, 10)
    ax2.set_ylim(0, 10)
    ax2.set_xlabel('log10 (Intensity)\nRep. 2', size=12)
    ax2.set_yticklabels('')
    ax3.scatter(Rep1, Rep3, color="gray", s=8, alpha=0.6)
    ax3.set_xlim(0, 10)
    ax3.set_ylim(0, 10)
    ax3.set_xlabel('log10 (Intensity)\nRep. 1', size=20)
    ax3.set_ylabel('Rep. 3\nlog10 (Intensity)', size=20)
    r1, p_value1 = st.pearsonr(Rep1, Rep2)
    r2, p_value2 = st.pearsonr(Rep2, Rep3)
    r3, p_value3 = st.pearsonr(Rep1, Rep3)
    ax1.text(5, 1.5, "$R^2=%s$" % str(r1**2)[0:6], size=13)
    ax2.text(5, 1.5, "$R^2=%s$" % str(r2**2)[0:6], size=13)
    ax3.text(5, 1.5, "$R^2=%s$" % str(r3**2)[0:6], size=13)
    plt.subplots_adjust(bottom=0.2, wspace=0.2, hspace=0.15)
    ftitle_rep = title + "_Reproducibility.png"
    fig.savefig(ftitle_rep)
    fig.clf()
    def set_color(cv):
        if cv <= 20:
            return "blue"
        else:
            return "gray"
    df_result = pd.DataFrame()
    lnAverage = [np.nanmean(df_lnInt.loc[i]) for i in df_lnInt.index]
    df_result['avr'] = [np.nanmean(df_Int.loc[i]) for i in df_Int.index]
    df_result['stdev'] = [np.std(df_Int.loc[i]) for i in df_Int.index]
    df_result['cv'] = df_result['stdev'] / df_result['avr'] * 100
    color_list = map(set_color, df_result.cv)
    plt.scatter(df_result['cv'], lnAverage, c=color_list, s=8, alpha=0.6)
    plt.xlim(0, 100)
    plt.ylim(2, 7)
    plt.axvline(x=20, ls='--', lw=0.8, color='black')
    plt.text(20, 6.2, "C.V. $\leqq$ 20\n(n = %d)" % len(df_result[df_result['cv'] <= 20]), size=20)
    plt.xlabel('C.V. (%)', size=20)
    plt.ylabel('Average log10(Intensity)', size=20)
    df_result.to_csv('csv.csv')
    plt.subplots_adjust(bottom=0.2, left=0.1)
    plt.title(fig_title, size=24)
    ftitle_cv = title + "_treePlot.png"
    plt.savefig(ftitle_cv)
    plt.clf()

#radar_data = pd.DataFrame()
#
#for k in keys:
#    v = func.RemoveDupli(dfs[k])
#    # Get average values and arrange these valeus for radar chart(0-4)
#    meanAcid = v['number_of_acidic_residue'].mean() * 4 / 10
#    meanBase = v['number_of_basic_residue'].mean() * 4 / 3
#    meanGRAVY = (v['GRAVY_score'].mean() + 1.6) * 10 * 4 / 10 
#    meanLength = (v['peptide_length'].mean() - 6) * 4 / 24
#    meanPhospho = v['number_of_phospho'].mean() * 4 / 2
#    radar_data.loc['No. of D, E\n(0, 10)', k] = meanAcid
#    radar_data.loc['No. of H, K, R\n(0, 3)', k] = meanBase
#    radar_data.loc['GRAVY\n(-1.6, -0.6)', k] = meanGRAVY
#    radar_data.loc['length\n(6, 30)', k] = meanLength
#    radar_data.loc['No. of phosphate groups\n(0, 2)', k] = meanPhospho
#
#print radar_data
#import pylab as pl
#
#
#pl.rcParams['xtick.major.size'] = 6 
#pl.rcParams['xtick.major.width'] = 1
#pl.rcParams['ytick.major.size'] = 6 
#pl.rcParams['ytick.major.width'] = 1
#
#class Radar(object):
#
#    def __init__(self, fig, titles, labels, rect=None):
#        if rect is None:
#            rect = [0.05, 0.05, 0.95, 0.95]
#        fig.patch.set_color('white')
#
#        self.n = len(titles)
#        self.angles = np.arange(90, 90+360, 360.0/self.n)
#        self.axes = [fig.add_subplot(111, projection="polar", label="axes%d" % i) 
#                         for i in range(self.n)]
#
#        for i in range(self.n):
#            self.ax = self.axes[i]
#            self.ax.set_thetagrids(self.angles, labels=titles, fontsize=20, va='center', ha='center')
#            self.ax.xaxis.grid(ls='-', lw=1.5, color='0.6')
#            self.ax.set_yticks([1,2,3,4])
#            self.ax.yaxis.grid(lw=1.0, alpha=0.5)
#
#        for ax in self.axes[1:]:
#            ax.patch.set_visible(False)
#            ax.grid("off")
#            ax.xaxis.set_visible(False)
#
#        for ax, angle, label in zip(self.axes, self.angles, labels):
#            ax.set_rgrids(range(1, 6), angle=angle, labels=label)
#            ax.spines["polar"].set_visible(False)
#            ax.set_ylim(0, 4)
#
#    def plot(self, values, *args, **kw):
#        angle = np.deg2rad(np.r_[self.angles, self.angles[0]])
#        values = np.r_[values, values[0]]
#        self.ax.plot(angle, values, *args, **kw)
#
#titles = radar_data.index
#
#labels = [
#    ['','','','',''],
#    ['','','','',''],
#    ['','','','',''],
#    ['','','','',''],
#    ['','','','','']
#]
#
#fig = pl.figure(figsize=(30/2.54, 30/2.54))
#radar_data = radar_data * 5 / 4
#
#radar = Radar(fig, titles, labels)
#radar.plot(radar_data.loc[:,'TMT-TiO$_2$'].values, "-", lw=2, color="black", label=radar_data.columns[0])
#radar.plot(radar_data.loc[:,'TiO$_2$-TMT'].values, "--", lw=2, color="black", label=radar_data.columns[1])
#
#radar.ax.legend(bbox_to_anchor=(1.05, 0.2), loc='center', frameon=False, fontsize=24)
#radar.ax.xaxis.grid(ls='-', lw=1.0)
#radar.ax.set_yticks([0,5])
#radar.ax.yaxis.grid(lw=1.5)
#pl.subplots_adjust(top=0.92, bottom=0.02, left=0.12, right=0.82)
#pl.show()
#
#
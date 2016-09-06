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

fnamepath = re.compile(u"D.+ManualInputFile\\\\|D.+ManualInputFile/")
fnameext = re.compile(u".raw")
filename_1 = '160403ko-tipTMT-full.txt'
filename_2 = '160403ko-tipTMT-full_2.txt'
data_1 = func.ReadFile(filename_1)
data_2 = func.ReadFile(filename_2)
data_HFBA = func.ReadFile('160319ko-solHFBA.txt')

data_HFBA[u'RawdataFile'] = data_HFBA[u'RawdataFile'].apply(
    lambda x: fnamepath.sub("", x))
data_HFBA[u'RawdataFile'] = data_HFBA[u'RawdataFile'].apply(
    lambda x: fnameext.sub("", x))
data_HFBA['RetTime'] = data_HFBA['RetTime'] + 120 * (data_HFBA['RawdataFile'].str[-1:].astype('int') - 1)

data = data_HFBA
#pd.concat([data_1, data_2], ignore_index=True)
data[u'RawdataFile'] = data[u'RawdataFile'].apply(
    lambda x: fnamepath.sub("", x))
data[u'RawdataFile'] = data[u'RawdataFile'].apply(
    lambda x: fnameext.sub("", x))
#data['RetTime'] = data['RetTime'] + 135 * (data['RawdataFile'].str[-1:].astype('int') - 1)
#data.loc[data['RawdataFile'].str[-1:].astype('int') != 1, 'RetTime']\
#    = data.loc[data['RawdataFile'].str[-1:].astype('int') != 1, 'RetTime'] - 15
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

fig = plt.figure()
fig.patch.set_color('white')
ax1 = fig.add_subplot(111)
ax1.hist(df['RetTime'], bins=240, range=[40,520], color='0.4', edgecolor='0.4')
ax1.set_xlabel('Retention time')
ax1.set_ylabel('Number of phosphopeptides')
plt.subplots_adjust(bottom=0.2, left=0.2, right=0.95)
plt.savefig('hist_rettime_HFBA.png')
plt.clf()

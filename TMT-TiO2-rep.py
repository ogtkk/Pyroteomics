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


filename = '160319ko-TMT-TiO2-fix.txt'
data = func.ReadFile(filename)

time = dt.now()
f_time = time.strftime('%Y%m%d%H%M%S')
f_name = 'For_paper'
if not os.path.isdir(f_name):
    os.mkdir(f_name)
os.chdir(f_name)
homedir = os.getcwd()

data = func.getSeqProperties(data)
#data = func.TMTlabelEfficiency(data)
#data = func.GRAVYScoreCalculator(data)
#data = func.pICalculator(data)
data = func.extractPhospho(data)
dfs_temp = func.DivideMergeFile(data)
dfs_temp = {k: func.RemoveDupli(v) for k, v in dfs_temp.items()}
df = pd.merge(dfs_temp.values(), how='inner')
df.to_csv('test_inner.csv')
#dfs['TMT-TiO$_2$'] = func.MergeFilesNoRemove([dfs_temp['160319ko02-TMT-TiO2.raw'], dfs_temp['160319ko03-TMT-TiO2.raw'], dfs_temp['160319ko04-TMT-TiO2.raw']])
#dfs['TiO$_2$-TMT'] = func.MergeFilesNoRemove([dfs_temp['160319ko05-TiO2-TMT.raw'], dfs_temp['160319ko06-TiO2-TMT.raw'], dfs_temp['160319ko07-TiO2-TMT.raw']])

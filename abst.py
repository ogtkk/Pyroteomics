import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st
import os
from collections import Counter
import re
import matplotlib as mpl

plt.rcParams['font.family'] = 'sans-serif' 
plt.rcParams['font.size'] = 32
plt.rcParams['xtick.major.pad'] = 20
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

data = pd.read_csv('abst_data_2.csv', index_col=0, header=0)
print data

fig, ax1 = plt.subplots()
fig.patch.set_color('white')
b0 = ax1.bar([0.2,1.2], data.iloc[0,:], width=0.4, color='white', edgecolor='black', hatch='//')
plt.ylim(3,5)
plt.ylabel(data.index[0], fontsize=40)
ax1.yaxis.set_label_coords(-0.18, 0.5)
plt.axvline(0.4, 0.42, 0.94, color='black', linestyle='--', linewidth='2.5')
plt.axvline(1.4, 0.89, 0.94, color='black', linestyle='--', linewidth='2.5')
plt.axhline(4.89, 0.182, 0.636, color='black', linestyle='--', linewidth='2.5')
plt.text(0.45, 4.3, '8.7-fold', fontsize=40)
ax2 = ax1.twinx()
b1 = ax2.bar([0.6,1.6], data.iloc[1,:], width=0.4, color='black')
plt.ylabel(data.index[1], fontsize=40)
ax2.yaxis.set_label_coords(1.24, 0.5)

plt.grid(False)
legend = plt.legend((b0, b1), list(data.index), bbox_to_anchor=(-0.2, 1.04, .8, .098), loc=3,
           mode="expand", borderaxespad=0., fontsize=36)
legend.get_frame().set_edgecolor('white')
legend.get_frame().set_alpha(0)
plt.subplots_adjust(bottom=0.1, top=0.80, left=0.2, right=0.8)
plt.xlim(0,2.2)
plt.xticks([0.6, 1.6], data.columns, ha='center', fontsize=40)
plt.show()
plt.clf()

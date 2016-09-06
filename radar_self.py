import Functions as func
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
from datetime import datetime as dt
import scipy.stats as st
import matplotlib as mpl



data = pd.read_csv('radar_chart.csv', header=0, index_col=0)

time = dt.now()
f_time = time.strftime('%Y%m%d%H%M%S')
f_name = 'For_paper'
if not os.path.isdir(f_name):
    os.mkdir(f_name)
os.chdir(f_name)
homedir = os.getcwd()


fig = plt.figure()
fig.patch.set_color('white')
ax = [fig.add_subplot(111, projection="polar", label="axes%d" % i) for i in range(5)])
angles = np.arange(90, 90+360, 360.0/5)

radar.plot(data.loc[:,'TMT-TiO2'].values, "-", lw=2, color="black", label=data.columns[0])
radar.plot(data.loc[:,'TiO2-TMT'].values, "--", lw=2, color="black", label=data.columns[1])
radar.ax.legend(bbox_to_anchor=(1.2, 0.2), loc='center')
pl.subplots_adjust(top=0.8, bottom=0.2, left=0.2, right=0.8)
pl.show()
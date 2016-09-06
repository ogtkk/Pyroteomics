import numpy as np
import pandas as pd
import pylab as pl

radar_data = pd.read_csv('../Analysis_csv/radar_chart.csv', header=0, index_col=0)
radar_data.columns = ['TMT-TiO$_2$', 'TiO$_2$-TMT']
radar_data = radar_data / 2
titles = radar_data.index





pl.rcParams['xtick.major.size'] = 6 
pl.rcParams['xtick.major.width'] = 1
pl.rcParams['ytick.major.size'] = 6 
pl.rcParams['ytick.major.width'] = 1
pl.rcParams['mathtext.default'] = 'regular'

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

radar = Radar(fig, titles, labels)
radar.plot(radar_data.loc[:,'TMT-TiO$_2$'].values, "-", lw=2, color="black", label=radar_data.columns[0])
radar.plot(radar_data.loc[:,'TiO$_2$-TMT'].values, "--", lw=2, color="black", label=radar_data.columns[1])

radar.ax.legend(bbox_to_anchor=(1.05, 0.2), loc='center', frameon=False, fontsize=24)
radar.ax.xaxis.grid(ls='-', lw=1.0)
radar.ax.set_yticks([0,5])
radar.ax.yaxis.grid(lw=1.5)
pl.subplots_adjust(top=0.92, bottom=0.02, left=0.12, right=0.82)
pl.show()

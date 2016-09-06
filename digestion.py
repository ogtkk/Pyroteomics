import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import Functions as func
import numpy as np
import matplotlib as mpl
from scipy import stats
import seaborn as sns
import matplotlib as mpl

plt.rcParams['font.size'] = 20
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['mathtext.default'] = 'regular'


data = pd.read_csv('digestion.txt', index_col=None, header=0, delimiter="\t")
data_1 = data[data.Dataname.str.endswith('_1')]
data_2 = data[data.Dataname.str.endswith('_2')]
data_TMT = data[data.Dataname.str.endswith('_TMT')]

i = 0
color = sns.color_palette("muted")
sns.set_style("whitegrid")
plt.grid(False)
plt.grid(axis="y")
K_0 = data_1['K_0']
K_1 = data_1['K_1']
K_2 = data_1['K_2']
K_3 = data_1['K_3']
peptID = K_0 + K_1 + K_2 + K_3
area0 = area0 / idn
area1 = area1 / idn
area2 = area2 / idn
area3 = area3 / idn
plt.bar(i + float(i)/2, area0, width=1, color=color[0],
        alpha=0.8, edgecolor=color[0])
plt.bar(i + float(i)/2, area1, bottom=area0, width=1,
        color=color[1], alpha=0.8, edgecolor=color[1])
plt.bar(i + float(i)/2, area2, bottom=area0+area1, width=1,
        color=color[2], alpha=0.8, edgecolor=color[2])
plt.bar(i + float(i)/2, area3, bottom=area0+area1+area2, width=1,
        color=color[3], alpha=0.8, edgecolor=color[2])
i += 1
plt.legend(["K*0", "K*1", "K*2", "K*3"],
           bbox_to_anchor=(1.05, 1), loc='upper left')
plt.subplots_adjust(bottom=0.3, left=0.15, right=0.65)
plt.xticks(map(lambda x: float(x) * 3/2 + 0.5, range(len(dfs.keys()))),
           sorted(dfs.keys()), ha='right', rotation=45)
plt.savefig("number_of_KID_bar.png")
plt.clf()

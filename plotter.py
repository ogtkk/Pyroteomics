import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import Functions as func
import numpy as np
import matplotlib as mpl
from scipy import stats
from beeswarm import *

plt.rcParams['font.size'] = 20
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['mathtext.default'] = 'regular'


def plot_all(dfs):
    TMTefficiency_barplot(dfs)
    chargeState_barplot(dfs)
    numberOfTMT_barplot(dfs)
    numberOfK_barplot(dfs)
    numberOfKID_barplot(dfs)
    multiphospho_peakarea(dfs)
    acidicresidue_peakarea(dfs)
    basicresidue_peakarea(dfs)
    pI_peakarea(dfs)
    RetTime_peakarea(dfs)
    multiphospho_hist(dfs)
    acidicresidue_hist(dfs)
    basicresidue_hist(dfs)
    pI_hist(dfs)
    GRAVYscore_hist(dfs)
    RetTime_hist(dfs)
    peptide_length_hist(dfs)
    numberofDhist(dfs)
    numberofEhist(dfs)
    numberofKhist(dfs)
    numberofRhist(dfs)
    #peptScoreKdeHistK(dfs)
    peptScoreHist(dfs)
    #peptScoreKdeHistK1(dfs)
    #peptScoreKdeHistR1(dfs)
    #peptScoreKdeHistCharge(dfs)
    #peptHalfWidthHist(dfs)
    #peptHalfWidthHistK(dfs)
    pepAreaBoxPlot(dfs)


def TMTefficiency_barplot(dfs):
    i = 0
    color = sns.color_palette("muted")
    sns.set_style("whitegrid")
    plt.grid(False)
    plt.grid(axis="y")

    for k, v in sorted(dfs.items()):
        farea = v[v[u'efficiency'].str.contains('Full')][u'PepArea'].sum()
        parea = v[v[u'efficiency'].str.contains('Partial')][u'PepArea'].sum()
        narea = v[v[u'efficiency'].str.contains('None')][u'PepArea'].sum()
        plt.bar(i + float(i)*3/2, farea, width=1.5, color=color[0],
                alpha=0.8, edgecolor=color[0])
        plt.bar(i + float(i)*3/2, parea, bottom=farea,
                width=1.5, color=color[1], alpha=0.8, edgecolor=color[1])
        plt.bar(i + float(i)*3/2, narea, bottom=farea+parea, width=1.5,
                color=color[2], alpha=0.8, edgecolor=color[2])
        i += 1

    plt.legend(["Full", "Partial", "None"],
               bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, left=0.2, right=0.75)
    plt.ylabel('Peak area')
    plt.xticks(map(lambda x: float(x) * 5/2 + 1.0, range(len(dfs.keys()))),
               sorted(dfs.keys()), ha='right', rotation=45)
    plt.savefig("TMTefficiency.png")
    plt.clf()


def Guanidination_efficiency_barplot(dfs):
    i = 0
    color = sns.color_palette("muted")
    sns.set_style("whitegrid")
    plt.grid(False)
    plt.grid(axis="y")

    for k, v in sorted(dfs.items()):
        farea = v[v[u'Guanidination_efficiency']
                  .str.contains('Full')][u'PepArea'].sum()
        parea = v[v[u'Guanidination_efficiency']
                  .str.contains('Partial')][u'PepArea'].sum()
        narea = v[v[u'Guanidination_efficiency']
                  .str.contains('None')][u'PepArea'].sum()
        oarea = v[v[u'Guanidination_efficiency']
                  .str.contains('Overlabel')][u'PepArea'].sum()
        plt.bar(i + float(i)*3/2, farea, width=1.5, color=color[0],
                alpha=0.8, edgecolor=color[0])
        plt.bar(i + float(i)*3/2, parea, bottom=farea,
                width=1.5, color=color[1], alpha=0.8, edgecolor=color[1])
        plt.bar(i + float(i)*3/2, narea, bottom=farea+parea, width=1.5,
                color=color[2], alpha=0.8, edgecolor=color[2])
        plt.bar(i + float(i)*3/2, oarea, bottom=farea+parea+narea, width=1.5,
                color=color[3], alpha=0.8, edgecolor=color[3])

        i += 1

    plt.legend(["Full", "Partial", "None", "Overlabel"],
               bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.4, left=0.2, right=0.75)
    plt.ylabel('Peak area')
    plt.xticks(map(lambda x: float(x) * 5/2 + 1.0, range(len(dfs.keys()))),
               sorted(dfs.keys()), ha='right', rotation=45)
    plt.savefig("Guanidination_efficiency.png")
    plt.clf()


def meanPeptideIDwithSD(dfs):
    i = 0

    color = sns.color_palette("muted")
    sns.set_style("whitegrid")
    plt.grid(False)
    plt.grid(axis="y")

    for k, v in sorted(dfs.items()):
        li = v[u'RawdataFile'].drop_duplicates().values
        idnumbers = [len(v[v[u'RawdataFile'] == fn]) for fn in li]
        meanID = np.nanmean(idnumbers)
        sd = np.std(idnumbers)
        plt.bar(i + float(i)*3/2, meanID, yerr=sd, ecolor="black",
                width=1.5, color=color[0], alpha=0.8, edgecolor=color[0])
        i += 1

    plt.subplots_adjust(bottom=0.3, left=0.2, right=0.75)
    plt.ylabel('Peptide ID')
    plt.xticks(map(lambda x: float(x) * 5/2 + 1.0, range(len(dfs.keys()))),
               sorted(dfs.keys()), ha='right', rotation=45)
    plt.savefig("peptideIDwithSD.png")
    plt.clf()


def chargeState_barplot(dfs):
    i = 0
    color = sns.color_palette("muted")
    sns.set_style("whitegrid")
    plt.grid(False)
    plt.grid(axis="y")

    for k, v in sorted(dfs.items()):
        area2 = v[v[u'Charge'] == 2][u'PepArea'].sum() / v[u'PepArea'].sum()
        area3 = v[v[u'Charge'] == 3][u'PepArea'].sum() / v[u'PepArea'].sum()
        area4 = v[v[u'Charge'] >= 4][u'PepArea'].sum() / v[u'PepArea'].sum()
        plt.bar(i + float(i)/2, area2, width=1, color=color[0],
                alpha=0.8, edgecolor=color[0])
        plt.bar(i + float(i)/2, area3, bottom=area2, width=1,
                color=color[1], alpha=0.8, edgecolor=color[1])
        plt.bar(i + float(i)/2, area4, bottom=area2+area3, width=1,
                color=color[2], alpha=0.8, edgecolor=color[2])
        i += 1

    plt.legend(["Charge 2", "Charge 3", "Charge 4+"],
               bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, left=0.15, right=0.65)
    plt.xticks(map(lambda x: float(x) * 3/2 + 0.5, range(len(dfs.keys()))),
               sorted(dfs.keys()), ha='right', rotation=45)
    plt.ylabel('Peak area ratio')
    plt.savefig("ChargeState.png")
    plt.clf()


def chargeStateID_barplot(dfs):
    i = 0
    color = sns.color_palette("muted")
    sns.set_style("whitegrid")
    plt.grid(False)
    plt.grid(axis="y")

    for k, v in sorted(dfs.items()):
        id2 = float(len(v[v[u'Charge'] == 2])) / float(len(v[u'Charge']))
        id3 = float(len(v[v[u'Charge'] == 3])) / float(len(v[u'Charge']))
        id4 = float(len(v[v[u'Charge'] >= 4])) / float(len(v[u'Charge']))
        plt.bar(i + float(i)/2, id2, width=1, color=color[0],
                alpha=0.8, edgecolor=color[0])
        plt.bar(i + float(i)/2, id3, bottom=id2, width=1,
                color=color[1], alpha=0.8, edgecolor=color[1])
        plt.bar(i + float(i)/2, id4, bottom=id2+id3, width=1,
                color=color[2], alpha=0.8, edgecolor=color[2])
        i += 1

    plt.legend(["Charge 2", "Charge 3", "Charge 4+"],
               bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, left=0.15, right=0.65)
    plt.xticks(map(lambda x: float(x) * 3/2 + 0.5, range(len(dfs.keys()))),
               sorted(dfs.keys()), ha='right', rotation=45)
    plt.ylabel('PSM')
    plt.savefig("ChargeStateID.png")
    plt.clf()


def numberOfTMT_barplot(dfs):
    i = 0
    color = sns.color_palette("muted")
    sns.set_style("whitegrid")
    plt.grid(False)
    plt.grid(axis="y")

    for k, v in sorted(dfs.items()):
        area0 = v[v[u'number_of_TMT'] == 0][u'PepArea'].sum()\
            / v[u'PepArea'].sum()
        area1 = v[v[u'number_of_TMT'] == 1][u'PepArea'].sum()\
            / v[u'PepArea'].sum()
        area2 = v[v[u'number_of_TMT'] == 2][u'PepArea'].sum()\
            / v[u'PepArea'].sum()
        area3 = v[v[u'number_of_TMT'] == 3][u'PepArea'].sum()\
            / v[u'PepArea'].sum()
        area4 = v[v[u'number_of_TMT'] == 4][u'PepArea'].sum()\
            / v[u'PepArea'].sum()
        plt.bar(i + float(i)/2, area0, width=1, color=color[0],
                alpha=0.8, edgecolor=color[0])
        plt.bar(i + float(i)/2, area1, bottom=area0, width=1,
                color=color[1], alpha=0.8, edgecolor=color[1])
        plt.bar(i + float(i)/2, area2, bottom=area0+area1, width=1,
                color=color[2], alpha=0.8, edgecolor=color[2])
        plt.bar(i + float(i)/2, area3, bottom=area0+area1+area2,
                width=1, color=color[3], alpha=0.8, edgecolor=color[2])
        plt.bar(i + float(i)/2, area4, bottom=area0+area1+area2+area3,
                width=1, color=color[4], alpha=0.8, edgecolor=color[2])

        i += 1

    plt.legend(["TMT*0", "TMT*1", "TMT*2", "TMT*3", "TMT*4"],
               bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, left=0.15, right=0.65)
    plt.xticks(map(lambda x: float(x) * 3/2 + 0.5, range(len(dfs.keys()))),
               sorted(dfs.keys()), ha='right', rotation=45)
    plt.savefig("number_of_TMT_bar.png")
    plt.clf()


def numberOfK_barplot(dfs):
    i = 0
    color = sns.color_palette("muted")
    sns.set_style("whitegrid")
    plt.grid(False)
    plt.grid(axis="y")

    for k, v in sorted(dfs.items()):
        area0 = v[v[u'number_of_k'] == 0][u'PepArea'].sum()\
            / v[u'PepArea'].sum()
        area1 = v[v[u'number_of_k'] == 1][u'PepArea'].sum()\
            / v[u'PepArea'].sum()
        area2 = v[v[u'number_of_k'] == 2][u'PepArea'].sum()\
            / v[u'PepArea'].sum()
        area3 = v[v[u'number_of_k'] == 3][u'PepArea'].sum()\
            / v[u'PepArea'].sum()
        plt.bar(i + float(i)/2, area0, width=1, color=color[0],
                alpha=0.8, edgecolor=color[0])
        plt.bar(i + float(i)/2, area1, bottom=area0, width=1,
                color=color[1], alpha=0.8, edgecolor=color[1])
        plt.bar(i + float(i)/2, area2, bottom=area0+area1, width=1,
                color=color[2], alpha=0.8, edgecolor=color[2])
        plt.bar(i + float(i)/2, area3, bottom=area0+area1+area2,
                width=1, color=color[3], alpha=0.8, edgecolor=color[2])

        i += 1

    plt.legend(["K*0", "K*1", "K*2", "K*3"],
               bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, left=0.15, right=0.65)
    plt.ylabel('Peak area ratio')
    plt.xticks(map(lambda x: float(x) * 3/2 + 0.5, range(len(dfs.keys()))),
               sorted(dfs.keys()), ha='right', rotation=45)
    plt.savefig("number_of_K_bar.png")
    plt.clf()


def numberOfKID_barplot(dfs):
    i = 0
    color = sns.color_palette("muted")
    sns.set_style("whitegrid")
    plt.grid(False)
    plt.grid(axis="y")
    plt.title("number of K")

    for k, v in sorted(dfs.items()):
        area0 = len(v[v[u'number_of_k'] == 0])
        area1 = len(v[v[u'number_of_k'] == 1])
        area2 = len(v[v[u'number_of_k'] == 2])
        area3 = len(v[v[u'number_of_k'] == 3])
        idn = area0 + area1 + area2 + area3
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


def multiphospho_peakarea(dfs):
    i = 0
    sns.set_style("whitegrid")
    plt.grid(False)
    plt.grid(axis="y")

    for k, v in sorted(dfs.items()):
        color = sns.color_palette("muted")[i % 6]
        x = v.groupby([u'number_of_phospho'])[u'PepArea'].sum()
        plt.bar(x.index + 0.1 + float(i)*2/5, x.values,
                width=0.4, color=color, alpha=0.8, edgecolor=color)
        i += 1

    plt.xlim([0, 5])
    plt.xticks(map(lambda x: float(x) + 0.5, range(5)),
               ["phospho_0", "phospho_1", "phospho_2",
                "phospho_3", "phospho_4"])
    plt.xlabel('number of phospho')
    plt.ylabel('Peak area')
    plt.legend(sorted(dfs.keys()), bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, right=0.65)
    plt.savefig("multiphospho_peakarea.png")
    plt.clf()


def acidicresidue_peakarea(dfs):
    i = 0
    sns.set_style("whitegrid")
    plt.grid(False)
    plt.grid(axis="y")

    for k, v in sorted(dfs.items()):
        x = v.groupby([u'number_of_acidic_residue'])[u'PepArea'].sum()
        color = sns.color_palette("muted")[i % 6]
        plt.bar(x.index + float(i)*3/10, x.values,
                width=0.3, color=color, alpha=0.8, edgecolor=color)
        i += 1

    plt.xlim([0, 20])
    plt.xlabel('number of D, E')
    plt.ylabel('Peak area')
    plt.legend(sorted(dfs.keys()), bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, left=0.1, right=0.6)
    plt.savefig("acidicresidue_peakarea.png")
    plt.clf()


def basicresidue_peakarea(dfs):
    i = 0
    sns.set_style("whitegrid")
    plt.grid(False)
    plt.grid(axis="y")
    plt.title("basic-residue peakarea")

    for k, v in sorted(dfs.items()):
        x = v.groupby([u'number_of_basic_residue'])[u'PepArea'].sum()
        color = sns.color_palette("muted")[i % 6]
        plt.bar(x.index + float(i)/10, x.values, width=0.1,
                color=color, alpha=0.8, edgecolor=color)
        i += 1

    plt.xlim([0, 10])
    plt.xticks(map(lambda x: float(x) + 0.5, range(10)), range(10))
    plt.legend(sorted(dfs.keys()), bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, right=0.65)
    plt.savefig("basicresidue_peakarea.png")
    plt.clf()


def pI_peakarea(dfs):
    i = 0
    bins = map(lambda x: float(x)/4, range(80))
    sns.set_style("whitegrid")
    plt.grid(False)
    plt.grid(axis="y")
    plt.title("pI peakarea")

    for k, v in sorted(dfs.items()):
        v[u'pI_bin'] = pd.cut(v.pI, bins)
        x = v.groupby([u'pI_bin'])[u'PepArea'].sum()
        x.index = map(lambda x: float(x)/4, range(79))
        color = sns.color_palette("muted")[i % 6]
        plt.bar(x.index + float(i)/10, x.values, width=0.1,
                color=color, alpha=0.8, edgecolor=color)
        i += 1

    plt.legend(sorted(dfs.keys()), bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, right=0.65)
    plt.savefig("pI_peakarea.png")
    plt.clf()


def RetTime_peakarea(dfs):
    i = 0
    bins = map(lambda x: float(x), range(110))
    sns.set_style("whitegrid")
    plt.grid(False)
    plt.grid(axis="y")
    plt.title("RetTime peakarea")

    for k, v in sorted(dfs.items()):
        v[u'RetTime_bin'] = pd.cut(v.RetTime, bins)
        x = v.groupby([u'RetTime_bin'])[u'PepArea'].sum()
        x.index = map(lambda x: float(x), range(109))
        color = sns.color_palette("muted")[i % 6]
        plt.bar(x.index + float(i)/2, x.values, width=0.4,
                color=color, alpha=1, edgecolor=color)
        i += 1

    plt.legend(sorted(dfs.keys()), bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, right=0.65)
    plt.savefig("RetTime_peakarea.png")
    plt.clf()


def jointplot(dfs):
    i = 0
    for k, v in sorted(dfs.items()):
        if k == "151218ko16-RP-C18.raw" or k == "151218ko17-sol.raw":

            x = v[u'number_of_phospho']
            y = v[u'pI']
            sns.set_style("whitegrid")
            sns.jointplot(x=x, y=y)
            i += 1
    plt.savefig("joint.png")
    plt.clf()


def multiphospho_hist(dfs):
    sns.set_style("white")
    x = []
    color = []
    i = 0
    for k, v in sorted(dfs.items()):
        x.append(v[u'number_of_phospho'])
        color.append(sns.color_palette("muted")[i % 6])
        i += 1
    plt.hist(x, color=color, histtype="bar",
             alpha=0.8, range=[0, 5], bins=5)
    plt.legend(sorted(dfs.keys()), bbox_to_anchor=(1.05, 1), loc='upper left', fontsize="20")
    plt.xticks(map(lambda x: float(x) + 0.5, range(5)), range(5))
    plt.xlabel('number of phospho')
    plt.ylabel('Peptide ID')
    plt.subplots_adjust(bottom=0.3, left=0.2, right=0.65)
    plt.savefig("multiphospho_hist.png")
    plt.clf()


def acidicresidue_hist(dfs):
    i = 0
    sns.set_style("white")

    for k, v in sorted(dfs.items()):
        color = sns.color_palette("muted")[i % 6]
        plt.hist(v[u'number_of_acidic_residue'].values, color=color,
                 histtype="step", linewidth=1.5, alpha=1,
                 range=[0, 10], bins=10)
        i += 1
    plt.xticks(map(lambda x: float(x) + 0.5, range(10)), range(10))
    plt.xlabel('number of D, E')
    plt.ylabel('Peptide ID')
    plt.legend(sorted(dfs.keys()), bbox_to_anchor=(1.05, 1), loc='upper left', fontsize="20")
    plt.subplots_adjust(bottom=0.3, right=0.65)
    plt.savefig("acidicresidue_hist.png")
    plt.clf()


def basicresidue_hist(dfs):
    i = 0
    sns.set_style("white")
    plt.title("basic-residue ID")

    for k, v in sorted(dfs.items()):
        color = sns.color_palette("muted")[i % 6]
        plt.hist(v[u'number_of_basic_residue'].values, color=color,
                 histtype="step", linewidth=1.5, alpha=1,
                 range=[0, 6], bins=6)
        i += 1
    plt.xticks(map(lambda x: float(x) + 0.5, range(6)), range(6))
    plt.legend(sorted(dfs.keys()), bbox_to_anchor=(1.05, 1), loc='upper left', fontsize="20")
    plt.subplots_adjust(bottom=0.3, right=0.65)
    plt.savefig("basicresidue_hist.png")
    plt.clf()


def peptide_length_hist(dfs):
    i = 0
    sns.set_style("white")
    for k, v in sorted(dfs.items()):
        color = sns.color_palette("muted")[i % 6]
        plt.hist(v[u'peptide_length'].values, color=color,
                 histtype="step", linewidth=1.5, alpha=1,
                 range=[5, 50], bins=45)
        i += 1
    plt.legend(sorted(dfs.keys()), bbox_to_anchor=(1.05, 1), loc='upper left', fontsize="20")
    plt.xlabel('Peptide length')
    plt.ylabel('Peptide ID')
    plt.subplots_adjust(bottom=0.3, right=0.65)
    plt.savefig("peptide_length_hist.png")
    plt.clf()


def GRAVYscore_hist(dfs):
    i = 0
    sns.set_style("white")

    for k, v in sorted(dfs.items()):
        color = sns.color_palette("muted")[i % 6]
        plt.hist(v[u'GRAVY_score'].values, color=color,
                 histtype="step", linewidth=1.5, alpha=1,
                 range=[-4, 4], bins=80)
        i += 1
    plt.legend(sorted(dfs.keys()), bbox_to_anchor=(1.05, 1), loc='upper left', fontsize="20")
    plt.xlabel('GRAVY score')
    plt.ylabel('Peptide ID')
    plt.subplots_adjust(bottom=0.3, right=0.65)
    plt.savefig("GRAVY_hist.png")
    plt.clf()


def pI_hist(dfs):
    i = 0
    sns.set_style("white")
    plt.title("pI")

    for k, v in sorted(dfs.items()):
        color = sns.color_palette("muted")[i % 6]
        plt.hist(v[u'pI'].values, color=color,
                 histtype="step", linewidth=1.5, alpha=1,
                 range=[0, 15], bins=60)
#        plt.axvline(np.median(v[u'pI']), color=color)
        i += 1

    plt.legend(sorted(dfs.keys()), bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, right=0.65)
    plt.savefig("pI_hist.png")
    plt.clf()


def RetTime_hist(dfs):
    i = 0
    sns.set_style("white")

    for k, v in sorted(dfs.items()):
        color = sns.color_palette("muted")[i % 6]
        plt.hist(v[u'RetTime'].values, color=color,
                 histtype="step", linewidth=1.5, alpha=1,
                 range=[0, 600], bins=300)
        i += 1

    plt.legend(sorted(dfs.keys()), bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xlabel('Ret time(min)')
    plt.xlim(0,600)
    plt.ylabel('Peptide ID')
    plt.subplots_adjust(bottom=0.3, right=0.65)
    plt.savefig("RetTime_hist.png")
    plt.clf()


def peakarea2Dhist(dfs):
    keys = sorted(dfs.keys())
    df = func.CombinePepArea([dfs[keys[0]], dfs[keys[1]]])
    df = df[df[u'Count'] == 2]
    plt.title("Peak area: n=" + str(len(df)))
    plt.xlabel(keys[0])
    plt.ylabel(keys[1])

    x = np.log2(df[u'PepArea_1'].values)
    y = np.log2(df[u'PepArea_2'].values)
    cmap = mpl.cm.jet
    cmap.set_under("w")
    plt.hist2d(x, y, bins=100, range=[[5, 40], [5, 40]], cmap=cmap, cmin=1)
    plt.colorbar()
    plt.show()


def numberofKhist(dfs):
    i = 0
    sns.set_style("white")
    plt.title("number of K")

    for k, v in sorted(dfs.items()):
        color = sns.color_palette("muted")[i % 6]
        plt.hist(v[u'number_of_k'].values, color=color,
                 histtype="step", linewidth=1.5, alpha=1,
                 range=[0, 4], bins=4)
        i += 1
    plt.xticks(map(lambda x: float(x) + 0.5, range(4)), range(4))
    plt.legend(sorted(dfs.keys()), bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, right=0.65)
    plt.savefig("numberK_hist.png")
    plt.clf()


def numberofRhist(dfs):
    i = 0
    sns.set_style("white")
    plt.title("number of R")

    for k, v in sorted(dfs.items()):
        color = sns.color_palette("muted")[i % 6]
        plt.hist(v[u'number_of_r'].values, color=color,
                 histtype="step", linewidth=1.5, alpha=1,
                 range=[0, 4], bins=4)
        i += 1
    plt.xticks(map(lambda x: float(x) + 0.5, range(4)), range(4))
    plt.legend(sorted(dfs.keys()), bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, right=0.65)
    plt.savefig("numberR_hist.png")
    plt.clf()


def numberofDhist(dfs):
    i = 0
    sns.set_style("white")
    plt.title("number of D")

    for k, v in sorted(dfs.items()):
        color = sns.color_palette("muted")[i % 6]
        plt.hist(v[u'number_of_d'].values, color=color,
                 histtype="step", linewidth=1.5, alpha=1,
                 range=[0, 10], bins=10)
        i += 1
    plt.xticks(map(lambda x: float(x) + 0.5, range(10)), range(10))
    plt.legend(sorted(dfs.keys()), bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, right=0.65)
    plt.savefig("numberD_hist.png")
    plt.clf()


def numberofEhist(dfs):
    i = 0
    sns.set_style("white")
    plt.title("number of E")

    for k, v in sorted(dfs.items()):
        color = sns.color_palette("muted")[i % 6]
        plt.hist(v[u'number_of_e'].values, color=color,
                 histtype="step", linewidth=1.5, alpha=1,
                 range=[0, 10], bins=10)
        i += 1
    plt.xticks(map(lambda x: float(x) + 0.5, range(10)), range(10))
    plt.legend(sorted(dfs.keys()), bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, right=0.65)
    plt.savefig("numberE_hist.png")
    plt.clf()


def peptScoreKdeHistK(dfs):
    sns.set_style("white")
    plt.title("Peptide Score")

    for k, v in sorted(dfs.items()):
        n0 = 'K=0 (' + str(len(v[v[u'number_of_k'] == 0])) + ')'
        n1 = 'K=1 (' + str(len(v[v[u'number_of_k'] == 1])) + ')'
        n2 = 'K>=2 (' + str(len(v[v[u'number_of_k'] >= 2])) + ')'
        color = sns.color_palette("muted")[0]
        sns.kdeplot(v[v[u'number_of_k'] == 0][u'PeptScore'].values,
                    color=color, label=n0, legend=False)
        color = sns.color_palette("muted")[1]
        sns.kdeplot(v[v[u'number_of_k'] == 1][u'PeptScore'].values,
                    color=color, label=n1, legend=False)
        color = sns.color_palette("muted")[2]
        sns.kdeplot(v[v[u'number_of_k'] >= 2][u'PeptScore'].values,
                    color=color, label=n2, legend=False)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.subplots_adjust(bottom=0.3, right=0.65)
        title = str(k) + "_PeptScore_kdehistk.png"
        plt.savefig(title)
        plt.clf()


def peptScoreKdeHistK1(dfs):
    sns.set_style("white")
    plt.title("Peptide Score")
    for k, v in sorted(dfs.items()):
        v = v[v[u'number_of_k'] == 1]
        nk = 'terminal K (' + str(len(v[v[u'Seq'].str.endswith('K')])) + ')'
        nr = 'terminal R (' + str(len(v[v[u'Seq'].str.endswith('R')])) + ')'
        color = sns.color_palette("muted")[0]
        sns.kdeplot(v[v[u'Seq'].str.endswith('K')][u'PeptScore'].values,
                    color=color, label=nk, legend=False)
        color = sns.color_palette("muted")[1]
        sns.kdeplot(v[v[u'Seq'].str.endswith('R')][u'PeptScore'].values,
                    color=color, label=nr, legend=False)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.subplots_adjust(bottom=0.3, right=0.65)
        title = str(k) + "_PeptScore_kdehistk1.png"
        plt.savefig(title)
        plt.clf()


def peptScoreHist(dfs):
    i = 0
    for k, v in sorted(dfs.items()):
        color = sns.color_palette("muted")[i % 6]
        plt.hist(v[u'PeptScore'].values, color=color,
                 histtype="step", linewidth=1.5, alpha=1,
                 range=[0, 150], bins=50)
        i += 1
    plt.legend(sorted(dfs.keys()), bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, right=0.65)
    plt.savefig("peptScore_hist.png")
    plt.clf()


def peptScoreKdeHistR1(dfs):
    sns.set_style("white")
    plt.title("Peptide Score")
    for k, v in sorted(dfs.items()):
        v = v[v[u'number_of_r'] == 1]
        nk = 'terminal K (' + str(len(v[v[u'Seq'].str.endswith('K')])) + ')'
        nr = 'terminal R (' + str(len(v[v[u'Seq'].str.endswith('R')])) + ')'
        color = sns.color_palette("muted")[0]
        sns.kdeplot(v[v[u'Seq'].str.endswith('K')][u'PeptScore'].values,
                    color=color, label=nk, legend=False)
        color = sns.color_palette("muted")[1]
        sns.kdeplot(v[v[u'Seq'].str.endswith('R')][u'PeptScore'].values,
                    color=color, label=nr, legend=False)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.subplots_adjust(bottom=0.3, right=0.65)
        title = str(k) + "_PeptScore_kdehistr1.png"
        plt.savefig(title)
        plt.clf()


def peptScoreKdeHistCharge(dfs):
    sns.set_style("white")
    for k, v in sorted(dfs.items()):
        c2 = 'Charge 2 (' + str(len(v[v[u'Charge'] == 2])) + ')'
        c3 = 'Charge 3 (' + str(len(v[v[u'Charge'] >= 3])) + ')'
        color = sns.color_palette("muted")[0]
        plt.title("Peptide Score")
        sns.kdeplot(v[v[u'Charge'] == 2][u'PeptScore'].values,
                    color=color, label=c2, legend=False)
        color = sns.color_palette("muted")[1]
        sns.kdeplot(v[v[u'Charge'] >= 3][u'PeptScore'].values,
                    color=color, label=c3, legend=False)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.subplots_adjust(bottom=0.3, right=0.65)
        title = str(k) + "_PeptScore_kdehistcharge.png"
        plt.savefig(title)
        plt.clf()


def peptHalfWidthHist(dfs):
    i = 0
    sns.set_style("white")
    for k, v in sorted(dfs.items()):
        color = sns.color_palette("muted")[i % 6]
        plt.title("HalfWidth")
        sns.kdeplot(v[u'PepHalfWidth'].values, color=color,
                    label=k, legend=False)
        i += 1
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, right=0.65)
    title = "HalfWidth.png"
    plt.savefig(title)
    plt.clf()


def peptHalfWidthHistK(dfs):
    sns.set_style("white")
    for k, v in sorted(dfs.items()):
        i = 0
        plt.title("HalfWidth")
        k0 = 'K=0 (' + str(len(v[v[u'number_of_k'] == 0])) + ')'
        k1 = 'K=1 (' + str(len(v[v[u'number_of_k'] == 1])) + ')'
        k2 = 'K=2 (' + str(len(v[v[u'number_of_k'] >= 2])) + ')'
        color = sns.color_palette("muted")[i % 6]
        sns.kdeplot(v[v[u'number_of_k'] == 0][u'PepHalfWidth'].values,
                    color=color, label=k0, legend=False)
        i += 1
        color = sns.color_palette("muted")[i % 6]
        sns.kdeplot(v[v[u'number_of_k'] == 1][u'PepHalfWidth'].values,
                    color=color, label=k1, legend=False)
        i += 1
        color = sns.color_palette("muted")[i % 6]
        sns.kdeplot(v[v[u'number_of_k'] >= 2][u'PepHalfWidth'].values,
                    color=color, label=k2, legend=False)

        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.subplots_adjust(bottom=0.3, right=0.65)
        title = str(k) + "_HalfWidthK.png"
        plt.savefig(title)
        plt.clf()


def pepAreaRank(dfs):
    data = func.CombinePepArea(dfs)
    keys = sorted(dfs.keys())
    data = data[data['Count'] == 2]
    data['AreaRatio'] = np.log2(data['PepArea_2'] / data['PepArea_1'])
    data = data.sort(columns='AreaRatio')
    plt.plot(data['AreaRatio'], marker="o", markersize=4.0, linestyle="none")
    plt.title(str(keys[1]) + '/' + str(keys[0]))
    plt.show()


def pepAreaBoxPlot(dfs):
    sns.set_style("white")
    plt.title("PepArea")
    for k, v in sorted(dfs.items()):
        missK = np.log2(v[v[u'Seq'].str.contains('.*K.+')][u'PepArea'].values)
        missR = np.log2(v[v[u'Seq'].str.contains('.*R.+')][u'PepArea'].values)
        plt.boxplot([missK, missR], labels=['missK', 'missR'])
        plt.grid(axis='y')
        title = str(k) + "_PeptArea_msclvg.png"
        plt.savefig(title)
        plt.clf()


def mzViolinPlot(data):
    sns.set_style("white")
    plt.title("m/z distribution")
    data['Charge2'] = 'Charge: 3+'
    data.ix[data[u'Charge'] == 2, 'Charge2'] = 'Charge: 2'
    sns.violinplot(x=data[u'RawdataFile'], y=data[u'ObsMz'],
                   pallete="muted",
                   split=True, inner="quartile")
    li = data[u'RawdataFile'].drop_duplicates().values
    leg = []
    for rawname in li:
        name = str(rawname) + ' ('\
            + str(len(data[data[u'RawdataFile'] == rawname])) + ')'
        leg.append(name)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.subplots_adjust(bottom=0.3, right=0.65)

    plt.savefig("MassChargedist.png")
    plt.clf()


def phosphoSiteLocalization(dfs):
    x = []
    for k, v in sorted(dfs.items()):
        s = pd.Series()
        s['phosphosite_1'] = v[u'ModDetail']\
            .str.extract('79.9663@[STY]:(\d{1,2})').astype(float)\
            / v[u'peptide_length']
        s['phosphosite_2'] = v[u'ModDetail']\
            .str.extract('79.9663@[STY]:\d{1,2}\,.?79.9663@[STY]:(\d{1,2})')\
            .astype(float) / v[u'peptide_length']
        s['phosphosite_3'] = v[u'ModDetail']\
            .str.extract('79.9663@[STY]:\d{1,2}\,.?79.9663@[STY]:\d{1,2}\,.?79.9663@[STY]:(\d{1,2})')\
            .astype(float) / v[u'peptide_length']
        s['phosphosite_4'] = v[u'ModDetail']\
            .str.extract('79.9663@[STY]:\d{1,2}\,.?79.9663@[STY]:\d{1,2}\,.?79.9663@[STY]:\d{1,2}\,.?79.9663@[STY]:(\d{1,2})')\
            .astype(float) / v[u'peptide_length']
        x.append(list(s['phosphosite_1'].dropna())
                 + list(s['phosphosite_2'].dropna())
                 + list(s['phosphosite_3'].dropna())
                 + list(s['phosphosite_4'].dropna()))
    sns.set_style("white")
    color = sns.color_palette("muted")

    plt.hist(x, bins=10, color=[color[0], color[1]],
             alpha=0.8, histtype='step')
    plt.show()


def lysineLocalization(dfs):
    x = []
    for k, v in sorted(dfs.items()):
        li = list(v[u'Seq'])
        k1 = [(seq.find('K') + 1) / len(seq) for seq in li]
        k2 = [(seq.find('K', seq.find('K') + 1) + 1)
              / len(seq) for seq in li]
        k3 = [(seq.find('K', seq.find('K', seq.find('K') + 1) + 1) + 1)
              / len(seq) for seq in li]
        posK = k1
        posK = [i for i in posK if i != 0]
        plt.hist(posK)
        plt.show()
        print [len(seq) for seq in li]
        print [seq.find('K') for seq in li]
    sns.set_style("white")
    color = sns.color_palette("muted")
    plt.hist(x, bins=10, color=[color[0], color[1]], alpha=0.8, histtype='step')
    plt.show()


def reporterionIntensity_boxplot(df, tagname='TMT10plex'):
    if tagname == 'TMT10plex':
        reporter_values = (
            u'iTRAQ:126.127725', u'iTRAQ:127.124760',
            u'iTRAQ:127.131079', u'iTRAQ:128.128114',
            u'iTRAQ:128.134433', u'iTRAQ:129.131468',
            u'iTRAQ:129.137787', u'iTRAQ:130.134822',
            u'iTRAQ:130.141141', u'iTRAQ:131.138176')
    elif tagname == 'TMT6plex':
        reporter_values = (
            u'iTRAQ:126.127725', u'iTRAQ:127.124760',
            u'iTRAQ:128.134433', u'iTRAQ:129.131468',
            u'iTRAQ:130.141141', u'iTRAQ:131.138176')
    reporter_intensities = [np.log10(df[cn].values)[~np.isinf(np.log10(df[cn].values))]
                            for cn in reporter_values]
    sns.set_style("ticks")
    sns.despine()
    plt.boxplot(reporter_intensities, labels=['126', '127N', '127C',
                                              '128N', '128C', '129N',
                                              '129C', '130N', '130C', '131'])
    plt.ylabel('log10 (Intensity)')
    title = "ReporterIntensity_boxplot.png"
    sns.despine()
    plt.savefig(title)
    plt.clf()
    beeswarm(reporter_intensities, method="hex", s=8)
    plt.show()
    plt.clf()


def reproducibilityCheck(df, tagname='TMT10plex'):
    if tagname == 'TMT10plex':
        reporter_values = (
            u'iTRAQ:126.127725', u'iTRAQ:127.124760',
            u'iTRAQ:127.131079', u'iTRAQ:128.128114',
            u'iTRAQ:128.134433', u'iTRAQ:129.131468',
            u'iTRAQ:129.137787', u'iTRAQ:130.134822',
            u'iTRAQ:130.141141', u'iTRAQ:131.138176')
    elif tagname == 'TMT6plex':
        reporter_values = (
            u'iTRAQ:126.127725', u'iTRAQ:127.124760',
            u'iTRAQ:128.134433', u'iTRAQ:129.131468',
            u'iTRAQ:130.141141', u'iTRAQ:131.138176')

    reporters = {#'RP-C18': [u'iTRAQ:129.137787', u'iTRAQ:130.134822', u'iTRAQ:130.141141']}
                 #'Solution': [u'iTRAQ:128.128114', u'iTRAQ:128.134433', u'iTRAQ:129.131468']}
                 'TMT-HM': [u'iTRAQ:126.127725', u'iTRAQ:127.124760', u'iTRAQ:127.131079']}

    df_t = pd.DataFrame()
    if bool(reporter_values) is True:
        for cn in reporter_values:
            series = df.ix[:, df.columns.map(
                lambda x: x.startswith(cn))]
            df_t = pd.concat([df_t, series], axis=1)
    df_t = df_t.join(df[u'SeqModModDetail'])
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
        sns.set_style("ticks")
        fig = plt.figure()
        fig.text(0.5, 0.85, fig_title, size=18, ha='center')
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
        ax3.set_xlabel('log10 (Intensity)\nRep. 1', size=12)
        ax3.set_ylabel('Rep. 3\nlog10 (Intensity)', size=12)
        r1, p_value1 = stats.pearsonr(Rep1, Rep2)
        r2, p_value2 = stats.pearsonr(Rep2, Rep3)
        r3, p_value3 = stats.pearsonr(Rep1, Rep3)
        ax1.text(5, 1.5, "$R^2=%s$" % str(r1**2)[0:6], size=13)
        ax2.text(5, 1.5, "$R^2=%s$" % str(r2**2)[0:6], size=13)
        ax3.text(5, 1.5, "$R^2=%s$" % str(r3**2)[0:6], size=13)
        plt.subplots_adjust(bottom=0.2, wspace=0.2, hspace=0.15)
        sns.despine()
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
        df_result.to_csv('cv.csv')
        color_list = map(set_color, df_result.cv)
        plt.scatter(df_result['cv'], lnAverage, c=color_list, s=8, alpha=0.6)
        plt.xlim(0, 100)
        plt.ylim(2, 7)
        plt.axvline(x=20, ls='--', lw=0.8, color='black')
        plt.text(20, 6.2, "C.V. $\leqq$ 20\n(n = %d)" % len(df_result[df_result['cv'] <= 20]), size=14)
        plt.xlabel('C.V. (%)', size=14)
        plt.ylabel('Average log10(Intensity)', size=14)
        plt.subplots_adjust(bottom=0.2, left=0.1)
        plt.title(fig_title, size=16)
        sns.despine()
        ftitle_cv = title + "_treePlot.png"
        plt.savefig(ftitle_cv)
        plt.clf()



def log2histID(df, tagname='TMT10plex'):
    if tagname == 'TMT10plex':
        reporter_values = (
            u'iTRAQ:126.127725', u'iTRAQ:127.124760',
            u'iTRAQ:127.131079', u'iTRAQ:128.128114',
            u'iTRAQ:128.134433', u'iTRAQ:129.131468',
            u'iTRAQ:129.137787', u'iTRAQ:130.134822',
            u'iTRAQ:130.141141', u'iTRAQ:131.138176')
    elif tagname == 'TMT6plex':
        reporter_values = (
            u'iTRAQ:126.127725', u'iTRAQ:127.124760',
            u'iTRAQ:128.134433', u'iTRAQ:129.131468',
            u'iTRAQ:130.141141', u'iTRAQ:131.138176')
    elif tagname == '160227':
        reporter_values = (
            u'iTRAQ:128.128114',
            u'iTRAQ:128.134433', u'iTRAQ:129.131468',
            u'iTRAQ:129.137787', u'iTRAQ:130.134822',
            u'iTRAQ:130.141141')

    reporters = {'RP-C18': [u'iTRAQ:129.137787', u'iTRAQ:130.134822', u'iTRAQ:130.141141'],
                 'Solution': [u'iTRAQ:128.128114', u'iTRAQ:128.134433', u'iTRAQ:129.131468']}

    df_temp = pd.DataFrame()
    if bool(reporter_values) is True:
        for cn in reporter_values:
            series = df.ix[:, df.columns.map(
                lambda x: x.startswith(cn))]
            df_temp = pd.concat([df_temp, series], axis=1)
    df_temp = df_temp.join(df[u'SeqModModDetail'])
    df_temp = df_temp.replace(0, np.nan)
    df_temp = df_temp.dropna()
    df_temp = df_temp.groupby([u'SeqModModDetail']).mean()
    for title, reporter in reporters.items():
        df_Int = pd.DataFrame()
        index = []
        average = []

        for cn in reporter:
            df_Int[cn] = df_temp[cn]

        for i in df_Int.index:
            index.append(i)
            average.append(np.nanmean(df_Int.loc[i]))
            df_temp[u'AVERAGE: ' + title] = pd.Series(average, index=index)


    df_temp['log2'] = np.log2(df_temp[u'AVERAGE: ' + 'RP-C18'] / df_temp[u'AVERAGE: ' + 'Solution'])
    median = np.median(df_temp['log2'])
    df_temp['SeqModModDetail'] = df_temp.index
    df = pd.merge(func.RemoveDupli(df), df_temp,
                  on='SeqModModDetail', how='left')
    color = sns.color_palette("muted")[0]
    sns.set_style("ticks")

    plt.hist(df['log2'], color=color,
             histtype="step", linewidth=1.5, alpha=1,
             range=[-4, 4], bins=60)
    plt.axvline(x=median, ls='--', lw=0.8, color='black')
    plt.text(median + 0.5, 160, 'median = %f' % median, size=12)
    plt.xlabel('log2 (RP-C18 / Solution)')
    plt.ylabel('Peptide ID')
    plt.subplots_adjust(bottom=0.2, left=0.1)
    sns.despine()
    plt.savefig('log2hist.png')
    plt.clf()
    plt.scatter(df['RetTime'], df['log2'], alpha=0.8, s=4)
    plt.axhline(y=median, ls='--', lw=0.6, color='black')
    plt.text(90, median-0.5, 'median = %f' % median, size=12)
    plt.xlabel('RetTime (min)')
    plt.ylabel('log2 (RP-C18 / Solution)')
    plt.xlim(0, 110)
    plt.ylim(-4, 4)
    plt.subplots_adjust(bottom=0.2, left=0.1)
    sns.despine()
    plt.savefig('RettimeVSlog2.png')
    plt.clf()


    def set_color(ratio):
        if ratio < -0.5:
            return "red"
        else:
            return "gray"


    color_list = map(set_color, df.log2)
    plt.scatter(df['RetTime'], df['log2'], c=color_list, alpha=0.8, s=4)
    plt.axhline(y=-0.5, ls='--', lw=0.6, color='black')
    plt.text(90, -1, 'ratio < -0.5\n(n = %d)' % len(df[df['log2'] < -0.5]), size=12)
    plt.text(95, 3, 'n = %d' % len(df_temp), size=14)
    plt.xlabel('RetTime (min)')
    plt.ylabel('log2 (RP-C18 / Solution)')
    plt.xlim(0, 110)
    plt.ylim(-4, 4)
    plt.subplots_adjust(bottom=0.2, left=0.1)
    sns.despine()
    plt.savefig('RettimeVSlog2-minus.png')
    plt.clf()
    bins = map(lambda x: float(x) * 10, range(11))
    df['RetTime_bin'] = pd.cut(df.RetTime, bins)
    df.boxplot(column='log2', by='RetTime_bin', grid=False, labels=[10,20,30,40,50,60,70,80,90,100])
    plt.grid(False)
    sns.despine()
    plt.show()

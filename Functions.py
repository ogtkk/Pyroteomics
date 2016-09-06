import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st
import os
from collections import Counter
import re

fnamepat = re.compile(u"D.+ManualInputFile\\\\|D.+ManualInputFile/")
ukey = u'SeqModModDetail'
reporter_values = (
    u'iTRAQ:126.127725', u'iTRAQ:127.124760',
    u'iTRAQ:127.131079', u'iTRAQ:128.128114',
    u'iTRAQ:128.134433', u'iTRAQ:129.131468',
    u'iTRAQ:129.137787', u'iTRAQ:130.134822',
    u'iTRAQ:130.141141', u'iTRAQ:131.138176')
control_1 = u'iTRAQ:126.127725'
control_2 = u'iTRAQ:131.138176'


def ReadFile(filename):
    data = pd.read_csv(filename, header=0, delimiter="\t")\
        .fillna(u'')
    data[ukey] = data[u'Seq'].astype('str') + data[u'Mod'].astype('str')\
        + data[u'ModDetail'].astype('str')
    return data


def DivideMergeFile(df):
    dfs = {}
    df[u'RawdataFile'] = df[u'RawdataFile'].apply(
        lambda x: fnamepat.sub("", x))
    li = df[u'RawdataFile'].drop_duplicates().values
    for fn in li:
        df_temp = df[df[u'RawdataFile'] == fn]
        dfs[fn] = df_temp
    return dfs


def RemoveDupli(df):
    df_unique = df.sort_values(by=u'PeptScore')\
        .drop_duplicates(ukey, keep='last')
    df_unique = df_unique.sort_index()
    return df_unique


def ppmPlotter(df):
    df[u'xy'] = (df[u'ObsMass'] - df[u'CalcMass']) / df[u'CalcMass'] * 1000000
    df = df.sort(columns=u'xy')
    xy = df[u'xy']
    plt.plot(xy, marker="o", markersize=4.0, linestyle="none")
    plt.show()


def RetTimeHist(df, ranges=(0, 110), bins=50, alpha=0.5, color='k'):
    x = df[u'RetTime'].values
    plt.hist(x, range=ranges, bins=bins, alpha=alpha, color=color)
    plt.xlabel('RetTime (min)')
    plt.ylabel('Peptide ID')
    plt.title('Elution Profile of Identified Peptides')
    plt.show()


def RatioCalculater(dfs, control=control_1):
    resultdfs = []

    for df in dfs.values():
        df_temp = pd.DataFrame()

        #   Remove peptides whose intensities of control reporter ions = 0
        df = df[df[control] != 0]

        #   Calculate log2 ratio of reporter ion intensities and
        #   normalize them by median values
        homedir = os.getcwd()
        if not os.path.isdir("hist"):
            os.mkdir("hist")
        os.chdir("hist")

        for cn in reporter_values:
            if cn != control:
                df_temp['RATIO: ' + cn + '/' + control] = np.log2(
                    df.loc[:, cn] / df.loc[:, control])\
                    - np.median(np.log2(df.loc[:, cn] / df.loc[:, control]))
                plt.hist(df_temp['RATIO: ' + cn + '/' + control].values,
                         bins=50, alpha=0.5)
                plt.savefig(cn + 'hist.png')

        os.chdir(homedir)
        #   Calculate average log2 ratio of each peptides
        df_temp = df_temp.join(df[ukey])
        df_temp = df_temp.replace([np.inf, -np.inf], np.nan)
        df_temp = df_temp.groupby([ukey]).mean()

        df = RemoveDupli(df)
        df = df.join(df_temp, on=ukey)

        resultdfs.append(df)

    return resultdfs


def MergeFiles(dfs):
    df_merge = RemoveDupli(pd.concat(dfs))
    df_merge = df_merge.reset_index(drop=True)
    return df_merge


def MergeFilesNoRemove(dfs):
    df_merge = pd.concat(dfs, ignore_index=True)
    df_merge = df_merge.reset_index(drop=True)
    return df_merge


def CombineQuantResults(dfs):
    dfs_temp = []
    i = 1
    #   Drop columns of reporter ion intensities and ratios
    for df in dfs:
        df_temp = df.ix[:, df.columns.map(
            lambda x:
            not(x.startswith('RATIO'))
            and not(x.startswith('iTRAQ')))]
        dfs_temp.append(df_temp)
    df_merge = MergeFiles(dfs_temp)

    #   Merge log2 ratio columns
    for df in dfs:
        df[u'Identified_' + str(i)] = pd.Series(
            True, index=df.index)
        df_merge = pd.merge(df_merge, df.ix[:, df.columns.map(
            lambda x: x.startswith('RATIO')
            or x.startswith(ukey)
            or x.startswith('Identified'))],
            on=ukey, how='left', suffixes=('_' + str(i-1), '_' + str(i)))
        i += 1
    df_merge = IDCounter(df_merge)
    return df_merge


def IDCounter(df):
    ind = []
    con = []
    s = df.ix[:, df.columns.map(lambda x: x.startswith('Identified'))]
    for i in range(len(s)):
        ind.append(i)
        con.append(s.ix[i].count())
    df[u'Count'] = pd.Series(con, index=ind)
    df = df.drop([x for x in df.columns if x.startswith('Identified')], axis=1)
    return df


def GetAPValue(df, count_criteria=3,
               control_1=control_1, control_2=control_2):
    #Get average and p-value of Welch's T-test
    df = df[df[u'Count'] >= count_criteria]
    df = df.reset_index(drop=True)
    c2_ratios = df.ix[:, df.columns.map(
        lambda x: x.startswith('RATIO: ' + control_2))]
    for cn in reporter_values:
        if cn != control_1:
            index = []
            average = []
            ttest = []
            tmp = df.ix[:, df.columns.map(
                lambda x: x.startswith('RATIO: ' + cn))]
            for i in df.index:
                index.append(i)
                average.append(np.nanmean(tmp.ix[i]))
                if cn != control_2:
                    ttest.append(st.ttest_ind(
                        tmp.ix[i][np.isfinite(tmp.ix[i])],
                        c2_ratios.ix[i][np.isfinite(c2_ratios.ix[i])],
                        equal_var=False)[1])
            df[u'AVERAGE: ' + cn] = pd.Series(average, index=index)
            if cn != control_2:
                df[u'TTEST: ' + cn] = pd.Series(ttest, index=index)
    return df


def ExtractUseColumns(df, necessaries=[]):
    basic_information = [u'AccNum', u'Seq',
                         u'Mod', u'ModDetail', u'SeqModModDetail']
    new_df = pd.DataFrame()
    b_series = df.ix[:, df.columns.isin(basic_information)]
    new_df = pd.concat([new_df, b_series], axis=1)
    if bool(necessaries) is True:
        for cn in necessaries:
            if not cn in basic_information:
                series = df.ix[:, df.columns.map(
                    lambda x: x.startswith(cn))]
                new_df = pd.concat([new_df, series], axis=1)
    return new_df


def VolcanoPlotter(df):
    for cn in reporter_values:
        if not cn in [control_1, control_2]:
            x = df.ix[:, df.columns.map(lambda x:
                      x.startswith(u'AVERAGE: ' + cn))].values
            y = -np.log10(df.ix[:, df.columns.map(lambda x:
                          x.startswith(u'TTEST: ' + cn))].values)
            plt.plot(x, y, marker="o", markersize=4.0, linestyle="none")

    return True


def CombinePepArea(dfs):
    dfs_temp = []
    i = 1
    #   Drop PepArea column
    for k, v in sorted(dfs.items()):
        df_temp = v.ix[:, v.columns.map(
            lambda x:
            not(x.startswith('PepArea')))]
        dfs_temp.append(df_temp)
    df_merge = MergeFiles(dfs_temp)

    #   Merge PepArea columns
    for k, v in sorted(dfs.items()):
        v[u'Identified_' + str(i)] = pd.Series(
            True, index=v.index)
        df_merge = pd.merge(df_merge, v.ix[:, v.columns.map(
            lambda x: x.startswith('PepArea')
            or x.startswith(ukey)
            or x.startswith('Identified'))],
            on=ukey, how='left', suffixes=('_' + str(i-1), '_' + str(i)))
        i += 1
    df_merge = IDCounter(df_merge)
    return df_merge


def Hist2D(df):
    df = df[df[u'Count'] == 2]
    stdline = range(10000000)
    x = np.log2(df[u'PepArea_1'].values)
    y = np.log2(df[u'PepArea_2'].values)
    plt.hist2d(x, y, bins=100, range=[[10, 30], [10, 30]])
    plt.plot(stdline, stdline, color='white')
    plt.show()


def TMTlabelEfficiency(df):
    df[u'number_of_TMT'] = 0
    df[u'number_of_S_TMT'] = 0
    df[u'number_of_T_TMT'] = 0
    df[u'number_of_Y_TMT'] = 0
    df[u'number_of_H_TMT'] = 0
    df[u'efficiency'] = "Partial"
    df.ix[df[u'Mod'].str.contains
          ('3\sTMT\s\(K\)|3\sTMT6plex\s\(K\)'), u'number_of_TMT'] = 2
    df.ix[df[u'Mod'].str.contains
          ('2\sTMT\s\(K\)|2\sTMT6plex\s\(K\)'), u'number_of_TMT'] = 1
    df.ix[df[u'Mod'].str.contains
          ('TMT\s\(K\)|TMT6plex\s\(K\)'), u'number_of_TMT'] += 1
    df.ix[df[u'Mod'].str.contains
          ('TMT\s\(N-term\)|TMT6plex\s\(N-term\)'), u'number_of_TMT'] += 1
    df.ix[df[u'Mod'].str.contains('TMT\s\(S\)'), u'number_of_S_TMT'] = 1
    df.ix[df[u'Mod'].str.contains('TMT\s\(T\)'), u'number_of_T_TMT'] = 1
    df.ix[df[u'Mod'].str.contains('TMT\s\(Y\)'), u'number_of_Y_TMT'] = 1
    df.ix[df[u'Mod'].str.contains('TMT\s\(H\)'), u'number_of_H_TMT'] = 1

    df.ix[df.apply(lambda x: x[u'number_of_k'] + 1 == x[u'number_of_TMT'],
          axis=1), u'efficiency'] = 'Full'
    df.ix[df.apply(lambda x: x[u'number_of_TMT'] == 0,
          axis=1), u'efficiency'] = 'None'

    return df


def GuanidinationEfficiency(df):
    df[u'number_of_Guanidination'] = 0
    df[u'number_of_K_Guanidination'] = 0
    df[u'number_of_Nterm_Guanidination'] = 0

    df[u'Guanidination_efficiency'] = "Partial"
    df.ix[df[u'Mod'].str.contains
          ('3\sGuanidinyl\s\(K\)'), u'number_of_K_Guanidination'] = 2
    df.ix[df[u'Mod'].str.contains
          ('2\sGuanidinyl\s\(K\)'), u'number_of_K_Guanidination'] = 1
    df.ix[df[u'Mod'].str.contains
          ('Guanidinyl\s\(K\)'), u'number_of_K_Guanidination'] += 1
    df.ix[df[u'Mod'].str.contains
          ('Guanidinyl\s\(N-term\)'), u'number_of_Nterm_Guanidination'] += 1
    df[u'number_of_Guanidination'] = df[u'number_of_K_Guanidination']\
        + df[u'number_of_Nterm_Guanidination']
    df.ix[df.apply(lambda x: x[u'number_of_K_Guanidination'] == 0,
          axis=1), u'Guanidination_efficiency'] = 'None'
    df.ix[df.apply(lambda x: x[u'number_of_k']
          == x[u'number_of_K_Guanidination'], axis=1),
          u'Guanidination_efficiency'] = 'Full'
    df.ix[df.apply(lambda x: x[u'number_of_Nterm_Guanidination'] != 0,
          axis=1), u'Guanidination_efficiency'] = 'Overlabel'

    return df


def getSeqProperties(df):
    df = AminoAcidCounter(df)
    df[u'number_of_phospho'] = df[u'ModDetail'].str.count('79.9663')
    df[u'number_of_acidic_residue'] = df[u'number_of_d'].values\
        + df[u'number_of_e'].values
    df[u'number_of_basic_residue'] = df[u'number_of_h'].values\
        + df[u'number_of_k'].values + df[u'number_of_r'].values
    df[u'peptide_length'] = df[u'Seq'].str.len()

    return df


def AminoAcidCounter(df):
    df[u'number_of_a'] = [Counter(x)['A'] for x in df[u'Seq'].values]
    df[u'number_of_r'] = [Counter(x)['R'] for x in df[u'Seq'].values]
    df[u'number_of_n'] = [Counter(x)['N'] for x in df[u'Seq'].values]
    df[u'number_of_d'] = [Counter(x)['D'] for x in df[u'Seq'].values]
    df[u'number_of_c'] = [Counter(x)['C'] for x in df[u'Seq'].values]
    df[u'number_of_q'] = [Counter(x)['Q'] for x in df[u'Seq'].values]
    df[u'number_of_e'] = [Counter(x)['E'] for x in df[u'Seq'].values]
    df[u'number_of_g'] = [Counter(x)['G'] for x in df[u'Seq'].values]
    df[u'number_of_h'] = [Counter(x)['H'] for x in df[u'Seq'].values]
    df[u'number_of_i'] = [Counter(x)['I'] for x in df[u'Seq'].values]
    df[u'number_of_l'] = [Counter(x)['L'] for x in df[u'Seq'].values]
    df[u'number_of_k'] = [Counter(x)['K'] for x in df[u'Seq'].values]
    df[u'number_of_m'] = [Counter(x)['M'] for x in df[u'Seq'].values]
    df[u'number_of_f'] = [Counter(x)['F'] for x in df[u'Seq'].values]
    df[u'number_of_p'] = [Counter(x)['P'] for x in df[u'Seq'].values]
    df[u'number_of_s'] = [Counter(x)['S'] for x in df[u'Seq'].values]
    df[u'number_of_t'] = [Counter(x)['T'] for x in df[u'Seq'].values]
    df[u'number_of_w'] = [Counter(x)['W'] for x in df[u'Seq'].values]
    df[u'number_of_y'] = [Counter(x)['Y'] for x in df[u'Seq'].values]
    df[u'number_of_v'] = [Counter(x)['V'] for x in df[u'Seq'].values]

    return df


def GRAVYScoreCalculator(df):
    HYDROPATHY_INDEX_OF_AA = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8,
                              'C': 2.5, 'M': 1.9, 'A': 1.8, 'G': -0.4,
                              'T': -0.7, 'W': -0.9, 'S': -0.8, 'Y': -1.3,
                              'P': -1.6, 'H': -3.2, 'E': -3.5, 'Q': -3.5,
                              'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5}
    df = AminoAcidCounter(df)
    df[u'Peptide_length'] = df[u'Seq'].str.len()
    df[u'GRAVY_score'] = (df[u'number_of_a'] * HYDROPATHY_INDEX_OF_AA['A'] +
                          df[u'number_of_r'] * HYDROPATHY_INDEX_OF_AA['R'] +
                          df[u'number_of_n'] * HYDROPATHY_INDEX_OF_AA['N'] +
                          df[u'number_of_d'] * HYDROPATHY_INDEX_OF_AA['D'] +
                          df[u'number_of_c'] * HYDROPATHY_INDEX_OF_AA['C'] +
                          df[u'number_of_q'] * HYDROPATHY_INDEX_OF_AA['Q'] +
                          df[u'number_of_e'] * HYDROPATHY_INDEX_OF_AA['E'] +
                          df[u'number_of_g'] * HYDROPATHY_INDEX_OF_AA['G'] +
                          df[u'number_of_h'] * HYDROPATHY_INDEX_OF_AA['H'] +
                          df[u'number_of_i'] * HYDROPATHY_INDEX_OF_AA['I'] +
                          df[u'number_of_l'] * HYDROPATHY_INDEX_OF_AA['L'] +
                          df[u'number_of_k'] * HYDROPATHY_INDEX_OF_AA['K'] +
                          df[u'number_of_m'] * HYDROPATHY_INDEX_OF_AA['M'] +
                          df[u'number_of_f'] * HYDROPATHY_INDEX_OF_AA['F'] +
                          df[u'number_of_p'] * HYDROPATHY_INDEX_OF_AA['P'] +
                          df[u'number_of_s'] * HYDROPATHY_INDEX_OF_AA['S'] +
                          df[u'number_of_t'] * HYDROPATHY_INDEX_OF_AA['T'] +
                          df[u'number_of_w'] * HYDROPATHY_INDEX_OF_AA['W'] +
                          df[u'number_of_y'] * HYDROPATHY_INDEX_OF_AA['Y'] +
                          df[u'number_of_v'] * HYDROPATHY_INDEX_OF_AA['V'])\
        / df[u'Peptide_length']

    return df


def pICalculator(df):
    pH_temp = []

    for seq in df[u'Seq']:
        pH = 6.5
        left_bound = 0.0
        right_bound = 14.0
        epsilon = 1
        charge = getCharge(seq)
        while epsilon > 0.01:
            if charge < 0:
                #drop it by half
                new_pH = (pH-left_bound)/2.0
                right_bound = pH
            else:
                #increase it by half
                new_pH = pH+(right_bound-pH)/2.0
                left_bound = pH
            epsilon = abs(pH-new_pH)
            charge = getCharge(seq, pH=new_pH)
            pH = new_pH
        pH_temp.append(pH)
    df[u'pI'] = pH_temp
    return df


def getCharge(seq, pH=7.0):
    PI_RESIDUE_CHARGE = {'D': (3.65, -1.0),
                         'E': (4.25, -1.0),
                         'C': (8.18, -1.0),
                         'Y': (10.07, -1.0),
                         'H': (6.00, 1.0),
                         'K': (10.53, 1.0),
                         'R': (12.48, 1.0)}

    d = PI_RESIDUE_CHARGE
    pH = float(pH)
    s = seq.upper()
    l = [(i, float(s.count(i))) for i in set(list(s))]
    #do n-term/c-term first
    charge = (1.0 / (1.0 + 10.0 ** (pH - 9.69)))
    charge += (-1.0 / (1.0 + 10.0 ** (-1.0 * (pH - 2.34))))
    for aa, count in l:
        v = d.get(aa, (0.0, 0.0))
        charge += (v[1] * count) / (1.0 + 10.0 ** (v[1] * (pH - v[0])))
    return charge


def extractPhospho(dfs):
    dfs = dfs[dfs[u'Mod'].str.contains('Phospho')]
    return dfs


def selectData(dfs, data):
    new_dfs = {}
    for dataname in data:
        new_dfs[dataname] = dfs[dataname]

    return new_dfs


def divideTMTnoTMT(dfs):
    new_dfs = {}
    for k, v in sorted(dfs.items()):
        TMT_key = 'TMT_' + str(k)
        noTMT_key = 'noTMT_' + str(k)
        TMT_df = v[v[u'efficiency'] != 'None']
        noTMT_df = v[v[u'efficiency'] == 'None']
        new_dfs[TMT_key] = TMT_df
        new_dfs[noTMT_key] = noTMT_df

    return new_dfs


def calculateCV(df, tagname='TMT10plex'):
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


def reporterIntCal(df):
    abRatio = np.matrix([[1,0,0.006,0,0,0,0,0,0,0],
                         [0,1,0,0.009,0,0,0,0,0,0],
                         [0.05,0,1,0,0.005,0,0,0,0,0],
                         [0,0.054,0,1,0,0.007,0,0,0,0],
                         [0,0,0.064,0,1,0,0.013,0,0,0],
                         [0,0,0,0.047,0,1,0,0.013,0,0],
                         [0,0,0,0,0.032,0,1,0,0.015,0],
                         [0,0,0,0.002,0,0.033,0,1,0,0.021],
                         [0,0,0,0,0,0,0.029,0,1,0],
                         [0,0,0,0,0,0,0,0.03,0,1]])
    reporter_values = (
        u'iTRAQ:126.127725', u'iTRAQ:127.124760',
        u'iTRAQ:127.131079', u'iTRAQ:128.128114',
        u'iTRAQ:128.134433', u'iTRAQ:129.131468',
        u'iTRAQ:129.137787', u'iTRAQ:130.134822',
        u'iTRAQ:130.141141', u'iTRAQ:131.138176')

    intensities = df.ix[:, df.columns.isin(reporter_values)]
    for i in intensities.index():
        X = np.linalg.solve(abRatio, np.matrix(intensities.i))
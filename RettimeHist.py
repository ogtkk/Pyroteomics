import Functions as func
import re

filename = '151009ko-SAS.txt'
data = func.ReadFile(filename)
data = func.TMTlabelEfficiency(data)
dfs = func.DivideMergeFile(data)
pattern = re.compile(u"D.+ManualInputFile\\\\")
dfs = dict([(pattern.sub("", k), v) for k, v in dfs.items()])

func.RetTimeHist(dfs['151009ko08-E2_200.raw'])


    dic['acidic_0~1'].append(len(v[v[u'number_of_d'] + v[u'number_of_e'] <= 1]))
    dic['acidic_2~3'].append(len(v[v[u'number_of_d'] + v[u'number_of_e'] >= 2
                             and v[u'number_of_d'] + v[u'number_of_e'] <= 3]))
    dic['acidic_4~5'].append(len(v[v[u'number_of_d'] + v[u'number_of_e'] == 2]))
    dic['acidic_6~7'].append(len(v[v[u'number_of_d'] + v[u'number_of_e'] == 3]))
    dic['acidic_8~9'].append(len(v[v[u'number_of_d'] + v[u'number_of_e'] == 4]))
    dic['acidic_>9'].append(len(v[v[u'number_of_d'] + v[u'number_of_e'] == 5]))
    dic['acidic_6'].append(len(v[v[u'number_of_d'] + v[u'number_of_e'] == 6]))
    dic['acidic_>6'].append(len(v[v[u'number_of_d'] + v[u'number_of_e'] > 6]))
    dic['acidic_0_area'].append(np.sum(v[v[u'number_of_d'] + v[u'number_of_e'] == 0][u'PepArea'].values))
    dic['acidic_1_area'].append(np.sum(v[v[u'number_of_d'] + v[u'number_of_e'] == 1][u'PepArea'].values))
    dic['acidic_2_area'].append(np.sum(v[v[u'number_of_d'] + v[u'number_of_e'] == 2][u'PepArea'].values))
    dic['acidic_3_area'].append(np.sum(v[v[u'number_of_d'] + v[u'number_of_e'] == 3][u'PepArea'].values))
    dic['acidic_4_area'].append(np.sum(v[v[u'number_of_d'] + v[u'number_of_e'] == 4][u'PepArea'].values))
    dic['acidic_5_area'].append(np.sum(v[v[u'number_of_d'] + v[u'number_of_e'] == 5][u'PepArea'].values))
    dic['acidic_6_area'].append(np.sum(v[v[u'number_of_d'] + v[u'number_of_e'] == 6][u'PepArea'].values))
    dic['acidic_>6_area'].append(np.sum(v[v[u'number_of_d'] + v[u'number_of_e'] > 6][u'PepArea'].values))
    func1 = lambda x: x * slope + intercept
    line1 = mpl.lines.Line2D([0, 10], [func1(0), func1(10)], color='r')
    ax1.add_line(line1)
    ax4 = fig.add_subplot(332, frame_on=False)
    ax5 = fig.add_subplot(336, frame_on=False)
    ax6 = fig.add_subplot(333, frame_on=False)
    ax4.set_xticks([])
    ax4.set_yticks([])
    ax5.set_xticks([])
    ax5.set_yticks([])
    ax6.set_xticks([])
    ax6.set_yticks([])

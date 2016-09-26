import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st
import os
from collections import Counter
import re

filename = '160319ko-TMT.txt'
os.chdir('../')
os.chdir('Python_script')


data = pd.read_csv(filename, header=0, delimiter="\t").fillna(u'')
fnamepat = re.compile(u"D.+ManualInputFile\\\\|D.+ManualInputFile/")
fnameext = re.compile(u".raw")
data['RawdataFile'] = data['RawdataFile'].apply(lambda x: fnamepat.sub("", x))
data['DatFile'] = data['RawdataFile'].apply(lambda x: fnameext.sub(".dat", x))
AAmss = {}
obsmz = {}
charge = {}
peakLs = {}

datfiles = data['DatFile'].drop_duplicates()
f = open(datfiles[0], 'r')
lines = f.readlines()
flagmss = 0
i = 1
flag_que = 0
for line in lines:
    if re.compile('IT_MODS\=.+').search(line):
        line = re.compile('IT_MODS\=').sub("", line)
        mods = line.split(',')
        for mod in mods:
            if re.compile('Phospho \(ST\)').search(mod):
                AAmss['delta_pST'] = 'delta' + str(i)
            elif re.compile('Phospho \(Y\)').search(mod):
                AAmss['delta_pY'] = 'delta' + str(i)
            elif re.compile('\(K\)').search(mod):
                AAmss['delta_K'] = 'delta' + str(i)
            elif re.compile('\(R\)').search(mod):
                AAmss['delta_R'] = 'delta' + str(i)
            elif re.compile('\(N-term\)').search(mod):
                AAmss['delta_Nterm'] = 'delta' + str(i)
            i += 1
    if re.compile('.+name\="unimod"').search(line):
        flagmss = 0
    if re.compile('name\="masses"').search(line):
        flagmss = 1
    elif flagmss == 1:
        temp_ar = re.match('([A-z0-9_]+)\=([0-9.]+)', line)
        if temp_ar:
            AAmss[temp_ar.groups()[0]] = float(temp_ar.groups()[1])

    temp_ar = re.search('qexp([0-9]+)\=([0-9.]+)\,([0-9]+)', line)
    if temp_ar:
        que = temp_ar.groups()[0]
        obsmz[que] = float(temp_ar.groups()[1])
        charge[que] = float(temp_ar.groups()[2])
    temp_ar = re.search('name\="query([0-9]+)"', line)
    if temp_ar:
        que = temp_ar.groups()[0]
        flag_que = 1
    elif flag_que == 1:
        temp_ar = re.search('Ions1\=([0-9e+.,:]+)', line)
        if temp_ar:
            peakLs[que] = temp_ar.groups()[0]
            flag_que = 0
AAmss['pS'] = AAmss['S'] + AAmss[AAmss['delta_pST']]
AAmss['pT'] = AAmss['T'] + AAmss[AAmss['delta_pST']]
AAmss['pY'] = AAmss['Y'] + AAmss[AAmss['delta_pY']]
AAmss['oM'] = AAmss['M'] + AAmss['Oxygen']
AAmss['tmtK'] = AAmss['K'] + AAmss[AAmss['delta_K']]
AAmss['tmtNterm'] = AAmss[AAmss['delta_Nterm']]

if AAmss['delta1'] == 0.997040:
    AAmss['S'] = AAmss['S'] + AAmss['delta1']
    AAmss['T'] = AAmss['T'] + AAmss['delta1']
    AAmss['Y'] = AAmss['Y'] + AAmss['delta1']

ic_b = float(AAmss['N_term']) - float(AAmss['Electron'])
ic_bb = float(AAmss['N_term']) / 2 + float(AAmss['Hydrogen']) / 2 - float(AAmss['Electron'])
ic_y = float(AAmss['C_term']) + float(AAmss['Hydrogen']) * 2 - float(AAmss['Electron'])
ic_yy = float(AAmss['C_term']) / 2 + float(AAmss['Hydrogen']) * 3 / 2 - float(AAmss['Electron'])


for i in data.index:
    flag_bNL = 0
    flag_yNL = 0
    match_bs = []
    match_ys = []
    delta_bs = {}
    delta_ys = {}
    int_bs = {}
    int_ys = {}
    o_ions = []
    o_ints = []
    dicionint = {}
    temp_ar = re.search('\&query\=([0-9]+)', data.ix[i, 'PeptsURL'])
    if temp_ar:
        dat = data.ix[i, 'DatFile']
        que = temp_ar.groups()[0]

    ponl = obsmz[que] - AAmss['NeutralLoss2'] / charge[que]
    seq = data.ix[i, 'Seq']
    mod = data.ix[i, 'ModDetail']
    peakL = peakLs[que]
    mzints = peakL.split(',')
    for mzint in mzints:
        ionint = mzint.split(':')
        o_ions.append(float(ionint[0]))
        o_ints.append(float(ionint[1]))
        dicionint[float(ionint[0])] = float(ionint[1])

#    o_ints = sort [a <=> b] o_ints
#    o_ints = reverse(o_ints)
#    for o_ino in o_ions:
#        if hashionint[_] >= o_ints[9]:
#            if abs(ponl - _) < 0.02:
#                hashponl[que] = '+'
    aa = list(seq)
    length = len(seq)
    moddetails = mod.split(',')
    for moddetail in moddetails:
        temp_ar = re.search('([0-9.\-]+)\@S\:|([0-9.\-]+)\@T\:|([0-9.\-]+)\@Y\:|([0-9.\-]+)\@M\:',
                            moddetail)
        if temp_ar:
            modmz = temp_ar.groups()[0]
        temp_arS = re.search('\@S\:([0-9]+)', moddetail)
        temp_arT = re.search('\@T\:([0-9]+)', moddetail)
        temp_arY = re.search('\@Y\:([0-9]+)', moddetail)
        temp_arM = re.search('\@M\:([0-9]+)', moddetail)

        if temp_arS:
            aa[int(temp_arS.group(1)) - 1] = 'pS'
        elif temp_arT:
            aa[int(temp_arT.group(1)) - 1] = 'pT'
        elif temp_arY:
            aa[int(temp_arY.group(1)) - 1] = 'pY'
        elif temp_arM:
            aa[int(temp_arM.group(1)) - 1] = 'oM'

    for k in range(length):
        if k == 0:
            mw_b = ic_b + AAmss[aa[0]]
            mw_bb = ic_bb + AAmss[aa[0]] / 2
            mw_y = ic_y + AAmss[aa[length - (k + 1)]]
            mw_yy = ic_yy + AAmss[aa[length - (k + 1)]] / 2
            if re.compile('TMT.+\(N-term\)').search(mod):
                mw_b = mw_b + AAmss['tmtNterm']
                mw_bb = mw_bb + AAmss['tmtNterm'] / 2
        else:
            mw_b += AAmss[aa[k]]
            mw_bb += AAmss[aa[k]] / 2
            mw_y += AAmss[aa[length - (k + 1)]]
            mw_yy += AAmss[aa[length - (k + 1)]] / 2
        if aa[k] == "pS" or aa[k] == "pT":
            flag_bNL += 1
        if (aa[length - (k + 1)] == "pS" or aa[length - (k + 1)] == "pT"):
            flag_yNL += 1

        for obsi, intensity in dicionint.items():
            sub_b = abs(mw_b - obsi)
            sub_bb = abs(mw_bb - obsi)
            sub_y = abs(mw_y - obsi)
            sub_yy = abs(mw_yy - obsi)
            if sub_b < 0.02:
                b_ion = 'b(%d)+' % (k + 1)
                match_bs.append(b_ion)
                delta_bs[b_ion] = int((mw_b - obsi) * 100) / 100
                int_bs[b_ion] = intensity
            if sub_bb < 0.02:
                bb_ion = 'b(%d)++' % (k + 1)
                match_bs.append(bb_ion)
                delta_bs[bb_ion] = int((mw_bb - obsi) * 100) / 100
                int_bs[bb_ion] = intensity
            if sub_y < 0.02:
                y_ion = 'y(%d)+' % (k + 1)
                match_ys.append(y_ion)
                delta_ys[y_ion] = int((mw_y - obsi) * 100) / 100
                int_ys[y_ion] = intensity
            if sub_yy < 0.02:
                yy_ion = 'y(%d)++' % (k + 1)
                match_ys.append(yy_ion)
                delta_ys[yy_ion] = int((mw_yy - obsi) * 100) / 100
                int_ys[yy_ion] = intensity
        data.ix[i, 'y_ions'] = ",".join(k+": "+str(v) for k, v in int_ys.items())
        data.ix[i, 'b_ions'] = ",".join(k+": "+str(v) for k, v in int_bs.items())
data.to_csv('test.csv')

            for i in range(flag_bNL):
                i += 1
                mz_nl = 98 * i
                if (abs (mw_b -AAmss[NeutralLoss2]*i- obsi) < 0.02)[
                    bnl_ion = "b\(k\)\-mz_nl"
                    push (match_b, bnl_ion)
                    hash_delta_b[bnl_ion] = int ((mw_b -AAmss[NeutralLoss2]*i- obsi) * 100) / 100
                ]
                if (abs (mw_bb -((AAmss[NeutralLoss2]*i)/2)- obsi) < 0.02)[
                    bbnl_ion = "b\(k\)\-mz_nl\+\+"
                    push (match_b, bbnl_ion)
                    hash_delta_b[bbnl_ion] = int ((mw_bb -((AAmss[NeutralLoss2]*i)/2)- obsi) * 100) / 100
                ]
            ]
            for (i = 1 i<=flag_yNL i++)[
                mz_nl = 98 * i
                if (abs (mw_y -AAmss[NeutralLoss2]*i- obsi) < 0.02)[
                    ynl_ion = "y\(k\)\-mz_nl"
                    push (match_y, ynl_ion)
                    hash_delta_y[ynl_ion] = int ((mw_y -AAmss[NeutralLoss2]*i- obsi) * 100) / 100
                ]
                if (abs (mw_yy -((AAmss[NeutralLoss2]*i)/2)- obsi) < 0.02)[
                    yynl_ion = "y\(k\)\-mz_nl\+\+"
                    push (match_y, yynl_ion)
                    hash_delta_y[yynl_ion] = int ((mw_yy -((AAmss[NeutralLoss2]*i)/2)- obsi) * 100) / 100
                ]
            ]
#    undef (%tmp)
#    match_b = grep(!tmp[_]++, match_b)
#    match_y = grep(!tmp[_]++, match_y)
#    foreach b_ion (match_b)[
#        delta_b = hash_delta_b[b_ion]
#        push (match_delta_b, delta_b)
#    ]
#    foreach y_ion (match_y)[
#        delta_y = hash_delta_y[y_ion]
#        push (match_delta_y, delta_y)
#    ]
#    hash_mbi[que] = join("\, ", match_b)
#    hash_myi[que] = join("\, ", match_y)
#    hash_dbi[que] = join("\, ", match_delta_b)
#    hash_dyi[que] = join("\, ", match_delta_y)
#]
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
################################################
#
#opendir TXTDIR, txtdir
#dlfiles = readdir TXTDIR
#
#foreach (dlfiles)[
#    if (/([0-9]+)/)[
#        entry = 1
#    ]
#    url = data[entry][url]
#    if (url =~ m|file\=\.\.(\/data\/[0-9]+\/F[0-9]+\.dat)\&query\=([0-9]+)|) [
#        dat = 1
#        dat = "//mascot2" . dat
#        que = 2
#    ]
#    data[entry][b] = hash_mbi[que]
#    data[entry][y] = hash_myi[que]
#    data[entry][db] = hash_dbi[que]
#    data[entry][dy] = hash_dyi[que]
#    data[entry][ponl] = hashponl[que]
#
#    undef(legend)
#    undef(tlegend)
#    flag_legend = 0
#    flag_tlegend = 0
#
#    open PV, "txtdir/_"
#
#    while(<PV>)[
#        if(/Fragmentation/)[
#            flag_legend = 1
#        ]elsif(/Data file/)[
#            flag_legend = 0
#        ]
#
#        if(/Monoisotopic/)[
#            flag_tlegend = 1
#        ]elsif(/TABLE BORDER/)[
#            flag_tlegend = 0
#        ]
#
#        if(flag_legend)[
#            legend .= _
#        ]
#        if(flag_tlegend)[
#            tlegend = tlegend . _ . "<Br>"
#            tlegend =~ s|\&nbsp\\&nbsp\\&nbsp\\(\<a href\=\"\.\.\/help\/results_help\.html\#PEP\" target\=_blank\>help\<\/a\>\)||
#        ]
#
#        if(/^MS\/MS Fragmentation of \<B\>\<FONT COLOR\=\#FF0000\>([^\<]+)\<\/FONT\>\<\/B\>\<BR\>/)[
#            seq[entry] = 1
#        ]
#    ]
#    close (PV)
#
#    data[entry][legend]  = legend
#    data[entry][tlegend] = tlegend
#
#]
#
#sortkeys = keys (%data)
#sortkeys = sort [a <=> b] sortkeys
#
#for(sortkeys)[
#    print REPORT "Peptide No._\n"
#    print REPORT "seq[_]\n"
#    print REPORT "Confirmed sites: data[_][conf]\n", "Ambiguous sites: data[_][amb]\n"
#    tlegend = data[_][tlegend]
#    tlegend =~ s/\<[^\>\<]+\>//gi
#    chomp(tlegend)
#    chomp(tlegend)
#    chomp(tlegend)
#    print REPORT "tlegend\n"
#    print REPORT "Matched b ions: data[_][b]\n"
#    print REPORT "Matched y ions: data[_][y]\n"
#    print REPORT "delta mass of b ions: data[_][db]\n"
#    print REPORT "delta mass of y ions: data[_][dy]\n"
#    print REPORT "\n"
#]
#
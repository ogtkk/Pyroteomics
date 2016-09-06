import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st
import os
from collections import Counter
import re

filename = '160319ko-TMT.txt'
data = pd.read_csv(filename, header=0, delimiter="\t").fillna(u'')
fnamepat = re.compile(u"D.+ManualInputFile\\\\|D.+ManualInputFile/")
fnameext = re.compile(u".raw")
data['RawdataFile'] = data['RawdataFile'].apply(lambda x: fnamepat.sub("", x))
data['DatFile'] = data['RawdataFile'].apply(lambda x: fnameext.sub(".dat", x))

datfiles = data['DatFile'].drop_duplicates()
f = open(datfiles[0], 'r')
lines = f.readlines()
for line in lines:
    flagmss = 0
    if re.compile('IT_MODS\=.+').match(line):
        line = re.compile('IT_MODS\=').sub("", line)
        mods = line.split(',')
        for mod in mods:
            if ('Phospho (ST)'):
                i = 1
            elif ('Phospho \(Y\)'):
                i = 1
            elif ('\(K\)'):
                i = 1
            elif ('\(R\)'):
                i = 1
    if re.compile('.+name\="masses"').match(line):
        flagmss = 1
    elif flagmss == 1:
        if re.compile('([A-Za-z0-9_]+)\=([0-9.\-]+)/)').match(line):
            AAmss[re.compile('=.+').sub('', line)] = re.compile('.+=').sub('', line)
        
    if re.compile('.+name\="header"').match(line):
        flagmss = 0;

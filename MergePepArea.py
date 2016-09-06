import Functions as func
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotter as plot
import os
import re
from datetime import datetime as dt

filename = '160403ko-tipTMT-full.txt'
p_AreaCol = ['PepIntensity', 'PepArea',
             'PepHalfWidth', 'PepSymmetryFacter',
             'PepTargetMz', 'PepRtimeStart',
             'PepRtimePeak', 'PepRtimeEnd',
             'PepSmoothCount', 'PepPeakCount',
             'PepTicRtime', 'PepTicInt',
             'PepTicIntRaw', 'PepPreInfoMz']

data = func.ReadFile(filename)
data['Key'] = data[u'no'].astype('str') + data[u'RawdataFile'].astype('str')
filelist = data.RawdataFile.drop_duplicates().values
os.chdir('160403ko-B1-XcaliburPepQuant1')
for fn in filelist:
    fn_temp = re.compile(u"D.+ManualInputFile\\\\|D.+ManualInputFile/").sub("", fn)
    fn_temp = re.compile(".raw").sub(".txt", fn_temp)
    peakareafile = func.ReadFile(fn_temp)
    peakareafile['Key'] = peakareafile[u'no'].astype('str') + peakareafile[u'RawdataFile'].astype('str')
    data.loc[data.Key.isin(peakareafile.Key), 'PepArea'] = peakareafile['PepArea']
data.to_csv('mergefile_A.csv')
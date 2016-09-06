
import Functions as func
import os
from datetime import datetime as dt

time = dt.now()
f_time = time.strftime('%Y%m%d%H%M%S')
filenames = ['Y1.txt', 'Y2.txt', 'Y3.txt', 'Y4.txt']
dfs = map(func.ReadFile, filenames)
dfs = dict(zip(filenames, dfs))
if not os.path.isdir(f_time + "_Quant"):
    os.mkdir(f_time + "_Quant")
os.chdir(f_time + "_Quant")
homedir = os.getcwd()
dfs = func.RatioCalculater(dfs)
dic = dict(zip(filenames, dfs))

if not os.path.isdir("Raw_RATIO"):
    os.mkdir("Raw_RATIO")
os.chdir("Raw_RATIO")
for k, v in dic.items():
    v.to_csv(k + '_RATIO.csv')
os.chdir(homedir)

df = func.CombineQuantResults(dfs)
df = func.GetAPValue(df)

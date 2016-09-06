import Functions as func
import re
filename = '151216ko-tipTMT.txt'
data = func.ReadFile(filename)
data = func.TMTlabelEfficiency(data)
dfs = func.DivideMergeFile(data)
pattern = re.compile(u"D.+ManualInputFile\\\\")

dfs = dict([(pattern.sub("", k), func.RemoveDupli(v)) for k, v in dfs.items()])
dfs = [dfs['151216ko08-sol.raw'], dfs['151216ko03-B.raw']]
df = func.CombinePepArea(dfs)
func.Hist2D(df)

import pandas as pd
import re
protPat = re.compile(u">sp\|.*\|")

f = open('ALBU-TRFE-LACB_BOVIN.fasta', 'r')
lines = f.readlines()
f.close()
prot = {}
aaCount = {}
aaArr = {}

for line in lines:
    if line.startswith(">"):
        n = line.find(" ")
        Accnum = line[:n]
        Accnum = protPat.sub("", Accnum)
        aaArr[Accnum] = []
    if not line.startswith(">"):
        aaArr[Accnum].extend(line[:-1])

data = pd.read_csv('160722ko-Set1-1.txt', header=0, delimiter="\t")\
    .fillna(u'')
def getResidueNum(df):
    for index in df.index:
        df.ix[index, "ResStartsAt"] = ''.join(aaArr[df.ix[index, 'AccNum']]).find(df.ix[index, 'Seq']) + 1

getResidueNum(data)
#data["ResStartsAt"] = ''.join(aaArr['ALBU_BOVIN']).find("ECCHGDLLECADDRADLAK") + 1
#data["ResEndsAt"] = data["ResStartsAt"] + len("ECCHGDLLECADDRADLAK")

data.to_csv('res.csv')


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st


data = pd.read_csv('summary_A-P.txt', header=0, delimiter="\t")
print data.ix[:, data.dtypes == np.float64]

import numpy as np
import pandas as pd
from pandas_plink import read_plink1_bin, write_plink1_bin
import random

dp = 'data/example_tiny/HumanOrigins249_tiny'
G = read_plink1_bin(dp+'.bed', verbose=False)

cols = ['iid', 'fid']
m = pd.DataFrame({x: list(G[x].values) for x in cols})
#m = m[m['iid'].isin(random.sample(list(m['iid']), 50))]
#m = m.sample(n=50, random_state=42, axis=0).sort_index()

# no, select random families, filter to them,
# then take sample of a slightly smaller size to ensure a varying number of inds per family
random.seed(42)
fams_all = list(np.unique(list(m['fid'])))
fam_sample = random.sample(fams_all, 10)
m = m[m['fid'].isin(fam_sample)]
m = m.sample(n=45, random_state=42, axis=0).sort_index()

print('Keeping %d families and %d individuals' % (len(fams_all), m.shape[0]))
G1 = G.sel(sample = list(m['iid']))
write_plink1_bin(G1, 'test/pqconv_testdata/HumanOrigins%d_tiny.bed' % m.shape[0])

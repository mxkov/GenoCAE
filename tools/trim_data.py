'''Trim the data in example_tiny (EIGENSTRAT only).

Usage:
  trim_data.py --dataprefix=<name> --outprefix=<name> [-m <num> -n <num> -k <num>]

Options:
  -h --help            show this screen
  --dataprefix=<name>  path and prefix of the data in the EIGENSTRAT format. if not absolute: assumed relative to GenoCAE/ directory
  --outprefix=<name>   path and prefix for the output data. if not absolute: assumed relative to GenoCAE/ directory
  -m <num>             number of SNPs to cut to. randomly sampled without replacement, order preserved. DEFAULT: all SNPs
  -n <num>             number of samples per population to cut to. DEFAULT: all samples
  -k <num>             if specified, cut the TOTAL number of samples further down to k. DEFAULT: not specified

'''

from docopt import docopt, DocoptExit
import os
import pandas as pd
from pathlib import Path
from random import seed, sample
from shutil import copy2
from tabulate import tabulate


def reportOutfile(filepath, margin = ''):
	print('{}File written: {}'.format(margin, filepath))


## Parsing args

try:
	arguments = docopt(__doc__, version='DataTrim 0.2.0')
except DocoptExit:
	print('Invalid command. Run `python3 trim_data.py --help` for more information.')
	exit(1)
for x in list(arguments.keys()):
	xnew = x.split('-')[-1]
	arg = arguments.pop(x)
	arguments[xnew] = arg

argtypes = {'dataprefix': 'path', 'outprefix': 'path',
            'm': 'num', 'n': 'num', 'k': 'num'}
margin1 = '\n'
for arg in argtypes.keys():
	if arguments[arg]:
		try:
			if argtypes[arg] == 'num':
				exec('{} = {}'.format(arg, int(arguments[arg])))
			elif argtypes[arg] == 'path':
				exec('{} = "{}"'.format(arg, arguments[arg]))
		except:
			print('{}Warning: Invalid argument {}'.format(margin1, arg))
			margin1 = ''
			exec('{} = None'.format(arg))
	else:
		exec('{} = None'.format(arg))

GCAE_DIR = Path(__file__).resolve().parents[1]
print('\nProject directory: ' + str(GCAE_DIR))
paths = [p if os.path.isabs(os.path.dirname(p)) else os.path.join(GCAE_DIR, p)
         for p in (dataprefix, outprefix)]
dataprefix, outprefix = paths[0], paths[1]

print('\nPARAMETERS:\n' +
      '\n'.join(['{}: {}'.format(arg, eval(arg)) for arg in arguments]) + '\n')

outdir = Path(outprefix).parent.absolute()
if not os.path.isdir(outdir):
	os.mkdir(outdir)
	print('\nCreated directory {}'.format(outdir))


## Processing the files

if m is None and n is None and k is None:
	print('\nWarning: m, n, k are not specified, copying files without changes')
	for ext in ('.snp', '.ind', '.eigenstratgeno'):
		copy2(dataprefix + ext, outprefix + ext)
		reportOutfile(outprefix + ext)
	exit()

# snips:
if m is None:
	print('\nWarning: m is not specified, using all SNPs')
	copy2(dataprefix + '.snp', outprefix + '.snp')
	margin1 = ''
else:
	snp = pd.read_table(dataprefix + '.snp', header = None, delimiter=r'\s+')
	seed(42)
	snips_sampled = sorted(sample(range(snp.shape[0]), m))
	snp = snp.iloc[snips_sampled]
	f = open(outprefix + '.snp', 'w')
	f.write(tabulate(snp.values.tolist(), tablefmt = 'plain'))
	f.close()
	margin1 = '\n'
reportOutfile(outprefix + '.snp', margin = margin1)

# individuals:
if n is None and k is None:
	print('\nWarning: neither n nor k are specified, using all individuals')
	copy2(dataprefix + '.ind', outprefix + '.ind')
	margin1 = ''
else:
	ind = pd.read_table(dataprefix + '.ind', header = None, delimiter=r'\s+')
	if n is None:
		ind1_list = ind.to_dict('records')
		ind1_indx = list(ind.index)
	else:
		grps = pd.unique(ind[2])
		ind1_list = []
		ind1_indx = []
		for g in grps:
			ind_g = ind[ind[2] == g][:n]
			ind1_list += ind_g.to_dict('records')
			ind1_indx += list(ind_g.index)
	if k is not None:
		ind1_list = ind1_list[:k]
		ind1_indx = ind1_indx[:k]
	ind1 = pd.DataFrame(ind1_list)
	ind1.to_csv(outprefix + '.ind',
	            header = None, index = False, sep = '\t')
	margin1 = '\n'
reportOutfile(outprefix + '.ind', margin = margin1)

# genomic:
f = open(dataprefix + '.eigenstratgeno', 'r')
esg = f.read().split('\n')[:-1]
f.close()
if m is not None:
	esg = [esg[i] for i in snips_sampled]
if n is not None or k is not None:
	for i in range(len(esg)):
		esg[i] = ''.join([esg[i][j] for j in ind1_indx])
f = open(outprefix + '.eigenstratgeno', 'w')
f.write('\n'.join(esg))
f.close()
reportOutfile(outprefix + '.eigenstratgeno', margin = '\n')


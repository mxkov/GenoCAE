'''Trim the data in example_tiny.

Usage:
  trim_data.py --dataprefix=<name> --outprefix=<name> [-m <num> -n <num> -k <num>]

Options:
  -h --help            show this screen
  --dataprefix=<name>  path and prefix of the data in EIGENSTRAT or PLINK format. if not absolute: assumed relative to GenoCAE/ directory
  --outprefix=<name>   path and prefix for the output data. if not absolute: assumed relative to GenoCAE/ directory
  -m <num>             number of SNPs to cut to. randomly sampled without replacement, order preserved. DEFAULT: all SNPs
  -n <num>             number of samples per population to cut to. DEFAULT: all samples
  -k <num>             if specified, cut the TOTAL number of samples further down to k. DEFAULT: not specified

'''

from docopt import docopt, DocoptExit
import os
import pandas as pd
from numpy import array
from pandas_plink import read_plink1_bin, write_plink1_bin
from pathlib import Path
from random import seed, sample
from shutil import copy2
from tabulate import tabulate


def reportOutfile(filepath, margin = ''):
	print('{}File written: {}'.format(margin, filepath))


def samplesnips(snp_data, m1, dataformat = ('esg', 'plink'), rseed = 42):
	'''
Keep m1 random SNPs from snp_data
	'''
	# sort it before sampling to ensure consistency
	# sorting by chromosome and position
	if dataformat == 'esg':
		snp_data = snp_data.sort_values([1, 3])
		snp_col = 0
	if dataformat == 'plink':
		snp_data = snp_data.sort_values([2, 3])
		snp_col = 1
	# now sample them
	seed(rseed)
	snp_sample_indx = sample(range(snp_data.shape[0]), m1)
	snp_data1 = snp_data.iloc[snp_sample_indx]
	# restore the original order
	snp_data1 = snp_data1.sort_index()
	# save the sample for later
	snp_sample = list(snp_data1[snp_col])
	return snp_data1, snp_sample


def trimIndFam(data, n1, k1, dataformat = ('esg', 'plink')):
	'''
Cut individual + family data to n1 members per family, then to k1 total members
	'''
	# identify the family column, sort data for order invariance
	if dataformat == 'esg':
		data = data.sort_values([2, 0])
		fid_col = 2
	if dataformat == 'plink':
		data = data.sort_values([0, 1])
		fid_col = 0
	if n1 is None:
		data1_list = data.to_dict('records')
		data1_indx = list(data.index)
	else:
		# reduce each family down to n1 individuals
		fams = pd.unique(data[fid_col])
		data1_list = []
		data1_indx = []
		for f in fams:
			data_f = data[data[fid_col] == f][:n1]
			data1_list += data_f.to_dict('records')
			data1_indx += list(data_f.index)
	if k1 is not None:
		# reduce the total number of individuals further down to k1
		data1_list = data1_list[:k1]
		data1_indx = data1_indx[:k1]
	# restore the original order
	data1 = pd.DataFrame(data1_list, index = data1_indx)
	data1 = data1.sort_index()
	return data1, data1.index



## Parsing args

try:
	arguments = docopt(__doc__, version='DataTrim 1.0.0')
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

# If eigenstrat files exist, process those.
# If plink files exist, process those.
# Process both if both exist.

dataformats = {'esg': ('.snp', '.ind', '.eigenstratgeno'),
               'plink': ('.bed', '.bim', '.fam')}
format_available = {frmt: all([os.path.isfile(dataprefix + ext) for ext in exts])
                    for (frmt, exts) in dataformats.items()}
if not all(format_available.values()):
	print('Error: none of the input files exist.')
	exit(1)

if m is None and n is None and k is None:
	print('\nWarning: m, n, k are not specified, copying files without changes')
	for frmt in dataformats:
		if format_available[frmt]:
			for ext in dataformats[frmt]:
				copy2(dataprefix + ext, outprefix + ext)
				reportOutfile(outprefix + ext)
	exit()

snips_sampled = None
ind1 = None


# EIGENSTRATGENO format
if format_available['esg']:
	current_format = 'esg'

	# snips:
	if m is None:
		print('\nWarning: m is not specified, using all SNPs')
		copy2(dataprefix + '.snp', outprefix + '.snp')
		margin1 = ''
	else:
		snp0 = pd.read_table(dataprefix + '.snp', header = None, delimiter=r'\s+')
		# roll m random snips
		snp, snips_sampled = samplesnips(snp0, m, dataformat = current_format)
		# save to file
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
		ind1, ind1_indx = trimIndFam(ind, n, k, dataformat = current_format)
		ind1.to_csv(outprefix + '.ind',
		            header = None, index = False, sep = '\t')
		margin1 = '\n'
	reportOutfile(outprefix + '.ind', margin = margin1)

	# genotypes:
	f = open(dataprefix + '.eigenstratgeno', 'r')
	esg = f.read().split('\n')[:-1]
	f.close()
	if m is not None:
		esg = [esg[i] for i in range(len(esg)) if snp0.iloc[i][0] in snips_sampled]
	if n is not None or k is not None:
		for i in range(len(esg)):
			esg[i] = ''.join([esg[i][j] for j in ind1_indx])
	f = open(outprefix + '.eigenstratgeno', 'w')
	f.write('\n'.join(esg))
	f.close()
	reportOutfile(outprefix + '.eigenstratgeno', margin = '\n')


# PLINK format
if format_available['plink']:
	current_format = 'plink'
	margin1 = '\n'

	G = read_plink1_bin(dataprefix + '.bed', verbose=False)

	# snips & genotypes:
	if m is None:
		print('{}Warning: m is not specified, using all SNPs'.format(margin1))
		margin1 = ''
	else:
		snp = pd.DataFrame(array([G.variant.values,
		                          G.snp.values,
		                          G.chrom.values,
		                          G.pos.values]).transpose())
		if snips_sampled is None:
			_, snips_sampled = samplesnips(snp, m, dataformat = current_format)
		vars_left = [snp.iloc[i][0] for i in range(snp.shape[0])
		             if snp.iloc[i][1] in snips_sampled]
		G = G.sel(variant = vars_left)

	# individuals & genotypes:
	if n is None and k is None:
		print('{}Warning: neither n nor k are specified, using all individuals'.format(margin1))
		margin1 = ''
	else:
		if ind1 is None:
			fam = pd.DataFrame(array([G.fid.values, G.iid.values]).transpose())
			fam, _ = trimIndFam(fam, n, k, dataformat = current_format)
			iids_left = list(fam[1])
		else:
			iids_left = list(ind1[0])
			# to preserve order:
			iids_left = [x for x in G.sample.values if x in iids_left]
		G = G.sel(sample = iids_left)

	# output
	outfiles = [outprefix + ext for ext in dataformats[current_format]]
	write_plink1_bin(G, *outfiles, verbose=False)
	for outfile in outfiles:
		reportOutfile(outfile, margin = '\n')
	


'''Trim a PLINK dataset.

Usage:
  trimdata_plink.py --outprefix=<name> --bed=<name> [--bim=<name> --fam=<name> --chroms=<name> -m <num> --traits=<name> -n <num> -k <num>]

Options:
  -h --help            show this screen
  --outprefix=<name>   path and prefix for the output data. if not absolute: assumed relative to GenoCAE/ directory
  --bed=<name>         path to the BED file. if not absolute: assumed relative to GenoCAE/ directory
  --bim=<name>         path to the BIM file. if not absolute: assumed relative to GenoCAE/ directory. DEFAULT: inferred from BED
  --fam=<name>         path to the FAM file. if not absolute: assumed relative to GenoCAE/ directory. DEFAULT: inferred from BED
  --chroms=<name>      chromosomes to keep, comma-separated with no spaces. DEFAULT: all chromosomes
  -m <num>             number of SNPs to keep AFTER discarding chromosomes. randomly sampled without replacement, order preserved. DEFAULT: all SNPs
  --traits=<name>      trait/phenotype/batch values to keep, comma-separated with no spaces. DEFAULT: all traits
  -n <num>             number of samples PER FAMILY to keep AFTER discarding traits. if larger than the size of a family, the entire family is kept. DEFAULT: all samples
  -k <num>             if specified, cut the TOTAL number of samples further down to k. DEFAULT: not specified

'''

from docopt import docopt, DocoptExit
import os
import pandas as pd
from numpy import array
from pandas_plink import read_plink1_bin, write_plink1_bin
from pathlib import Path
from random import seed, sample
from re import sub
from shutil import copy2


def reportOutfile(filepath, margin = ''):
	print('{}File written: {}'.format(margin, filepath))


def filterByValues(data, values, col = ('chrom', 'trait'), margin = ''):
	'''
Filter data by values in column col
	'''
	names = {'chrom': 'chromosome', 'trait': 'trait'}
	name = names[col]
	# all values present in the dataset
	values_all = list(pd.unique(data[col]))
	# values given in values but not present in the dataset
	values_res = set(values)-set(values_all)
	# if there are any invalid values in values:
	if len(values_res) > 0:
		if len(values_res) == 1:
			grammar = ['', 'is']
		else:
			grammar = ['s', 'are']
		# keep only valid values in values
		values = list(set(values)-set(values_res))
		if len(values) > 0:
			if len(values) == 1:
				msg = '{} {}'.format(name, values[0])
			else:
				msg = '{}s {}'.format(name, ', '.join(values))
		else:
			# if all values are invalid,
			# disregard the values parameter altogether
			# and use all values
			values = None
			msg = 'all {}s'.format(name)
		print('{}{}{} {} {} '.format(margin,
		                             name.capitalize(),
		                             grammar[0],
		                             ', '.join(values_res),
		                             grammar[1])
		      + 'not present in the dataset; using {}'.format(msg))
		margin = ''
	# if values is not empty, filter by values
	if values is not None:
		data = data[data[col].isin(values)]
	return data, values, margin


def samplesnips(data_v, m1, rseed = 42):
	'''
Keep m1 random SNPs from data_v
	'''
	# sort it before sampling to ensure consistency
	# sorting by chromosome and position
	data_v = data_v.sort_values(['chrom', 'pos'])
	# now sample them
	seed(rseed)
	snp_sample_indx = sample(range(data_v.shape[0]), m1)
	data_v = data_v.iloc[snp_sample_indx]
	# restore the original order
	data_v = data_v.sort_index()
	# save the sample for later
	snp_sample = list(data_v['snp'])
	return data_v, snp_sample


def trimIndFam(data_s, n1, k1):
	'''
Cut individual + family data to n1 members per family, then to k1 total members
	'''
	# sort data_s for order invariance
	data_s = data_s.sort_values(['fid', 'iid'])
	if n1 is None:
		data_s1_list = data_s.to_dict('records')
		data_s1_indx = list(data_s.index)
	else:
		# reduce each family down to n1 individuals
		fams = pd.unique(data_s['fid'])
		data_s1_list = []
		data_s1_indx = []
		for f in fams:
			data_s_f = data_s[data_s['fid'] == f][:n1]
			data_s1_list += data_s_f.to_dict('records')
			data_s1_indx += list(data_s_f.index)
	if k1 is not None:
		# reduce the total number of individuals further down to k1
		data_s1_list = data_s1_list[:k1]
		data_s1_indx = data_s1_indx[:k1]
	# restore the original order
	data_s1 = pd.DataFrame(data_s1_list, index = data_s1_indx)
	data_s1 = data_s1.sort_index()
	return data_s1, data_s1.index



if __name__ == '__main__':

	## Parsing args

	try:
		arguments = docopt(__doc__, version='DataTrim-Plink 1.0.0')
	except DocoptExit:
		print('\nInvalid command. Run `python3 trimdata_plink.py --help` for more information.\n')
		exit(1)
	for x in list(arguments.keys()):
		xnew = x.split('-')[-1]
		arg = arguments.pop(x)
		arguments[xnew] = arg

	argtypes = {'m': 'num', 'n': 'num', 'k': 'num',
	            'chroms': 'multistr', 'traits': 'multistr'}
	margin1 = '\n'
	for arg in argtypes.keys():
		if arguments[arg]:
			try:
				if argtypes[arg] == 'num':
					exec('{} = {}'.format(arg, int(arguments[arg])))
				elif argtypes[arg] == 'multistr':
					vals = arguments[arg].split(',')
					exec('{} = vals'.format(arg))
			except:
				print('{}Warning: Invalid argument {}'.format(margin1, arg))
				margin1 = ''
				exec('{} = None'.format(arg))
		else:
			exec('{} = None'.format(arg))

	GCAE_DIR = Path(__file__).resolve().parents[1]
	print('\nProject directory: ' + str(GCAE_DIR))
	for parg in ('outprefix', 'bed', 'bim', 'fam'):
		if not arguments[parg] and parg in ('bim', 'fam'):
			p = sub(r'.bed$', '.{}'.format(parg), bed)
			if not os.path.isfile(p):
				print('\nError: path to {} file is not supplied '.format(parg)
				      + 'and cannot be inferred from bed. Does the file exist?\n')
				exit(1)
			exec('{} = p'.format(parg))
			continue
		if os.path.isabs(os.path.dirname(arguments[parg])):
			exec('{} = "{}"'.format(parg, arguments[parg]))
		else:
			p = os.path.join(GCAE_DIR, arguments[parg])
			exec('{} = p'.format(parg))
		if parg == 'bed' and not os.path.isfile(bed):
			print('\nError: BED file {} does not exist.\n'.format(bed))
			exit(1)

	print('\nPARAMETERS:\n' +
	      '\n'.join(['{}: {}'.format(arg, eval(arg)) for arg in arguments]) + '\n')

	outdir = Path(outprefix).parent.absolute()
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
		print('\nCreated directory {}'.format(outdir))


	## Processing the files

	dataformats = {'plink': ('bed', 'bim', 'fam')}
	format_available = {frmt: all([os.path.isfile(eval(ext)) for ext in exts])
	                    for (frmt, exts) in dataformats.items()}

	if not all(format_available.values()):
		print('\nError: some of the input files do not exist.\n')
		exit(1)

	arg_undefined = [arg is None for arg in (chroms, m, traits, n, k)]
	if all(arg_undefined):
		print('\nchroms, m, traits, n, k are not specified, '
		      + 'copying files without changes')
		for frmt in dataformats:
			if format_available[frmt]:
				for ext in dataformats[frmt]:
					outfile = '.'.join((outprefix, ext))
					copy2(eval(ext), outfile)
					reportOutfile(outfile)
		exit()

	margin1 = '\n'

	# Read data
	G = read_plink1_bin(bed, bim=bim, fam=fam, verbose=False)
	# Extract metadata for samples and variants
	G_meta_s = {field: G[field].values
	            for field in ('sample', 'fid', 'iid', 'trait')}
	G_meta_s = pd.DataFrame(data = G_meta_s)
	G_meta_v = {field: G[field].values
	            for field in ('variant', 'snp', 'chrom', 'pos')}
	G_meta_v = pd.DataFrame(data = G_meta_v)
	# TODO: enable recording list of samples/variants to a file for future reuse?

	# chromosomes:
	if chroms is None:
		print('{}chroms are not specified, using all chromosomes'.format(margin1))
		margin1 = ''
	else:
		G_meta_v, chroms, margin1 = filterByValues(G_meta_v,
		                                           chroms,
		                                           col = 'chrom',
		                                           margin = margin1)
		G = G.sel(variant = list(G_meta_v['variant']))

	# snips:
	if m is None:
		print('{}m is not specified, using all SNPs'.format(margin1))
		margin1 = ''
	else:
		G_meta_v, snips_sampled = samplesnips(G_meta_v, m)
		G = G.sel(variant = list(G_meta_v['variant']))

	# traits / phenotypes / batches / whatever:
	if traits is None:
		print('{}traits are not specified, using all traits'.format(margin1))
		margin1 = ''
	else:
		G_meta_s, traits, margin1 = filterByValues(G_meta_s,
		                                           traits,
		                                           col = 'trait',
		                                           margin = margin1)
		G = G.sel(sample = list(G_meta_s['sample']))

	# individuals:
	if n is None and k is None:
		print('{}n, k are not specified, using all individuals'.format(margin1))
		margin1 = ''
	else:
		G_meta_s, _ = trimIndFam(G_meta_s, n, k)
		G = G.sel(sample = list(G_meta_s['sample']))

	# output
	outfiles = ['.'.join((outprefix, ext)) for ext in dataformats['plink']]
	write_plink1_bin(G, *outfiles, verbose=False)
	for outfile in outfiles:
		reportOutfile(outfile, margin = '\n')
	print('\n')


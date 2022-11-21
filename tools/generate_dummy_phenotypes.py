'''Generate dummy phenotypes from EIGENSTRAT data.

Usage:
  generate_dummy_phenotypes.py --dataprefix=<name> --mode=<num> --superpop_file=<name> [--rseed=<num> --sigma=<num> --margin_factor=<num> --padding_factor=<num> --generator=<name>]

Options:
  -h --help               show this screen
  --dataprefix=<name>     path and prefix of the .ind file. if not absolute: assumed relative to GenoCAE/ directory.
  --mode=<num>            1 for random, 2 for deterministic phenotypes.
  --superpop_file=<name>  path to the superpopulations file. if not absolute: assumed relative to GenoCAE/ directory.
  --rseed=<num>           for mode=1: random seed. DEFAULT: 888
  --sigma=<num>           for mode=1: st.dev of the normal distribution of dummy phenotype values for each population. DEFAULT: 1.0
  --margin_factor=<num>   for mode=1: the uniform distribution for superpopulation means starts at value (margin_factor * sigma * max_population_size). recommended >=3.0. DEFAULT: 3.0
  --padding_factor=<num>  for mode=1: how far away from each other superpopulation means should be expected to be. recommended >1.0. DEFAULT: 1.2
  --generator=<name>      REQUIRED for mode=2: path to a .py script containing two functions: get_markers() and generate(..). see generator.py for an example and more info. the path should be relative to the current working directory.

'''

from docopt import docopt, DocoptExit
import os
from pathlib import Path
from numpy.random import seed
from scipy.stats import uniform, norm
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import rgb2hex
from re import sub
import seaborn as sns


def read_sp_ind(data_prefix, superpopfile):
	'''Reads superpopulation and individual data.'''
	sp_data = pd.read_csv(superpopfile, header = None, names = ['PID', 'SP'])

	ind_data = pd.read_csv(data_prefix + '.ind', sep = None, engine = 'python',
	                       header = None, names = ['IID', None, 'PID'], usecols = [0,2])

	popsizes = ind_data['PID'].value_counts()
	popsizes = popsizes.sort_index()
	populations = list(popsizes.index)
	max_popsize = max(popsizes)
	superpopulations = sorted(set(sp_data['SP']))
	M_sp = len(superpopulations)
	print('Samples by superpopulation:\n{}\n'.format(sp_data['SP'].value_counts().to_string()))

	return sp_data, ind_data, popsizes, populations, max_popsize, superpopulations, M_sp


def printParams(args):
	print('\nProject directory: ' + str(GCAE_DIR))
	print('\nPARAMETERS:\n' +
	      '\n'.join(['{}: {}'.format(arg, eval(arg))
	                 for arg in args if args[arg] is not None]) + '\n\n')



## Parsing the common arguments

print('\n')

try:
	arguments = docopt(__doc__, version='PhenoGen 1.0.0')
except DocoptExit:
	print('Invalid command. Run `python3 generate_dummy_phenotypes.py -h` for more info.\n')
	exit(1)
for k in list(arguments.keys()):
	knew = k.split('--')[-1]
	arg = arguments.pop(k)
	arguments[knew] = arg

GCAE_DIR = Path(__file__).resolve().parents[1]

dataprefix = arguments['dataprefix']
if not os.path.isabs(os.path.dirname(dataprefix)):
	dataprefix = os.path.join(GCAE_DIR, dataprefix)

if arguments['mode'] not in ('1', '2'):
	print('Invalid value for mode: {}. '.format(arguments['mode'])
	      + 'Try `python3 generate_dummy_phenotypes.py -h` for more info.\n')
	exit(1)
mode = int(arguments['mode'])

superpop_file = arguments['superpop_file']
if not os.path.isabs(os.path.dirname(superpop_file)):
	superpop_file = os.path.join(GCAE_DIR, os.path.dirname(superpop_file),
	                             Path(superpop_file).name)


#
## MODE 1: Random
#

if mode == 1:

	## Processing the mode-specific arguments
	arg_defaults = {'rseed': 888, 'sigma': 1.0, 'margin_factor': 3.0, 'padding_factor': 1.2}
	for arg in arg_defaults.keys():
		if arguments[arg]:
			try:
				if type(arg_defaults[arg]) is int:
					exec('{} = {}'.format(arg, int(arguments[arg])))
				else:
					exec('{} = {}'.format(arg, float(arguments[arg])))
				continue
			except:
				print('Invalid argument {}, using default {} = {}'.format(arg, arg,
				                                                          arg_defaults[arg]))
		exec('{} = {}'.format(arg, arg_defaults[arg]))
		arguments[arg] = arg_defaults[arg]
	printParams(arguments)

	## Reading superpop and ind-pop data, getting stats
	(spdata, inddata, pop_sizes, pops,
	 max_pop_size, superpops, M) = read_sp_ind(dataprefix, superpop_file)

	## Generating

	seed(rseed)

	# Superpopulation means
	max_sigma = sigma * max_pop_size
	loc_superpops = margin_factor * max_sigma
	scale_superpops = padding_factor * max_sigma * (M-1)
	means_superpops = pd.Series(abs(uniform.rvs(loc = loc_superpops,
	                                            scale = scale_superpops, size = M)),
	                            index = superpops)
	print('Generated superpopulation means:\n{}\n'.format(means_superpops.to_string()))
	#print('Max sigma {}, mean {}, sigma {}\n'.format(max_sigma, loc_superpops, scale_superpops))

	# Population means
	means_pops_df = pd.DataFrame({'PID': pop_sizes.index, 'N': pop_sizes.values,
	                              'SP': pop_sizes.index.map(lambda x: spdata[spdata['PID'] == x].values[0,1])})
	means_pops_df.insert(3, 'mean',
	                     means_pops_df.apply(lambda x: abs(norm.rvs(loc = means_superpops[x['SP']],
	                                                                scale = x['N']*sigma)),
	                                         axis=1))
	means_pops = pd.Series(means_pops_df['mean'].values, index = means_pops_df['PID'])

	# Individual phenotypes
	phenos = inddata.copy(deep = True)
	phenos.insert(2, 'mean', phenos['PID'].map(means_pops))
	phenos.insert(3, 'value', phenos.apply(lambda x: abs(norm.rvs(loc = x['mean'], scale = sigma)),
	                                       axis=1))
	phenos = phenos.iloc[:, [1,0,3]]


#
## MODE 2: Deterministic
#

if mode == 2:

	## Trying to import the generator functions
	if arguments['generator']:
		if not os.path.isfile(arguments['generator']):
			print('Generator file {} not found.\n'.format(arguments['generator']))
			exit(1)
		gen_module = sub(r'.py$', '', arguments['generator'])
		if '.' in gen_module:
			print('Invalid path to generator (should not contain dots).\n')
			exit(1)
		gen_module_init = os.path.join('/'.join(gen_module.split('/')[:-1]), '__init__.py')
		if not os.path.isfile(gen_module_init):
			f = open(gen_module_init, 'w')
			f.close()
		gen_module = gen_module.replace('/', '.')
		try:
			exec('from {} import get_markers, generate'.format(gen_module))
		except:
			print('Failed to import module `{}`.\n'.format(gen_module))
			exit(1)
		generator = arguments['generator']
	else:
		print('--generator is not specified but is required for mode = 2.\n')
		exit(1)
	printParams(arguments)

	## Reading superpop and ind-pop data, getting stats
	(spdata, inddata, pop_sizes, pops,
	 max_pop_size, superpops, M) = read_sp_ind(dataprefix, superpop_file)

	## Reading snp data
	f = open(dataprefix + '.snp', 'rt')
	snps = [line.strip().split()[0] for line in f.read().strip().split('\n')]
	f.close()

	## Reading genotype data & attaching them to ind-pop data
	f = open(dataprefix + '.eigenstratgeno', 'rt')
	genolines = f.read().strip().split('\n')
	f.close()
	if len(genolines) < len(snps):
		print('Warning: snp file has more snips than eigenstratgeno file')
		snps = snps[:len(genolines)]
	elif len(genolines) > len(snps):
		print('Warning: eigenstratgeno file has more snips than snp file')
		genolines = genolines[:len(snps)]
	genodata = inddata.copy(deep = True)
	genodata.insert(2, 'value', 0.0)
	genocols = list(genodata.columns)
	genodata = pd.concat([genodata] + [pd.Series([int(x) for x in line]) for line in genolines],
	                     axis = 1)
	genodata.columns = genocols + snps

	## Generating
	markers = get_markers()
	genodata['value'] = genodata.apply(lambda x: generate(*[x[mr] for mr in markers]), axis=1)
	phenos = genodata.iloc[:, [1,0,2]]



## Output

outfile = dataprefix + '.phe'
phenos.to_csv(outfile, sep = '\t', index = False)
print('\nPhenotypes written to file {}\n'.format(outfile))


## Plotting
phenos.insert(3, 'SP', phenos['PID'].map(pd.Series(spdata['SP'].values, index = spdata['PID'])))

pos = 1
N_pops_total = len(set(phenos['PID']))
dpi_val = 300.0
for i in range(M):
	sp = superpops[i]
	N_pops = len(set(phenos[phenos['SP'] == sp]['PID']))
	if N_pops == 0:
		continue
	if i == 0:
		ax0 = plt.subplot(1, N_pops_total+M, (pos, pos+N_pops))
		ax = ax0
	else:
		ax = plt.subplot(1, N_pops_total+M, (pos, pos+N_pops), sharey = ax0)
	pos += N_pops+1
	sns.stripplot(data = phenos[phenos['SP'] == sp], x = 'PID', y = 'value',
				  jitter = 0.2, size = 5, edgecolor = 'black', linewidth = 0.5)
	spacepos = sp.find(' ')
	if spacepos >= 0:
		sp_lab = sp[:spacepos] + '\n' + sp[spacepos+1:]
	else:
		sp_lab = sp
	ax.set_xlabel(sp_lab, fontsize = 8)
	ax.xaxis.set_label_position('top')
	ax.set_ylabel(None)
	if i != 0:
		for tick in ax.yaxis.get_major_ticks():
			tick.tick1line.set_visible(False)
			tick.tick2line.set_visible(False)
			tick.label1.set_visible(False)
			tick.label2.set_visible(False)
	plt.tick_params('x', labelrotation = 90.0, labelsize = 5)
	plt.grid(True)
	plt.subplots_adjust(hspace = 0, wspace=0)

plotfile = dataprefix + '.png'
plt.gcf().set_size_inches(3000.0/dpi_val, 1500.0/dpi_val)
plt.savefig(plotfile, dpi = dpi_val, bbox_inches = 'tight')
plt.close()
print('\nPlot written to file {}\n'.format(plotfile))



'''Generate dummy phenotypes from EIGENSTRAT data.

Usage:
  generate_dummy_phenotypes.py --dataprefix=<name> [--superpop_file=<name> --rseed=<num> --sigma=<num> --margin_factor=<num> --padding_factor=<num>]

Options:
  -h --help               show this screen
  --dataprefix=<name>     path and prefix of the .ind file. if not absolute: assumed relative to GenoCAE/ directory.
  --superpop_file=<name>  path to the superpopulations file. if not absolute: assumed relative to GenoCAE/ directory. DEFAULT: data/HO_superpopulations
  --rseed=<num>           random seed. DEFAULT: 888
  --sigma=<num>           st.dev of the normal distribution of dummy phenotype values for each population. DEFAULT: 1.0
  --margin_factor=<num>   the uniform distribution for superpopulation means starts at value (margin_factor * sigma * max_population_size). Recommended >=3.0. DEFAULT: 3.0
  --padding_factor=<num>  how far away from each other superpopulation means should be expected to be. Recommended >1.0. DEFAULT: 1.2

'''

from docopt import docopt, DocoptExit
import os
from pathlib import Path
from numpy.random import seed
from scipy.stats import uniform, norm
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import rgb2hex
import seaborn as sns


## Parsing arguments

try:
	arguments = docopt(__doc__, version='PhenoGen 0.1.0')
except DocoptExit:
	print('Invalid command. Run `python3 generate_dummy_phenotypes.py --help` for more information.')
	exit(1)
for k in list(arguments.keys()):
	knew = k.split('--')[-1]
	arg = arguments.pop(k)
	arguments[knew] = arg

GCAE_DIR = Path(__file__).resolve().parents[1]
print('\nProject directory: ' + str(GCAE_DIR))

if arguments['dataprefix']:
	dataprefix = arguments['dataprefix']
	if not os.path.isabs(os.path.dirname(dataprefix)):
		dataprefix = os.path.join(GCAE_DIR, dataprefix)

if arguments['superpop_file']:
	superpop_file = arguments['superpop_file']
	if not os.path.isabs(os.path.dirname(superpop_file)):
		superpop_file = os.path.join(GCAE_DIR,
									 os.path.dirname(superpop_file),
									 Path(superpop_file).name)
else:
	superpop_file = os.path.join(GCAE_DIR, 'data', 'HO_superpopulations')

if arguments['rseed']:
	rseed = int(arguments['rseed'])
else:
	rseed = 888

arg_defaults = {'sigma': 1.0, 'margin_factor': 3.0, 'padding_factor': 1.2}
for arg in arg_defaults.keys():
	if arguments[arg]:
		try:
			exec('{} = {}'.format(arg, float(arguments[arg])))
			continue
		except:
			print('Invalid argument {}, using default {} = {}'.format(arg, arg, arg_defaults[arg]))
	exec('{} = {}'.format(arg, arg_defaults[arg]))

print('\nPARAMETERS:\n' +
	  '\n'.join(['{}: {}'.format(arg, eval(arg)) for arg in arguments]) + '\n\n')


## Reading superpop and ind-pop data, getting stats

spdata = pd.read_csv(superpop_file, header = None, names = ['PID', 'SP'])

inddata = pd.read_csv(dataprefix + '.ind', sep = None, engine = 'python',
					  header = None, names = ['IID', None, 'PID'], usecols = [0,2])

pop_sizes = inddata['PID'].value_counts()
pop_sizes = pop_sizes.sort_index()
pops = list(pop_sizes.index)
max_pop_size = max(pop_sizes)
superpops = sorted(set(spdata['SP']))
M = len(superpops)
print('Samples by superpopulation:\n{}\n'.format(spdata['SP'].value_counts().to_string()))


## Generating

seed(rseed)

# Superpopulation means
max_sigma = sigma * max_pop_size
loc_superpops = margin_factor * max_sigma
scale_superpops = padding_factor * max_sigma * (M-1)
means_superpops = pd.Series(abs(uniform.rvs(loc = loc_superpops, scale = scale_superpops, size = M)),
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



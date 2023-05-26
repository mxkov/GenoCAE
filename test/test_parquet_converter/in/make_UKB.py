import numpy as np
import os
import pandas as pd
from pandas_plink import read_plink1_bin, write_plink1_bin
from pathlib import Path
import random
import re



def make_filelist(infile, outfile):
	ukbdir = '/mnt/ukb/'
	df0 = pd.read_table(infile, header = None, sep = r'\s+')
	all_files = list(df0[8])

	yeet = ['s488147', 'MT', 'X', 'Y']
	for item in yeet:
		all_files = [x for x in all_files if item not in x]

	df = pd.DataFrame(data = {ext: [ukbdir+x for x in all_files
	                                if x.endswith(ext)]
	                          for ext in ('bed', 'bim', 'fam')})
	df.to_csv(outfile, header = False, index = False)



def take_samples(fdf, N, M, outdir,
                 filecsv = 'files_UKB.csv',
                 random_state = 42):
	if not os.path.isdir(outdir):
		os.makedirs(outdir)
	nfiles = fdf.shape[0]
	numlen = len(str(nfiles))
	EXTS = ('bed', 'bim', 'fam')
	filecsv_df = {ext: [] for ext in EXTS}
	inds_i = None
	random.seed(random_state)
	for i in range(nfiles):
		arg_dict = {ext: fdf[k][i] for k, ext in enumerate(EXTS)}
		arg_dict['verbose'] = False
		G = read_plink1_bin(**arg_dict)
		# extract metadata
		cols = ['variant', 'sample']
		meta = {x: list(G[x].values) for x in cols}
		num_snps = len(meta['variant'])
		num_inds = len(meta['sample'])
		# take N random samples and M random variants;
		# samples should be the same every iteration!
		if inds_i is None:
			inds_i = sorted(random.sample(list(range(num_inds)), N))
		snps_i = sorted(random.sample(list(range(num_snps)), M))
		G1 = G.sel(variant = [meta['variant'][i1] for i1 in snps_i],
		           sample  = [meta['sample' ][i2] for i2 in inds_i])
		outfile = 'UKB-seed{}-{}'.format(random_state,
		                                 str(i+1).zfill(numlen))
		outfile = os.path.join(outdir, outfile)
		write_plink1_bin(G1, outfile+'.bed', verbose = False)
		print('{} samples and {} variants '.format(N, M) +
		      'written to files {}'.format(outfile))
		for ext in EXTS:
			filecsv_df[ext].append(outfile+'.'+ext)
	filecsv = os.path.join(outdir, filecsv)
	filecsv_df = pd.DataFrame(data = filecsv_df)
	filecsv_df.to_csv(filecsv, header = False, index = False)
	print('\nFile list written to {}'.format(filecsv))



if __name__ == '__main__':

	tdd = Path(__file__).resolve().parents[0]
	grd = Path(__file__).resolve().parents[3]

	lf0 = os.path.join(grd, 'data', 'UKBfiles.txt')
	lf  = re.sub('.txt$', '.csv', lf0)
	#make_filelist(lf0, lf)
	# i called this already, then sorted the lines by hand

	files_df = pd.read_csv(lf, header = None)
	take_samples(files_df, 5000, 100, os.path.join(tdd, 'UKB'))

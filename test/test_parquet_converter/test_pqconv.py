from math import ceil
import numpy as np
import os
import pandas as pd
from pandas_plink import read_plink1_bin, write_plink1_bin
import pyarrow.parquet as pq
import pytest
import re
import shutil
from subprocess import Popen, DEVNULL, STDOUT
import xarray as xr



def cleanup(filepaths):
	"""Remove files given in filepaths"""
	for f in filepaths:
		if os.path.isfile(f):
			os.remove(f)
		if os.path.isdir(f):
			shutil.rmtree(f)


def make_sessions(files_df, testoutdir, ssn_descr):
	"""Split CSVs with input files into sessions, return new CSV paths"""

	# Interpret session description
	n_fs = files_df.shape[0]  # Number of FileSets
	if ssn_descr is None:
		sessions = []
	elif ssn_descr.startswith('even'):
		# we'll need to split all input file sets into K roughly even sessions
		K = int(ssn_descr.replace('even', ''))
		sessions = sorted(list(set([i*n_fs // K for i in range(1, K)])))
	elif ssn_descr == 'edges':
		if n_fs > 2:
			sessions = [1, n_fs-1]
		else:
			sessions = []
	elif ssn_descr == 'all':
		sessions = list(range(1, n_fs))
	else:
		raise ValueError('Invalid ssn_descr value')

	ssn_files = np.split(files_df, sessions)  # list of dfs!
	lenlen = len(str(len(ssn_files)-1))
	file_csvs = []
	for i, df in enumerate(ssn_files):
		sid = str(i).zfill(lenlen)
		f = os.path.join(testoutdir, f'testdata_{sid}.csv')
		df.to_csv(f, header=False, index=False)
		file_csvs.append(f)
	return file_csvs


def split_local(bedfile, outpref, testoutdir,
                ssn_id = None, data_split = 1):
	"""Split local test example into parts by snips"""

	data_split_len = len(str(data_split-1))

	G = read_plink1_bin(bedfile, verbose = False)
	n_inds, n_snps = G.shape
	s = G['sample'].values
	v = G['variant'].values

	testdata_files = {}
	EXTS = ('bed', 'bim', 'fam')
	for ext in EXTS:
		testdata_files[ext] = []

	chunk_size = int(ceil(n_snps/data_split))
	start = 0
	chunk_num = -1
	while start < n_snps:
		chunk_snps = v[start:(start+chunk_size)]
		start += chunk_size
		chunk_num += 1
		G1 = G.sel(variant = chunk_snps)
		outfile = outpref + '_' + str(chunk_num).zfill(data_split_len)
		outfile = os.path.join(testoutdir, outfile)
		write_plink1_bin(G1, outfile + '.bed', verbose = False)
		for ext in EXTS:
			testdata_files[ext].append(outfile + '.' + ext)

	testdata_df = pd.DataFrame(data = testdata_files)
	filecsvs = make_sessions(testdata_df, testoutdir, ssn_id)

	return filecsvs


def call_pqconv(testoutdir, filecsvs, outpref, opts = None):
	"""Call converter for test example"""
	# opts: chunk_size, n_batches, batch_col, drop_inds, drop_snps

	outpref = os.path.join(testoutdir, outpref)

	lenlen = len(str(len(filecsvs)-1))
	for i, filecsv in enumerate(filecsvs):
		cmd  = 'python3 tools/parquet_converter.py convert'
		cmd += ' --filecsv={}'.format(filecsv)
		cmd += ' --out={}'.format(outpref)
		if opts is not None:
			for arg in opts:
				if opts[arg] is not None:
					cmd += ' --{}={}'.format(arg, opts[arg])
		lf = os.path.join(testoutdir,
		                  'convert{}.log'.format(str(i).zfill(lenlen)))
		cmd += ' --logfile='+lf
		#print(cmd)
		Popen(cmd, shell = True, stderr = STDOUT).wait()

	cmd  = 'python3 tools/parquet_converter.py merge'
	cmd += ' --out={}'.format(outpref)
	cmd += ' --logfile={}'.format(os.path.join(testoutdir, 'merge.log'))
	cmd += ' --overwrite'
	#cmd += ' --remove_tmpfiles'
	#print(cmd)
	Popen(cmd, shell = True, stderr = STDOUT).wait()


def check_result(testoutdir, files_df, outpref, opts = None):
	"""Check if converter reconstructed original data from test example"""

	outpref = os.path.join(testoutdir, outpref)
	resdict = {}

	## 1) Check if the combined BIM and FAM dfs are equal to the originals
	df_tst = {}
	df_ref = {}
	drop_opts = {'bim': 'drop_snps', 'fam': 'drop_inds'}
	for ext in ('bim', 'fam'):
		# Test output
		if ext == 'fam':
			if 'n_batches' in opts and opts['n_batches'] is not None:
				n_batches = max(opts['n_batches'], 1)
			else:
				n_batches = 10
			L = len(str(n_batches-1))
			df_tst[ext] = None
			samples_by_batch = []
			for i in range(n_batches):
				f = '{}_{}.{}'.format(outpref,str(i).zfill(L), ext)
				df_ = pd.read_table(f, header = None, delimiter=r'\s+',
				                    dtype = 'object')
				samples_by_batch.append(list(df_[1]))
				if df_tst[ext] is None:
					df_tst[ext] = df_
				else:
					df_tst[ext] = pd.concat((df_tst[ext], df_),
					                        ignore_index = True)
					df_tst[ext].drop_duplicates(inplace = True)
		else:
			df_tst[ext] = pd.read_table(outpref + '.' + ext, header = None,
			                            delimiter=r'\s+', dtype = 'object')

		# Original (reference) data
		df_ref[ext] = None
		extcol = 2 if ext == 'fam' else 1
		for i in range(files_df.shape[0]):
			df_ = pd.read_table(files_df[extcol][i], header = None,
			                    delimiter=r'\s+', dtype = 'object')
			if df_ref[ext] is None:
				df_ref[ext] = df_
			else:
				df_ref[ext] = pd.concat((df_ref[ext], df_),
				                        ignore_index = True)
				df_ref[ext].drop_duplicates(inplace = True)
		# drop items from reference data
		dopt = drop_opts[ext]
		if (opts is not None and drop_opts[ext] in opts
		                     and opts[drop_opts[ext]] is not None):
			ddf = pd.read_csv(opts[drop_opts[ext]], header = None,
			                  dtype = 'object')
			dropdata = list(ddf.itertuples(index = False, name = None))
			for col, val in dropdata:
				df_ref[ext] = df_ref[ext].loc[df_ref[ext][int(col)] != val]

		# need to ensure the same item order in both test and ref data
		df_ref[ext].sort_values(by = [1], inplace = True,
		                        key = lambda x: x.map({item: order
		                            for order, item in
		                            enumerate(list(df_tst[ext][1]))}))
		df_ref[ext].reset_index(drop = True, inplace = True)

		# need to compare them now.
		# but the floats could be messed up due to floating-point errors.
		# thus need to round them.
		# this mostly concerns column 2 in BIM files;
		# gotta convert it to numeric first
		if ext == 'bim':
			df_ref[ext][2] = pd.to_numeric(df_ref[ext][2]).round(6)
			df_tst[ext][2] = pd.to_numeric(df_tst[ext][2]).round(6)

		key = ext.upper() + ' restored correctly'
		resdict[key] = df_tst[ext].equals(df_ref[ext])

	## 2) Check genotypes
	key = 'Genotypes restored correctly'

	# Genodata 1: original data
	genodata1 = None
	for i in range(files_df.shape[0]):
		g1_G = read_plink1_bin(bed = files_df[0][i],
		                       bim = files_df[1][i],
		                       fam = files_df[2][i],
		                       verbose = False)
		# Apply the same sorting and dropping, fill NAs
		g1_G = g1_G.assign_coords(variant = g1_G['snp'].values)
		snps_here = list(g1_G['variant'].values)
		inds_here = list(g1_G['sample'].values)
		snps_sort = [x for x in list(df_ref['bim'][1]) if x in snps_here]
		inds_sort = [x for x in list(df_ref['fam'][1]) if x in inds_here]
		g1_G = g1_G.sel(variant = snps_sort, sample  = inds_sort)
		g1_G = g1_G.fillna(9.0)
		if genodata1 is not None:
			genodata1 = xr.merge([genodata1, g1_G]
			   ).to_array().squeeze().drop_vars('variable').rename('genotype')
		else:
			genodata1 = g1_G.copy(deep = True)
	# Sort once again
	snps_here = list(genodata1['variant'].values)
	inds_here = list(genodata1['sample'].values)
	snps_sort = [x for x in list(df_ref['bim'][1]) if x in snps_here]
	inds_sort = [x for x in list(df_ref['fam'][1]) if x in inds_here]
	genodata1 = genodata1.sel(variant = snps_sort, sample  = inds_sort)
	# Extract
	genodata1 = genodata1.values

	# Genodata 2: combined (reconstructed) data
	genodata2 = None
	iids_sch = []
	for i in range(n_batches):
		f = '{}_{}.parquet'.format(outpref,str(i).zfill(L))
		g2_pq = pq.read_table(f)
		iids_sch.extend([g2_pq.schema.field(j).name
		                 for j in range(g2_pq.shape[1])])
		g2_df = g2_pq.to_pandas(ignore_metadata = True)
		g2 = g2_df.to_numpy().T
		sch_fields = [field.name for field in g2_pq.schema]
		cond1 = ( sch_fields != samples_by_batch[i])
		cond2 = (g2.shape[0] != len(samples_by_batch[i]))
		if cond1 or cond2:
			resdict[key] = False
			return resdict
		if genodata2 is None:
			genodata2 = g2
		else:
			genodata2 = np.concatenate((genodata2, g2), axis=0)
	# Compare
	if genodata1.shape == genodata2.shape:
		resdict[key] = np.array_equiv(genodata1, genodata2)
	else:
		resdict[key] = False

	## 3) Check table schema as well
	key = 'Table schema consistent with FAM'
	iids_fam = list(df_tst['fam'][1])
	resdict[key] = (iids_sch == iids_fam)

	return resdict


def tsting_pipeline(where = 'local',
                    out_files = 'pqtest-result',
                    test_indir  = 'test/test_parquet_converter/in',
                    test_outdir = 'test/test_parquet_converter/out',
                    datasplit = 12,
                    sessions_id = None,
                    options = None):
	"""All testing steps combined for a single converter test"""

	if not os.path.isdir(test_outdir):
		os.makedirs(test_outdir)

	# Set up input data
	file_csv = os.path.join(test_indir, f'files_{where}.csv')
	fdf = pd.read_csv(file_csv, header = None, dtype = 'object')
	if where == 'local':
		# Split input data into many files and group them into sessions
		file_csvs = split_local(fdf[0][0], 'pqtest', test_outdir,
		                        ssn_id = sessions_id,
		                        data_split = datasplit)
	if where == 'UKB':
		# Split input files into sessions
		file_csvs = make_sessions(fdf, test_outdir, sessions_id)

	call_pqconv(test_outdir, file_csvs, out_files, opts = options)

	results = check_result(test_outdir, fdf, out_files, opts = options)
	return list(results.values())


### The test function itself ###

def test_pqconv(f_where, f_datasplit, f_sessions,
                f_n_batches, f_batch_col, f_chunk_size,
                f_drop_inds, f_drop_snps):
	"""Top-level test function for tools/parquet_converter.py"""

	td  = 'test/test_parquet_converter'
	tid = os.path.join(td, 'in')
	tod = os.path.join(td, 'out')

	reslist = tsting_pipeline(where       = f_where,
	                          out_files   = 'pqtest-result',
	                          test_indir  = tid,
	                          test_outdir = tod,
	                          datasplit   = f_datasplit,
	                          sessions_id = f_sessions,
	                          options     = {'n_batches':  f_n_batches,
	                                         'batch_col':  f_batch_col,
	                                         'chunk_size': f_chunk_size,
	                                         'drop_inds':  f_drop_inds,
	                                         'drop_snps':  f_drop_snps})
	testfiles = [os.path.join(tod, f) for f in os.listdir(tod)]
	#cleanup([f for f in testfiles if not f.endswith('.log')])
	cleanup(testfiles)
	assert all(reslist)


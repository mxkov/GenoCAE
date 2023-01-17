from math import ceil
import numpy as np
import os
import pandas as pd
from pandas_plink import read_plink1_bin, write_plink1_bin
import pyarrow.parquet as pq
from subprocess import Popen, DEVNULL, STDOUT


def cleanup(filepaths):
	'''cleanup(filepaths): remove files given in filepaths'''
	for f in filepaths:
		if os.path.isfile(f):
			os.remove(f)


def drop_from_tiny(dataprefix, dataprefix_out, dropfile):
	'''Drop certain samples from example_tiny to test dropping'''

	drop = pd.read_csv(dropfile, header = None)
	drop = list(drop.itertuples(index = False, name = None))

	fam = pd.read_table(dataprefix + '.fam',
	                    header = None, delimiter=r'\s+')
	dropped = []
	for col, val in drop:
		dropped1 = list(fam.loc[fam[col] == val, 1])
		dropped.extend(dropped1)
	dropped.sort()

	G = read_plink1_bin(dataprefix + '.bed', verbose = False)
	s = G['sample'].values
	G = G.sel(sample = [x for x in s if x not in dropped])
	write_plink1_bin(G, dataprefix_out + '.bed', verbose = False)
	newfiles = [dataprefix_out + ext for ext in ('.bed', '.bim', '.fam')]

	return newfiles, dropped


def split_tiny(dataprefix, testfilename, testdatadir):
	'''Create a test example by splitting example_tiny (PLINK) into 12 chunks of snips'''

	G = read_plink1_bin(dataprefix + '.bed', verbose = False)
	n_inds, n_snps = G.shape
	s = G['sample'].values
	v = G['variant'].values

	if not os.path.isdir(testdatadir):
		os.makedirs(testdatadir)
	testexample_files = []
	EXTS = ('bed', 'bim', 'fam')

	n_chunks = 12
	n_chunks_len = 2
	chunk_size = int(ceil(n_snps/n_chunks))
	start = 0
	chunk_num = -1
	while start < n_snps:
		chunk_snps = v[start:(start+chunk_size)]
		start += chunk_size
		chunk_num += 1
		G1 = G.sel(variant = chunk_snps)
		outfile = testfilename.replace('[0-9]?[0-9]?',
		                               str(chunk_num).zfill(n_chunks_len))
		outfile = os.path.join(testdatadir, outfile)
		write_plink1_bin(G1, outfile + '.bed', verbose = False)
		for ext in EXTS:
			testexample_files.append(outfile + '.' + ext)

	return testexample_files


def call_pqconv(testdatadir, testfilename, outfilename,
                chunk_size = None, drop = None):
	'''Call parquet_converter.py for the test example'''

	testfilepaths = os.path.join(testdatadir, testfilename)
	outprefix     = os.path.join(testdatadir, outfilename)

	cmd = 'python3 tools/parquet_converter.py'
	EXTS = ('bed', 'bim', 'fam')
	for ext in EXTS:
		cmd += ' --{}={}'.format(ext, testfilepaths + '.' + ext)
	cmd += ' --out={}'.format(outprefix)
	if drop is not None:
		cmd += ' --drop={}'.format(drop)
	if chunk_size is not None:
		cmd += ' --chunksize={}'.format(chunk_size)

	Popen(cmd, shell = True, stdout = DEVNULL, stderr = STDOUT).wait()
	#os.system(cmd)


def check_result(dataprefix, testdatadir, outfilename):
	'''Check if parquet_converter.py managed to reconstruct the original data from the test example'''

	outprefix = os.path.join(testdatadir, outfilename)
	resdict = {}

	# Check if the combined BIM and FAM dataframes are equal to the originals
	for ext in ('bim', 'fam'):
		df_ref = pd.read_table(dataprefix + '.' + ext,
		                       header = None, delimiter=r'\s+')
		df_tst = pd.read_table(outprefix  + '.' + ext,
		                       header = None, delimiter=r'\s+')
		# need to compare them now.
		# but the floats could be messed up due to floating-point errors.
		# thus need to round them.
		df_ref = df_ref.round(6)
		df_tst = df_tst.round(6)
		key = ext.upper() + ' restored correctly'
		resdict[key] = df_tst.equals(df_ref)

	# Genodata 1: original data
	genodata1_G = read_plink1_bin(dataprefix + '.bed', verbose = False)
	genodata1_G = genodata1_G.fillna(9.0)
	genodata1 = genodata1_G.values
	# Genodata 2: combined (reconstructed) data
	genodata2_pq = pq.read_table(outprefix + '.parquet')
	genodata2_df = genodata2_pq.to_pandas(ignore_metadata = True)
	genodata2 = genodata2_df.to_numpy().T
	# Compare
	key = 'Genotypes restored correctly'
	if genodata1.shape == genodata2.shape:
		resdict[key] = np.array_equiv(genodata1, genodata2)
	else:
		resdict[key] = False

	return resdict


def tsting_pipeline(out_files = 'pqtest_result',
                    chunk = None, drop_file = None):
	'''All testing steps combined for a single test'''

	dp  = 'data/example_tiny/HumanOrigins249_tiny'
	tfn = 'pqtest_[0-9]?[0-9]?_file'
	tdd = 'test/test_parquet_converter/'
	dp1 = os.path.join(tdd, 'HumanOrigins249_tiny_dropped')

	filelist = split_tiny(dp, tfn, tdd)
	call_pqconv(tdd, tfn, out_files,
	            chunk_size = chunk, drop = drop_file)
	cleanup(filelist)
	if drop_file is None:
		results = check_result(dp, tdd, out_files)
	else:
		filelist1, _ = drop_from_tiny(dp, dp1, drop_file)
		results = check_result(dp1, tdd, out_files)
		cleanup(filelist1)

	msg = '\n\nTEST RESULTS'
	msg_add = []
	if chunk is not None:
		msg_add.append(' chunk size {}'.format(chunk))
	if drop_file is not None:
		msg_add.append(' drop file {}'.format(drop_file))
	if len(msg_add) > 0:
		msg += ' for' + ','.join(msg_add)
	msg += ':\n'
	print(msg)
	for tst in results:
		print('{}: {}'.format(tst, results[tst]))

	results_list = list(results.values())
	return results_list


def test_pqconv_1():
	'''Test 1: single chunk of snips'''
	reslist = tsting_pipeline(out_files = 'pqtest_result_1',
	                          chunk = None)
	assert all(reslist)

def test_pqconv_2():
	'''Test 2: several chunk of snips, size 800'''
	reslist = tsting_pipeline(out_files = 'pqtest_result_2',
	                          chunk = 800)
	assert all(reslist)

def test_pqconv_3():
	'''Test 3: single chunk of snips, with dropping'''
	reslist = tsting_pipeline(out_files = 'pqtest_result_3',
	                          chunk = None,
	                          drop_file = 'test/drop_pqconvtest.csv')
	assert all(reslist)

def test_pqconv_4():
	'''Test 4: several chunk of snips, size 1234, with dropping'''
	reslist = tsting_pipeline(out_files = 'pqtest_result_4',
	                          chunk = 1234,
	                          drop_file = 'test/drop_pqconvtest.csv')
	assert all(reslist)


if __name__ == '__main__':
	test_pqconv_1()
	test_pqconv_2()
	test_pqconv_3()
	test_pqconv_4()
	print('')


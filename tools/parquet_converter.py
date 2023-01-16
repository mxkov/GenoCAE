'''Convert multiple PLINK files into a single Apache Parquet dataset.
All path arguments are assumed relative to GenoCAE/ directory if not absolute.
Regex are only accepted for file names, not directories. Hence, all files of a specific format should be in the same directory.

Usage:
  parquet_converter.py --bed=<name> --bim=<name> --fam=<name> --out=<name> [--normalize=<name> --chunksize=<num>]

Options:
  -h --help           show this screen
  --bed=<name>        path + filename regex for the BED files
  --bim=<name>        path + filename regex for the BIM files
  --fam=<name>        path + filename regex for the FAM files
  --out=<name>        path and prefix for the output data (BIM + FAM + Parquet)
  --normalize=<name>  normalization method: genotypewise01, smartPCAstyle, standard, or no normalization if not specified
  --chunksize=<num>   size of SNP chunks to use when building the dataset, mostly for testing/debugging purposes. DEFAULT: None (auto)
'''
# TODO: implement normalization (don't forget nan imputation, see data_handler)

# TODO and to reflect: eigenstrat format is inverted relative to plink.
# that is, each genotype value is the number of ref alleles in eigenstrat
# and the number of alt alleles in plink.
# But that's just interpretation, that the model is probably agnostic of.
# Does it matter? Probably not. But we need to be sure.
# Train the model on the same dataset in the two formats and compare the predictions.

from docopt import docopt, DocoptExit
from math import ceil
import os
import pandas as pd
from pandas_plink import read_plink1_bin
from pathlib import Path
from psutil import virtual_memory
import pyarrow as pa
import pyarrow.parquet as pq
import random
import re
from shutil import rmtree
from string import ascii_lowercase, digits
import traceback
import warnings
import xarray as xr


class GenoMeta:
	'''A Genomic Metadata structure (snips from BIM files or individuals from FAM files).

	Attributes:
	  fileindex: dict
	      each key is item ID (item = snip or individual), and the value is a set of i`s of all BED files that contain genotypes for this item
	  filelist: list
	      BIM or FAM filepaths corresponding to the i`s
	  bedlist: list
	      BED filepaths corresponding to the i`s
	  df: pandas.DataFrame
	      full BIM or FAM dataframe
	  _start_index: int
	      tip of the filelist, in case we decide to append more files to an existing GenoMeta object

	Methods:
	  append(filelist_add, bedlist_add, drop = [], idcol = 1):
	      add data from more BIM or FAM files

See GenoMeta.__init__.__doc__ for more info on how to build a GenoMeta object.
	'''

	def __init__(self, filelist_init, bedlist_init,
	                   drop_init = [], make_index = False):
		'''Constructs a GenoMeta object.

	Parameters:
	  filelist_init: list
	      BIM or FAM filepaths to 1) initialize self.filelist with and 2) to read self.df from
	  bedlist_init: list
	      BED filepaths to initialize self.bedlist with; assumed to be aligned with filelist_init
	  drop_init: list
	      list of (col, val) tuples; for each pair, values equal to val will be dropped from column col in self.df. DEFAULT: empty list
	  make_index: bool
	      whether to create self.fileindex. DEFAULT: False
		'''
		if make_index:
			self.fileindex = {}
		else:
			self.fileindex = None
		self.filelist = []
		self.bedlist = []
		self.df = None
		self._start_index = 0
		self.append(filelist_init, bedlist_init, drop = drop_init)

	def append(self, filelist_add, bedlist_add, drop = [], idcol = 1):
		'''Adds more data to a GenoMeta object.

	Parameters:
	  filelist_add: list
	      BIM or FAM filepaths to add to self.filelist and self.df
	  bedlist_add: list
	      BED filepaths to add to self.bedlist; assumed to be aligned with filelist_add
	  drop: list
	      list of (col, val) tuples; for each pair, values equal to val will be dropped from column col in the dataframe constructed from filelist_add. DEFAULT: empty list
	  idcol: int
	      index of the ID column in the BIM/FAM files, starting with 0. DEFAULT: 1 (snip ID and IID in BIM and FAM formats)
		'''
		# TODO: !!!!! DROP SHOULD BE DONE IN A DIFFERENT PLACE AS WELL
		# WHEN LOADING Gs
		for i, filepath in enumerate(filelist_add):
			df_ = pd.read_table(filepath,
			                    header = None, delimiter=r'\s+')
			for col, val in drop:
				df_ = df_.loc[df_[col] != val]
			if self.df is None:
				self.df = df_.copy(deep = True)
			else:
				self.df = pd.concat((self.df, df_),
				                    ignore_index = True)
				self.df.drop_duplicates(subset = idcol, inplace = True)
			if self.fileindex is not None:
				# Index of the current file:
				i_current = self._start_index + i
				# Look at each sample/snip id,
				# update the set of files where it's present
				for id_ in df_[idcol]:
					if id_ in self.fileindex:
						self.fileindex[id_].add(i_current)
					else:
						self.fileindex[id_] = {i_current}
		self.filelist.extend(filelist_add)
		self.bedlist.extend(bedlist_add)
		self._start_index += i+1
		# TODO: figure out how to handle duplicates


def MissingPairsWarning(indent = '\n'):
	'''Raised when the provided BIM and FAM files include variant+sample pairs with genotypes not given in the BED files provided.
This does NOT include genotypes explicitly specified as missing in a BED file; rather, it's when a variant+sample pair is missing completely.
For example, this can happen when n-1 BED files each have a different subset of variants for all samples, but one sample is missing from the last n-th file.

These genotypes will then be processed in the same way as regular missing genotypes, but the warning needs to be raised.
	'''
	msg = ('{}Warning: Genotypes for some variant+sample '.format(indent)
	     + ' pairs are not present in the provided BED files')
	_ = warnings.warn(msg)
	# TODO: how about specifying which pairs are missing huh


def Sorter(col, ref):
	'''Sorter(col, ref): sort col in the same order as ref, return the result.'''
	return col.map({item: order for order, item in enumerate(ref)})


def normalize_genos(genotypes, method = None):
	'''Apply per-variant normalization to genotypes.

	Parameters:
	  genotypes: numpy.ndarray
	      Genotype array of shape (variants, samples)
	  method: str
	      Normalization method to use

	Returns:
	  numpy.ndarray of shape (variants, samples): normalized genotypes
	'''
	if method is None:
		return genotypes
	else:
		raise NotImplementedError    # TODO



if __name__ == '__main__':
	# just some flags for debugging
	debug_filepaths = False
	debug_prechunks = False

	# Parse arguments
	try:
		arguments = docopt(__doc__, version='Parquet Converter 1.0.0')
	except DocoptExit:
		print('\nInvalid command. Run `python3 parquet_converter.py --help` for more information.\n')
		exit(1)
	for x in list(arguments.keys()):
		xnew = x.split('-')[-1]
		arg = arguments.pop(x)
		arguments[xnew] = arg

	EXTS = ('bed', 'bim', 'fam')

	# Handle the paths
	GCAE_DIR = Path(__file__).resolve().parents[1]
	print('\nProject directory: ' + str(GCAE_DIR))
	for parg in ('bed', 'bim', 'fam', 'out'):
		if not os.path.isabs(os.path.dirname(arguments[parg])):  # note: works with regex dirnames
			arguments[parg] = os.path.join(GCAE_DIR, arguments[parg])

	# Handle chunksize
	int_args = ('chunksize',)
	for arg in int_args:
		if arguments[arg]:
			try:
				arguments[arg] = int(arguments[arg])
			except:
				print('\nInvalid value of argument ' +
				      '{}: {}'.format(arg, arguments[arg]))
				exit(1)

	# Print parameters
	print('\nPARAMETERS:\n' +
	      '\n'.join(['{}: {}'.format(arg, arguments[arg])
	                 for arg in arguments]) +
	      '\n')

	# Create the output directory if it doesn't exist
	outdir = Path(arguments['out']).parent.absolute()
	if not os.path.isdir(outdir):
		os.makedirs(outdir)    # TODO: should be put in all tools instead of mkdir!
		print('\nCreated directory {}'.format(outdir))

	# Time to do something meaningful

	# Get lists of BIM, FAM, BED files
	files = {}
	for ext in EXTS:
		path_pattern = arguments[ext]
		filedir = str(Path(path_pattern).parent.absolute())
		filename_pattern = str(Path(path_pattern).name)
		filename_regex = re.compile(filename_pattern)
		files[ext] = [os.path.join(filedir, f) for f in os.listdir(filedir)
		              if filename_regex.match(f)]
		files[ext] = sorted(files[ext])
	# if files from the sorted lists don't match each other,
	# then it's on the user, sorry

	# Report unmatched files
	file_counts = {ext: len(files[ext]) for ext in EXTS}
	n_full_sets = min(file_counts.values())
	for ext in EXTS:
		if file_counts[ext] > n_full_sets:
			extra_files = [f for i,f in enumerate(files[ext])
			               if i >= n_full_sets]
			print('\nWarning: no match for the following ' +
			      '{} files\n{}'.format(ext.upper(),
			                            '\n'.join(extra_files)))

	if debug_filepaths:
		test = pd.DataFrame(data = {ext: files[ext] for ext in EXTS})
		testfile = os.path.join(outdir, 'files.csv')
		test.to_csv(testfile,
		            index = False, sep = '\t')
		print('\nFile written: {}\n'.format(testfile))
		exit()

	# Read all BIMs & FAMs, combine them,
	# make an index of BIM files for snips
	print('\nReading BIM and FAM data...')
	gmeta = {}
	gmeta['bim'] = GenoMeta(files['bim'], files['bed'],
	                        make_index = True)
	gmeta['fam'] = GenoMeta(files['fam'], files['bed'],
	                        drop_init = [(5, 'redacted3')])
	# Write the combined BIM and FAM
	for ext in ('bim', 'fam'):
		outfile = arguments['out'] + '.' + ext
		gmeta[ext].df.to_csv(outfile,
		                     header = None, index = False, sep = '\t')
		print('Combined {} data written to file {}'.format(ext.upper(),
		                                                   outfile))

	# Find snip chunk size
	n_inds = gmeta['fam'].df.shape[0]
	n_snps = gmeta['bim'].df.shape[0]
	if arguments['chunksize']:
		chunk_size = arguments['chunksize']
	else:
		k = 0.4    # the free ram fraction we're willing to sacrifice for a chunk
		# TODO: wait, no, see how much pandasplink actually takes up!
		chunk_size = min(round(virtual_memory()[1]*k / (4*n_inds)), n_snps)
	n_chunks = int(ceil(n_snps/chunk_size))
	print('\nUsing {} chunks of size {} '.format(n_chunks, chunk_size) +
	      '(last chunk might not be full)')

	if debug_prechunks:
		exit()

	# If we have multiple chunks, we have to write each chunk
	# to its own parquet file before combining them.
	# If that's the case, create a temp dir for those smaller files.
	if n_chunks > 1:
		outprefix_base = os.path.basename(arguments['out'])
		tmpdir = os.path.join(outdir, '.' + outprefix_base + '_tmp')
		while os.path.isdir(tmpdir):
			# if a dir exists with the same name,
			# make up random names until we come up with something unique
			insert = ''.join(random.sample(ascii_lowercase + digits, 5))
			tmpdir = os.path.join(outdir,
			                      '.' + '_'.join((outprefix_base,
			                                      insert, 'tmp')))
		os.mkdir(tmpdir)
		# And also get the number of digits, for zfilling later
		n_chunks_len = len(str(n_chunks))

	# Now iterate over snip chunks;
	# in each iteration, find and read all snips in the chunk
	# from ALL SAMPLES, normalize, save.
	print('\nPreparing chunks...')
	iids_all = gmeta['fam'].df[1]
	start = 0
	chunk_num = -1
	while start < n_snps:
		# Within this loop, we're working with one chunk of snips.
		#
		# Create chunk:
		chunk_snps = gmeta['bim'].df[1][start:(start+chunk_size)]
		start += chunk_size
		chunk_num += 1
		print('Processing chunk #{}'.format(chunk_num+1), end='\r')
		# Now, gotta combine file sets for all snips in this chunk
		# to find which BED files we need to read
		# to get all these snips for all samples.
		chunk_dict = {s: gmeta['bim'].fileindex[s] for s in chunk_snps}
		chunk_files = set.union(*list(chunk_dict.values()))
		chunk_files = sorted(list(chunk_files))
		chunk_data = None
		done_with_files = False
		try:
			# Read all genotypes for this chunk and combine them
			for i in chunk_files:
				# Read PLINK
				G = read_plink1_bin(files['bed'][i],
				                    files['bim'][i],
				                    files['fam'][i],
				                    verbose = False)
				# Extract metadata for variants and samples
				meta_v = pd.DataFrame(data = {x: G[x].values
				                              for x in ('variant','snp')})
				meta_s = pd.DataFrame(data = {x: G[x].values
				                              for x in ('sample', 'iid')})
				# Assign snip ids to the variant field.
				# The "variant{i} pattern is too ambiguous and messes everything up.
				G = G.assign_coords(variant = G['snp'].values)
				# Select chunk snips, filter G, replace nans
				chunk_snps_present = list(meta_v['snp'].loc[meta_v['snp'].isin(chunk_snps)])
				G = G.sel(variant = chunk_snps_present)
				G = G.fillna(9.0)
				# Combine
				if chunk_data is not None:
					chunk_data = xr.merge([chunk_data, G]#, compat='override'
					   ).to_array().squeeze().drop(labels='variable').rename('genotype')
					# Fun fact: this automatically picks a reasonable chunksize for the resulting dask array.
					# (To be clear, that's a different kind of chunk.)
					# TODO: gotta memory profile this thing...
				else:
					chunk_data = G.copy(deep = True)
			done_with_files = True
			# The result has to be checked for nans!
			# Those would be gaps in the merged array:
			# pairs (sample, variant) not present in the input files,
			# not even as explicitly specified missing genotypes.
			# This is also why we fill missing genotypes before merging.
			if xr.DataArray.isnull(chunk_data).any():
				MissingPairsWarning()
				chunk_data = chunk_data.fillna(9.0)
				# TODO: this might eat memory.
				# Maybe we should just shut up and fill all missing values with 9s.

			# Sort the chunk to align it with combined BIM/FAM data
			# (just in case)
			meta_v = pd.DataFrame(data = {field: chunk_data[field].values
			                              for field in ('variant','snp')})
			meta_v.sort_values(by = 'snp', inplace = True,
			                   key = lambda x: Sorter(x, chunk_snps))
			meta_s = pd.DataFrame(data = {field: chunk_data[field].values
			                              for field in ('sample', 'iid')})
			meta_s.sort_values(by = 'iid', inplace = True,
			                   key = lambda x: Sorter(x, iids_all))
			chunk_data = chunk_data.sel(variant = list(meta_v['variant']),
			                            sample  = list(meta_s['sample']))

			# Extract genotypes from the chunk
			# (and hope this doesn't eat all your RAM)
			chunk_data = chunk_data.values.T
			# The transposed array is numpy.ndarray of shape (variants, samples).

			# After getting this chunk of genotypes for all samples,
			# this is where we normalize
			if arguments['normalize']:
				chunk_data = normalize_genos(chunk_data,
				                             method=arguments['normalize'])

			# Convert into a pyarrow table
			chunk_data = pa.Table.from_pandas(pd.DataFrame(chunk_data))

			# Write this chunk to a smaller parquet file
			# (if we have multiple chunks)
			if n_chunks > 1:
				outfilename = (outprefix_base
				               + str(chunk_num).zfill(n_chunks_len)
				               + '.parquet')
				outfile = os.path.join(tmpdir, outfilename)
				pq.write_table(chunk_data, outfile)

		except Exception as e:
			if done_with_files:
				insert = ''
			else:
				insert = ' and file {}'.format(files['bed'][i])
			print('\nWhile processing ' +
			      'chunk #{}{}, '.format(chunk_num+1, insert) +
			      'a following exception occured:')
			print(traceback.format_exc())
			exit(1)

	# Finally, we combine all the chunks into one parquet file.
	outfile = arguments['out'] + '.parquet'
	if n_chunks > 1:
		print('\n\nCombining chunks...')
		# get the schema from the first temp file
		tmpfiles = sorted(os.listdir(tmpdir))
		tfile = os.path.join(tmpdir, tmpfiles[0])
		sch = pq.read_table(tfile).schema
		# iterate over temp files, combine them
		with pq.ParquetWriter(outfile, sch,
		                      version = '2.6',
		                      compression = 'gzip',
		                      use_dictionary = True,
		                      data_page_size = 2097152,  # 2MB
		                      write_statistics = True) as writer:
			for chunk_num, tfilename in enumerate(tmpfiles):
				print('Writing chunk #{}'.format(chunk_num+1), end='\r')
				tfile = os.path.join(tmpdir, tfilename)
				writer.write_table(pq.read_table(tfile))
		print('\nGenotype data from all chunks written to file '
		      + str(outfile))
		# remove the temp dir
		print('Deleting the temporary directory...')
		rmtree(tmpdir)
	else:
		# write our single chunk into the final file
		pq.write_table(chunk_data, outfile)
		print('\n\nGenotype data written to file ' + str(outfile))

	print('\nDone!\n')


"""Convert multiple PLINK filesets into an Apache Parquet dataset.

The input dataset may consist of any number of PLINK filesets (BED + BIM + FAM), provided as paths in a CSV file.
The output is a set of FAM + Parquet file pairs, with a single BIM file.
All path arguments are assumed relative to GenoCAE/ directory if not absolute.

Regardless of how samples and variants are distributed among the input filesets, the final Parquet dataset will be split into FAM + Parquet filesets by samples - with n_batches filesets, or "batches", in total. Each batch will thus contain all variants, in a single final BIM file.


*** HOW TO RUN ***

If the number of variants is really large, doing such a conversion in one step can be very time- and memory-inefficient. Thus, the pipeline is split into two stages: "convert" and "merge".

** STEP 1) First, you run as many "convert" stages as you want, as long as a) the input PLINK filesets do not overlap between the runs and b) the batching options --n_batches and --batch_col remain the same for all runs.
Batching is applied at every "convert" run. If --chunk_size is given and is less than the number of variants in the input, or if all variants don't fit into RAM, each batch is further split into chunks by variants. Larger chunk size means better compression in the final Parquet files.

Note that the output of "convert" stages is not written to the path given by --out, but rather to a temporary subdirectory within the same directory, with '_chunks' appended to its name. This temporary subdirectory is later processed in the "merge" stage.

** STEP 2) Then you do a single "merge" run, with the same --out option as your "convert" runs. This will merge all variant chunks within each sample batch, producing a final dataset, which is actually placed to the --out path this time.

Because of the temporary subdirectory, you have to run the "merge" stage even if you only had one "convert" run and one chunk (in which case the "merge" run would simply rearrange the files).


Usage:
  parquet_converter.py convert --filecsv=<name> --out=<name> [--n_batches=<num> --batch_col=<num> --chunk_size=<num> --drop_inds=<name> --drop_snps=<name> --memprof --logfile=<name>]
  parquet_converter.py merge --out=<name> [--memprof --logfile=<name> --overwrite --remove_tmpfiles]

Options:
  -h --help           show this screen
  --filecsv=<name>    path to a CSV file with paths to input BED, BIM and FAM files (as 3 columns, in this order, with no header)
  --out=<name>        path and prefix for the output data (BIM + FAM + Parquet)
  --n_batches=<num>   number of sample batches to split the dataset into (not to be confused with training batches!). DEFAULT: 10
  --batch_col=<num>   index of a column (between 0 and 5) in FAM files, to ensure that samples with equal batch_col values end up in the same batch (if possible). DEFAULT: None (split into equal batches instead)
  --chunk_size=<num>  size of variant chunks to use when building the dataset, larger = better. DEFAULT: None (as many as fits into RAM)
  --drop_inds=<name>  path to a text file with lines in format 'i,var'; for every such line, all samples with value 'var' in column i (starting from 0) will be dropped from the FAM files
  --drop_snps=<name>  same as drop_inds, but for variants and BIM files
  --memprof           do memory profiling
  --logfile=<name>    attempt to write stdout to this file in real time (doesn't include memory profiling if --memprof is given)
  --overwrite         overwrite the final output files if they already exist; if not specified and the files exist, an error is raised
  --remove_tmpfiles   remove the temporary directory created by `parquet_converter.py convert` after merging the chunks in it, along with its contents
"""


from docopt import docopt, DocoptExit
from math import ceil, floor
from memory_profiler import profile
import numpy as np
import os
import pandas as pd
from pandas_plink import read_plink1_bin
from parquet_converter_exceptions import *
from pathlib import Path
from psutil import virtual_memory
import pyarrow as pa
import pyarrow.parquet as pq
import re
import shutil
import sys
import traceback
import xarray as xr


class InstantPrinter:
	"""Copy stdout to a log file

	When running on some clusters, you don\'t get the program output
	until the job is finished.
	For long and risky jobs, this is annoying.

	Attributes
	----------
	logfile: str
	    path to the log file
	printed_lines: list
	    list of strings printed so far

	Methods
	-------
	write(msg):
	    print msg to stdout and write it to the log file
	"""
	def __init__(self, logfile = 'pqconv.log'):
		"""Make an InstantPrinter

		Arguments
		---------
		logfile: str, optional
		    path to the log file.
		    DEFAULT: pqconv.log
		"""
		self.logfile = logfile
		self.printed_lines = []
		if os.path.isfile(self.logfile):
			os.remove(self.logfile)
	def write(self, msg, end = None):
		"""Print msg to stdout and write it to the log file"""
		if end is not None:
			print(msg, end = end)
		else:
			print(msg)
		f = open(self.logfile, 'a')
		f.write(str(msg) + '\n')
		f.close()
		self.printed_lines.append(str(msg))


class GenoMeta:
	"""Describes variants from BIM files or samples from FAM files

	GenoMeta = Genomic Metadata

	Attributes
	----------
	fileindex: dict
	    each key is item ID (item = variant or sample), and the value
	    is a set of indexes of all BED files that contain genotypes
	    for this item
	filelist: list
	    BIM or FAM filepaths corresponding to the indexes of BED files
	bedlist: list
	    BED filepaths corresponding to the indexes of BED files
	df: pandas.DataFrame
	    full BIM or FAM dataframe (from all files combined)
	split: list or 1-D numpy.ndarray
	    sequence of sample or variant indexes to split the genotype data
	    into batches of samples or chunks of variants, respectively
	    (depending on whether the current instance is FAM or BIM).
	    self.split[1:-1] is supposed to be fed to numpy.split(..),
	    self.split[0] is always 0 and self.split[-1] is always
	    the total number of items (samples or variants)
	dropped: list
	    list of item IDs dropped from self.df at user's request
	_start_index: int
	    last index of the filelist, in case we decide to append more
	    files to an existing GenoMeta object

	Methods
	-------
	append(filelist_add, bedlist_add, dropdata = None, idcol = 1):
	    add data from more BIM or FAM files
	set_batches(batching_column = None, nbatches = 1):
	    create self.split (for GenoMeta instances with FAM data)
	set_chunks(chunksize = None):
	    create self.split (for GenoMeta instances with BIM data)
	write(output_files, writer1 = None):
	    write self.df to one or many BIM or FAM files,
	    depending on self.split
	"""

	def __init__(self, filelist_init, bedlist_init,
	                   drop_init = None, make_index = False):
		"""Construct a GenoMeta object

		Arguments
		---------
		filelist_init: list
		    BIM or FAM filepaths to: 1) initialize self.filelist with
		    and 2) read self.df from
		bedlist_init: list
		    BED filepaths to initialize self.bedlist with;
		    assumed to be already aligned with filelist_init
		drop_init: pandas.DataFrame, optional
		    dataframe with rows in format (col, val); for each pair,
		    values equal to val will be dropped from column col
		    (starting from 0) in self.df.
		    DEFAULT: None
		make_index: bool, optional
		    whether to create self.fileindex.
		    DEFAULT: False
		"""
		if make_index:
			self.fileindex = {}
		else:
			self.fileindex = None
		self.filelist = []
		self.bedlist = []
		self.df = None
		self.split = None
		self.dropped = []
		self._start_index = 0
		self.append(filelist_init, bedlist_init, dropdata = drop_init)

	def append(self, filelist_add, bedlist_add, dropdata = None, idcol = 1):
		"""Add more data to a GenoMeta object

		Arguments
		---------
		filelist_add: list
		    same as filelist_init in self.__init__(..)
		bedlist_add: list
		    same as bedlist_init in self.__init__(..)
		dropdata: pandas.DataFrame, optional
		    same as drop_init in self.__init__(..)
		idcol: int, optional
		    index of the ID column in the BIM/FAM files, starts from 0.
		    DEFAULT: 1 (variant ID and IID in BIM and FAM formats,
		                respectively)
		"""
		if dropdata is None:
			dropdata = []
		else:
			dropdata = list(dropdata.itertuples(index=False, name=None))

		for i, filepath in enumerate(filelist_add):
			df_ = pd.read_table(filepath, header = None,
			                    delimiter = r'\s+', dtype = 'object')
			for col, val in dropdata:
				col1 = int(col)
				self.dropped.extend(list(df_.loc[df_[col1] == val, idcol]))
				df_ = df_.loc[df_[col1] != val]
			self.dropped.sort()
			if self.df is None:
				self.df = df_.copy(deep = True)
			else:
				self.df = pd.concat((self.df, df_),
				                    ignore_index = True)
				self.df.drop_duplicates(subset = idcol, inplace = True)
			# If init was called with make_index == True:
			if self.fileindex is not None:
				# Index of the current file:
				i_current = self._start_index + i
				# Look at each sample/variant id,
				# update the set of files where it's present
				for id_ in df_[idcol]:
					if id_ in self.fileindex:
						self.fileindex[id_].add(i_current)
					else:
						self.fileindex[id_] = {i_current}
		self.filelist.extend(filelist_add)
		self.bedlist.extend(bedlist_add)
		self._start_index += i+1
		# (Re)initialize self.split if a proper split hasn't been made yet
		if self.split is None or len(self.split) <= 2:
			self.split = np.array([0, self.df.shape[0]]).astype('int32')
		else:
			self.split = np.append(self.split, self.df.shape[0])

	def set_batches(self, batching_column = None, nbatches = 1):
		"""Create a splitting into sample batches (FAM instances only)

		In other words, generate self.split.

		Arguments
		---------
		batching_column: int, optional
		    index of the column in the FAM files (between 0 and 5)
		    to determine the split (so that samples with the same value
		    in batching_column end up in the same batch, if possible).
		    DEFAULT: None (split into equally-sized batches instead)
		nbatches: int, optional
		    how many batches to split into.
		    DEFAULT: 1 (no split)
		"""
		if nbatches in (0, 1):
			return
		if batching_column is None:
			self.split = np.linspace(0, self.df.shape[0],
			                         num = nbatches+1,
			                         dtype = 'int32')
			assert len(self.split) == nbatches+1
			return

		bcol = self.df.columns[batching_column]
		# Sort self.df by batch column!
		# (AND by index, to avoid ties getting randomly shuffled)
		self.df[self.df.shape[1]] = self.df.index.astype(str)
		icol = self.df.columns[self.df.shape[1]-1]
		self.df[icol] = [x.zfill(len(str(self.df.shape[0]-1)))
		                 for x in self.df[icol]]
		self.df.sort_values(by = [bcol, icol], inplace = True,
		                    key = lambda x: x.str.lower())
		self.df = self.df.iloc[:, :(self.df.shape[1]-1)]
		self.df.reset_index(drop = True, inplace = True)

		cats = self.df[bcol].unique()  # CATegorieS, order preserved
		if len(cats) >= nbatches:
			catsplit = np.linspace(0, len(cats), num = nbatches+1,
			                       dtype = 'int32')
			cat_batches = np.split(cats, catsplit[1:-1])
			# This is the split to be fed to numpy.split:
			#self.split = np.array([self.df[self.df[bcol].isin(cb)].index[-1]+1
			#                       for cb in cat_batches])[:-1].astype('int32')
			# And the self.split below has to be taken with [1:-1]
			self.split = np.array([0] +
			                      [self.df[self.df[bcol].isin(cb)].index[-1]+1
			                       for cb in cat_batches]
			                      ).astype('int32')
		else:
			# Disregard batching_column
			self.split = np.linspace(0, self.df.shape[0],
			                         num = nbatches+1,
			                         dtype = 'int32')
		assert len(self.split) == nbatches+1

	def set_chunks(self, chunksize = None):
		"""Create a splitting into variant chunks (BIM instances only)

		In other words, generate self.split.

		Arguments
		---------
		chunksize: int, optional
		    number of variants per chunk (last chunk might not be full).
		    DEFAULT: None (not set, no splitting is performed)
		"""
		if chunksize is None:
			return
		nchunks = int(ceil(self.df.shape[0]/chunksize))
		self.split = np.array([i*chunksize for i in range(nchunks)] +
		                      [self.df.shape[0]]).astype('int32')

	def write(self, output_files, writer1 = None):
		"""Write the metadata from self.df to one or more files

		(depending on self.split)

		Arguments
		---------
		output_files: list
		    list of paths for the output files to be written to
		    (with extension); should be the same length as the number
		    of batches/chunks in the current instance
		writer1: object, optional
		    a printer object with a write(..) method.
		    DEFAULT: None (print everything to stdout instead)

		Raises
		------
		ValueError
		    if len(output_files) does not match len(self.split)-1
		"""
		if writer1 is None:
			printer = print
		else:
			printer = writer1.write

		if len(output_files) != len(self.split)-1:
			raise ValueError('Too few output filepaths ' +
			                 'supplied to GenoMeta.write()')
		out_dfs = np.array_split(self.df, self.split[1:-1])
		assert len(output_files) == len(out_dfs)

		nfiles = 0
		for i in range(len(out_dfs)):
			out_dfs[i].to_csv(output_files[i], header = None,
			                  index = False, sep = '\t')
			nfiles += 1
		if nfiles == 0:
			printer('No data; nothing to write')
		else:
			suffx = '' if nfiles == 1 else 's'
			printer('Data written to {} file{}'.format(nfiles, suffx))
		return nfiles


def GetTempDir(dirname, filename_prefix):
	"""Generate temporary directory path in a consistent manner"""
	return os.path.join(dirname, filename_prefix + '_chunks')


class OutFilename:
	"""Output filename generator

	Ensures that filenames are generated and parsed in a consistent way
	by all functions that use them.
	For default attribute values, see OutFilename.__init__(..) docs.

	Attributes
	----------
	prefix: str
	    first part of the filename, possibly with path
	delim: str
	    delimiter between filename parts
	ext: str
	    file format (no dot): 'parquet', 'bim' or 'fam'
	id: dict
	    self.id['batch'] is batch id and self.id['chunk'] is chunk id
	id_pos: dict
	    positions of batch id and self id in the filename;
	    same dict keys as self.id
	zfills: dict
	    number of characters to force the ids to in the filename
	    by applying zfill(..) method to them;
	    same dict keys as self.id
	id_pos_defaults: dict
	    dict with ('parquet', 'bim', 'fam') as keys
	    and default self.id_pos values as values
	name: str
	    the filename itself

	Methods
	-------
	replace_delim():
	    make sure self.delim is not contained in self.prefix
	makename():
	    generate self.name, return self
	set(name1):
	    set self.name to name1, parse it according to self.id_pos,
	    and set all the relevant attibutes accordingly; return self
	"""

	def __init__(self, prefix = 'out', delim = '_', ext = 'parquet',
	                   batch_id = None, chunk_id = None,
	                   batch_zfill = 0, chunk_zfill = 0,
	                   batch_id_pos = None, chunk_id_pos = None):
		"""Construct an OutFilename object

		(set the attributes and call self.makename())

		Arguments
		---------
		prefix: str, optional
		    set self.prefix.
		    DEFAULT: 'out'
		delim: str, optional
		    set self.delim.
		    DEFAULT: '_'
		ext: str, optional
		    set self.ext; if other than 'parquet', 'bim' or 'fam',
		    defaults to 'parquet'.
		    DEFAULT: 'parquet'
		batch_id: int, optional
		    set self.id['batch'].
		    DEFAULT: None (skip batch id)
		chunk_id: int, optional
		    set self.id['chunk'].
		    DEFAULT: None (skip chunk id)
		batch_zfill: int, optional
		    set self.zfills['batch'].
		    DEFAULT: 0
		chunk_zfill: int, optional
		    set self.zfills['chunk'].
		    DEFAULT: 0
		batch_id_pos: int, optional
		    set self.id_pos['batch'].
		    if None, set according to self.id_pos_defaults (hardcoded)
		chunk_id_pos: int, optional
		    set self.id_pos['chunk'].
		    if None, set according to self.id_pos_defaults (hardcoded)
		"""
		if ext not in ('parquet', 'bim', 'fam'):
			ext = 'parquet'
		self.delim = delim
		self.ext = ext

		self.prefix = prefix
		self.replace_delim()

		self.id = {'batch': batch_id, 'chunk': chunk_id}
		self.id_pos = {}
		self.zfills = {'batch': batch_zfill, 'chunk': chunk_zfill}

		# Position attributes (self.id_pos) are most important.
		# Basically the whole point of this class.
		#
		# First, set default values if not specified
		self.id_pos_defaults = {'parquet': {'batch': 1, 'chunk': 2},
		                        'bim': {'batch': None, 'chunk': 1},
		                        'fam': {'batch': 1, 'chunk': None}}
		if batch_id_pos is None:
			self.id_pos['batch'] = self.id_pos_defaults[self.ext]['batch']
		else:
			self.id_pos['batch'] = batch_id_pos
		if chunk_id_pos is None:
			self.id_pos['chunk'] = self.id_pos_defaults[self.ext]['chunk']
		else:
			self.id_pos['chunk'] = chunk_id_pos
		# Make sure the positions of different ids
		# have different values >0 (if defined at all)
		if self.ext in ('parquet', 'bim'):
			if self.id_pos['chunk'] <= 0:
				self.id_pos['chunk'] = 1
		if self.ext in ('parquet', 'fam'):
			if self.id_pos['batch'] <= 0:
				self.id_pos['batch'] = 1
		if None not in self.id_pos.values():
			if  self.id_pos['chunk'] == self.id_pos['batch']:
				self.id_pos['chunk'] =  self.id_pos['batch'] + 1

		_ = self.makename()

	def replace_delim(self):
		"""Make sure the delimiter is not contained in the prefix

		(otherwise everything falls apart)

		Arguments
		---------
		none
		"""
		prefix_name = os.path.basename(self.prefix)
		prefix_path = Path(self.prefix).parent
		repl = '-' if self.delim != '-' else '_'
		self.prefix = os.path.join(prefix_path,
		                           prefix_name.replace(self.delim, repl))

	def makename(self):
		"""Generate self.name

		Arguments
		---------
		none

		Returns
		-------
		self
		"""
		posits = [self.id_pos['batch'], self.id_pos['chunk']]
		posits = [x for x in posits if x is not None]
		nameparts = ['' for i in range(max(posits)+1)]
		nameparts[0] = self.prefix
		for t in ('batch', 'chunk'):
			if None not in (self.id_pos[t], self.id[t]):
				nameparts[self.id_pos[t]]=str(self.id[t]).zfill(self.zfills[t])

		self.name = self.delim.join(nameparts) + '.' + self.ext
		return self

	def set(self, name1):
		"""Set self.name to a given name and parse attributes from it

		Attributes parsed are: prefix, id, zfills.
		Attributes delim, ext, id_pos are NOT parsed from the new name
		and are instead assumed to have been set by the class constructor
		or by hand.

		Arguments
		---------
		name1: str
		    set self.name

		Returns
		-------
		self
		"""
		self.name = name1
		# self.ext needs to be set in init,
		# because that determines initialization of attributes
		#
		# Parse the filename only, without the path!
		# Then attach the path back
		name2   = os.path.basename(name1)
		namedir = Path(name1).parent
		nameparts = re.sub('.'+self.ext+'$', '', name2).split(self.delim)
		# Get & set the ids permitted by self.ext (see init)
		prefix_end = 9999
		for t in ('batch', 'chunk'):
			if self.id_pos[t] is not None:
				try:
					x = nameparts[self.id_pos[t]]
					self.zfills[t] = max(self.zfills[t], len(x))
					self.id[t] = int(x)
				except:
					self.zfills[t] = 0
					self.id[t] = None
				prefix_end = min(prefix_end, self.id_pos[t])
		nameparts = nameparts[:prefix_end]
		self.prefix = os.path.join(namedir, self.delim.join(nameparts))
		self.replace_delim()
		return self


class ExistingOutfiles:
	"""Keeps track of existing output files from ParquetConverter(..)

	Used for verification and merging.

	Attributes
	----------
	temp_dir: str
	    path to the temporary directory where the temporary output
	    files are stored (except FAM files)
	out_dir: str
	    path to the final directory where the final output is written
	    (FAM files are written here directly)
	prefix_base: str
	    output filename prefix without path
	tmpfiles: list
	    list of temporary Parquet files as OutFilename objects
	    (with full paths)
	tmpfiles_bim: list
	    list of temporary BIM files as OutFilename objects
	    (with full paths)
	chunks_zfill: int
	    number of digits in chunk IDs in the names of existing files
	chunks_found: int
	    how many variant chunks existing files represent
	last_chunk: int
	    ID of the last chunk found in existing files
	batches_found: int
	    how many sample batches existing files represent
	last_batch: int
	    ID of the last batch found in existing files

	Methods
	-------
	check():
	    check file existence and set attributes
	check_bim(current_gm_bim):
	    check the contents of existing BIM files, make sure they don't
	    share any variants with current_gm_bim
	check_fam(current_gm_fam):
	    check the contents of existing FAM files, make sure they have
	    the same samples and batches as current_gm_fam
	rename(nchunks):
	    expand the number of digits in the chunk IDs in the names of
	    existing files, in case self.last_chunk+nchunks is longer than
	    the current limit
	"""

	def __init__(self, temp_dir, out_dir, prefix_base):
		"""Set paths to directories and filename prefix

		Arguments
		---------
		temp_dir: str
		    set self.temp_dir
		out_dir: str
		    set self.out_dir
		prefix_base: str
		    set self.prefix_base
		"""
		self.temp_dir = temp_dir
		self.out_dir = out_dir
		self.prefix_base = prefix_base

	def check(self):
		"""Examine files in the directory set by self.temp_dir

		Make sure each Parquet file there has a corresponding BIM file,
		and set attributes.

		Attributes set by this method:
		  tmpfiles, tmpfiles_bim (lists of filepaths),
		  chunks_found,  last_chunk, chunks_zfill,
		  batches_found, last_batch.

		Arguments
		---------
		none

		Raises
		------
		FileNotFoundError
		    if a Parquet file does not have a corresponding BIM file
		"""
		tmpfiles1 = sorted(os.listdir(self.temp_dir))
		tmpfiles1 = [f for f in tmpfiles1
		             if (os.path.isfile(os.path.join(self.temp_dir, f)) and
		                 f.startswith(self.prefix_base) and
		                 f.endswith('.parquet'))]
		tmpfiles1_bim = []
		chunks_found_set   = set()
		batches_found_set  = set()
		self.tmpfiles      = []
		self.tmpfiles_bim  = []
		self.chunks_zfill  =  0
		self.chunks_found  =  0
		self.last_chunk    = -1
		self.batches_found =  0
		self.last_batch    = -1
		for f in tmpfiles1:
			fn = OutFilename().set(os.path.join(self.temp_dir, f))
			c = fn.id['chunk']
			b = fn.id['batch']
			if None in (c, b):
				continue
			self.tmpfiles.append(fn)
			# Check if the BIM file is there
			fn_bim = OutFilename(prefix = os.path.join(self.temp_dir,
			                                           self.prefix_base),
			                     ext = 'bim', chunk_id = c,
			                     chunk_zfill  = fn.zfills['chunk'])
			if fn_bim.name not in tmpfiles1_bim:
			# ^ to make sure we don't end up with a bunch of duplicated
			# BIM files in self.tmpfiles_bim
				if not os.path.isfile(fn_bim.name):
					err_msg = ('File {} '.format(fn.name) +
					           'does not have a corresponding BIM file: ' +
					           'file {} does not exist'.format(fn_bim.name))
					raise FileNotFoundError(err_msg)
				tmpfiles1_bim.append(fn_bim.name)
				self.tmpfiles_bim.append(fn_bim)
			# Add the chunk & batch
			self.chunks_zfill = max(self.chunks_zfill, fn.zfills['chunk'])
			chunks_found_set.add(c)
			batches_found_set.add(b)
			self.chunks_found  = len(chunks_found_set)
			self.batches_found = len(batches_found_set)
			self.last_chunk = max(self.last_chunk, c)
			self.last_batch = max(self.last_batch, b)

	def check_bim(self, current_gm_bim):
		"""Examine the contents of existing BIM files in self.temp_dir

		Compare them to the new BIM data about to be written to that
		directory (current_gm_bim) to make sure they don't share any
		variants (won't be able to merge the chunks otherwise).

		Arguments
		---------
		current_gm_bim: GenoMeta
		    GenoMeta object with the new BIM data to compare existing
		    BIM files to

		Raises
		------
		DuplicateVariantsException
		    if some variants from the new BIM data are already present
		    in the existing BIM files
		"""
		df_bim = GenoMeta([fn.name for fn in self.tmpfiles_bim], []).df
		dupes = current_gm_bim.df.merge(df_bim, how = 'inner',
		                                indicator = False).shape[0]
		if dupes > 0:
			err_msg = ('Found {} common variants '.format(dupes) +
			           'between the new data and existing BIM files')
			raise DuplicateVariantsException(err_msg)

	def check_fam(self, current_gm_fam):
		"""Examine the contents of existing FAM files in self.temp_dir

		Make sure they match the new FAM data supplied in current_gm_fam
		(won't be able to merge the chunks otherwise).

		Arguments
		---------
		current_gm_fam: GenoMeta
		    GenoMeta object with the new FAM data to compare existing
		    FAM files to

		Raises
		------
		BadBatchesException
		    if the sample batches in the new data do not match
		    the batches in the existing FAM files
		"""
		tmpfiles_fam1 = sorted(os.listdir(self.out_dir))
		tmpfiles_fam1 = [f for f in tmpfiles_fam1
		                 if (f.startswith(self.prefix_base) and
		                     f.endswith('.fam'))]
		tmpfiles_fam1 = [os.path.join(self.out_dir, f)
		                 for f in tmpfiles_fam1]
		tmpfiles_fam1 = [f for f in tmpfiles_fam1 if os.path.isfile(f)]
		err_msg = ('The current sample batches ' +
		           'are inconsistent with existing files')
		if len(tmpfiles_fam1) != len(current_gm_fam.split)-1:
			raise BadBatchesException(err_msg)
		current_fams = np.array_split(current_gm_fam.df,
		                              current_gm_fam.split[1:-1])
		if len(tmpfiles_fam1) != len(current_fams):
			raise BadBatchesException(err_msg)
		for i, f in enumerate(tmpfiles_fam1):
			df_fam = GenoMeta([f], []).df.reset_index(drop = True)
			current_fams[i].reset_index(drop = True, inplace = True)
			if not df_fam.equals(current_fams[i]):
				raise BadBatchesException(err_msg)

	def rename(self, nchunks):
		"""Update number of digits in chunk numbers of existing filenames

		If the new number of the last variant chunk (after new data
		have been written) will have more digits than the old one,
		rename existing files to update the number of digits.

		Arguments
		---------
		nchunks: int
		    number of variant chunks in the new data
		"""
		if self.last_chunk == -1:
			return
		maxlen_new = len(str(self.last_chunk + nchunks))
		if maxlen_new > self.chunks_zfill:
			self.chunks_zfill = maxlen_new
			tmpfiles_all = {'parquet': self.tmpfiles,
			                'bim':     self.tmpfiles_bim}
			for ext1 in ('parquet', 'bim'):
				for fn in tmpfiles_all[ext1]:
					fpath_old = fn.name
					fn.zfills['chunk'] = self.chunks_zfill
					fn = fn.makename()
					fpath_new = fn.name
					shutil.move(fpath_old, fpath_new)
					assert os.path.isfile(fpath_new)
		# Final checks
		zfs  = [fn.zfills['chunk'] for fn in self.tmpfiles]
		zfs += [fn.zfills['chunk'] for fn in self.tmpfiles_bim]
		assert all([x == self.chunks_zfill for x in zfs])


def Sorter(col, ref):
	"""Sort col in the same order as ref"""
	return col.map({item: order for order, item in enumerate(ref)})



def GetGenotypes(gm_bim, gm_fam, variants = None, samples = None,
                 g_dtype = 'int32', writer1 = None):
	"""Extract and transpose genotypes using a GenoMeta object

	Arguments
	---------
	gm_bim: GenoMeta
	    object with BIM data and BIM + BED filepaths for extraction
	gm_fam: GenoMeta
	    object with FAM data and FAM + BED filepaths for extraction
	variants: list, optional
	    list of variant IDs to extract.
	    DEFAULT: None (all variants)
	samples: list, optional
	    list of sample IDs to extract.
	    DEFAULT: None (all samples)
	g_dtype: str, optional
	    numpy dtype to convert the extracted genotypes to.
	    DEFAULT: 'int32'
	writer1: object, optional
	    a printer object with a write(..) method.
	    DEFAULT: None (print everything to stdout instead)

	Returns
	-------
	genotypes: numpy.ndarray
	    transposed genotypes, shape (variants, samples);
	    None if an exception occurred during extraction of genotypes.
	bad_file: str
	    path to a BED filename if an exception occurred while processing
	    data read from that file, None otherwise
	e: Exception
	    an exception if it occurred during extraction of genotypes,
	    None otherwise

	Raises
	------
	BEDMismatchException
	    if gm_bim and gm_fam do not have the exact same list of BEDs
	"""
	if writer1 is None:
		printer = print
	else:
		printer = writer1.write
	if variants is None:
		variants = list(gm_bim.df[1])
	if samples  is None:
		samples  = list(gm_fam.df[1])

	if gm_bim.bedlist == gm_fam.bedlist:
		bednums_all = gm_bim.bedlist
	else:
		raise BEDMismatchException('Inconsistent BED file lists')

	# Combine file sets for all selected variants and samples
	# to find which of those BED files we need to read
	bednums_bim = [gm_bim.fileindex[v] for v in variants]
	bednums_fam = [gm_fam.fileindex[s] for s in samples]
	bednums = set.union(*(bednums_bim + bednums_fam))
	bednums = sorted(list(bednums))

	genotypes = None
	done_with_files = False
	try:
		for i in bednums:
			G = read_plink1_bin(bednums_all[i],
			                    gm_bim.filelist[i],
			                    gm_fam.filelist[i],
			                    verbose = False)
			G_snps    = list(G['snp'].values)
			G_samples = list(G['sample'].values)
			# Assigning snip ids to the variant field.
			# The 'variant{i}' pattern is too ambiguous
			# and messes everything up.
			G = G.assign_coords(variant = G_snps)
			G = G.sel(variant = [x for x in G_snps
			                     if (x in variants and
			                         x not in gm_bim.dropped)])
			G = G.sel(sample  = [x for x in G_samples
			                     if (x in samples and
			                         x not in gm_fam.dropped)])
			G = G.fillna(9.0)
			if genotypes is not None:
				genotypes = xr.merge([genotypes, G]
				 #, compat='override' <- this actually breaks everything
				   ).to_array().squeeze().drop_vars('variable').rename('genotype')
				# Fun fact: this automatically picks a reasonable
				# chunksize for the resulting dask array.
				# (To be clear, that's a different kind of chunk.)
			else:
				genotypes = G.copy(deep = True)

		done_with_files = True
		# The result has to be checked for nans!
		# Those would be gaps in the merged array:
		# pairs (sample, variant) not present in the input files,
		# not even as explicitly specified missing genotypes.
		# This is also why we fill missing genotypes BEFORE merging.
		if xr.DataArray.isnull(genotypes).any():
			printer('\nWarning: Genotypes for some variant+sample ' +
			        'pairs are not present in the provided BED files')
			genotypes = genotypes.fillna(9.0)
			# Could this eat memory?..
			# Maybe we should just shut up and fill all missing values with 9s.

		genotypes = genotypes.astype(g_dtype)
		# (can't do this before filling all missing values with G.fillna(),
		#  since nans are floats and get replaced with zeros)

		# Sort the chunk to align it with combined BIM/FAM data
		# (just in case)
		variants = [x for x in variants
		            if x in genotypes['variant'].values]
		samples  = [x for x in samples
		            if x in genotypes['sample'].values]
		meta_v = pd.DataFrame(data = {field: genotypes[field].values
		                              for field in ('variant','snp')})
		meta_v.sort_values(by = 'snp', inplace = True,
		                   key = lambda x: Sorter(x, variants))
		meta_s = pd.DataFrame(data = {field: genotypes[field].values
		                              for field in ('sample', 'iid')})
		meta_s.sort_values(by = 'iid', inplace = True,
		                   key = lambda x: Sorter(x, samples))
		genotypes = genotypes.sel(variant = list(meta_v['variant']),
		                          sample  = list(meta_s['sample']))

		genotypes = genotypes.values.T  # Might eat your RAM
		assert genotypes.shape == (len(variants), len(samples))
		return genotypes, None, None

	except Exception as e:
		if done_with_files:
			return None, None, e
		else:
			return None, bednums_all[i], e



def ParquetConvert(filecsv = None, out = 'out',
                   n_batches = 10, batch_col = None, chunk_size = None,
                   drop_inds = None, drop_snps = None,
                   writer = None):
	"""Convert a PLINK dataset to Parquet

	The input dataset may consist of any number of BED files,
	each accompanied by a BIM file and a FAM file.

	The output dataset will be written to one or multiple Parquet files,
	each likewise accompanied by a BIM file and a FAM file.
	If n_batches > 1, the dataset will be split into that many sample
	batches.
	If chunk_size is given and is less than the number of variants
	in the input data, or if all variants don't fit into available RAM,
	each batch will be further split into chunks of variants that can
	later be merged with ParquetMerge(..).

	If the input dataset is split into multiple BED + BIM + FAM
	filesets by variants, the function can be called multiple times
	for different subset of such a dataset, so that they can be
	processed in several runs rather than all at once.
	With that in mind, the output files are actually written to
	a temporary subdirectory within the output directory pointed to by
	the out argument, rather than directly to that output directory.
	ParquetMerge(..) can then be used to convert those temporary files
	into proper output files, so it's good to always call that function
	after ParquetConvert(..) has done its job, even if every batch has
	only one chunk.

	Arguments
	---------
	filecsv: str
	    path to a CSV file with paths to input BED, BIM and FAM files
	    as three columns in this order, with no header
	    (required)
	out: str, optional
	    intended path and prefix for the final output files;
	    preliminary output files made by this function will be written
	    to a temporary subdirectory within the same directory instead.
	    DEFAULT: 'out'
	n_batches: int, optional
	    number of sample batches to split the dataset into
	    (not to be confused with training batches!).
	    DEFAULT: 10
	batch_col: int, optional
	    index of the column (between 0 and 5) in FAM files, to ensure
	    that samples with equal batch_col values end up in the same
	    batch (if possible).
	    DEFAULT: None (split into equal batches instead)
	chunk_size: int, optional
	    size of variant chunks to use when building the dataset,
	    the larger the better.
	    DEFAULT: None (as many as fits into RAM)
	drop_inds: str, optional
	    path to a text file with lines in format 'i,var';
	    for every such line, all samples with value 'var' in column i
	    (starting from 0) will be dropped from the FAM files.
	    DEFAULT: None (keep everything)
	drop_snps: str, optional
	    same as drop_inds, but for variants and BIM files.
	    DEFAULT: None (keep everything)
	writer: object, optional
	    a printer object with a write(..) method.
	    DEFAULT: None (print everything to stdout instead)

	Raises
	------
	TypeError
	    if required arguments are missing from the function call
	ValueError
	    if the data provided have incorrect format
	NoItemsException
	    if the BIM or FAM files provided are empty
	    (i.e. no variants or samples could be found)
	InsufficientRAMException
	    if not even one variant can fit into available RAM
	BadBatchesException
	    if the sample batches in the new data do not match the batches
	    in the existing FAM files
	"""
	if writer is None:
		print1 = print
	else:
		print1 = writer.write
	if filecsv is None:
		raise TypeError('CSV file with paths to input files is required')

	EXTS = ('bed', 'bim', 'fam')

	outprefix_base = os.path.basename(out)
	outdir = Path(out).parent
	if not os.path.isdir(outdir):
		os.makedirs(outdir)
		print1('\nCreated directory {}'.format(outdir))

	files = {}
	files_df = pd.read_csv(filecsv, header = None, dtype = 'object')
	if files_df.shape[0] == 0:
		print1('\nInput file table is empty, nothing to convert.')
		return
	if files_df.shape[1] != len(EXTS):
		raise ValueError('CSV file with input filepaths ' +
		                 'should have {} columns'.format(len(EXTS)))
	for i in range(len(EXTS)):
		files[EXTS[i]] = list(files_df[i])
		files[EXTS[i]].sort()

	drop_dfs = []
	for drop_arg in (drop_inds, drop_snps):
		if drop_arg:
			drop_dfs.append(pd.read_csv(drop_arg, header = None,
			                            dtype = 'object'))
		else:
			drop_dfs.append(None)

	print1('\nReading BIM and FAM data...')
	gmeta = {}
	gmeta['bim'] = GenoMeta(files['bim'], files['bed'],
	                        drop_init = drop_dfs[1],
	                        make_index = True)
	gmeta['fam'] = GenoMeta(files['fam'], files['bed'],
	                        drop_init = drop_dfs[0],
	                        make_index = True)
	for ext in ('bim', 'fam'):
		assert gmeta[ext].filelist == files[ext]
		assert gmeta[ext].bedlist  == files['bed']
	n_snps = gmeta['bim'].df.shape[0]
	n_inds = gmeta['fam'].df.shape[0]
	if n_snps == 0:
		raise NoItemsException('No variants found in the BIM files')
	if n_inds == 0:
		raise NoItemsException('No samples found in the FAM files')

	if n_batches <= 1:
		n_batches = 1
		print1('All samples will be kept in a single batch')
	else:
		print1('The dataset will be split into ' +
		       '{} sample batches'.format(n_batches))
	# And also get the number of digits, for zfilling later
	batch_numlen = len(str(n_batches-1))
	gmeta['fam'].set_batches(batching_column = batch_col,
	                         nbatches = n_batches)
	# Get max sample batch size
	batchsizes = [gmeta['fam'].split[i+1] - gmeta['fam'].split[i]
	              for i in range(gmeta['fam'].split.shape[0]-1)]
	assert 0 not in batchsizes
	assert sum(batchsizes) == n_inds
	assert len(batchsizes) == n_batches
	max_batch_size = max(batchsizes)

	# Find variant chunk size
	#
	# Reducing dtype size saves both RAM and time
	genotypes_dtype = 'int8'
	nbytes = np.dtype(genotypes_dtype).itemsize
	if chunk_size is None:
		# To maximize Parquet conversion, we want to maximize chunk size.
		# So we take as many variants as we can comfortably fit into RAM.
		#
		# Using a rough equation derived by measuring RAM consumption.
		# Assumes that no other big jobs are running on the machine.
		# The 0.9 multiplier added for extra security.
		# If this isn't enough, god help you.
		###
		max_ram = 0.9*virtual_memory().available
		chunk_size = floor(max_ram / (nbytes*(3.0*max_batch_size)))
		###
		chunk_size = min(chunk_size, n_snps)
		if chunk_size == 0:
			raise InsufficientRAMException('Insufficient RAM! ' +
			                               'Try more batches')
	n_chunks = int(ceil(n_snps/chunk_size))
	assert n_chunks > 0
	print1('\nUsing {} chunks of size {} '.format(n_chunks, chunk_size) +
	       '(last chunk might not be full)')
	gmeta['bim'].set_chunks(chunksize = chunk_size)

	chunk_num = -1
	chunk_numlen = 0
	# When starting from chunk #0, chunk_num starts from -1 rather than 0,
	# because I want to increment it at the start of the while-loop
	# so it doesn't get lost at the bottom of it.
	#
	# In case we have multiple chunks, we have to write each chunk
	# to its own parquet file before combining them.
	# If that's the case, we need a temp dir for those smaller files.
	#
	# If we're doing several convert runs with a merge run in the end,
	# then the temp dir already exists,
	# and we need to be careful not to overwrite anything there.
	# So we need to check if there are any existing chunks in the tmpdir,
	# and maybe we should add that to chunk_num
	# and resume the count from chunk_num+1.
	tmpdir = GetTempDir(outdir, outprefix_base)
	if os.path.isdir(tmpdir):
		print1('Temporary directory {} already exists.'.format(tmpdir)+
		       ' Checking...')
		eo = ExistingOutfiles(tmpdir, outdir, outprefix_base)
		eo.check()
		if eo.batches_found != n_batches or eo.last_batch != n_batches-1:
			err_msg = ('The current number of batches ' +
			           'is inconsistent with existing files')
			raise BadBatchesException(err_msg)
		eo.check_bim(gmeta['bim'])
		eo.check_fam(gmeta['fam'])
		eo.rename(n_chunks)
		if eo.last_chunk == -1:
			print1('No existing chunks found, starting from chunk 0')
		else:
			print1('{} existing chunks found'.format(eo.chunks_found))
		chunk_num = eo.last_chunk
		chunk_numlen = eo.chunks_zfill
	else:
		os.mkdir(tmpdir)
	if chunk_numlen == 0:
		chunk_numlen = len(str(n_chunks-1))

	# Write the BIM(s) and FAM(s).
	# BIMs are gonna be chunk BIMs: written to tmpdir
	# and later merged along with parquets.
	# FAMs are gonna be batch FAMs: written in batches to the final path.
	outfiles = {}
	outfiles['fam'] = [OutFilename(prefix = out, ext = 'fam',
	                               batch_id = i,
	                               batch_zfill = batch_numlen).name
	                   for i in range(n_batches)]
	outfiles['bim'] = [OutFilename(prefix = os.path.join(tmpdir,
	                                                     outprefix_base),
	                               ext = 'bim', chunk_id = chunk_num+i+1,
	                               chunk_zfill = chunk_numlen).name
	                   for i in range(n_chunks)]
	fams_exist = [os.path.isfile(f) for f in outfiles['fam']]
	bims_exist = [os.path.isfile(f) for f in outfiles['bim']]
	assert (not any(fams_exist)) or all(fams_exist)
	assert  not any(bims_exist)
	for ext in ('bim', 'fam'):
		if ext == 'fam' and all(fams_exist):
			continue
		print1('\nWriting combined {} data...'.format(ext.upper()))
		_ = gmeta[ext].write(outfiles[ext], writer1 = writer)

	#
	# This is where I'd put an if-else for normalization,
	# IF I HAD ONE
	#
	# Since there's no normalization,
	# iterate over BATCHES first, then over chunks.
	print1('\nConverting genotypes to Parquet...')
	chunk_num_0 = chunk_num
	for i in range(n_batches):
		print1('\nConverting batch #{}'.format(i))
		ind1 = gmeta['fam'].split[i]
		ind2 = gmeta['fam'].split[i+1]
		iids_batch = list(gmeta['fam'].df[1][ind1:ind2])
		#iids_batch = [str(iid) for iid in iids_batch]

		# Now iterate over chunks
		chunk_num = chunk_num_0
		start = 0
		while start < n_snps:
			chunk_snps = list(gmeta['bim'].df[1][start:(start+chunk_size)])
			start += chunk_size
			chunk_num += 1
			print1('Processing chunk #{}'.format(chunk_num), end='\r')

			genos, badfile, exc = GetGenotypes(gmeta['bim'], gmeta['fam'],
			                                   variants = chunk_snps,
			                                   samples  = iids_batch,
			                                   g_dtype  = genotypes_dtype,
			                                   writer1  = writer)
			if exc is not None:
				if badfile is None:
					insert1 = ' and'
					insert2 = ''
				else:
					insert1 = ','
					insert2 = ' and file {}'.format(badfile)
				print1('\n\nWhile converting' +
				       ' batch #{}{}'.format(i, insert1) +
				       ' chunk #{}{}'.format(chunk_num, insert2) +
				       ', a following exception occurred:')
				print1(''.join(traceback.format_exception(exc)))
				return
			sch = pa.schema([pa.field(iid, pa.from_numpy_dtype(genos.dtype))
			                 for iid in iids_batch])
			genos = pd.DataFrame(genos, columns = iids_batch)
			genos = pa.Table.from_pandas(genos, schema = sch)
			assert genos.shape == (len(chunk_snps), len(iids_batch))
			ofn = OutFilename(prefix = outprefix_base,
			                  batch_id = i, chunk_id = chunk_num,
			                  batch_zfill = batch_numlen,
			                  chunk_zfill = chunk_numlen)
			outfile = os.path.join(tmpdir, ofn.name)
			assert not os.path.isfile(outfile)
			pq.write_table(genos, outfile)

	print1('\nConversion done!\n')



def ParquetMerge(out = None, remove_tmpfiles = False,
                 overwrite = False, writer = None):
	"""Convert temporary output files to proper output files

	Inspect the temporary directory made by ParquetConvert(..),
	merge variant chunks within each sample batch if necessary,
	and write the result to final output files
	(a Parquet + BIM + FAM set for every batch).

	Arguments
	---------
	out: str, optional
	    path and prefix for the final output files, for real this time
	    (required); this path should also contain the temporary
	    directory previously made by ParquetConvert(..)
	remove_tmpfiles: bool, optional
	    if True, the temporary files will be removed as they are
	    processed, along with the temporary directory itself in the end.
	    DEFAULT: False
	overwrite: bool, optional
	    if True, overwrite existing final output files;
	    otherwise, throw an error if the files exist.
	    DEFAULT: False
	writer: object, optional
	    a printer object with a write(..) method
	    DEFAULT: None (print everything to stdout instead)

	Raises
	------
	TypeError
	    if required arguments are missing from the function call
	FileNotFoundError
	    if the temporary directory does not exist
	BadBatchesException
	    if the sample batches in the new data do not match the batches
	    in the existing FAM files
	FileExistsException
	    if a file with the specified output path already exists
	    and overwrite is False
	"""
	if writer is None:
		print1 = print
	else:
		print1 = writer.write
	if out is None:
		raise TypeError('Output prefix is required')

	outprefix_base = os.path.basename(out)
	outdir = Path(out).parent

	tmpdir = GetTempDir(outdir, outprefix_base)
	if not os.path.isdir(tmpdir):
		raise FileNotFoundError('Temporary directory not found')
	print1('Checking the temporary directory: {}'.format(tmpdir))
	eo = ExistingOutfiles(tmpdir, outdir, outprefix_base)
	eo.check()
	if len(eo.tmpfiles) == 0:
		print1('No relevant files in the temporary directory, ' +
		       'nothing to merge.\n')
		return
	if eo.batches_found == 0:
		raise BadBatchesException('No batches found in {}'.format(tmpdir))
	print1('{} batches found'.format(eo.batches_found))

	# Might need to set the highest possible thrift limits
	# for reading our parquet files
	int32_t_MAX = 2**31-1

	batches_seen = []
	for fn in eo.tmpfiles:
		b = fn.id['batch']
		if b in batches_seen:
			continue
		batches_seen.append(b)
		print1('\nProcessing batch ' + str(b).zfill(fn.zfills['batch']))

		# Set outfile name, check if the file exists.
		# Initializing filename with fam ext
		# to ensure appropriate filename format,
		# then swapping ext to parquet.
		ofn = OutFilename(prefix = out, ext = 'fam',
		                  batch_id = fn.id['batch'],
		                  batch_zfill = fn.zfills['batch'])
		ofn.ext = 'parquet'
		outfile = ofn.makename().name
		if os.path.isfile(outfile) and not overwrite:
			raise FileExistsException('File {} already '.format(outfile) +
			                       'exists and overwrite is set to False')

		tmpfiles_batch = [fn1 for fn1 in eo.tmpfiles
		                  if fn1.id['batch'] == fn.id['batch']]
		assert len(tmpfiles_batch) > 0

		tfile = tmpfiles_batch[0].name

		# Special case of a single chunk
		if len(tmpfiles_batch) == 1:
			if remove_tmpfiles:
				dest = shutil.move(tfile, outfile)
			else:
				dest = shutil.copyfile(tfile, outfile)
			assert os.path.isfile(dest)
			print1('Batch written to file ' + str(dest))
			continue

		first_table = pq.read_table(tfile,
		                    thrift_string_size_limit = int32_t_MAX,
		                    thrift_container_size_limit = int32_t_MAX)
		sch0 = first_table.schema
		with pq.ParquetWriter(outfile, sch0,
		                      version = '2.6',
		                      compression = 'gzip',
		                      use_dictionary = True,
		                      data_page_size = 2097152,  # 2 MiB
		                      write_statistics = True) as pqwriter:
			for fn1 in tmpfiles_batch:
				print1('Writing chunk #{}'.format(fn1.id['chunk']),
				       end='\r')
				if first_table is None:
					pqt = pq.read_table(fn1.name,
					            thrift_string_size_limit = int32_t_MAX,
					            thrift_container_size_limit = int32_t_MAX)
					assert pqt.schema == sch0
					pqwriter.write_table(pqt)
				else:
					pqwriter.write_table(first_table)
					first_table = None
				if remove_tmpfiles:
					os.remove(fn1.name)
		print1('\nBatch written to file ' + str(outfile))

	print1('\nWriting combined BIM data...')
	gmeta_bim = GenoMeta([fn_bim.name for fn_bim in eo.tmpfiles_bim], [])
	_ = gmeta_bim.write([out+'.bim'], writer1 = writer)

	if remove_tmpfiles:
		print1('\nDeleting the temporary directory {}...'.format(tmpdir))
		shutil.rmtree(tmpdir)

	print1('\nMerging done!\n')



if __name__ == '__main__':

	try:
		arguments = docopt(__doc__, version='Parquet Converter 2.0.0')
	except DocoptExit:
		print('\nInvalid command. ' +
		      'Run \'python3 parquet_converter.py --help\' ' +
		      'for more information.\n')
		sys.exit(1)
	for x in list(arguments.keys()):
		xnew = x.split('-')[-1]
		arg = arguments.pop(x)
		arguments[xnew] = arg

	GCAE_DIR = Path(__file__).resolve().parents[1]

	logfile = arguments.pop('logfile')
	if logfile:
		if not os.path.isabs(logfile):
			logfile = os.path.join(GCAE_DIR, logfile)
		w = InstantPrinter(logfile = logfile)
		print1 = w.write
	else:
		w = None
		print1 = print

	memprof = arguments.pop('memprof')
	if memprof:
		ParquetConvert = profile()(ParquetConvert)
		ParquetMerge   = profile()(ParquetMerge)

	# Mode and mode-specific arguments
	convert = arguments.pop('convert')
	merge   = arguments.pop('merge')
	assert convert or merge  # should be taken care of by docopt but idk
	if convert and merge:
		print1('\nBoth convert and merge are given, which one?\n')
		sys.exit(1)
	args_all = list(arguments.keys())
	args_common  = ['out']
	args_merge   = ['overwrite', 'remove_tmpfiles']
	if convert:
		# convert mode
		mode_args = [arg for arg in args_all
		             if arg in args_common or arg not in args_merge]
	else:
		# merge mode
		mode_args = args_common + args_merge
	# Discard args irrelevant to the mode
	for arg in args_all:
		if arg not in mode_args:
			_ = arguments.pop(arg)

	# Handle path arguments
	print1('\nProject directory: ' + str(GCAE_DIR))
	path_args = ['filecsv', 'out', 'drop_inds', 'drop_snps']
	path_args = [arg for arg in path_args if arg in mode_args]
	for arg in path_args:
		if arguments[arg]:
			if not os.path.isabs(arguments[arg]):
			# note: ^ works with regex dirnames
				arguments[arg] = os.path.join(GCAE_DIR, arguments[arg])
	if not arguments['out']:
		arguments['out'] = 'out'

	# Handle integer arguments
	if convert:
		int_args = ('n_batches', 'chunk_size', 'batch_col')
		for arg in int_args:
			if arguments[arg]:
				try:
					arguments[arg] = int(arguments[arg])
				except:
					print1('\nInvalid value of argument ' +
					       '{}: {}\n'.format(arg, arguments[arg]))
					sys.exit(1)
		if not arguments['n_batches']:
			arguments['n_batches'] = 10
		arg = 'batch_col'
		if arguments[arg]:
			if arguments[arg] not in range(6):
				print1('\nInvalid value of argument ' +
				       '{}: {} '.format(arg, arguments[arg]) +
				       '(should be between 0 and 5)\n')
				sys.exit(1)

	print1('\nPARAMETERS:\n' +
	       'mode: {}\n'.format('convert' if convert else 'merge') +
	       '\n'.join(['{}: {}'.format(arg, arguments[arg])
	                  for arg in arguments]) + '\n' +
	       'memprof: {}\n'.format(memprof) +
	       'logfile: {}\n'.format(logfile))

	if convert:
		ParquetConvert(**arguments, writer = w)

	if merge:
		ParquetMerge(**arguments, writer = w)


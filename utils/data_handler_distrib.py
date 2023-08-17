import glob
import numpy as np
import os
from pathlib import Path
import pyarrow.parquet as pq
import tensorflow as tf
import utils.normalization as normalization


def _bytes_feature(value):
	if isinstance(value, type(tf.constant(0))):
		value = value.numpy()
	return tf.train.Feature(bytes_list=tf.train.BytesList(value=[value]))

def _int64_feature(value):
	return tf.train.Feature(int64_list=tf.train.Int64List(value=[value]))

def make_example(genos_, indpop_):
	features_ = {
		'len': _int64_feature(genos_.shape[0]),
		'snp': _bytes_feature(tf.io.serialize_tensor(genos_)),
		'ind': _bytes_feature(tf.io.serialize_tensor(indpop_[0])),
		'pop': _bytes_feature(tf.io.serialize_tensor(indpop_[1]))
	}
	return tf.train.Example(features=tf.train.Features(feature=features_))

def decode_example(x, geno_dtype=np.float32):
	features_ = {
		'len': tf.io.FixedLenFeature([], tf.int64),
		'snp': tf.io.FixedLenFeature([], tf.string),
		'ind': tf.io.FixedLenFeature([], tf.string),
		'pop': tf.io.FixedLenFeature([], tf.string)
	}
	example = tf.io.parse_single_example(x, features_)
	genos_ = tf.io.parse_tensor(example['snp'],
	                            out_type=tf.as_dtype(geno_dtype))
	ind_ = tf.io.parse_tensor(example['ind'], out_type=tf.string)
	pop_ = tf.io.parse_tensor(example['pop'], out_type=tf.string)
	return genos_, [ind_,pop_]


class data_generator_distrib:
	"""docstring"""
	# TODO: go over all attrs used, make sure they exist & are set

	def __init__(self, filebase,
	             normalization_mode="genotypewise01",
	             normalization_options={"flip": False,
	                                    "missing_val": 0.0},
	             impute_missing=True):
		self.filebase = filebase
		self.missing_val = normalization_options["missing_val"]
		self.normalization_mode = normalization_mode
		self.normalization_options = normalization_options
		self.impute_missing = impute_missing
		# maybe take & set more attrs here
		self.tf_dataset = None

		self.get_ind_pop_list()
		self.get_n_markers()
		self._define_samples()
		# no loading here

	def get_ind_pop_list(self):
		self.ind_pop_list = np.empty(shape=(0,2), dtype=str)
		fam_files = sorted(glob.glob(self.filebase+"*.fam"))
		for file in fam_files:
			indpop = np.genfromtxt(file, usecols=(1,0), dtype=str)
			self.ind_pop_list = np.concatenate((self.ind_pop_list, indpop), axis=0)

	def get_n_markers(self):
		self.n_markers = 0
		bim_files = sorted(glob.glob(self.filebase+"*.bim"))
		for file in bim_files:
			self.n_markers += len(np.genfromtxt(file, usecols=(1), dtype=str))

	def _define_samples(self):
		self.n_total_samples = len(self.ind_pop_list)
		self.n_train_samples = len(self.ind_pop_list)
		self.n_valid_samples = 0

		self.sample_idx_all   = np.arange(self.n_total_samples)
		self.sample_idx_train = np.arange(self.n_train_samples)
		self.sample_idx_valid = np.arange(self.n_valid_samples)

	def define_validation_set(self, validation_split=0.2, random_state=None):
		# TODO: should be stratified by population in the general case
		self.n_valid_samples = np.floor(self.n_total_samples)
		self.n_train_samples = self.n_total_samples - self.n_valid_samples

		np.random.seed(random_state)
		self.sample_idx_valid = np.random.choice(self.sample_idx_all,
		                                         size=self.n_valid_samples,
		                                         replace=False)
		self.sample_idx_valid = np.sort(self.sample_idx_valid)
		train_idx = np.in1d(self.sample_idx_train, self.sample_idx_valid,
		                    invert=True)
		self.sample_idx_train = self.sample_idx_train[train_idx]


	# The key part (supposedly)
	def generator(self, pref_chunk_size_, training_=True, shuffle_=True):
		# handle data loading
		# https://www.tensorflow.org/guide/data
		# https://stackoverflow.com/q/68164440

		if training_:
			n_samples = self.n_train_samples
			cur_sample_idx = self.sample_idx_train[np.arange(0, n_samples)]
		else:
			n_samples = self.n_valid_samples
			cur_sample_idx = self.sample_idx_valid[np.arange(0, n_samples)]
		if shuffle_:
			cur_sample_idx = tf.random.shuffle(cur_sample_idx)
		cur_sample_idx = tf.cast(cur_sample_idx, tf.int32)

		# TODO: make sure to properly support multiple files, everywhere
		pq_paths = sorted(glob.glob(self.filebase+"*.parquet"))
		int32_t_MAX = 2**31-1
		pqds = pq.ParquetDataset(path_or_paths = pq_paths,
		                         thrift_string_size_limit = int32_t_MAX,
		                         thrift_container_size_limit = int32_t_MAX,
		                         use_legacy_dataset = False)
		# OBS! might not preserve column order. rely on schema instead.
		pqds_schema = pqds.schema
		inds_sch = [entry.name for entry in pqds_schema]
		inds_fam = list(self.ind_pop_list[:,0])
		if inds_sch != inds_fam:
			raise ValueError("Parquet schema inconsistent with FAM files")

		chunk_size = pref_chunk_size_ - pref_chunk_size_ % self.batch_size
		num_chunks = np.ceil(n_samples / chunk_size)

		chunks_read = 0
		while chunks_read < num_chunks:

			start = chunk_size * chunks_read
			end   = chunk_size *(chunks_read+1)
			chunk_idx = cur_sample_idx[start:end,:]
			batches_per_chunk = np.ceil(len(chunk_idx)/self.batch_size)

			inds_to_read = list(self.ind_pop_list[chunk_idx,0])
			chunk = pqds.read(columns = inds_to_read,
			                  use_threads = True,  # TODO: try without
			                  use_pandas_metadata = False)
			sch   = chunk.schema
			chunk = chunk.to_pandas(self_destruct=True).to_numpy(dtype=self.geno_dtype)
			# TODO: if you use float16, other scripts should support that
			chunk = chunk.T
			assert chunk.shape[0] == batches_per_chunk*self.batch_size
			assert chunk.shape[1] == self.n_markers

			chunks_read += 1

			batches_read = 0
			last_batch = False
			while batches_read < batches_per_chunk:

				start_ = self.batch_size * batches_read
				end_   = self.batch_size *(batches_read+1)
				batch = chunk[start_:end_]
				if end_ >= chunk.shape[0]:
					last_batch = True

				batch_idx    = cur_sample_idx[start_:end_]
				batch_indpop = self.ind_pop_list[batch_idx,:]

				batches_read += 1

				yield batch, batch_indpop, last_batch


	def _normalize(self, x, indpop, last_batch):

		missing = tf.where(x == 9)
		a = tf.ones(shape=tf.shape(missing)[0], dtype=x.dtype)
		a = tf.sparse.SparseTensor(indices=missing, values=a,
		                           dense_shape=x.shape)

		if self.impute_missing:
			# TODO: get most_common_genos
			b = tf.gather(self.most_common_genos, indices=missing[:,1])-9
			b = tf.sparse.SparseTensor(indices=missing, values=b,
			                           dense_shape=x.shape)
			x = tf.sparse.add(x, b)

		if self.normalization_mode == "genotypewise01":
			if self.normalization_options["flip"]:
				x = -(x-2)/2
				if not self.impute_missing:
					x = tf.sparse.add(x, a*(3.5-self.missing_val))
			else:
				x = x/2
				if not self.impute_missing:
					x = tf.sparse.add(x, a*(self.missing_val-4.5))

		elif self.normalization_mode in ("standard", "smartPCAstyle"):
			# TODO
			raise NotImplementedError("Only genotypewise01 normalization "+
			                          "method supported for now")

		return x, indpop, last_batch


	def _sparsify(self, x, indpop):
		fraction = 0.0
		# TODO: get fraction from data opts sparsifies
		if not self.missing_mask_input:
			inputs = tf.expand_dims(x, axis=-1)
			return inputs, x, indpop

		mask = tf.experimental.numpy.full(shape=x.shape, fill_value=1,
		                                  dtype=x.dtype)

		probs = tf.random.uniform(shape=x.shape, minval=0, maxval=1)
		where_sparse = tf.where(probs < fraction)
		b = tf.repeat(-1, tf.shape(where_sparse)[0])
		b = tf.sparse.SparseTensor(indices=where_sparse, values=b,
		                           dense_shape=mask.shape)
		mask = tf.sparse.add(mask, b)

		inputs = tf.math.add(tf.math.multiply(x, mask),
		                     -1*self.missing_val*(mask-1))
		inputs = tf.stack([inputs, mask], axis=-1)

		return inputs, x, indpop
		# TODO: do we really need to store both inputs and x? inefficient.
		#       would need to change that in the calling scope first.


	def write_TFRecords(self, dataset, outprefix_,
	                    num_workers_=1, training_=True):
		if training_:
			mode = "train"
			n_samples = self.n_train_samples
		else:
			mode = "valid"
			n_samples = self.n_valid_samples

		outdir = Path(outprefix_).parent
		if not os.path.isdir(outdir):
			os.makedirs(outdir)

		# At least 10 MB in each shard, but at most 10 shards per worker
		# TODO: rethink this. we might end up with A LOT more shards.
		#       how about you remove this upper limit for now
		#       and reinstate it later if needed
		genos_size = n_samples * self.n_markers * self.geno_dtype.itemsize
		total_shards = min(10*num_workers_, np.ceil(genos_size*1e-7))

		batches_per_shard = np.ceil(n_samples/(self.batch_size*total_shards))

		batch_count = 0
		shard_count = 0
		for batch, batch_indpop, last_batch in dataset:
			if batch_count == 0:
				genos    = batch
				indpop   = batch_indpop
			else:
				genos    = tf.concat([genos, batch], axis=0)
				indpop   = tf.concat([indpop, batch_indpop], axis=0)
			batch_count += 1

			if batch_count != batches_per_shard and not last_batch:
				continue

			shard_id = str(shard_count).zfill(len(str(total_shards-1)))
			shard_filename = f"{outprefix_}_{mode}_{shard_id}.tfrecords"
			writer = tf.io.TFRecordWriter(shard_filename)
			# TODO: compression maybe?

			for i in range(genos.shape[0]):
				cur_geno = tf.cast(genos[i,:], tf.as_dtype(self.geno_dtype))
				example = make_example(cur_geno, indpop[i,:])
				writer.write(example.SerializeToString())

			shard_count += 1
			writer.close()


	# Another key part
	def create_tf_dataset(self, pref_chunk_size, outprefix,
	                      num_workers=1, geno_dtype=np.float32,
	                      training=True, shuffle=True, overwrite=False):
		# feed self.generator to tf.data.Dataset.from_generator()
		# possibly with TFRecord: https://stackoverflow.com/q/59458298
		# see also https://www.tensorflow.org/guide/data_performance

		# TODO: how is validation handled with tf.data, exactly? the training arg

		existing_tfr_files = sorted(glob.glob(outprefix+"*.tfrecords"))
		if len(existing_tfr_files) == 0 or overwrite:

			self.geno_dtype = geno_dtype
			gen_outshapes = (
				tf.TensorSpec(shape=(None, self.n_markers),
				              dtype=tf.as_dtype(self.geno_dtype)),
				tf.TensorSpec(shape=(None, 2), dtype=tf.string),
				tf.TensorSpec(shape=(1, 1), dtype=tf.bool)
			)
			gen_args = (pref_chunk_size, training, shuffle)
			# TODO: can't you, like... shard it right away? instead of batching?
			#       note that shards cannot be larger than chunks though
			#       (assuming no dtype shenanigans)
			ds = tf.data.Dataset.from_generator(self.generator,
			                                    output_signature=gen_outshapes,
			                                    args=gen_args)

			# TODO: if training, prep norm scaler
			# (if norm methods other than genotype-wise are applicable)

			ds = ds.map(self._normalize, num_parallel_calls=tf.data.AUTOTUNE)

			# Apparently can't write TFRecords in parallel.
			# TODO: make sure this works with multiprocessing (check in calling scope)
			if "SLURM_PROCID" in os.environ and int(os.environ["SLURM_PROCID"]) == 0:
				self.write_TFRecords(ds, outprefix,
				                     num_workers_=num_workers,
				                     training_=training)

		current_tfr_files = sorted(glob.glob(outprefix+"*.tfrecords"))
		ds = tf.data.from_tensor_slices(filenames=current_tfr_files)

		worker_id = int(os.environ["SLURM_PROCID"])
		ds = ds.shard(num_workers, worker_id)
		ds = ds.interleave(tf.data.TFRecordDataset,
		                   num_parallel_calls=tf.data.AUTOTUNE,
		                   cycle_length=num_workers, block_length=1)
		ds = ds.map(lambda d: decode_example(d, geno_dtype=self.geno_dtype),
		            num_parallel_calls=tf.data.AUTOTUNE)

		ds = ds.prefetch(tf.data.AUTOTUNE)
		ds = ds.batch(self.batch_size//num_workers)
		# Sparsifying changes dataset structure!
		ds = ds.map(self._sparsify,  num_parallel_calls=tf.data.AUTOTUNE)

		self.tf_dataset = ds

	# distribute dataset: https://stackoverflow.com/q/59185729


import glob
import numpy as np
import pyarrow.parquet as pq
import tensorflow as tf
import utils.normalization as normalization


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

		self._get_dims()
		self._define_samples()
		# no loading here

	def _define_samples(self):
		# TODO
		return

	def _get_dims(self):
		"""Count markers and samples"""
		# TODO
		return

	def _normalize(self):
		# TODO, as mapping
		return

	def _sparsify(self, mask, keep_fraction):
		# TODO, as mapping
		return

	def define_validation_set(self, validation_split=0.2):
		# TODO
		return

	# The key part (supposedly)
	def generator(self, pref_chunk_size_, geno_dtype_=np.float16,
	              training_=True, shuffle_=True):
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
		sch0 = pqds.schema
		# TODO: assert consistency with self.ind_pop_list

		chunk_size = pref_chunk_size_ - pref_chunk_size_ % self.batch_size
		num_chunks = np.ceil(n_samples / chunk_size)

		chunks_read = 0
		while chunks_read < num_chunks:

			start = chunk_size * chunks_read
			end   = chunk_size *(chunks_read+1)
			chunk_idx = cur_sample_idx[start:end]
			batches_per_chunk = np.ceil(len(chunk_idx)/self.batch_size)

			# TODO: store ind_pop list as an attr?
			inds_to_read = list(self.ind_pop_list[chunk_idx,0])
			chunk = pqds.read(columns = inds_to_read,
			                  use_threads = True,  # TODO: try without
			                  use_pandas_metadata = False)
			sch   = chunk.schema
			chunk = chunk.to_pandas(self_destruct=True).to_numpy(dtype=geno_dtype_)
			# TODO: if you use float16, other scripts should support that
			chunk = chunk.T
			assert chunk.shape[0] == batches_per_chunk*self.batch_size
			assert chunk.shape[1] == self.n_markers

			chunks_read += 1

			batches_read = 0
			while batches_read < batches_per_chunk:

				start_ = self.batch_size * batches_read
				end_   = self.batch_size *(batches_read+1)
				batch = chunk[start_:end_]

				batch_idx  = cur_sample_idx[start_:end_]
				batch_inds = self.ind_pop_list[batch_idx,:]

				batches_read += 1

				yield batch, batch_inds


	# Another key part
	def create_tf_dataset(self, pref_chunk_size, geno_dtype=np.float16,
	                      training=True, shuffle=True):
		# feed self.generator to tf.data.Dataset.from_generator()
		# possibly with TFRecord: https://stackoverflow.com/q/59458298
		# see also https://www.tensorflow.org/guide/data_performance

		gen_outshapes = (
			tf.TensorSpec(shape=(None, self.n_markers),
			              dtype=tf.as_dtype(geno_dtype)),
			tf.TensorSpec(shape=(None, 2), dtype=tf.string)
		)
		gen_args = (pref_chunk_size, geno_dtype, training, shuffle)
		ds = tf.data.Dataset.from_generator(self.generator,
		                                    output_signature = gen_outshapes,
		                                    args = gen_args)

		# TODO: if training, prep norm scaler
		# (if norm methods other than genotype-wise are applicable)

		ds = ds.map(self._normalize, num_parallel_calls=tf.data.AUTOTUNE)

		# TODO: stash as TFRecord?
		# TODO: prefetch, mb batching

		self.tf_dataset = ds

		# TODO: apparently interleave is good for concurrent reading & mapping?
		# can also be used to preprocess many input files? (as per docstring)

	# distribute dataset: https://stackoverflow.com/q/59185729

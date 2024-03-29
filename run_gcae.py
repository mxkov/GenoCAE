"""GenoCAE.

Usage:
  run_gcae.py train --datadir=<name> --data=<name> --model_id=<name> --train_opts_id=<name> --data_opts_id=<name> --epochs=<num> [--resume_from=<num> --trainedmodeldir=<name> --patience=<num> --save_interval=<num> --start_saving_from=<num> --pheno_model_id=<name> --phenotype_index=<num>]
  run_gcae.py project --datadir=<name>   [ --data=<name> --model_id=<name>  --train_opts_id=<name> --data_opts_id=<name> --superpops=<name> --epoch=<num> --trainedmodeldir=<name>   --pdata=<name> --trainedmodelname=<name> --pheno_model_id=<name> --phenotype_index=<num>]
  run_gcae.py plot --datadir=<name> [  --data=<name>  --model_id=<name> --train_opts_id=<name> --data_opts_id=<name>  --superpops=<name> --epoch=<num> --trainedmodeldir=<name>  --pdata=<name> --trainedmodelname=<name> --pheno_model_id=<name> --phenotype_index=<num>]
  run_gcae.py animate --datadir=<name>   [ --data=<name>   --model_id=<name> --train_opts_id=<name> --data_opts_id=<name>  --superpops=<name> --epoch=<num> --trainedmodeldir=<name> --pdata=<name> --trainedmodelname=<name> --pheno_model_id=<name> --phenotype_index=<num>]
  run_gcae.py evaluate --datadir=<name> --metrics=<name>  [  --data=<name>  --model_id=<name> --train_opts_id=<name> --data_opts_id=<name>  --superpops=<name> --epoch=<num> --trainedmodeldir=<name>  --pdata=<name> --trainedmodelname=<name> --pheno_model_id=<name> --phenotype_index=<num>]

Options:
  -h --help                  show this screen
  --datadir=<name>           directory where sample data is stored. if not absolute: assumed relative to GenoCAE/ directory. DEFAULT: data/
  --data=<name>              file prefix, not including path, of the data files (EIGENSTRAT of PLINK format)
  --trainedmodeldir=<name>   base path where to save model training directories. if not absolute: assumed relative to GenoCAE/ directory. DEFAULT: ae_out/
  --model_id=<name>          model id, corresponding to a file models/{model_id}.json
  --train_opts_id=<name>     train options id, corresponding to a file train_opts/{train_opts_id}.json
  --data_opts_id=<name>      data options id, corresponding to a file data_opts/{data_opts_id}.json
  --epochs<num>              number of epochs to train
  --resume_from<num>         saved epoch to resume training from. set to -1 for latest saved epoch. DEFAULT: None (don't resume)
  --save_interval<num>       epoch intervals at which to save state of model. DEFAULT: None (don't save)
  --start_saving_from<num>   number of epochs to train before starting to save model state. DEFAULT: 0.
  --trainedmodelname=<name>  name of the model training directory to fetch saved model state from when project/plot/evaluating
  --pdata=<name>             file prefix, not including path, of data to project/plot/evaluate. if not specified, assumed to be the same the model was trained on.
  --epoch<num>               epoch at which to project/plot/evaluate data. DEFAULT: all saved epochs
  --superpops<name>          path+filename of file mapping populations to superpopulations. used to color populations of the same superpopulation in similar colors in plotting. if not absolute path: assumed relative to GenoCAE/ directory.
  --metrics=<name>           the metric(s) to evaluate, e.g. hull_error of f1 score. can pass a list with multiple metrics, e.g. "f1_score_3,f1_score_5". DEFAULT: f1_score_3
  --patience=<num>           stop training after this number of epochs without improving lowest validation. DEFAULT: None
  --pheno_model_id=<name>    phenotype model id, e.g. "p1", corresponding to a file models/{pheno_model_id}.json
  --phenotype_index=<num>    number of the phenotype of interest in file {data}.phe, 0 for the first phenotype after the ID columns. DEFAULT: 0.

"""

from docopt import docopt, DocoptExit
import tensorflow as tf
from tensorflow.keras import Model, layers
from datetime import datetime
from utils.data_handler import get_saved_epochs, get_projected_epochs, write_h5, read_h5, get_coords_by_pop, data_generator_ae, data_generator_pheno, convex_hull_error, f1_score_kNN, plot_genotype_hist, to_genotypes_sigmoid_round, to_genotypes_invscale_round, GenotypeConcordance, get_pops_with_k, get_ind_pop_list_from_map, get_baseline_gc, write_metric_per_epoch_to_csv
from utils.visualization import plot_coords_by_superpop, plot_clusters_by_superpop, plot_coords, plot_coords_by_pop, make_animation, write_f1_scores_to_csv
import utils.visualization
import utils.layers
import json
import numpy as np
import time
import os
import glob
import math
import matplotlib.pyplot as plt
import csv
import copy
import h5py
import matplotlib.animation as animation
from pathlib import Path

GCAE_DIR = Path(__file__).resolve().parent
class Autoencoder(Model):

	def __init__(self, model_architecture, n_markers, noise_std, regularizer):
		'''

		Initiate the autoencoder with the specified options.
		All variables of the model are defined here.

		:param model_architecture: dict containing a list of layer representations
		:param n_markers: number of markers / SNPs in the data
		:param noise_std: standard deviation of noise to add to encoding layer during training. False if no noise.
		:param regularizer: dict containing regularizer info. False if no regularizer.
		'''
		super(Autoencoder, self).__init__()
		self.all_layers = []
		self.n_markers = n_markers
		self.noise_std = noise_std
		self.residuals = dict()
		self.marker_spec_var = False

		print("\n______________________________ Building model ______________________________")
		# variable that keeps track of the size of layers in encoder, to be used when constructing decoder.
		ns=[]
		ns.append(n_markers)

		first_layer_def = model_architecture["layers"][0]
		layer_module = getattr(eval(first_layer_def["module"]), first_layer_def["class"])
		layer_args = first_layer_def["args"]
		try:
			activation = getattr(tf.nn, layer_args["activation"])
			layer_args.pop("activation")
			first_layer = layer_module(activation=activation, **layer_args)

		except KeyError:
			first_layer = layer_module(**layer_args)
			activation = None

		self.all_layers.append(first_layer)
		print("Adding layer: " + str(layer_module.__name__) + ": " + str(layer_args))

		if first_layer_def["class"] == "conv1d" and "strides" in layer_args.keys() and layer_args["strides"] > 1:
			ns.append(int(first_layer.shape[1]))

		# add all layers except first
		for layer_def in model_architecture["layers"][1:]:
			layer_module = getattr(eval(layer_def["module"]), layer_def["class"])
			layer_args = layer_def["args"]

			for arg in ["size", "layers", "units", "shape", "target_shape", "output_shape", "kernel_size", "strides"]:

				if arg in layer_args.keys():
					layer_args[arg] = eval(str(layer_args[arg]))

			if layer_def["class"] == "MaxPool1D":
				ns.append(int(math.ceil(float(ns[-1]) / layer_args["strides"])))

			if layer_def["class"] == "Conv1D" and "strides" in layer_def.keys() and layer_def["strides"] > 1:
				raise NotImplementedError

			print("Adding layer: " + str(layer_module.__name__) + ": " + str(layer_args))

			if "name" in layer_args and (layer_args["name"] == "i_msvar" or layer_args["name"] == "nms"):
				self.marker_spec_var = True

			if "activation" in layer_args.keys():
				activation = getattr(tf.nn, layer_args["activation"])
				layer_args.pop("activation")
				this_layer = layer_module(activation=activation, **layer_args)
			else:
				this_layer = layer_module(**layer_args)

			self.all_layers.append(this_layer)

		if noise_std:
			self.noise_layer = tf.keras.layers.GaussianNoise(noise_std)

		self.ns = ns
		self.regularizer = regularizer

		if self.marker_spec_var:
			random_uniform = tf.random_uniform_initializer()
			self.ms_variable = tf.Variable(random_uniform(shape = (1, n_markers), dtype=tf.float32), name="marker_spec_var")
			self.nms_variable = tf.Variable(random_uniform(shape = (1, n_markers), dtype=tf.float32), name="nmarker_spec_var")
		else:
			print("No marker specific variable.")


	def call(self, input_data, is_training = True, verbose = False):
		'''
		The forward pass of the model. Given inputs, calculate the output of the model.

		:param input_data: input data
		:param is_training: if called during training
		:param verbose: print the layers and their shapes
		:return: output of the model (last layer) and latent representation (encoding layer)

		'''

		# if we're adding a marker specific variables as an additional channel
		if self.marker_spec_var:
			# Tiling it to handle the batch dimension

			ms_tiled = tf.tile(self.ms_variable, (tf.shape(input_data)[0], 1))
			ms_tiled = tf.expand_dims(ms_tiled, 2)
			nms_tiled = tf.tile(self.nms_variable, (tf.shape(input_data)[0], 1))
			nms_tiled = tf.expand_dims(nms_tiled, 2)
			concatted_input = tf.concat([input_data, ms_tiled], 2)
			input_data = concatted_input

		if verbose:
			print("inputs shape " + str(input_data.shape))

		first_layer = self.all_layers[0]
		counter = 1

		if verbose:
			print("layer {0}".format(counter))
			print("--- type: {0}".format(type(first_layer)))

		x = first_layer(inputs=input_data)

		if "Residual" in first_layer.name:
			out = self.handle_residual_layer(first_layer.name, x, verbose=verbose)
			if not out == None:
				x = out
		if verbose:
			print("--- shape: {0}".format(x.shape))

		# indicator if were doing genetic clustering (ADMIXTURE-style) or not
		have_encoded_raw = False

		# initializing encoded data
		encoded_data = None

		# do all layers except first
		for layer_def in self.all_layers[1:]:
			try:
				layer_name = layer_def.cname
			except:
				layer_name = layer_def.name

			counter += 1

			if verbose:
				print("layer {0}: {1} ({2}) ".format(counter, layer_name, type(layer_def)))

			if layer_name == "dropout":
				x = layer_def(x, training = is_training)
			else:
				# some future changes might require specifying the `training` arg here as well
				x = layer_def(x)

			# If this is a clustering model then we add noise to the layer first in this step
			# and the next layer, which is sigmoid, is the actual encoding.
			if layer_name == "encoded_raw":
				have_encoded_raw = True
				if self.noise_std:
					x = self.noise_layer(x, training = is_training)
				encoded_data_raw = x

			# If this is the encoding layer, we add noise if we are training
			if layer_name == "encoded":
				if self.noise_std and not have_encoded_raw:
					x = self.noise_layer(x, training = is_training)
				encoded_data = x

			if "Residual" in layer_name:
				out = self.handle_residual_layer(layer_name, x, verbose=verbose)
				if not out == None:
					x = out

			# inject marker-specific variable by concatenation
			if "i_msvar" in layer_name and self.marker_spec_var:
				x = self.injectms(verbose, x, layer_name, ms_tiled, self.ms_variable)

			if "nms" in layer_name and self.marker_spec_var:
				x = self.injectms(verbose, x, layer_name, nms_tiled, self.nms_variable)

			if verbose:
				print("--- shape: {0}".format(x.shape))

		if self.regularizer and encoded_data is not None:
			reg_module = eval(self.regularizer["module"])
			reg_name = getattr(reg_module, self.regularizer["class"])
			reg_func = reg_name(float(self.regularizer["reg_factor"]))

			# if this is a clustering model then the regularization is added to the raw encoding, not the softmaxed one
			if have_encoded_raw:
				reg_loss = reg_func(encoded_data_raw)
			else:
				reg_loss = reg_func(encoded_data)
			self.add_loss(reg_loss)

		return x, encoded_data


	def handle_residual_layer(self, layer_name, input, verbose=False):
		suffix = layer_name.split("Residual_")[-1].split("_")[0]
		res_number = suffix[0:-1]
		if suffix.endswith("a"):
			if verbose:
				print("encoder-to-decoder residual: saving residual {}".format(res_number))
			self.residuals[res_number] = input
			return None
		if suffix.endswith("b"):
			if verbose:
				print("encoder-to-decoder residual: adding residual {}".format(res_number))
			residual_tensor = self.residuals[res_number]
			res_length = residual_tensor.shape[1]
			if len(residual_tensor.shape) == 3:
				x = tf.keras.layers.Add()([input[:,0:res_length,:], residual_tensor])
			if len(residual_tensor.shape) == 2:
				x = tf.keras.layers.Add()([input[:,0:res_length], residual_tensor])

			return x

	def injectms(self, verbose, x, layer_name, ms_tiled, ms_variable):
		if verbose:
			print("----- injecting marker-specific variable")

		# if we need to reshape ms_variable before concatting it
		if not self.n_markers == x.shape[1]:
			d = int(math.ceil(float(self.n_markers) / int(x.shape[1])))
			diff = d*int(x.shape[1]) - self.n_markers
			ms_var = tf.reshape(tf.pad(ms_variable,[[0,0],[0,diff]]), (-1, x.shape[1],d))
			# Tiling it to handle the batch dimension
			ms_tiled = tf.tile(ms_var, (tf.shape(x)[0],1,1))

		else:
			# Tiling it to handle the batch dimension
			ms_tiled = tf.tile(ms_variable, (x.shape[0],1))
			ms_tiled = tf.expand_dims(ms_tiled, 2)

		if "_sg" in layer_name:
			if verbose:
				print("----- stopping gradient for marker-specific variable")
			ms_tiled = tf.stop_gradient(ms_tiled)

		if verbose:
			print("ms var {}".format(ms_variable.shape))
			print("ms tiled {}".format(ms_tiled.shape))
			print("concatting: {0} {1}".format(x.shape, ms_tiled.shape))

		x = tf.concat([x, ms_tiled], 2)

		return x

@tf.function
def run_optimization(model, optimizer, loss_function, input, targets,
                     phenomodel=None, phenotargets=None, pheno_loss_function=None):
	'''
	Run one step of optimization process based on the given data.

	:param model: a tf.keras.Model
	:param optimizer: a tf.keras.optimizers
	:param loss_function: a loss function
	:param input: input genotype data
	:param targets: target genotype data
	:param phenomodel: a tf.keras.Model
	:param phenotargets: target phenotype data
	:return: loss function values for model and phenomodel

	'''
	def combine_gradients(grad1, grad2):
		alpha_nom = tf.constant(0.)
		alpha_denom = tf.constant(1.0e-30)
		for g1, g2 in zip(grad1, grad2):
			if g1 is not None and g2 is not None:
				gdiff = g2 - g1
				alpha_nom += tf.math.reduce_sum(gdiff * g2)
				alpha_denom += tf.math.reduce_sum(gdiff * gdiff)
		alpha = alpha_nom / alpha_denom
		alpha_capped = tf.clip_by_value(alpha, 0., 1.)
		grad = []
		for g1, g2 in zip(grad1, grad2):
			if g1 is None:
				grad.append(g2)
			elif g2 is None:
				grad.append(g1)
			else:
				grad.append(g1*(alpha_capped) + g2*(1.0-alpha_capped))
		return grad, alpha
	
	if pheno_loss_function is None:
		print("Warning: pheno loss function not specified, using default")
		def pheno_loss_function(y_pred, y_true):
			return tf.math.reduce_sum(tf.square(y_pred - y_true)) * 1e-2

	allvars = model.trainable_variables + (phenomodel.trainable_variables if phenomodel is not None else [])

	with tf.GradientTape() as g:
		output, _ = model(input, is_training=True)
		loss_value = loss_function(y_pred = output, y_true = targets)
		loss_value += sum(model.losses)
	gradients = g.gradient(loss_value, allvars)

	if phenomodel is not None:
		with tf.GradientTape() as g2:
			_, encoded_data = model(input, is_training=True)
			phenoutput, _ = phenomodel(encoded_data, is_training=True)
			pheno_loss_value = pheno_loss_function(phenoutput, phenotargets)
		phenogradients = g2.gradient(pheno_loss_value, allvars)
		gradients, _ = combine_gradients(gradients, phenogradients)
	else:
		pheno_loss_value = None

	optimizer.apply_gradients(zip(gradients, allvars))
	return loss_value, pheno_loss_value


def get_batches(n_samples, batch_size):
	n_batches = n_samples // batch_size

	n_samples_last_batch = n_samples % batch_size
	if n_samples_last_batch > 0:
		n_batches += 1
	else:
		n_samples_last_batch = batch_size

	return n_batches, n_samples_last_batch

def alfreqvector(y_pred):
	'''
	Get a probability distribution over genotypes from y_pred.
	Assumes y_pred is raw model output, one scalar value per genotype.

	Scales this to (0,1) and interprets this as a allele frequency, uses formula
	for Hardy-Weinberg equilibrium to get probabilities for genotypes [0,1,2].

	:param y_pred: (n_samples x n_markers) tensor of raw network output for each sample and site
	:return: (n_samples x n_markers x 3 tensor) of genotype probabilities for each sample and site
	'''
	if len(y_pred.shape) == 2:
		alfreq = tf.keras.activations.sigmoid(y_pred)
		alfreq = tf.expand_dims(alfreq, -1)
		return tf.concat(((1-alfreq) ** 2, 2 * alfreq * (1 - alfreq), alfreq ** 2), axis=-1)
	else:
		return tf.nn.softmax(y_pred)

def save_model_weights(epoch, train_dir, weights_dir, model, prefix=""):
	if model is None:		# happens to phenomodel sometimes
		return
	weights_file_prefix = os.path.join(train_dir, weights_dir, "{}{}".format(prefix, epoch))
	startTime = datetime.now()
	model.save_weights(weights_file_prefix, save_format ="tf")
	save_time = (datetime.now() - startTime).total_seconds()
	print("-------- Saving weights: {0} time: {1}".format(weights_file_prefix, save_time))




if __name__ == "__main__":
	print("tensorflow version {0}".format(tf.__version__))
	tf.keras.backend.set_floatx('float32')

	try:
		arguments = docopt(__doc__, version='GenoCAE 1.1.0')
	except DocoptExit:
		print("Invalid command. Run 'python run_gcae.py --help' for more information.")
		exit(1)

	for k in list(arguments.keys()):
		knew = k.split('--')[-1]
		arg=arguments.pop(k)
		arguments[knew]=arg

	if arguments["trainedmodeldir"]:
		trainedmodeldir = arguments["trainedmodeldir"]
		if not os.path.isabs(trainedmodeldir):
			trainedmodeldir = os.path.join(GCAE_DIR, trainedmodeldir)

	else:
		trainedmodeldir = os.path.join(GCAE_DIR, "ae_out")

	if arguments["datadir"]:
		datadir = arguments["datadir"]
		if not os.path.isabs(datadir):
			datadir = os.path.join(GCAE_DIR, datadir)

	else:
		datadir = os.path.join(GCAE_DIR, "data")

	if arguments["trainedmodelname"]:
		trainedmodelname = arguments["trainedmodelname"]
		train_directory = os.path.join(trainedmodeldir, trainedmodelname)

		namesplit = trainedmodelname.split(".")
		data_opts_id = namesplit[4]
		train_opts_id = namesplit[3]
		model_id = namesplit[1]
		data = namesplit[5]
		pheno_model_id = namesplit[2]
		if pheno_model_id == '_':
			pheno_model_id = None

	else:
		data = arguments["data"]
		data_opts_id = arguments["data_opts_id"]
		train_opts_id = arguments["train_opts_id"]
		model_id = arguments["model_id"]
		pheno_model_id = arguments["pheno_model_id"]
		train_directory = False

	with open(os.path.join(GCAE_DIR, "data_opts", data_opts_id+".json")) as data_opts_def_file:
		data_opts = json.load(data_opts_def_file)

	with open(os.path.join(GCAE_DIR, "train_opts", train_opts_id+".json")) as train_opts_def_file:
		train_opts = json.load(train_opts_def_file)

	with open(os.path.join(GCAE_DIR, "models", model_id+".json")) as model_def_file:
		model_architecture = json.load(model_def_file)

	if pheno_model_id is not None:
		with open(os.path.join(GCAE_DIR, "models", pheno_model_id+".json")) as pheno_model_def_file:
			pheno_model_architecture = json.load(pheno_model_def_file)
	else:
		pheno_model_architecture = None

	for layer_def in model_architecture["layers"]:
		if "args" in layer_def.keys() and "name" in layer_def["args"].keys() and "encoded" in layer_def["args"]["name"] and "units" in layer_def["args"].keys():
			n_latent_dim = layer_def["args"]["units"]

	# indicator of whether this is a genetic clustering or dimensionality reduction model
	doing_clustering = False
	for layer_def in model_architecture["layers"][1:-1]:
		if "encoding_raw" in layer_def.keys():
			doing_clustering = True

	print("\n______________________________ arguments ______________________________")
	for k in arguments.keys():
		print(k + " : " + str(arguments[k]))
	print("\n______________________________ data opts ______________________________")
	for k in data_opts.keys():
		print(k + " : " + str(data_opts[k]))
	print("\n______________________________ train opts ______________________________")
	for k in train_opts:
		print(k + " : " + str(train_opts[k]))
	print("______________________________")


	batch_size = train_opts["batch_size"]
	learning_rate = train_opts["learning_rate"]
	regularizer = train_opts["regularizer"]

	superpopulations_file = arguments['superpops']
	if superpopulations_file and not os.path.isabs(os.path.dirname(superpopulations_file)):
		superpopulations_file = os.path.join(GCAE_DIR,
		                                     os.path.dirname(superpopulations_file),
		                                     Path(superpopulations_file).name)

	norm_opts = data_opts["norm_opts"]
	norm_mode = data_opts["norm_mode"]
	validation_split = data_opts["validation_split"]

	if "sparsifies" in data_opts.keys():
		sparsify_input = True
		missing_mask_input = True
		n_input_channels = 2
		sparsifies = data_opts["sparsifies"]

	else:
		sparsify_input = False
		missing_mask_input = False
		n_input_channels = 1

	if "impute_missing" in data_opts.keys():
		fill_missing = data_opts["impute_missing"]

	else:
		fill_missing = False

	if fill_missing:
		print("Imputing originally missing genotypes to most common value.")
	else:
		print("Keeping originally missing genotypes.")
		missing_mask_input = True
		n_input_channels = 2

	if not train_directory:
		phm_id = pheno_model_id if pheno_model_id is not None else "_"
		train_directory = os.path.join(trainedmodeldir,
		                               ".".join(("ae", model_id, phm_id,
		                                         train_opts_id, data_opts_id, data)))

	if arguments["pdata"]:
		pdata = arguments["pdata"]
	else:
		pdata = data

	data_prefix = os.path.join(datadir, pdata)
	results_directory = os.path.join(train_directory, pdata)
	try:
		os.mkdir(results_directory)
	except OSError:
		pass

	encoded_data_file = os.path.join(train_directory, pdata, "encoded_data.h5")

	if "noise_std" in train_opts.keys():
		noise_std = train_opts["noise_std"]
	else:
		noise_std = False

	if (arguments['evaluate'] or arguments['animate'] or arguments['plot']):

		if os.path.isfile(encoded_data_file):
			encoded_data = h5py.File(encoded_data_file, 'r')
		else:
			print("------------------------------------------------------------------------")
			print("Error: File {0} not found.".format(encoded_data_file))
			print("------------------------------------------------------------------------")
			exit(1)

		epochs = get_projected_epochs(encoded_data_file)

		if arguments['epoch']:
			epoch = int(arguments['epoch'])
			if epoch in epochs:
				epochs = [epoch]
			else:
				print("------------------------------------------------------------------------")
				print("Error: Epoch {0} not found in {1}.".format(epoch, encoded_data_file))
				print("------------------------------------------------------------------------")
				exit(1)

		if doing_clustering:
			if arguments['animate']:
				print("------------------------------------------------------------------------")
				print("Error: Animate not supported for genetic clustering model.")
				print("------------------------------------------------------------------------")
				exit(1)


			if arguments['plot'] and not superpopulations_file:
				print("------------------------------------------------------------------------")
				print("Error: Plotting of genetic clustering results requires a superpopulations file.")
				print("------------------------------------------------------------------------")
				exit(1)

	else:
		dg = data_generator_ae(data_prefix,
		                       normalization_mode = norm_mode,
		                       normalization_options = norm_opts,
		                       impute_missing = fill_missing)
		n_markers = copy.deepcopy(dg.n_markers)

		try:
			phenotype_index = int(arguments["phenotype_index"])
		except:
			phenotype_index = 0
		dg_ph = data_generator_pheno(data_prefix + ".phe",
		                             pt_index = phenotype_index,
		                             phenomodel_defined = (pheno_model_architecture is not None))


		loss_def = train_opts["loss"]
		loss_class = getattr(eval(loss_def["module"]), loss_def["class"])
		if "args" in loss_def.keys():
			loss_args = loss_def["args"]
		else:
			loss_args = dict()
		loss_obj = loss_class(**loss_args)

		def get_originally_nonmissing_mask(genos):
			'''
			Get a boolean mask representing missing values in the data.
			Missing value is represented by float(norm_opts["missing_val"]).

			Uses the presence of missing_val in the true genotypes as indicator, missing_val should not be set to
			something that can exist in the data set after normalization!!!!

			:param genos: (n_samples x n_markers) genotypes
			:return: boolean mask of the same shape as genos
			'''
			orig_nonmissing_mask = tf.not_equal(genos, float(norm_opts["missing_val"]))

			return orig_nonmissing_mask

		if loss_class == tf.keras.losses.CategoricalCrossentropy or loss_class == tf.keras.losses.KLDivergence:

			def loss_func(y_pred, y_true):
				y_pred = y_pred[:, 0:n_markers]

				if not fill_missing:
					orig_nonmissing_mask = get_originally_nonmissing_mask(y_true)

				y_pred = alfreqvector(y_pred)
				y_true = tf.one_hot(tf.cast(y_true * 2, tf.uint8), 3)*0.9997 + 0.0001

				if not fill_missing:
					y_pred = y_pred[orig_nonmissing_mask]
					y_true = y_true[orig_nonmissing_mask]

				return loss_obj(y_pred = y_pred, y_true = y_true)


		else:
			def loss_func(y_pred, y_true):

				y_pred = y_pred[:, 0:n_markers]
				y_true = tf.convert_to_tensor(y_true)

				if not fill_missing:
					orig_nonmissing_mask = get_originally_nonmissing_mask(y_true)
					y_pred = y_pred[orig_nonmissing_mask]
					y_true = y_true[orig_nonmissing_mask]

				return loss_obj(y_pred = y_pred, y_true = y_true)

		# TODO: perhaps enable customizable pheno loss function
		def pheno_loss_func(y_pred, y_true):
			return tf.math.reduce_mean(tf.square(y_pred - y_true)) * 1e-2


	# defining some constants before responding to the 'app' arguments
	ae_weights_dir = "weights"
	pheno_weights_dir = "weights_pheno"


	if arguments['train']:

		epochs = int(arguments["epochs"])

		try:
			save_interval = int(arguments["save_interval"])
		except:
			save_interval = epochs

		try:
			start_saving_from = int(arguments["start_saving_from"])
		except:
			start_saving_from = 0

		try:
			patience = int(arguments["patience"])
		except:
			patience = epochs

		try:
			resume_from = int(arguments["resume_from"])
			if resume_from < 1:
				saved_epochs = get_saved_epochs(train_directory, weights_directory=ae_weights_dir)
				resume_from = saved_epochs[-1]
		except:
			resume_from = False

		dg.define_validation_set(validation_split = validation_split)
		input_valid, targets_valid, ind_pop_list_valid  = dg.get_valid_set(0.0)
		phenotargets_valid = dg_ph.generate(ind_pop_list_valid)

		# if we do not have missing mask input, remeove that dimension/channel from the input that data generator returns
		if not missing_mask_input:
			input_valid = input_valid[:,:,0, np.newaxis]

		n_unique_train_samples = copy.deepcopy(dg.n_train_samples)
		n_valid_samples = copy.deepcopy(dg.n_valid_samples)

		assert n_valid_samples == len(input_valid)
		assert n_valid_samples == len(targets_valid)

		if "n_samples" in train_opts.keys() and int(train_opts["n_samples"]) > 0:
			n_train_samples = int(train_opts["n_samples"])
		else:
			n_train_samples = n_unique_train_samples

		batch_size_valid = batch_size
		n_train_batches, n_train_samples_last_batch = get_batches(n_train_samples, batch_size)
		n_valid_batches, n_valid_samples_last_batch = get_batches(n_valid_samples, batch_size_valid)

		train_times = []
		train_epochs = []
		save_epochs = []

		# get one sample (two samples?..) to run through optimization
		# to reload model weights and optimizer variables
		#
		# NOTE: this piece is placed here in Richel's code,
		# but I'm moving it under `if resume_from` as per Kristiina's code.
		# Hence I'm commenting it out here. Might bring it back later.
		#
		#input_init, targets_init, ind_pop_list = dg.get_train_batch(0.0, 1)
		#phenotargets_init = dg_ph.generate(ind_pop_list)
		#dg.reset_batch_index()
		#if not missing_mask_input:
		#	input_init = input_init[:,:,0, np.newaxis]

		############### setup learning rate schedule ##############
		step_counter = resume_from * n_train_batches
		if "lr_scheme" in train_opts.keys():
			schedule_module = getattr(eval(train_opts["lr_scheme"]["module"]), train_opts["lr_scheme"]["class"])
			schedule_args = train_opts["lr_scheme"]["args"]

			if "decay_every" in schedule_args:
				decay_every = int(schedule_args.pop("decay_every"))
				decay_steps = n_train_batches * decay_every
				schedule_args["decay_steps"] = decay_steps

			lr_schedule = schedule_module(learning_rate, **schedule_args)

			# use the schedule to calculate what the lr was at the epoch were resuming from
			updated_lr = lr_schedule(step_counter)
			lr_schedule = schedule_module(updated_lr, **schedule_args)

			print("Using learning rate schedule {0}.{1} with {2}".format(train_opts["lr_scheme"]["module"],
			                                                             train_opts["lr_scheme"]["class"],
			                                                             schedule_args))
		else:
			lr_schedule = False

		print("\n______________________________ Data ______________________________")
		print("N unique train samples: {0}".format(n_unique_train_samples))
		print("--- training on : {0}".format(n_train_samples))
		print("N valid samples: {0}".format(n_valid_samples))
		print("N markers: {0}".format(n_markers))
		print("")

		autoencoder = Autoencoder(model_architecture, n_markers, noise_std, regularizer)
		if pheno_model_architecture is not None:
			pheno_model = Autoencoder(pheno_model_architecture, 2, noise_std, regularizer)
		else:
			pheno_model = None
		optimizer = tf.optimizers.Adam(learning_rate = lr_schedule)

		if resume_from:
			print("\n______________________________ Resuming training from epoch {0} ______________________________".format(resume_from))
			weights_file_prefix = os.path.join(train_directory, ae_weights_dir, str(resume_from))
			print("Reading weights from {0}".format(weights_file_prefix))
			if pheno_model is not None:
				pheno_weights_file_prefix = os.path.join(train_directory, pheno_weights_dir, str(resume_from))
				print("Reading phenomodel weights from {0}".format(pheno_weights_file_prefix))

			# get a single sample to run through optimization to reload weights and optimizer variables
			input_init, targets_init, ind_pop_list_init = dg.get_train_batch(0.0, 1)
			phenotargets_init = dg_ph.generate(ind_pop_list_init)
			dg.reset_batch_index()
			if not missing_mask_input:
				input_init = input_init[:,:,0, np.newaxis]

			# This initializes the variables used by the optimizers,
			# as well as any stateful metric variables
			run_optimization(autoencoder, optimizer, loss_func, input_init, targets_init,
			                 phenomodel=pheno_model, phenotargets=phenotargets_init,
			                 pheno_loss_function=pheno_loss_func)
			autoencoder.load_weights(weights_file_prefix)
			if pheno_model is not None:
				pheno_model.load_weights(pheno_weights_file_prefix)

		print("\n______________________________ Train ______________________________")

		# a small run-through of the model with just 2 samples for printing the dimensions of the layers (verbose=True)
		print("Model layers and dimensions:")
		print("-----------------------------")

		input_test, targets_test, _ = dg.get_train_set(0.0)
		if not missing_mask_input:
			input_test = input_test[:,:,0, np.newaxis]
		output_test, encoded_data_test = autoencoder(input_test[0:2], is_training = False, verbose = True)
		if pheno_model is not None:
			phenoutput_test, _ = pheno_model(encoded_data_test, is_training = False, verbose = True)

		######### Create objects for tensorboard summary ###############################

		train_writer = tf.summary.create_file_writer(os.path.join(train_directory, 'train'))
		valid_writer = tf.summary.create_file_writer(os.path.join(train_directory, 'valid'))
		if pheno_model is not None:
			train_pheno_writer = tf.summary.create_file_writer(os.path.join(train_directory, 'train_pheno'))
			valid_pheno_writer = tf.summary.create_file_writer(os.path.join(train_directory, 'valid_pheno'))

		######################################################

		# train losses per epoch
		losses_t = []
		# valid losses per epoch
		losses_v = []
		# min loss stats
		min_valid_loss = np.inf
		min_valid_loss_epoch = None

		if pheno_model is not None:
			# pheno train losses per epoch
			pheno_losses_t = []
			# pheno valid losses per epoch
			pheno_losses_v = []
			# min loss stats
			min_valid_pheno_loss = np.inf
			min_valid_pheno_loss_epoch = None

		for e in range(1,epochs+1):
			startTime = datetime.now()
			dg.shuffle_train_samples()
			effective_epoch = e + resume_from
			losses_t_batches = []
			losses_v_batches = []
			if pheno_model is not None:
				pheno_losses_t_batches = []
				pheno_losses_v_batches = []

			for ii in range(n_train_batches):
				step_counter += 1

				if sparsify_input:
					sparsify_fraction = sparsifies[step_counter % len(sparsifies)]
				else:
					sparsify_fraction = 0.0

				# Generating batches
				# last batch is probably not full
				if ii == n_train_batches - 1:
					batch_input, batch_target, batch_ind_pop_list = dg.get_train_batch(sparsify_fraction, n_train_samples_last_batch)
				else:
					batch_input, batch_target, batch_ind_pop_list = dg.get_train_batch(sparsify_fraction, batch_size)
				phenotargets_batch = dg_ph.generate(batch_ind_pop_list)

				# TODO temporary solution: should fix data generator so it doesnt bother with the mask if not needed
				if not missing_mask_input:
					batch_input = batch_input[:,:,0,np.newaxis]

				train_batch_loss, train_batch_pheno_loss = \
					run_optimization(autoencoder, optimizer,
					                 loss_func, batch_input, batch_target,
					                 phenomodel=pheno_model,
					                 phenotargets=phenotargets_batch,
					                 pheno_loss_function = pheno_loss_func)
				losses_t_batches.append(train_batch_loss)
				if pheno_model is not None:
					pheno_losses_t_batches.append(train_batch_pheno_loss)

			train_loss_this_epoch = np.average(losses_t_batches)
			if pheno_model is not None:
				train_pheno_loss_this_epoch = np.average(pheno_losses_t_batches)
			
			with train_writer.as_default():
				tf.summary.scalar('loss', train_loss_this_epoch, step = step_counter)
				if lr_schedule:
					tf.summary.scalar("learning_rate", optimizer._decayed_lr(var_dtype=tf.float32), step = step_counter)
				else:
					tf.summary.scalar("learning_rate", learning_rate, step = step_counter)

			if pheno_model is not None:
				with train_pheno_writer.as_default():
					tf.summary.scalar('pheno loss', train_pheno_loss_this_epoch, step = step_counter)



			train_time = (datetime.now() - startTime).total_seconds()
			train_times.append(train_time)
			train_epochs.append(effective_epoch)
			losses_t.append(train_loss_this_epoch)
			if pheno_model is not None:
				pheno_losses_t.append(train_pheno_loss_this_epoch)

			print("")
			print("Epoch: {}/{}...".format(effective_epoch, epochs+resume_from))
			print("--- Train loss: {:.4f}  time: {}".format(train_loss_this_epoch, train_time))
			if pheno_model is not None:
				print("--- Train pheno loss: {:.4f}  time: {}".format(train_pheno_loss_this_epoch, train_time))


			if n_valid_samples > 0:

				startTime = datetime.now()

				for jj in range(n_valid_batches):
					start = jj*batch_size_valid
					if jj == n_valid_batches - 1:
						input_valid_batch = input_valid[start:]
						targets_valid_batch = targets_valid[start:]
						phenotargets_valid_batch = (phenotargets_valid[start:]
						                            if phenotargets_valid is not None else None)
					else:
						input_valid_batch = input_valid[start:start+batch_size_valid]
						targets_valid_batch = targets_valid[start:start+batch_size_valid]
						phenotargets_valid_batch = (phenotargets_valid[start:start+batch_size_valid]
						                            if phenotargets_valid is not None else None)

					output_valid_batch, encoded_data_valid_batch = autoencoder(input_valid_batch, is_training = False)
					valid_loss_batch = loss_func(y_pred = output_valid_batch, y_true = targets_valid_batch)
					valid_loss_batch += sum(autoencoder.losses)
					losses_v_batches.append(valid_loss_batch)

					if pheno_model is not None:
						phenoutput_valid_batch, _ = pheno_model(encoded_data_valid_batch,
						                                        is_training = False)
						valid_pheno_loss_batch = pheno_loss_func(y_pred = phenoutput_valid_batch,
						                                         y_true = phenotargets_valid_batch)
						pheno_losses_v_batches.append(valid_pheno_loss_batch)

				valid_loss_this_epoch = np.average(losses_v_batches)
				losses_v.append(valid_loss_this_epoch)
				with valid_writer.as_default():
					tf.summary.scalar('loss', valid_loss_this_epoch, step=step_counter)

				if pheno_model is not None:
					valid_pheno_loss_this_epoch = np.average(pheno_losses_v_batches)
					pheno_losses_v.append(valid_pheno_loss_this_epoch)
					with valid_pheno_writer.as_default():
						tf.summary.scalar('pheno loss', valid_pheno_loss_this_epoch,
						                  step=step_counter)

				valid_time = (datetime.now() - startTime).total_seconds()

				if valid_loss_this_epoch <= min_valid_loss:
					min_valid_loss = valid_loss_this_epoch
					prev_min_val_loss_epoch = min_valid_loss_epoch
					min_valid_loss_epoch = effective_epoch

					# Recording phenoloss corresponding to min AE loss,
					# NOT min phenoloss. This is likely temporary.
					if pheno_model is not None:
						min_valid_pheno_loss = valid_pheno_loss_this_epoch
						prev_min_val_pheno_loss_epoch = min_valid_pheno_loss_epoch
						min_valid_pheno_loss_epoch = effective_epoch

					if e > start_saving_from:
						for dirname in (ae_weights_dir, pheno_weights_dir):
							for f in glob.glob(os.path.join(train_directory, dirname, "min_valid.{}.*".format(prev_min_val_loss_epoch))):
								os.remove(f)
						# Again, saving phenomodel weights corresponding to min AE loss,
						# NOT min phenoloss. For now.
						save_model_weights(effective_epoch, train_directory,
						                   ae_weights_dir, autoencoder,
						                   prefix = "min_valid.")
						save_model_weights(effective_epoch, train_directory,
						                   pheno_weights_dir, pheno_model,
						                   prefix = "min_valid.")

				evals_since_min_valid_loss = effective_epoch - min_valid_loss_epoch
				print("--- Valid loss: {:.4f}  time: {} min loss: {:.4f} epochs since: {}".format(
										valid_loss_this_epoch, valid_time, min_valid_loss, evals_since_min_valid_loss))

				if pheno_model is not None:
					print("--- Valid pheno loss: {:.4f}  time: {} min loss: {:.4f} epochs since: {}".format(
					      valid_pheno_loss_this_epoch, valid_time, min_valid_pheno_loss,
					      effective_epoch - min_valid_pheno_loss_epoch))

				if evals_since_min_valid_loss >= patience:
					break

			if e % save_interval == 0 and e > start_saving_from :
				save_model_weights(effective_epoch, train_directory,
				                   ae_weights_dir, autoencoder)
				save_model_weights(effective_epoch, train_directory,
				                   pheno_weights_dir, pheno_model)




		save_model_weights(effective_epoch, train_directory,
		                   ae_weights_dir, autoencoder)
		save_model_weights(effective_epoch, train_directory,
		                   pheno_weights_dir, pheno_model)

		outfilename = os.path.join(train_directory, "train_times.csv")
		write_metric_per_epoch_to_csv(outfilename, train_times, train_epochs)

		# recording and plotting train losses
		outfilename = os.path.join(train_directory, "losses_from_train_t.csv")
		epochs_t_combined, losses_t_combined = write_metric_per_epoch_to_csv(outfilename, losses_t, train_epochs)
		fig, ax = plt.subplots()
		plt.plot(epochs_t_combined, losses_t_combined, label="train", c="orange")

		# recording and plotting valid losses (same plot)
		if n_valid_samples > 0:
			outfilename = os.path.join(train_directory, "losses_from_train_v.csv")
			epochs_v_combined, losses_v_combined = write_metric_per_epoch_to_csv(outfilename, losses_v, train_epochs)
			plt.plot(epochs_v_combined, losses_v_combined, label="valid", c="blue")
			min_valid_loss_epoch = epochs_v_combined[np.argmin(losses_v_combined)]
			plt.axvline(min_valid_loss_epoch, color="black")
			plt.text(min_valid_loss_epoch + 0.1, 0.5,'min valid loss at epoch {}'.format(int(min_valid_loss_epoch)),
					 rotation=90,
					 transform=ax.get_xaxis_text1_transform(0)[0])

		# saving this plot
		plt.xlabel("Epoch")
		plt.ylabel("Loss function value")
		plt.legend()
		# I don't want PDF :<
		plt.savefig(os.path.join(train_directory, "losses_from_train.png"), dpi=300)
		plt.close()

		if pheno_model is not None:
			# recording and plotting train pheno losses
			# (separate plot, bc different loss function)
			outfilename = os.path.join(train_directory, "losses_from_train_t_pheno.csv")
			epochs_t_combined, pheno_losses_t_combined = \
				write_metric_per_epoch_to_csv(outfilename,
				                              pheno_losses_t,
				                              train_epochs)
			fig, ax = plt.subplots()
			plt.plot(epochs_t_combined, pheno_losses_t_combined,
			         label="train pheno", c="magenta")
			# marking min valid autoencoder loss
			if n_valid_samples > 0:
				plt.axvline(min_valid_loss_epoch, color="black")
				plt.text(min_valid_loss_epoch + 0.1, 0.5,
				         "min valid AE loss at epoch {}".format(int(min_valid_loss_epoch)),
				         rotation=90, transform=ax.get_xaxis_text1_transform(0)[0])
			# recording and plotting valid pheno losses
			if n_valid_samples > 0:
				outfilename = os.path.join(train_directory, "losses_from_train_v_pheno.csv")
				epochs_v_combined, pheno_losses_v_combined = \
					write_metric_per_epoch_to_csv(outfilename,
					                              pheno_losses_v,
					                              train_epochs)
				plt.plot(epochs_v_combined, pheno_losses_v_combined,
				         label="valid pheno", c="green")
			plt.xlabel("Epoch")
			plt.ylabel("Loss function value")
			plt.legend()
			plt.savefig(os.path.join(train_directory, "losses_from_train_pheno.png"), dpi=300)
			plt.close()

		print("Done training. Wrote to {0}".format(train_directory))

	if arguments['project']:

		projected_epochs = get_projected_epochs(encoded_data_file)

		if arguments['epoch']:
			epoch = int(arguments['epoch'])
			epochs = [epoch]

		else:
			epochs = get_saved_epochs(train_directory, weights_directory=ae_weights_dir)

		for projected_epoch in projected_epochs:
			try:
				epochs.remove(projected_epoch)
			except:
				continue

		print("Projecting epochs: {0}".format(epochs))
		print("Already projected: {0}".format(projected_epochs))

		batch_size_project = 50
		sparsify_fraction = 0.0

		_, _, ind_pop_list_train_reference = dg.get_train_set(sparsify_fraction)

		write_h5(encoded_data_file, "ind_pop_list_train", np.array(ind_pop_list_train_reference, dtype='S'))

		n_unique_train_samples = copy.deepcopy(dg.n_train_samples)

		autoencoder = Autoencoder(model_architecture, n_markers, noise_std, regularizer)
		if pheno_model_architecture is not None:
			pheno_model = Autoencoder(pheno_model_architecture, 2, noise_std, regularizer)
		else:
			pheno_model = None
		optimizer = tf.optimizers.Adam(learning_rate = learning_rate)

		# loss function of the train set per epoch
		losses_train = []
		# phenomodel loss function of the train set per epoch
		if pheno_model is not None:
			pheno_losses_train = []
		# genotype concordance of the train set per epoch
		genotype_concs_train = []

		genotype_concordance_metric = GenotypeConcordance()

		scatter_points_per_epoch = []
		colors_per_epoch = []
		markers_per_epoch = []
		edgecolors_per_epoch = []

		for epoch in epochs:
			print("########################### epoch {0} ###########################".format(epoch))
			weights_file_prefix = os.path.join(train_directory, ae_weights_dir, str(epoch))
			print("Reading weights from {0}".format(weights_file_prefix))
			if pheno_model is not None:
				pheno_weights_file_prefix = os.path.join(train_directory, pheno_weights_dir, str(epoch))
				print("Reading phenomodel weights from {0}".format(pheno_weights_file_prefix))

			input, targets, _= dg.get_train_batch(sparsify_fraction, 1)
			if not missing_mask_input:
				input = input[:,:,0, np.newaxis]

			# This initializes the variables used by the optimizers,
			# as well as any stateful metric variables
			# run_optimization(autoencoder, optimizer, loss_func, input, targets)
			autoencoder.load_weights(weights_file_prefix)
			if pheno_model is not None:
				pheno_model.load_weights(pheno_weights_file_prefix)

			if batch_size_project:
				dg.reset_batch_index()

				n_train_batches = (n_unique_train_samples // batch_size_project) + 1
				n_train_samples_last_batch = n_unique_train_samples % batch_size_project


				ind_pop_list_train = np.empty((0,2))
				encoded_train = np.empty((0, n_latent_dim))
				decoded_train = None
				targets_train = np.empty((0, n_markers))
				phenoutput_train = None
				phenotargets_train = None

				loss_value_per_train_batch = []
				if pheno_model is not None:
					pheno_loss_value_per_train_batch = []
				genotype_conc_per_train_batch = []

				for b in range(n_train_batches):

					if b == n_train_batches - 1:
						input_train_batch, targets_train_batch, ind_pop_list_train_batch = dg.get_train_batch(sparsify_fraction, n_train_samples_last_batch)
					else:
						input_train_batch, targets_train_batch, ind_pop_list_train_batch = dg.get_train_batch(sparsify_fraction, batch_size_project)
					phenotargets_train_batch = dg_ph.generate(ind_pop_list_train_batch)

					if not missing_mask_input:
						input_train_batch = input_train_batch[:,:,0, np.newaxis]

					decoded_train_batch, encoded_train_batch = autoencoder(input_train_batch, is_training = False)
					loss_train_batch = loss_func(y_pred = decoded_train_batch, y_true = targets_train_batch)
					loss_train_batch += sum(autoencoder.losses)
					loss_value_per_train_batch.append(loss_train_batch)

					ind_pop_list_train = np.concatenate((ind_pop_list_train, ind_pop_list_train_batch), axis=0)
					encoded_train = np.concatenate((encoded_train, encoded_train_batch), axis=0)
					if decoded_train is None:
						decoded_train = np.copy(decoded_train_batch[:,0:n_markers])
					else:
						decoded_train = np.concatenate((decoded_train, decoded_train_batch[:,0:n_markers]), axis=0)
					targets_train = np.concatenate((targets_train, targets_train_batch[:,0:n_markers]), axis=0)

					if pheno_model is not None:
						phenoutput_train_batch, _ = pheno_model(encoded_train_batch,
						                                        is_training = False)
						pheno_loss_train_batch = pheno_loss_func(y_pred = phenoutput_train_batch,
						                                         y_true = phenotargets_train_batch)
						pheno_loss_value_per_train_batch.append(pheno_loss_train_batch)

						if phenoutput_train is None:
							phenoutput_train = np.copy(phenoutput_train_batch)
						else:
							phenoutput_train = np.concatenate((phenoutput_train,
							                                   phenoutput_train_batch), axis=0)
						if phenotargets_train is None:
							phenotargets_train = np.copy(phenotargets_train_batch)
						else:
							phenotargets_train = np.concatenate((phenotargets_train,
							                                     phenotargets_train_batch), axis=0)


				ind_pop_list_train = np.array(ind_pop_list_train)
				encoded_train = np.array(encoded_train)
				loss_value = np.average(loss_value_per_train_batch)
				if pheno_model is not None:
					pheno_loss_value = np.average(pheno_loss_value_per_train_batch)

				if epoch == epochs[0]:
					assert len(ind_pop_list_train) == dg.n_train_samples, "{0} vs {1}".format(len(ind_pop_list_train), dg.n_train_samples)
					assert len(encoded_train) == dg.n_train_samples, "{0} vs {1}".format(len(encoded_train), dg.n_train_samples)
					assert list(ind_pop_list_train[:,0]) == list(ind_pop_list_train_reference[:,0])
					assert list(ind_pop_list_train[:,1]) == list(ind_pop_list_train_reference[:,1])
			else:
				input_train, targets_train, ind_pop_list_train = dg.get_train_set(sparsify_fraction)
				phenotargets_train = dg_ph.generate(ind_pop_list_train)

				if not missing_mask_input:
					input_train = input_train[:,:,0, np.newaxis]

				decoded_train, encoded_train = autoencoder(input_train, is_training = False)
				loss_value = loss_func(y_pred = decoded_train, y_true = targets_train)
				loss_value += sum(autoencoder.losses)

				if pheno_model is not None:
					phenoutput_train, _ = pheno_model(encoded_train, is_training = False)
					pheno_loss_value = pheno_loss_func(y_pred = phenoutput_train,
					                                   y_true = phenotargets_train)


			genotype_concordance_metric.reset_states()

			if not fill_missing:
				orig_nonmissing_mask = get_originally_nonmissing_mask(targets_train)
			else:
				orig_nonmissing_mask = np.full(targets_train.shape, True)

			if train_opts["loss"]["class"] == "MeanSquaredError" and (data_opts["norm_mode"] == "smartPCAstyle" or data_opts["norm_mode"] == "standard"):
				try:
					scaler = dg.scaler
				except:
					print("Could not calculate predicted genotypes and genotype concordance. No scaler available in data handler.")
					genotypes_output = np.array([])
					true_genotypes = np.array([])

				genotypes_output = to_genotypes_invscale_round(decoded_train[:, 0:n_markers], scaler_vals = scaler)
				true_genotypes = to_genotypes_invscale_round(targets_train, scaler_vals = scaler)
				genotype_concordance_metric.update_state(y_pred = genotypes_output[orig_nonmissing_mask],
														 y_true = true_genotypes[orig_nonmissing_mask])


			elif train_opts["loss"]["class"] == "BinaryCrossentropy" and data_opts["norm_mode"] == "genotypewise01":
				genotypes_output = to_genotypes_sigmoid_round(decoded_train[:, 0:n_markers])
				true_genotypes = targets_train
				genotype_concordance_metric.update_state(y_pred = genotypes_output[orig_nonmissing_mask], y_true = true_genotypes[orig_nonmissing_mask])

			elif train_opts["loss"]["class"] in ["CategoricalCrossentropy", "KLDivergence"] and data_opts["norm_mode"] == "genotypewise01":
				genotypes_output = tf.cast(tf.argmax(alfreqvector(decoded_train[:, 0:n_markers]), axis = -1), tf.float16) * 0.5
				true_genotypes = targets_train
				genotype_concordance_metric.update_state(y_pred = genotypes_output[orig_nonmissing_mask], y_true = true_genotypes[orig_nonmissing_mask])

			else:
				print("Could not calculate predicted genotypes and genotype concordance. Not implemented for loss {0} and normalization {1}.".format(train_opts["loss"]["class"], data_opts["norm_mode"]))
				genotypes_output = np.array([])
				true_genotypes = np.array([])

			genotype_concordance_value = genotype_concordance_metric.result()

			losses_train.append(loss_value)
			if pheno_model is not None:
				pheno_losses_train.append(pheno_loss_value)
			genotype_concs_train.append(genotype_concordance_value)

			if superpopulations_file:
				coords_by_pop = get_coords_by_pop(data_prefix, encoded_train, ind_pop_list = ind_pop_list_train)

				if doing_clustering:
					plot_clusters_by_superpop(coords_by_pop,
					                          os.path.join(results_directory,
					                                       "clusters_e_{}".format(epoch)),
					                          superpopulations_file,
					                          write_legend = epoch == epochs[0])
				else:
					scatter_points, colors, markers, edgecolors = \
						plot_coords_by_superpop(coords_by_pop,
						                        os.path.join(results_directory,
						                                     "dimred_e_{}_by_superpop".format(epoch)),
						                        superpopulations_file,
						                        plot_legend = epoch == epochs[0])

					scatter_points_per_epoch.append(scatter_points)
					colors_per_epoch.append(colors)
					markers_per_epoch.append(markers)
					edgecolors_per_epoch.append(edgecolors)

			else:
				try:
					coords_by_pop = get_coords_by_pop(data_prefix, encoded_train, ind_pop_list = ind_pop_list_train)
					plot_coords_by_pop(coords_by_pop,
					                   os.path.join(results_directory,
					                                "dimred_e_{}_by_pop".format(epoch)))
				except:
					plot_coords(encoded_train,
					            os.path.join(results_directory,
					                         "dimred_e_{}".format(epoch)))


			write_h5(encoded_data_file, "{}_encoded_train".format(epoch), encoded_train)
			if phenoutput_train is not None:
				dg_ph.write(os.path.join(results_directory, "pheno_e_{}.phe".format(epoch)),
				            ind_pop_list_train, phenoutput_train, include_stored = True)

		try:
			plot_genotype_hist(np.array(genotypes_output),
			                   os.path.join(results_directory,
			                                "output_as_genotypes_e{}".format(epoch)))
			plot_genotype_hist(np.array(true_genotypes),
			                   os.path.join(results_directory, "true_genotypes"))
		except:
			pass

		############################### losses ##############################

		## Autoencoder losses

		outfilename = os.path.join(results_directory, "losses_from_project.csv")
		epochs_combined, losses_train_combined = write_metric_per_epoch_to_csv(outfilename, losses_train, epochs)

		plt.plot(epochs_combined, losses_train_combined,
		         label="all data", c="red")
		plt.xlabel("Epoch")
		plt.ylabel("Loss function value")
		plt.legend()
		plt.savefig(os.path.join(results_directory, "losses_from_project.png"), dpi=300)
		plt.close()

		## Phenomodel losses

		if pheno_model is not None:

			outfilename = os.path.join(results_directory, "losses_from_project_pheno.csv")
			epochs_combined, pheno_losses_train_combined = \
				write_metric_per_epoch_to_csv(outfilename, pheno_losses_train, epochs)

			plt.plot(epochs_combined, pheno_losses_train_combined,
			         label="all pheno data", c="blue")
			plt.xlabel("Epoch")
			plt.ylabel("Loss function value")
			plt.legend()
			plt.savefig(os.path.join(results_directory, "losses_from_project_pheno.png"), dpi=300)
			plt.close()


		############################### gconc ###############################
		try:
			baseline_genotype_concordance = get_baseline_gc(true_genotypes)
		except:
			baseline_genotype_concordance = None

		outfilename = os.path.join(results_directory, "genotype_concordances.csv")
		epochs_combined, genotype_concs_combined = write_metric_per_epoch_to_csv(outfilename, genotype_concs_train, epochs)

		plt.plot(epochs_combined, genotype_concs_combined, label="train", c="orange")
		if baseline_genotype_concordance:
			plt.plot([epochs_combined[0], epochs_combined[-1]], [baseline_genotype_concordance, baseline_genotype_concordance], label="baseline", c="black")

		plt.xlabel("Epoch")
		plt.ylabel("Genotype concordance")
		plt.savefig(os.path.join(results_directory, "genotype_concordances.png"), dpi=300)
		plt.close()



	if arguments['animate']:

		print("Animating epochs {}".format(epochs))

		FFMpegWriter = animation.writers['ffmpeg']
		scatter_points_per_epoch = []
		colors_per_epoch = []
		markers_per_epoch = []
		edgecolors_per_epoch = []

		ind_pop_list_train = read_h5(encoded_data_file, "ind_pop_list_train")

		for epoch in epochs:
			print("########################### epoch {0} ###########################".format(epoch))

			encoded_train = read_h5(encoded_data_file, "{0}_encoded_train".format(epoch))

			coords_by_pop = get_coords_by_pop(data_prefix, encoded_train, ind_pop_list = ind_pop_list_train)
			name = ""

			if superpopulations_file:
				scatter_points, colors, markers, edgecolors = \
					plot_coords_by_superpop(coords_by_pop, name, superpopulations_file, plot_legend=False, savefig=False)
				suffix = "_by_superpop"
			else:
				try:
					scatter_points, colors, markers, edgecolors = plot_coords_by_pop(coords_by_pop, name, savefig=False)
					suffix = "_by_pop"
				except:
					scatter_points, colors, markers, edgecolors = plot_coords(encoded_train, name, savefig=False)
					suffix = ""

			scatter_points_per_epoch.append(scatter_points)
			colors_per_epoch.append(colors)
			markers_per_epoch.append(markers)
			edgecolors_per_epoch.append(edgecolors)

		make_animation(epochs, scatter_points_per_epoch, colors_per_epoch,
		               markers_per_epoch, edgecolors_per_epoch,
		               os.path.join(results_directory,
		                            "{0}{1}".format("dimred_animation", suffix)))

	if arguments['evaluate']:

		print("Evaluating epochs {}".format(epochs))

		# all metrics assumed to have a single value per epoch
		if arguments['metrics']:
			metric_names = arguments['metrics'].split(",")
		else:
			metric_names = ["f1_score_3"]

		metrics = dict()

		for m in metric_names:
			metrics[m] = []

		ind_pop_list_train = read_h5(encoded_data_file, "ind_pop_list_train")
		pop_list = []

		for pop in ind_pop_list_train[:, 1]:
			try:
				pop_list.append(pop.decode("utf-8"))
			except:
				pass

		for epoch in epochs:
			print("########################### epoch {0} ###########################".format(epoch))

			encoded_train = read_h5(encoded_data_file, "{0}_encoded_train".format(epoch))

			coords_by_pop = get_coords_by_pop(data_prefix, encoded_train, ind_pop_list = ind_pop_list_train)

			### count how many f1 scores were doing
			f1_score_order = []
			num_f1_scores = 0
			for m in metric_names:
				if m.startswith("f1_score"):
					num_f1_scores += 1
					f1_score_order.append(m)

			f1_scores_by_pop = {}
			f1_scores_by_pop["order"] = f1_score_order

			for pop in coords_by_pop.keys():
				f1_scores_by_pop[pop] = ["-" for i in range(num_f1_scores)]
			f1_scores_by_pop["avg"] = ["-" for i in range(num_f1_scores)]

			for m in metric_names:

				if m == "hull_error":
					coords_by_pop = get_coords_by_pop(data_prefix, encoded_train, ind_pop_list = ind_pop_list_train)
					n_latent_dim = encoded_train.shape[1]
					if n_latent_dim == 2:
						min_points_required = 3
					else:
						min_points_required = n_latent_dim + 2
					hull_error = convex_hull_error(coords_by_pop, plot=False, min_points_required= min_points_required)
					print("------ hull error : {}".format(hull_error))

					metrics[m].append(hull_error)

				elif m.startswith("f1_score"):
					this_f1_score_index = f1_score_order.index(m)

					k = int(m.split("_")[-1])
					# num_samples_required = np.ceil(k/2.0) + 1 + (k+1) % 2
					num_samples_required = 1

					pops_to_use = get_pops_with_k(num_samples_required, coords_by_pop)

					if len(pops_to_use) > 0 and "{0}_{1}".format(m, pops_to_use[0]) not in metrics.keys():
						for pop in pops_to_use:
							try:
								pop = pop.decode("utf-8")
							except:
								pass
							metric_name_this_pop = "{0}_{1}".format(m, pop)
							metrics[metric_name_this_pop] = []


					f1_score_avg, f1_score_per_pop = f1_score_kNN(encoded_train, pop_list, pops_to_use, k = k)
					print("------ f1 score with {0}NN :{1}".format(k, f1_score_avg))
					metrics[m].append(f1_score_avg)
					assert len(f1_score_per_pop) == len(pops_to_use)
					f1_scores_by_pop["avg"][this_f1_score_index] =  "{:.4f}".format(f1_score_avg)

					for p in range(len(pops_to_use)):
						try:
							pop = pops_to_use[p].decode("utf-8")
						except:
							pop = pops_to_use[p]

						metric_name_this_pop = "{0}_{1}".format(m, pop)
						metrics[metric_name_this_pop].append(f1_score_per_pop[p])
						f1_scores_by_pop[pops_to_use[p]][this_f1_score_index] =  "{:.4f}".format(f1_score_per_pop[p])

				else:
					print("------------------------------------------------------------------------")
					print("Error: Metric {0} is not implemented.".format(m))
					print("------------------------------------------------------------------------")

			write_f1_scores_to_csv(results_directory, "epoch_{0}".format(epoch), superpopulations_file, f1_scores_by_pop, coords_by_pop)

		for m in metric_names:

			plt.plot(epochs, metrics[m], label="train", c="orange")
			plt.xlabel("Epoch")
			plt.ylabel(m)
			plt.savefig(os.path.join(results_directory, m+".png"), dpi=300)
			plt.close()

			outfilename = os.path.join(results_directory, m+".csv")
			with open(outfilename, mode='w') as res_file:
				res_writer = csv.writer(res_file, delimiter=',')
				res_writer.writerow(epochs)
				res_writer.writerow(metrics[m])

	if arguments['plot']:

		print("Plotting epochs {}".format(epochs))

		ind_pop_list_train = read_h5(encoded_data_file, "ind_pop_list_train")
		pop_list = []

		for pop in ind_pop_list_train[:, 1]:
			try:
				pop_list.append(pop.decode("utf-8"))
			except:
				pass

		for epoch in epochs:
			print("########################### epoch {0} ###########################".format(epoch))

			encoded_train = read_h5(encoded_data_file, "{0}_encoded_train".format(epoch))

			coords_by_pop = get_coords_by_pop(data_prefix, encoded_train, ind_pop_list = ind_pop_list_train)

			if superpopulations_file:

				coords_by_pop = get_coords_by_pop(data_prefix, encoded_train, ind_pop_list = ind_pop_list_train)

				if doing_clustering:
					plot_clusters_by_superpop(coords_by_pop,
					                          os.path.join(results_directory,
					                                       "clusters_e_{}".format(epoch)),
					                          superpopulations_file,
					                          write_legend = epoch == epochs[0])
				else:
					scatter_points, colors, markers, edgecolors = \
						plot_coords_by_superpop(coords_by_pop,
						                        os.path.join(results_directory,
						                                     "dimred_e_{}_by_superpop".format(epoch)),
						                        superpopulations_file,
						                        plot_legend = epoch == epochs[0])

			else:
				try:
					plot_coords_by_pop(coords_by_pop,
					                   os.path.join(results_directory,
					                                "dimred_e_{}_by_pop".format(epoch)))
				except:
					plot_coords(encoded_train,
					            os.path.join(results_directory,
					                         "dimred_e_{}".format(epoch)))



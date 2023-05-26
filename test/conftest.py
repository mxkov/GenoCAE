import os
import pytest


## CLI OPTIONS

def pytest_addoption(parser):
	# Since pytest is too stupid to be able to read conftest.py files
	# from subdirectories, I have to use a single conftest.py for ALL tests.
	# To run tests for a specific module, we'll have to use --what.
	parser.addoption('--what',  action = 'store')
	parser.addoption('--where', action = 'store')


## FIXTURES

@pytest.fixture
def f_what(request):
	return request.config.getoption('--what')

@pytest.fixture
def f_where(request):
	return request.config.getoption('--where')

@pytest.fixture
def f_datasplit(request):
	return request.param

# with fewer batches than categories, equal, and more
@pytest.fixture(params=[1, None, 13])
def f_n_batches(request):
	return request.param

@pytest.fixture
def f_batch_col(request):
	return request.param

@pytest.fixture(params=[None, 845])
def f_chunk_size(request):
	return request.param

@pytest.fixture
def f_drop_inds(request):
	return request.param

@pytest.fixture
def f_drop_snps(request):
	return request.param

@pytest.fixture(params=[None, 'even3', 'edges', 'all'])
def f_sessions(request):
	return request.param


## FIXTURE PARAMETRIZATION

def pytest_generate_tests(metafunc):
	what  = metafunc.config.getoption('what')
	if what not in ('pqconv'):
		raise ValueError('what must refer to an existing test')
	where = metafunc.config.getoption('where')
	if where not in (None, 'local', 'UKB'):
		raise ValueError('where must be None, "local" or "UKB"')

	## test for tools/parquet_converter.py
	if what == 'pqconv' and where is not None:
		pq_indir = os.path.join('test', 'test_parquet_converter', 'in')
		# Define Fixture Params depending on where
		fp = {}
		fp['local'] = {
		   'f_datasplit': [1, 10, 16],
		   'f_batch_col': [None, 0],
		   'f_drop_inds': [None,
		                   os.path.join(pq_indir, 'drop_inds_local_1.csv')
		                  ],
		   'f_drop_snps': [None,
		                   os.path.join(pq_indir, 'drop_snps_local_1.csv'),
		                   os.path.join(pq_indir, 'drop_snps_local_2.csv')
		                  ]
		              }
		fp['UKB']   = {
		   'f_datasplit': [1],
		   'f_batch_col': [None, 5],
		   'f_drop_inds': [None,
		                   os.path.join(pq_indir, 'drop_inds_UKB_1.csv')
		                  ],
		   'f_drop_snps': [None,
		                   os.path.join(pq_indir, 'drop_snps_UKB_1.csv')
		                  ]
		              }
		# Set the params
		for fixt in fp[where]:
			if fixt in metafunc.fixturenames:
				metafunc.parametrize(fixt, fp[where][fixt])


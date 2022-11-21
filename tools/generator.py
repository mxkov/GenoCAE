def get_markers():
	'''
	Returns a list of SNP IDs used to generate phenotypes.
	The SNP IDs are taken from the .snp file in the input dataset.
	'''
	return ['rs6082956', 'rs2164169']


def generate(x1, x2):
	'''
	Specifies the exact operations performed on the SNPs to generate phenotypes.
	The number of arguments should be equal to the number of SNP IDs returned by get_markers().
	'''
	a = 2.5
	b = 1.5
	return a*x1 + b*x2


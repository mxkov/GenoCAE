"""Exceptions for parquet_converter.py"""

class CustomException(Exception):
	def __init__(self, msg):
		super().__init__(msg)

class DuplicateVariantsException(Exception):
	def __init__(self, msg):
		super().__init__(msg)

class BadBatchesException(Exception):
	def __init__(self, msg):
		super().__init__(msg)

class BEDMismatchException(Exception):
	def __init__(self, msg):
		super().__init__(msg)

class NoItemsException(Exception):
	def __init__(self, msg):
		super().__init__(msg)

class InsufficientRAMException(Exception):
	def __init__(self, msg):
		super().__init__(msg)

class FileExistsException(Exception):
	def __init__(self, msg):
		super().__init__(msg)


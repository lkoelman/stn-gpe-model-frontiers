"""
Logging utilities.

CREDITS

MultiProcessingHandler was copied from https://github.com/jruere/multiprocessing-logging

"""

from __future__ import absolute_import, division, unicode_literals

import logging
import multiprocessing
import sys
import threading
import traceback

DEFAULT_FORMAT = '[%(levelname)s::%(name)s@%(filename)s:%(lineno)s] %(message)s'

################################################################################
# Custom logging functions
################################################################################

# Create a very verbose log level named "ANAL"

# see https://docs.python.org/2/library/logging.html#logging-levels
DEBUG_ANAL_LVL = 5 
logging.addLevelName(DEBUG_ANAL_LVL, "ANAL")

def anal_logfun(self, message, *args, **kws):
	if self.isEnabledFor(DEBUG_ANAL_LVL):
		# Yes, logger takes its '*args' as 'args'.
		self._log(DEBUG_ANAL_LVL, message, args, **kws) 

logging.Logger.anal = anal_logfun


################################################################################
# Logging helpers
################################################################################

getLogger = logging.getLogger

def configure(logger, level=logging.WARNING, format=DEFAULT_FORMAT,
			  stream=sys.stderr, propagate=1):
	"""
	Configure logger.
	"""
	pass

def setFormat(logger=None, format=DEFAULT_FORMAT):
	"""
	Set logging format for the root logger.

	Call this function before importing any modules that configure the root logger.
	"""
	if logger is None:
		logger = logging.getLogger() # None -> root logger
	if len(logger.handlers) == 0:
		fmt = logging.Formatter(format)
		handler = logging.StreamHandler(stream=sys.stderr)
		handler.setLevel(logger.level)
		handler.setFormatter(fmt)
		logger.addHandler(handler)
	else:
		fmt = logging.Formatter(format)
		for h in logger.handlers:
			h.setFormatter(fmt)


def setConfig(level=logging.WARNING, format=DEFAULT_FORMAT):
	"""
	Configure level and format for all loggers.
	"""
	logging.basicConfig(level=level, format=format)


def getBasicLogger(name=None, level=logging.WARNING, format=DEFAULT_FORMAT,
				   stream=sys.stderr, propagate=1, propagate_config=False):
	"""
	Similar to logging.basicConfig() except this works on a new logger rather
	than on the root logger (logging.root).

	@param	copy_handlers_from		Logger or list(Logger) to copy handlers from

	@note	see https://github.com/python/cpython/blob/master/Lib/logging/__init__.py
	"""
	if isinstance(level, str):
		level = getLogLevel(level)

	# root = logging.getLogger()
	# if len(root.handlers) == 0:
	# 	# FIXME: this changes properties of all loggers
	# 	logging.basicConfig(level=logging.DEBUG, format=format)

	# Creates logger if not present, returns RootLogger if name is None
	# By default messages will propagate to parent (root logger) handlers
	logger = logging.getLogger(name)
	logger.setLevel(level)
	logger.propagate = propagate

	# Get handlers
	# handlers = list(logger.handlers) # new list
	# if copy_handlers_from is not None:
	# 	if isinstance(copy_handlers_from, logging.Logger):
	# 		src_loggers = [copy_handlers_from]
	# 	else:
	# 		src_loggers = copy_handlers_from
	# 	for src in src_loggers:
	# 		if src is not logger:
	# 			handlers.extend(src.handlers)
	
	if not logger.propagate and (len(logger.handlers) == 0):
		# Create new handler 
		fmt = logging.Formatter(format)
		handler = logging.StreamHandler(stream=stream)
		handler.setLevel(level)
		handler.setFormatter(fmt)
		logger.addHandler(handler)

	if propagate_config:
		# This configures the root logger
		logging.basicConfig(level=level, format=format)


	# for h in handlers:
	# 	# Only set formatter if explicitly specified or none present
	# 	# if (h.formatter is None) or (usr_format is not None):
	# 	# 	h.setFormatter(fmt)
		
	# 	# Add handler if handler with same stream not present
	# 	if h in logger.handlers:
	# 		continue
	# 	elif any((h.stream==hdlr.stream for hdlr in logger.handlers)):
	# 		# don't log to the same stream twice
	# 		continue
	# 	else:
	# 		logger.addHandler(h)

	return logger


def getLogLevel(level):
	"""
	Get log level from string
	"""
	if level in logging._nameToLevel:
		return level
	elif level.lower() in ['verbose', 'debug']: # DEBUG equivalents
		return logging.DEBUG
	elif level.lower() in ['silent', 'quiet', 'warning']: # WARNING equivalents
		return logging.WARNING
	elif level.lower() in ['anal', 'highly_verbose']: # ANAL equivalents
		return DEBUG_ANAL_LVL
	else:
		raise ValueError(level)


def setLogLevel(level, logger_names=None):
	"""
	Set log level of all loggers with given names to level.

	@param level	any log level accepted by logging.setLevel() or one of following:
					'verbose', 'quiet', 'silent', 'anal'
	"""
	level = getLogLevel(level)

	if logger_names is None:
		logger_names = logging.Logger.manager.loggerDict.keys()

	for logname in logger_names:
		if logname in logging.Logger.manager.loggerDict:
			slogger = logging.getLogger(logname) # creates new if not in loggerDict
			slogger.setLevel(level)


def install_mp_handler(logger=None):
	"""
	Wraps the handlers in the given Logger with an MultiProcessingHandler.

	:param logger: whose handlers to wrap. By default, the root logger.

	@note	this does not seem to work with ipyparallel, since it does not
			use the native multiprocessing module

	USAGE

		import logutils
		logutils.install_mp_handler()
	
	"""
	if logger is None:
		logger = logging.getLogger()

	for i, orig_handler in enumerate(list(logger.handlers)):
		handler = MultiProcessingHandler(
			'mp-handler-{0}'.format(i), sub_handler=orig_handler)

		logger.removeHandler(orig_handler)
		logger.addHandler(handler)


class MultiProcessingHandler(logging.Handler):
	"""
	Usage: see install_mp_handler()
	"""

	def __init__(self, name, sub_handler=None):
		super(MultiProcessingHandler, self).__init__()

		if sub_handler is None:
			sub_handler = logging.StreamHandler()

		self.sub_handler = sub_handler
		self.queue = multiprocessing.Queue(-1)
		self.setLevel(self.sub_handler.level)
		self.setFormatter(self.sub_handler.formatter)
		# The thread handles receiving records asynchronously.
		t = threading.Thread(target=self.receive, name=name)
		t.daemon = True
		t.start()

	def setFormatter(self, fmt):
		logging.Handler.setFormatter(self, fmt)
		self.sub_handler.setFormatter(fmt)

	def receive(self):
		while True:
			try:
				record = self.queue.get() # get messages put on queue from other threads
				self.sub_handler.emit(record) # put them on the local thread stream handler
			except (KeyboardInterrupt, SystemExit):
				raise
			except EOFError:
				break
			except:
				traceback.print_exc(file=sys.stderr)

	def send(self, s):
		self.queue.put_nowait(s)

	def _format_record(self, record):
		# ensure that exc_info and args
		# have been stringified. Removes any chance of
		# unpickleable things inside and possibly reduces
		# message size sent over the pipe.
		if record.args:
			record.msg = record.msg % record.args
			record.args = None
		if record.exc_info:
			self.format(record)
			record.exc_info = None

		return record

	def emit(self, record):
		try:
			s = self._format_record(record)
			self.send(s)
		except (KeyboardInterrupt, SystemExit):
			raise
		except:
			self.handleError(record)

	def close(self):
		self.sub_handler.close()
		logging.Handler.close(self)


class ZeroMQSocketHandler(object):
	"""
	TODO: implement ZMQ stream handler for loggers on ipyparallel engines

	http://olympiad.cs.uct.ac.za/docs/python-docs-3.2/howto/logging-cookbook.html#subclassing-queuehandler-a-zeromq-example

	http://ipyparallel.readthedocs.io/en/latest/development/connections.html#iopub

		Example of how to use this to subscribe to engines' stdout/stderr: https://github.com/ipython/ipyparallel/blob/master/examples/iopubwatcher.py

		Strategy: look up code where stdout/stderr is published to ZMQ iopub socket, and
		use a similar approach in this handler class to send/receive logging data.
	"""
	pass

"""
Utilities for working with Jupyter and the Juputer notebook.

@author	Lucas Koelman

@date	30/01/2018
"""

import matplotlib.pyplot as plt
from IPython.core.ultratb import AutoFormattedTB

def notebook_show_figs_after_exception():
	"""
	Insert notebook cell exception handler that fixes matplotlib notebook issue 
	where figures don't show up after exception is thrown anywhere in notebook
	"""

	# initialize the formatter for making the tracebacks into strings
	itb = AutoFormattedTB(mode = 'Plain', tb_offset = 1)

	# this function will be called on exceptions in any cell
	def custom_exc(shell, etype, evalue, tb, tb_offset=None):

	    # still show the error within the notebook, don't just swallow it
	    shell.showtraceback((etype, evalue, tb), tb_offset=tb_offset)

	    # grab the traceback and make it into a list of strings
	    stb = itb.structured_traceback(etype, evalue, tb)
	    sstb = itb.stb2text(stb)
	    
	    # Our custom code
	    for fig_num in plt.get_fignums():
	        plt.figure(fig_num).canvas.draw()

	# this registers a custom exception handler for the whole current notebook
	get_ipython().set_custom_exc((Exception,), custom_exc)
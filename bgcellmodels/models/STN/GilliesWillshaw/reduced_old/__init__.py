# override modules to export
# __all__ = ["module_a", "module_b", "module_c"]

# make classes available at package level
# from mymodule import myclass

# Make sure top-level modules are on Python path
# import sys, os.path
# scriptdir, scriptfile = os.path.split(__file__)
# repo_root = os.path.normpath(os.path.join(scriptdir, '..', '..'))
# sys.path.append(repo_root)

# Create new logging level to follow low-level operations (extremely verbose)
import logging
DEBUG_ANAL_LVL = 5 # see https://docs.python.org/2/library/logging.html#logging-levels
logging.addLevelName(DEBUG_ANAL_LVL, "ANAL")

def logfun(self, message, *args, **kws):
    if self.isEnabledFor(DEBUG_ANAL_LVL):
    	# Yes, logger takes its '*args' as 'args'.
        self._log(DEBUG_ANAL_LVL, message, args, **kws) 

logging.Logger.anal = logfun
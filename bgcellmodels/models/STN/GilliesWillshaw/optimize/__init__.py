# override modules to export
# __all__ = ["module_a", "module_b", "module_c"]

# make classes available at package level
# from mymodule import myclass

# Make sure top-level modules are on Python path
import sys, os.path
scriptdir, scriptfile = os.path.split(__file__)
repo_root = os.path.normpath(os.path.join(scriptdir, '..', '..'))
sys.path.append(repo_root)
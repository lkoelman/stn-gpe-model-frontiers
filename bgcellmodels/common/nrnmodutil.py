"""
Utilities for dealing with NEURON MOD files

@author Lucas Koelman
"""

import re

def get_mod_name(hobj):
    """
    Get NEURON mechanism name of given synapse object

    @param  hobj        HocObject: synapse POINT_PROCESS
    """
    if hasattr(hobj, 'htype'):
        hoc_name = hobj.htype.hname() # for wrapped HocObject, mechanism name is in htype attribute
    else:
        hoc_name = hobj.hname()
    match_mod = re.search(r'^[a-zA-Z0-9]+', hoc_name)
    modname = match_mod.group()
    return modname


def get_mod_block(mod_text, block_name):
    """
    Read contents of named MOD file block.
    """
    re_block = re.compile("^" + block_name + r"\s+{(.*?)}", re.DOTALL | re.MULTILINE)
    match = re_block.search(mod_text)
    if match is None:
        raise ValueError(block_name + " block not found in mod file.")
    if len(match.groups()) != 1:
        raise ValueError("Found {} occurrences of {} block".format(
                         len(match.groups()), block_name))
    return match.group(0)


def get_mod_mechname(mod_file, filename=True):
    """
    Read the NEURON mechanism name from contents of .mod file
    """
    if filename:
        with open(mod_file, 'r') as modfile:
            modtext = modfile.read()
    else:
        modtext = mod_file

    neuron_block = get_mod_block(modtext, 'NEURON')
    re_mechname = re.compile(r"(ARTIFICIAL_CELL|POINT_PROCESS|SUFFIX)[ \t]+(\w+)")
    match = re_mechname.search(neuron_block)
    if match is None or len(match.groups()) != 2:
        raise ValueError("No singular match for mechanism name in block:\n" + neuron_block)
    return match.group(2) # 0 is full match


def get_mod_parameters(mod_file, filename=True):
    """
    Get contents of PARAMETERS block in mod file.

    @param  modfile : str
            Full path to MOD file.

    @param  filename : bool
            True if mod_file is a filename, False if it's the contents

    @return parameters : list(tuple(str, float, str, str))
            A list of parameters. Each parameter is specified as a tuple
            (name, value, units, comment).
    """
    if filename:
        with open(mod_file, 'r') as modfile:
            modtext = modfile.read()
    else:
        modtext = mod_file

    str_param_block = get_mod_block(modtext, 'PARAMETER')

    # regex tested using https://regex101.com/ and GLUsyn.mod file
    re_param_defs = re.compile(
        r"^[ \t]*(\w+)[ \t]*=[ \t]*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)[ \t]*\(?(\w+)?\)?[ \t]*:?(.*?\n)?",
        re.MULTILINE)

    # group index and converter for name, value, units, comment
    tokens = ((1, str), (2, float), (4, str), (5, lambda t: t.strip()))
    param_defs = [tuple(token[1](match.group(token[0])) for token in tokens) 
                        for match in re_param_defs.finditer(str_param_block)]
    return param_defs


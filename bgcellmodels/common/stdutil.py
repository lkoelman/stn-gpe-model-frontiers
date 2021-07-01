"""
General Python utility functions, additions to the standard library.

@author     Lucas Koelman
"""

from enum import IntEnum


class dotdict(dict):
    """
    dot.notation access to dictionary attributes.
    """
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


class Bunch(object):
    """
    Bunch or struct-like object for data storage using dot syntax
    """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


class IntEnumDescriptor(IntEnum):

    def to_str(self):
        """
        Convert Enum instance to string
        """
        return self.name.lower()

    @classmethod
    def from_str(cls, descr):
        """
        Get Enum instance from string.
        """
        return cls._member_map_[descr.upper()]

    @classmethod
    def get(cls, descr):
        """
        Get Enum instance from description.
        """
        if isinstance(descr, cls):
            return descr
        elif isinstance(descr, str):
            return cls.from_str(descr)
        else:
            raise ValueError("Cannot convert {} to class {}".format(
                    descr, cls))

    to_descr = to_str
    from_descr = from_str


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    """
    same as Python >= 3.5 math.isclose
    """
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def getdictvals(d, *args, **kwargs):
    """
    Get dictionary values at given keys as tuple.

    @param      *args : *list(str)
                Expanded list of keys of arbitrary length

    @param      as_dict: bool (keyword argument)
                Return result as dict instead of tuple
    
    @return     vals : tuple
                Value associated with given keys
    """
    if kwargs.get("as_dict", False):
        return {k: d[k] for k in args}
    else:
        return tuple((d[k] for k in args))


def eval_context(**context):
    """
    Evaluate statement in given context.

    @param      statement : str
                First argument to eval()

    @param      do_format : bool
                If true, format(**locals) will be called on statement
                before passing it to eval().

    @param      globals : dict
                Second argument to eval(). If None, use 'caller_globals' if it
                is given, else use empty dict {}.

    @param      locals : dict
                Third argument to eval(). If None, use 'caller_locals' if it
                is given, else use empty dict {}.

    @param      caller_locals : dict OR list(dic)
                Dicts containing local variables to be used as third argument
                to 'eval()'. Dicts entries be updated (with override) in order
                of appearance and finally with 'locals' argument if not None. 

    @return     result of eval() of statement with given locals and globals.
    """

    stmt = context["statement"]

    # None means copy 
    context_locals = context.get("locals", None)
    
    # 'caller_locals' can be either a dict or a list of candidate dicts
    caller_locals = context.get("caller_locals", {})
    if isinstance(caller_locals, dict):
        locals_dict = dict(caller_locals)
    else: # accumulate passed dicts, with override!
        locals_dict = {}
        for ldict in caller_locals:
            locals_dict.update(ldict)

    # Final update must be context locals to ensure that keys are not overridden
    if context_locals is not None:
        locals_dict.update(context_locals)


    globals_dict = context.get("globals", None)
    if not isinstance(globals_dict, dict):
        globals_dict = context.get("caller_globals", {}) # if default is None it will use the ones from this module


    if context.get("do_format", False):
        stmt = stmt.format(**locals_dict)

    return eval(stmt, globals_dict, locals_dict)

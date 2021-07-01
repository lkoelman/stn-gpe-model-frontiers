"""
NEURON unit interopability using 'Pint' or 'Quantities' package.

@author Lucas Koelman

@see    http://pint.readthedocs.io
"""

import re
from neuron import h


def set_units_module(module_name='pint'):
    global QUANTITY_TYPE, GET_UNITS_FUNC, CONVERT_UNITS_FUNCNAME
    global DIMENSIONALITY_ERROR, UNITS_MODULE_NAME, PINT_IMPORTED, QUANTITIES_IMPORTED

    if module_name == 'pint':

        import pint
        PINT_IMPORTED = True

        ureg = pint.UnitRegistry()
        QUANTITY_TYPE= ureg.Quantity
        # QUANTITY_BASE = pint.quantity._Quantity # Base class
        GET_UNITS_FUNC = ureg
        CONVERT_UNITS_FUNCNAME = 'to'
        DIMENSIONALITY_ERROR = pint.errors.DimensionalityError

        # Define NEURON units that are not in pint's default units
        # TODO: add full list of NMODL/modlunit units
        ureg.define('Ohm = ohm')
        ureg.define('mho = 1/ohm')
        ureg.define('cm2 = cm^2')
        ureg.define('m2 = m^2')

    elif module_name == 'quantities':

        import quantities as pq

        QUANTITIES_IMPORTED = True
        QUANTITY_TYPE = pq.Quantity
        GET_UNITS_FUNC = pq.unit_registry.__getitem__
        CONVERT_UNITS_FUNCNAME = 'rescale'
        DIMENSIONALITY_ERROR = ValueError

        pq.UnitQuantity('nanovolt', pq.V * 1e-9, symbol='nV')
        pq.UnitQuantity('nanoAmpere', pq.A * 1e-9, symbol='nA')
        pq.UnitQuantity('mho', 1 / pq.Ohm, symbol='mho')
        pq.UnitQuantity('cm2', pq.cm ** 2, symbol='cm2')
        pq.UnitQuantity('m2', pq.m ** 2, symbol='m2')

    else:
        raise ValueError('Units module "{}" not supported'.format(module_name))

    UNITS_MODULE_NAME = module_name
    print("Using units module '{}'".format(UNITS_MODULE_NAME))


def Quantity(*args, **kwargs):
    try:
        make_quantity = QUANTITY_TYPE
    except NameError:
        set_units_module()
        make_quantity = QUANTITY_TYPE
    return make_quantity(*args, **kwargs)


def get_nrn_units(nrn_obj, attr, hoc_classname=None):
    """
    Get units of NEURON variable as a pint.Quantity object.

    @param  nrn_obj : nrn.HocObject
            NEURON object

    @param  attr : str
            Variable name of NEURON object

    @param  hoc_classname : str
            if 'attr' is not a mechanism, Section, or global Hoc variable name,
            specify the Hoc classname here. E.g. to set 'Exp2Syn.tau1',
            use attr='tau1' and hoc_classname='Exp2Syn'

    @return q : pint.Quantity
            Units of given NEURON variable

    @throws ValueError
            Variable name is not found by Hoc.
    """
    # TODO: see h.units() documentation: extract classname from nrn_obj
    #       so we can pass it to h.units('classname.attr'), e.g. if nrn_obj
    #       is h.Exp2Syn
    if hoc_classname is None:
        full_attr = attr
    else:
        full_attr = '{}.{}'.format(hoc_classname, attr)
    try:
        # NOTE: h.units() can return '' -> treated as dimensionless by Pint
        # NOTE: h.units() does not accept unicode strings
        nrn_units = h.units(str(full_attr))
    except RuntimeError as e:
        if e.args[0] == 'hoc error':
            raise ValueError('Attribute {} not recognized by hoc'.format(attr))

    units_nodash = nrn_units.replace('-', '*')
    units_exponents = re.sub(r'([a-zA-Z]+)(\d)',r'\1^\2', units_nodash)
    
    try:
        make_units_func = GET_UNITS_FUNC
    except NameError:
        set_units_module()
        make_units_func = GET_UNITS_FUNC

    target_units = make_units_func(units_exponents)
    return target_units


def to_nrn_units(quantity, nrn_obj, attr, hoc_classname=None):
    """
    Convert quantity to same units as NEURON variable.

    ARGUMENTS
    ---------

    @param  nrn_obj : nrn.HocObject
            Top-level Hoc interpreter

    @param  attr : str
            Variable name of NEURON object

    @param  quantity : pint.Quantity
            Quantity object consisting of value and units

    @param  hoc_classname : str
            if 'attr' is not a mechanism, Section, or global Hoc variable name,
            specify the Hoc classname here. E.g. to set 'Exp2Syn.tau1',
            use attr='tau1' and hoc_classname='Exp2Syn'


    EXCEPTIONS
    ----------

    @throws err : pint.errors.DimensionalityError
            Error thrown in case of dimensionality (not units) mismatch.

    @return q : pint.Quantity
            Original quantity converted to units of the NEURON object.


    USAGE
    -----
    
        > import neuron, units
        > quantity = units.Quantity(value, param_spec['units'])
        > converted_quantity = units.to_nrn_units(quantity, neuron.h, 'gnabar_hh')
        > value = converted_quantity.magnitude
    """
    target_units = get_nrn_units(nrn_obj, attr, hoc_classname)
    return getattr(quantity, CONVERT_UNITS_FUNCNAME)(target_units)


def compatible_units(nrn_obj, attr, quantity, hoc_classname=None):
    """
    Check if units of given quantity are compatible with those of NEURON variable.

    @see    get_nrn_units() for description of arguments

    @param  quantity : pint.Quantity
            Quantity object consisting of value and units

    @return compatible : bool
            True if units are compatible.
    """
    try:
        to_nrn_units(quantity, nrn_obj, attr, hoc_classname)
    except DIMENSIONALITY_ERROR:
        return True
    else:
        return False


def same_units(nrn_obj, attr, quantity, hoc_classname=None):
    """
    Check if units of given quantity are the same as those of NEURON variable.

    @see    get_nrn_units() for description of arguments

    @param  quantity : pint.Quantity
            Quantity object consisting of value and units

    @return same : bool
            True if units are the same
    """
    nrn_quantity = get_nrn_units(nrn_obj, attr, hoc_classname)
    return nrn_quantity.units == quantity.units


def set_nrn_quantity(nrn_obj, attr, quantity, hoc_classname=None):
    """
    Set attribute of NEURON object using a pint.Quantity object that includes
    a value and units.

    @see    get_nrn_units() for description of arguments

    @param  quantity : pint.Quantity
            Quantity object consisting of value and units

    @throws err : pint.errors.DimensionalityError
            Error thrown in case of dimensionality (not units) mismatch.
    """
    val_converted = to_nrn_units(quantity, nrn_obj, attr, hoc_classname)
    setattr(nrn_obj, attr, val_converted.magnitude)


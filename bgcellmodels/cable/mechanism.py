"""
Cell mechanism
"""

from __future__ import division
import inspect, functools
import sympy
import numpy as np

import bgcellmodels.common.units as pint_units

from symtypes import QuantitativeExpr
# NOTE:
# Pint quantity has attributes dimensionality, units, magnitude
# Sympy quantity has attributes dimensionality, scale_factor

# See http://docs.sympy.org/latest/modules/numeric-computation.html
# and backends with examples in https://github.com/sympy/sympy/blob/master/sympy/printing
# from sympy.utilities.lambdify import lambdify # sympy.lambdify
# from sympy.utilities.autowrap import ufuncify
# from sympy.printing.theanocode import theano_function
# from sympy.printing.llvmjitcode import llvm_callable


# Symbolic computation with Sympy
exp = sympy.exp
log = sympy.log
# pow is built-in function

# Shorthands for units
QuantityType = pint_units.QuantityType
ureg = pint_units.ureg # so we can use variable name 'units'
mV  = ureg('mV')
mA  = ureg('mA')
ms  = ureg('ms')
S   = ureg('S')
cm2 = ureg('cm^2')
m2  = ureg('m^2')
nodim = ureg('dimensionless')

# Symbolic variable for membrane voltage
v = sympy.symbols('v')
# Vm = v * mV

from bgcellmodels.common import logutils
logger = logutils.getBasicLogger('mechanism')

# TODO: find good way to combine symbolic expressions and units.
#       You cannot expr.subs(quantity), but you can multiply
#       a quantity with a symbolic expression and retrieve it via
#       '.magnitude' property. However, we want to ensure that when we
#       substitute quantities into a symbolic expression, all the units
#       match up.

# SEE:  http://docs.sympy.org/latest/modules/physics/units/index.html
#       https://github.com/sympy/sympy/tree/master/sympy/physics/units/tests
#       Possibly use sympy built-in unit system.

# METHOD 1
#       Make wrapper object that keeps for each symbolic variable, the
#       Symbol and units object. Can do this by subclassing Symbol/Quantity,
#       depending on which behaviour you want to keep. Or by creating true
#       wrapper class ~Neo.spiketrain. Then to check units you can do one of
#       the following: lambdify the expression, pass quantities and see
#       that units match (no error) and final units are correct. You can
#       lambdify by finding intersection with expr.free_symbols and providing
#       that as arguments list, then give params in that order.
#       This method allows you to use sympy or external units/quantities module.

# TODO: Use method 1, and make everything a quantitative expression. When doing
#       operation, simultaneously combine the associated units and return object
#       with combined units. The problem is that both sympy and Pint try to
#       'hijack' numerical operationts (addition,multiplication,...) so you need
#       to update one of the classes to play well with the other.

# METHOD 2
#       Use sympy.physics.units.Quantity object vor all variables.
#       See http://docs.sympy.org/latest/modules/physics/units/examples.html
#       E.g. first you can subs() q.dimensions and see that final dimensionality
#       is correct, then subs() q.scale_factor and check that final scale is correct.


class MechanismType(type):
    """
    Metaclass for cell mechanism.

    @see    https://docs.python.org/3/reference/datamodel.html#metaclasses

    @note   metaclass returns a new type for classes that use it
    """

    def __new__(mcls, name, bases, namespace, **kwds):
        """
        Create a class object for a class definition that has MechanismType
        as its metaclass.

        @param      mcls : type
                    the metaclass object (prototype describing type MechanismType)

        @return     cls : MechanismType
                    a new Class object (an instance of type MechanismType)

        @note       This method is called once every time a class is defined with 
                    MechanismType as its metaclass. Whereas hat class' __new__ 
                    method is called every time an instance is created.
        """
        # Process mechanism variables declared in class definition
        newattrs = MechanismType.group_namespace_vars(namespace)
        cls = super(MechanismType, mcls).__new__(mcls, name, bases, newattrs)
        logger.debug("{} :: called {} and returned {}".format(mcls, '__new__', cls))
        return cls


    def __init__(cls, name, bases, namespace):
        """
        @param      cls : type
                    the class object returned by __new__

        @note       This method is called once every time a class is defined with 
                    MechanismType as its metaclass. Whereas hat class' __new__ 
                    method is called every time an instance is created.
        """
        obj = super(MechanismType, cls).__init__(name, bases, namespace)
        # Common parameters and dynamics for mechanism class
        # NOTE: only call class method, not static method in MechanismBase
        cls.PARAMETERS()
        cls.DYNAMICS()
        logger.debug("{} :: called {} and returned {}".format(cls, '__init__', obj))
        return obj



    @staticmethod
    def define_parameter(name, default_value, units):
        # param = SymbolicQuantity(name, default_value=default_value, units=units)
        param = QuantitativeExpr(name, default_val=default_value, units=units)
        param._param_name = name
        param._mech_attr_type = 'parameter'
        return param


    @staticmethod
    def define_state(name, power, units='dimensionless'):
        state = QuantitativeExpr(name, default_val=1.0, units=units)
        state._mech_attr_type = 'statevar'
        state._state_name = name
        state.power = power
        return state


    @staticmethod
    def state_steadystate(expr, state):
        # TODO: replace exprelr by sympy Piecewise
        expr = QuantitativeExpr(expr)
        expr._mech_attr_type = 'statevar_steadystate'
        expr._state_name = state
        return expr


    @staticmethod
    def state_timeconst(expr, state):
        expr = QuantitativeExpr(expr)
        expr._mech_attr_type = 'statevar_timeconst'
        expr._state_name = state
        return expr


    @staticmethod
    def state_derivative(expr, state):
        expr = QuantitativeExpr(expr)
        expr._mech_attr_type = 'statevar_derivative'
        expr._state_name = state
        return expr


    @staticmethod
    def define_current(expr, ion):
        expr = QuantitativeExpr(expr)
        expr._mech_attr_type = 'current'
        expr._ion_species = ion
        return expr


    @staticmethod
    def group_namespace_vars(namespace):
        """
        Extract mechanism variables declares as variables in given namespace.

        @param  namespace : dict
                The namespace dict, e.g. a result from call to locals()
        """

        # Store definitions of states etc. in structured form
        newattrs = {}
        newattrs['_MECH_PARAMS'] = {}
        newattrs['_MECH_STATE_VARS'] = {}
        newattrs['_MECH_STATE_DERIV'] = {} # expression for derivative of state variable
        newattrs['_MECH_STATE_INF'] = {} # (optional) expression for steady state value of state
        newattrs['_MECH_STATE_TAU'] = {} # (optional) expression for time constant of state
        newattrs['_MECH_MEMB_CURRENTS'] = {} # (optional) expression for time constant of state

        for attr_name, attr_val in namespace.iteritems():
            # Process special attributes created in class scope
            attr_type = getattr(attr_val, '_mech_attr_type', None)
            
            if attr_type == 'parameter':
                newattrs['_MECH_PARAMS'][attr_val._param_name] = attr_val

            elif (attr_type is None
                    and isinstance(attr_val, QuantityType)
                    and isinstance(attr_val.m, (int, float))):
                
                newattrs['_MECH_PARAMS'][attr_name] = MechanismBase.define_parameter(
                                                    attr_name,
                                                    attr_val.magnitude,
                                                    str(attr_val.units))
            
            elif attr_type == 'statevar':
                newattrs['_MECH_STATE_VARS'][attr_val._state_name] = attr_val
            
            elif attr_type == 'statevar_steadystate':
                newattrs['_MECH_STATE_INF'][attr_val._state_name] = attr_val
            
            elif attr_type == 'statevar_timeconst':
                newattrs['_MECH_STATE_TAU'][attr_val._state_name] = attr_val
            
            elif attr_type == 'statevar_derivative':
                newattrs['_MECH_STATE_DERIV'][attr_val._state_name] = attr_val

            elif attr_type == 'current':
                ion_currents = newattrs['_MECH_MEMB_CURRENTS'].setdefault(attr_val._ion_species, {})
                if attr_name in ion_currents.keys():
                    raise ValueError('Duplicate current name \'{}\' for ion \'{}\'').format(
                        attr_name, attr_val._ion_species)
            
            else:
                # NOTE: NEED TO KEEP BUILTIN PYTHON NAMESPACE VARIABLES
                newattrs[attr_name] = attr_val

        return newattrs


def check_inf_nan(expr):
    """
    Check for inf and NaN in Sympy expression

    @return     bool
                True if ok

    @throws     ValueError
                Raises ValueError if any of the atoms of expr is inf or NaN
    """
    # Test for inf and nan
    atoms = expr.atoms(sympy.Number) # leaf nodes of expression tree
    for a in atoms:
        if not a.is_finite or a==sympy.nan:
            raise Exception('Found infinite or NaN in expression {}'.format(expr))


class MechanismBase:
    """
    Base class for Cell Mechanism
    """
    __metaclass__ = MechanismType


#     def __new__(cls, *args, **kwargs):
#         """
#         @param      cls : type
#                   the class object created by metaclass
# 
#         @return     self : MechanismBase
#                   new instance of type cls
#         """
#         logger.debug("{} :: called {}".format(cls, '__init__'))
#         obj = super(MechanismBase, cls).__new__(cls, *args, **kwargs)
#         return obj


    def __init__(self, *args, **kwargs):
        """
        Initialize instantiated mechanism.
        """
        super(MechanismBase, self).__init__(*args, **kwargs)
        # Copy class variables representing parameters so they can be modified
        self._IMECH_PARAMS = {
            p : MechanismBase.define_parameter(
                                v._param_name,
                                v.val,
                                v.units) 
                for p,v in self._MECH_PARAMS.iteritems()
        }


    @classmethod
    def PARAMETERS(subcls):
        """
        Called by subclasses to define parameters.

        @effect     updates instantiated parameters in self._IMECH_PARAMS
                    with parameters declared in local scope of caller.

        USAGE
        -----

            >>> class MyChannel(MechanismBase):
            >>>     @classmethod
            >>>     def PARAMETERS(cls):
            >>>         param1 = val1
            >>>         super(cls, cls).PARAMETERS()
        
        """
        caller_locals_dict = inspect.currentframe().f_back.f_locals
        def_mech_vars = MechanismType.group_namespace_vars(caller_locals_dict)
        subcls._MECH_PARAMS.update(def_mech_vars['_MECH_PARAMS'])
        for pname, pval in def_mech_vars['_MECH_PARAMS'].iteritems():
            setattr(subcls, pname, pval)
    

    @classmethod
    def DYNAMICS(subcls):
        """
        Called by subclasses to define dynamics.

        @effect     injects instantiated parameters in local scope of caller.
        """
        caller_locals_dict = inspect.currentframe().f_back.f_locals
        mech_vars = MechanismType.group_namespace_vars(caller_locals_dict)
        for vars_name in ['_MECH_STATE_VARS', '_MECH_STATE_TAU', '_MECH_STATE_INF',
                            '_MECH_STATE_DERIV', '_MECH_MEMB_CURRENTS']:
            getattr(subcls, vars_name).update(mech_vars[vars_name])


    def get_param(self, name):
        """
        Get instantiated parameter.

        @return     QuantitativeExpr
        """
        return self._IMECH_PARAMS[name]


    def make_v_func(self, expr):
        """
        Make ufunc by substituting mechanism parameters and leaving other symbols
        as free symbols
        """
        raw_expr = expr.as_expr()

        # Make expression with 'v' as only free symbol
        subs_params = [self.get_param(p.name) for p in raw_expr.free_symbols if p!=v]
        substitutions = {p.as_expr():p.val for p in subs_params}
        v_expr = raw_expr.subs(substitutions)

        # Test for inf and nan
        check_inf_nan(v_expr)

        # Debug lambda expression
        logger.debug('Attempting to lambdify expression {}'.format(v_expr))
        if logger.level <= 20:
            from sympy.utilities.lambdify import lambdastr
            from sympy.printing.lambdarepr import NumPyPrinter as Printer
            printer = Printer()
            lstr = lambdastr([v], v_expr, printer=printer)
            logger.debug('Lambda string passed to eval() is: {}'.format(lstr))

        func = sympy.lambdify([v], v_expr, 'numpy') # func = ufuncify(v, raw_expr)
        return func


    def make_full_func(self, expr):
        """
        Turn expression into function with keyword arguments
        equal to its free symbols.

        @return     callable
        """
        
        raw_expr = expr.as_expr()

        # Test for inf and nan
        check_inf_nan(raw_expr)

        # Lambdify
        args = [p for p in raw_expr.free_symbols]
        argnames = [p.name for p in args]
        func = sympy.lambdify(args, raw_expr, 'numpy') # func = ufuncify(v, raw_expr)

        @functools.wraps(func)
        def wrapped(**_kwargs):
            """
            Put kwargs in same order as initial args list passed to lambdify
            """
            args_ordered = [_kwargs[k] for k in argnames]
            return func(*args_ordered)
        
        return wrapped


    def check_units(self, expr):
        """
        Check units by converting expression to a function, and executing it with
        mechanism parameters as Quantities.

        @return     bool
                    True if expression can be evaluated with parameters as quantities

        @throws     pint.
        """
        func = self.make_func(expr)
        # - get all equations
        # for all equations:
        # - make full func
        # - gather parameters
        # - substitute quantity (param.val * param.units)


    def compile_mod():
        # TODO: compile Python representation of channel to .mod file
        pass


    def compile_nrn_c():
        # TODO: compile Python representation directly to C file. See mod2c
        pass


    def compile_nrn_lib():
        # TODO: Compile to library that can be linked into libnrnmech.so, see nrnivmodl script
        #       This will require making an adapted nrnivmodl script.
        pass


    def register_nonvint():
        """
        Register with NEURON's nonvint_block_supervisor

        @see    example at
                https://github.com/nrnhines/nrntest/blob/master/nrniv/rxd/kchanarray.py
                or see RxD module
        """
        # TODO: base on kchanarray example, except pass compiled C-functions as callbacks
        #       rather than slow Python functions
        pass


    def plot_steadystate_gating(self):
        """
        Plot steady state values and time constants of gating vars
        """
        # TODO: collect gating functions and plot
        import matplotlib.pyplot as plt
        for state_name in self._MECH_STATE_VARS.keys():

            # Get steady state and time constant
            s_inf = self._MECH_STATE_INF.get(state_name, None)
            s_tau = self._MECH_STATE_TAU.get(state_name, None)
            if s_inf is None and s_tau is None:
                continue

            # Convert quantity to expression
            fn_inf = self.make_v_func(s_inf)
            fn_tau = self.make_v_func(s_tau)
            v_range = np.linspace(-100, 100, 400)
            sst_range = fn_inf(v_range)
            tau_range = fn_tau(v_range)

            fig, ax = plt.subplots()
            ax.set_title('State variable {}'.format(state_name))
            ax.plot(v_range, sst_range, label=r'$x_{\inf}$')
            ax.plot(v_range, tau_range, label=r'$\tau_x$')
            ax.set_xlabel('V (mV)')
            ax.legend()

            # fig, axes = plt.subplots(2,1)
            # axes[0].plot(v_range, fn_inf(v_range))
            # axes[0].text(.5, .5, '${}$'.format(sympy.latex(s_inf.expr)),
            #                             horizontalalignment='right',
            #                             verticalalignment='top')
            # axes[1].plot(v_range, fn_tau(v_range))
            # axes[1].text(.5, .5, '${}$'.format(sympy.latex(s_tau.expr)),
            #                             horizontalalignment='right',
            #                             verticalalignment='top')

        plt.show(block=False)


# Alias static methods
mech = MechanismType
Parameter = mech.define_parameter
State = mech.define_state
Current = mech.define_current
state_steadystate = mech.state_steadystate
state_timeconst = mech.state_timeconst
state_derivative = mech.state_derivative


class HHNaChannel(MechanismBase):
    """
    Hodgkin-Huxley Na channel mechanism.
    """

    # PARAMETER block
    gnabar = Parameter('gnabar', 0.12, 'S/cm^2')
    ENa = Parameter('gnabar', 50.0, 'mV')

    # STATE block
    m = State('m', power=3)
    h = State('h', power=1)

    # RHS_EXPRESSIONS
    ## m gate
    alpha_m = 0.1 * -(v+40) / (exp(-(v+40)/(10)) - 1)
    beta_m = 4.0 * exp(-(v+65)/(18))
    minf = state_steadystate(alpha_m / (alpha_m + beta_m), state='m')
    tau_m = state_timeconst(1.0 / (alpha_m + beta_m), state='m')

    ## h gate
    alpha_h = 0.07 * exp(-(v+65)/(20))
    beta_h = 1.0 / (exp(-(v+35)/(10)) + 1)
    hinf = state_steadystate(alpha_h / (alpha_h + beta_h), state='h')
    tau_h = state_timeconst(1.0 / (alpha_h + beta_h), state='h')

    # DERIVATIVE block
    dm = state_derivative((minf - m) / tau_m, state='m')
    dh = state_derivative((hinf - h) / tau_h, state='h')

    # BREAKPOINT block
    ina = Current(gnabar * m**m.power * h**h.power * (v-ENa), ion='na')


if __name__ == '__main__':
    # Create channel
    chan = HHNaChannel()
    chan.plot_steadystate_gating()
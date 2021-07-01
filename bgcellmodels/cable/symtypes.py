"""
Sympy expressions with units
"""

from __future__ import division
import sympy
import sympy.physics.units as sympy_units

Expr = sympy.Expr

from bgcellmodels.common.units import QuantityType
from bgcellmodels.common import logutils
logger = logutils.getBasicLogger('symtypes')

import bgcellmodels.common.units as pint_units


# class SymbolicQuantity(sympy.Symbol):
#     """
#     Wrapper for sympy.Symbol that also stores units.
#     """
#     def __new__(*args, **kwargs):
#         """
#         Symbol name is passed via __new__ so need to override.
#         """
#         cls = args[0]
#         return super(SymbolicQuantity, cls).__new__(*args, **kwargs)

#     def __init__(self, name, default_value=1.0, units=None):

#         super(SymbolicQuantity, self).__init__()
#         self._default_value = default_value

#         # Create Pint units
#         unit_str = units
#         self._pint_units = pint_units.Quantity(default_value, unit_str)
#         base_units = self._pint_units.to_base_units()
#         dims_exponents = base_units.dimensionality.items()
        
#         # Create SymPy units
#         if len(dims_exponents) == 0:
#             sp_dims = 1 # '1' is dimensionless in sympy
#         else:
#             dim0, exp0 = dims_exponents[0]
#             sp_dims = getattr(sympy_units.dimensions, dim0.strip('[]')) ** exp0
#         for dim_str, dim_exp in dims_exponents[1:]:
#             sp_dims *= getattr(sympy_units.dimensions, dim_str.strip('[]')) ** dim_exp
        
#         sp_scale = base_units.magnitude
#         self._sympy_units = sympy_units.Quantity(name, sp_dims, sp_scale)


def as_raw_expr(expr):
    """
    Recurse through expression three and rebuild raw Sympy expression.

    @see        http://docs.sympy.org/latest/tutorial/manipulation.html
    """
    if expr.is_Atom:
        return expr
    elif hasattr(expr, 'expr') and expr.expr.is_Atom:
        return expr.expr
    else:
        return expr.func(*[as_raw_expr(a) for a in expr.args])


class QuantitativeExpr(Expr):
    """
    Wrapper around sympy expressions so we can store custom attributes.
    """

    __slots__ = [
        'q', 'val', 'expr', 
        '_udict', 
        '_mech_attr_type',
        '_param_name',
        '_state_name',
        '_ion_species',
        'power'
    ]

    def __new__(cls, expr, *args, **kwargs):
        """
        Symbol name is passed via __new__ so need to override.

        @param  expr : <str/sympy.Expr/QuantitativeExpr>

        @param  units: QuantityType
                A quantity object indicating units

        @param  default_val : float/int
                Default value. Only
        """

        quant = kwargs.pop('units', 'dimensionless')
        default_val = kwargs.pop('default_val', None)

        if isinstance(quant, (str, unicode)):
            quant = pint_units.Quantity(1.0, quant)

        if isinstance(expr, (str, unicode)):
            expr = sympy.symbols(expr)
        elif default_val is not None:
            raise ValueError('Only provide default values for parameters (symbols)')
        
        obj = Expr.__new__(cls, expr, *args, **kwargs) # type QuantitativeExpr
        obj.q = quant
        obj.val = default_val
        obj.expr = expr.as_expr()

        if isinstance(expr, QuantitativeExpr):
            obj._udict = expr._udict.copy()
        else:
            obj._udict = {}
        
        return obj


    # def __init__(self, expr):
    #     raw_expr = expr.as_expr() # getattr(expr, 'expr', expr)
    #     self.expr = raw_expr # Original wrapped expression
        
    #     return super(QuantitativeExpr, self).__init__() # just calles object.__init__

    @property
    def units(self):
        return self.q
    
    @property
    def func(self):
        # prevent returning self.__class__
        return self.expr.func
    

    def as_expr(self):
        """
        Return as raw sympy Expression (for printing etc.)

        @note   A SymPy expression e is e.func(*e.args)

        @see    http://docs.sympy.org/latest/tutorial/manipulation.html#recursing-through-an-expression-tree
        """
        # expr = self.expr
        # args = expr.args
        # if len(args) == 0: # or: expr.is_Atom/is_Symbol
        #     args = (str(expr),) # for Symbol
        # else:
        #     args = [a.as_expr() for a in args]
        # # return expr.func(*[a.as_expr() for a in expr.args])
        # # TODO: make diagram to understand recursive call and see if overriding arithmethic operations is necessary to prevent hijacking (debug +/- operations)
        # return expr.func(*args)
        return as_raw_expr(self)


    # def __repr__(self):
    #     return repr(self.as_expr())
    #     # in sympy.core.basic.py :
    #     # from sympy.printing import sstr
    #     # return sstr(self.expr, order=None)

    # def __str__(self):
    #     return str(self.as_expr())



# class QuantExpr(Expr):
#     """
#     Wrapper around sympy expressions.
#     """
#     __slots__ = ['q', '_udict'] # Quantity and user data
#     raise_units_exceptions = True

#     def __new__(cls, expr, quant):
#         """
#         Symbol name is passed via __new__ so need to override.

#         TODO: make sure the top-level atoms in tree are raw symbol/expression
#         so that as_expr works. E.g. when creating parameters, they need to be passed
#         as a raw symbol with units.
#         """
#         if isinstance(expr, (str, unicode)):
#             expr = sympy.symbols(expr)
        
#         obj = sympy.Expr.__new__(cls, expr)
#         obj.q = quant
#         if isinstance(expr, QuantExpr):
#             obj._udict = expr._udict.copy()
#         else:
#             obj._udict = {}
#         return obj


#     def as_expr(self):
#         """
#         Return as sympy Expression

#         @note   A SymPy expression e is e.func(*e.args)

#         @see    http://docs.sympy.org/latest/tutorial/manipulation.html#recursing-through-an-expression-tree
#         """
#         expr = self
#         return expr.func(*[a.as_expr() for a in expr.args])


#     ############################################################################
#     # All mathematical operators should also check units
#     ############################################################################

#     # def __abs__(self):
#     #     return abs(self.expr)

#     # def __neg__(self):
#     #     return -self.expr

#     def __div__(self, x):
#         """ Division operator """
#         if isinstance(x, QuantExpr):
#             return QuantExpr(self.as_expr() / x.as_expr(), self.q / x.q)
#         elif isinstance(x, QuantityType):
#             return QuantExpr(self.as_expr() / x.magnitude, self.q / (1.0 * x.units))
#         else: # x has no quantity
#             msg = '{} divided by unitless {}: assuming dimensionless'.format(self, x)
#             if self.raise_units_exceptions:
#                 raise ValueError(msg)
#             else:
#                 logger.warn(msg)
#                 return QuantExpr(self.as_expr() / x, self.q)


#     def __rdiv__(self, x):
#         """ Reverse division """
#         if isinstance(x, QuantExpr):
#             return QuantExpr(x.as_expr() / self.as_expr(), x.q / self.q)
#         elif isinstance(x, QuantityType):
#             return QuantExpr(x.magnitude / self.as_expr(), (1.0 * x.units) / self.q)
#         else:
#             # no quanity associated with other
#             logger.warn('unitless {} divided by {}: assuming dimensionless'.format(x, self))
#             return QuantExpr(x / self.as_expr(), self.q)


#     def __truediv__(self, x):
#         """ Float division """
#         return self.__div__(x)


#     def __rtruediv__(self, x):
#         """ Reverse float division """
#         return self.__rdiv__(x)


#     def __mul__(self, x):
#         if isinstance(x, QuantExpr):
#             return QuantExpr(self.as_expr() * x.as_expr(), self.q * x.q)
#         elif isinstance(x, QuantityType):
#             return QuantExpr(self.as_expr() * x.magnitude, self.q * (1.0 * x.units))
#         else:
#             # no quantity associated with other
#             logger.warn('{} multiplied by unitless {}: assuming dimensionless'.format(self, x))
#             return QuantExpr(self.as_expr() * x, self.q)


#     def __rmul__(self, x):
#         if isinstance(x, QuantExpr):
#             return QuantExpr(x.as_expr() * self.as_expr(), x.q * self.q)
#         elif isinstance(x, QuantityType):
#             return QuantExpr(x.magnitude * self.as_expr(), (1.0 * x.units) * self.q)
#         else:
#             # no quanity associated with other
#             logger.warn('unitless {} mutliplied by {}: assuming dimensionless'.format(x, self))
#             return QuantExpr(x * self.as_expr(), self.q)
    

#     def __add__(self, x):
#         # convert to our own units
#         if isinstance(x, QuantExpr):
#             scaled_q = x.q.to(self.q)
#             scale = scaled_q.m / x.q.m
#             return QuantExpr(self.as_expr() + (scale*x.as_expr()), self.q)
#         elif isinstance(x, QuantityType):
#             scaled_q = x.to(self.q)
#             scale = scaled_q.m / x.m
#             return QuantExpr(self.as_expr() + scaled_q.m, self.q)
#         else:
#             # no quantity associated with other
#             logger.warn('{} multiplied by unitless {}: assuming same units'.format(self, x))
#             return QuantExpr(self.as_expr() + x, self.q)

    # def __radd__(self, x):
    #     TODO

    # def __sub__(self, x):
    #     TODO

    # def __rsub__(self, x):
    #     TODO

    # def __pow__(self, x):
    #     return self.expr ** getattr(x, 'expr', x)

    # def __or__(self, x):
    #     return self.expr or getattr(x, 'expr', x)

    # def __eq__(self, x):
    #     return self.expr == getattr(x, 'expr', x)

    # def __ne__(self, x):
    #     return self.expr != getattr(x, 'expr', x)

    # def __gt__(self, x):
    #     return self.expr > getattr(x, 'expr', x)

    # def __ge__(self, x):
    #     return self.expr >= getattr(x, 'expr', x)

    # def __lt__(self, x):
    #     return self.expr < getattr(x, 'expr', x)

    # def __le__(self, x):
    #     return self.expr <= getattr(x, 'expr', x)

"""
PyNN utilities.
"""

import numpy

def connection_plot(projection, positive='O', zero='.', empty=' ', spacer=''):
    """
    Fixed version of same method in pyNN/utility/__init__.py

    @return     utf_plot, matrix : tuple(unicode, np.array)
                Tuple containing the connection matrix in UTF-8 format and as
                raw connection weights.
    """
    connection_array = numpy.array(projection.get('weight', format='array',
                                gather='all', multiple_synapses='sum'))
    # connection_array = np.nan_to_num(connection_array) # also replaces inf and -inf
    nan_mask = numpy.isnan(connection_array)
    connection_array[nan_mask] = -1.0

    image = numpy.zeros_like(connection_array, dtype=unicode)
    # old_settings = numpy.seterr(invalid='ignore')  # ignore the complaint that x > 0 is invalid for NaN
    image[connection_array > 0] = positive
    image[connection_array == 0] = zero
    # numpy.seterr(**old_settings)  # restore original floating point error settings
    image[nan_mask] = empty
    utf_plot = '\n'.join([spacer.join(row) for row in image])
    return utf_plot, connection_array
"""
Utility functions for plotting using matplotlib.
"""
import os.path
from datetime import datetime

import matplotlib.pyplot as plt


def save_figure(fname, fig=None, dir=None, timestamp=False, **kwargs):
    """
    Save given or current figure.
    For LaTeX embedding, use extension pdf/pgf/eps in the figure name.

    kwargs: see https://matplotlib.org/api/_as_gen/matplotlib.pyplot.savefig.html
    """
    kwargs.setdefault('bbox_inches', 'tight') # prevents cropping
    kwargs.setdefault('transparent', True) # transparent background, see also 'frameon'
    if timestamp:
        fname += '_' + datetime.now().strftime('%Y.%m.%d-%H.%M.%S')
    fname += '.' + kwargs.setdefault('format', 'pdf')
    
    fig_filepath = os.path.join(dir, fname)
    if fig is None:
        plt.savefig(fig_filepath, **kwargs) # save current figure
    else:
        fig.savefig(fig_filepath, **kwargs) # save specific figure

    print("Figure saved to file {}".format(fig_filepath))
    return fig_filepath


def offset_show_twin_yax(ax, offset=1.2):
    """
    Having been created by twinx, ax has its frame off, so the line of its
    detached spine is invisible.  First, activate the frame but make the patch
    and spines invisible. Then make the detached spine visible.
    """
    ax.spines["right"].set_position(("axes", offset))
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

    ax.spines["right"].set_visible(True)


def hide_axis(ax, frame=True, x=True, y=True):
    """
    Hide axis lines and frame box of figure.
    """
    ax.set_frame_on(not frame)
    ax.get_xaxis().set_visible(not x)
    ax.get_yaxis().set_visible(not y)


def set_axes_size(w, h, ax=None):
    """
    Set the size of the axes (plotting area) instead of the whole figure
    which is the default.

    w, h: width, height in inches
    """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)


def get_style_colors():
    """
    Get colors of the currently active style.
    """
    return plt.rcParams['axes.prop_cycle'].by_key()['color']
# Copyright (C) 2016 Miriam Cabero Mueller
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""Module to generate figures with scatter plots and histograms.
"""

import numpy
import itertools
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot
import matplotlib.gridspec as gridspec
from pycbc.results import str_utils
#pyplot.rcParams.update({'text.usetex': True})

def create_corner_axes(parameters, labels=None):
    """Given a list of parameters, creates a figure with an axis for
    every possible combination of the parameters.

    Parameters
    ----------
    parameters : list
        Names of the variables to be plotted.
    labels : {None, list}, optional
        A list of names for the parameters.

    Returns
    -------
    fig : pyplot.figure
        The figure that was created.
    axes : array
        A 2D array of the axes that were created.
    axis_dict : dict
        A dictionary mapping the parameter combinations to the axis and their
        location in the subplots grid; i.e., the key, values are:
        `{('param1', 'param2'): (pyplot.axes, row index, column index)}`
    """
    if labels is None:
        labels = parameters
    elif len(labels) != len(parameters):
        raise ValueError("labels and parameters must be same length")
    # Create figure with adequate size for number of parameters.
    ndim = len(parameters)
    if ndim < 3:
        fsize = (8, 7)
    else:
        fsize = (ndim*3 - 1, ndim*3 - 2)
    fig = pyplot.figure(figsize=fsize)
    gs = gridspec.GridSpec(ndim, ndim)
    # create grid of axis numbers to easily create axes in the right locations
    axes = numpy.arange(ndim**2).reshape((ndim, ndim))

    # Select possible combinations of plots and establish rows and columns.
    combos =  list(itertools.combinations(range(ndim),2))
    # add the diagonals
    combos += [(ii, ii) for ii in range(ndim)]

    # create the mapping between parameter combos and axes
    axis_dict = {}
    # cycle over all the axes, setting thing as needed
    for nrow in range(ndim):
        for ncolumn in range(ndim):
            ax = pyplot.subplot(gs[axes[nrow, ncolumn]])
            # map to a parameter index
            px = ncolumn
            py = nrow
            if (px, py) in combos:
                axis_dict[parameters[px], parameters[py]] = (ax, nrow, ncolumn)
                # x labels only on bottom
                if nrow + 1 == ndim:
                    ax.set_xlabel('{}'.format(labels[px]))
                else:
                    pyplot.setp(ax.get_xticklabels(), visible=False)
                # y labels only on left and non-diagonal
                if ncolumn == 0 and nrow != 0:
                    ax.set_ylabel('{}'.format(labels[py]))
                else:
                    pyplot.setp(ax.get_yticklabels(), visible=False)
            else:
                # make non-used axes invisible
                ax.axis('off')
    return fig, axis_dict

def get_scale_fac(fig, fiducial_width=8, fiducial_height=7):
    """Gets a factor to scale fonts by for the given figure. The scale
    factor is relative to a figure with dimensions
    (`fiducial_width`, `fiducial_height`).
    """
    width, height = fig.get_size_inches()
    return (width*height/(fiducial_width*fiducial_height))**0.5

def scatter_histogram(parameters, data, zvals, labels=None, cbar_label=None,
                      vmin=None, vmax=None, mins=None, maxs=None,
                      cmap=None):
    """Generate a figure with several scatter plots and histograms.

    Parameters
    ----------
    parameters: list
        Names of the variables to be plotted.
    data: recarray or dict
        Record array or dictionary containing the data to be plotted in the x-
        and y-axes for all variables in parameters.
    zvals: array
        Data to be plotted on the z-axis (color for the scatter plots).
    labels: {None, list}, optional
        A list of names for the parameters.
    cbar_label : {None, str}
        Specify a label to add to the colorbar.
    vmin: {None, float}, optional
        Minimum value for the colorbar. If None, will use the minimum of zvals.
    vmax: {None, float}, optional
        Maximum value for the colorbar. If None, will use the maxmimum of zvals.
    mins: {None, dict}, optional
        Minimum value for the axis of each variable in parameters.
        If None, it will use the minimum of the corresponding variable in data.
    maxs: {None, dict}, optional
        Maximum value for the axis of each variable in parameters.
        If None, it will use the maximum of the corresponding variable in data.
    cmap : {None, pyplot.cmap}
        The color map to use for color points. If None, will use the
        matplotlib default.
    """

    if len(parameters) == 1:
        raise ValueError('You need at least two parameters.')
    if labels is None:
        labels = parameters

    # Sort zvals to get higher values on top in scatter plots
    sort_indices = zvals.argsort()
    zvals = zvals[sort_indices]
    data = data[sort_indices]

    # set up the figure with the corner axes
    fig, axis_dict = create_corner_axes(parameters, labels=labels)

    # Plot histograms in diagonal
    for param in parameters:
        ax, nrow, ncol = axis_dict[param, param]
        values = data[param]
        ax.hist(values, bins=50, color='navy', histtype='step')
        # 90th percentile
        values5 = numpy.percentile(values,5)
        values95 = numpy.percentile(values,95)
        ax.axvline(x=values5, ls='dashed')
        ax.axvline(x=values95, ls='dashed')
        # Median value
        valuesM = numpy.median(values)
        ax.axvline(x=valuesM, ls='dashed')
        negerror = valuesM - values5
        poserror = values95 - valuesM
        fmt = '$' + str_utils.format_value(valuesM, negerror,
              plus_error=poserror, ndecs=2) + '$'
        ax.set_title('{} = {}'.format(labels[nrow], fmt))
        # Remove y-ticks
        ax.set_yticks([])

    # Arguments for scatter plots
    if mins is None:
        mins = {p:data[p].min() for p in parameters}
    if maxs is None:
        maxs = {p:data[p].max() for p in parameters}

    # Fill figure with scatter plots
    for px, py in axis_dict:
        if px == py:
            continue
        ax, _, _ = axis_dict[px, py]
        plt = ax.scatter(x=data[px], y=data[py], c=zvals, s=5, 
                    edgecolors='none', vmin=vmin, vmax=vmax, cmap=cmap)
        ax.set_xlim(mins[px], maxs[px])
        ax.set_ylim(mins[py], maxs[py])

    # adjust tick number
    for px, py in axis_dict:
        ax, _, _ = axis_dict[px, py]
        # keep number of ticks <= 3 for ndim > 3
        if len(parameters) > 3:
            ax.set_xticks(reduce_ticks(ax, 'x', maxticks=3))
            ax.set_yticks(reduce_ticks(ax, 'y', maxticks=3))

    # compute font size based on fig size
    scale_fac = get_scale_fac(fig)
    
    fig.subplots_adjust(right=0.85, wspace=0.03)
    cbar_ax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    cb = fig.colorbar(plt, cax=cbar_ax)
    if cbar_label is not None:
        cb.set_label(cbar_label, fontsize=12*scale_fac)
    cb.ax.tick_params(labelsize=8*scale_fac)

    return fig, axis_dict

def reduce_ticks(ax, which, maxticks=3):
    """Given a pyplot axis, resamples its `which`-axis ticks such that are at most
    `maxticks` left.

    Parameters
    ----------
    ax : axis
        The axis to adjust.
    which : {'x' | 'y'}
        Which axis to adjust.
    maxticks : {3, int}
        Maximum number of ticks to use.

    Returns
    -------
    array
        An array of the selected ticks.
    """
    ticks = getattr(ax, 'get_{}ticks'.format(which))()
    if len(ticks) > maxticks:
        # make sure the left/right value is not at the edge
        minax, maxax = getattr(ax, 'get_{}lim'.format(which))()
        dw = abs(maxax-minax)/10.
        start_idx, end_idx = 0, len(ticks)
        if ticks[0] < minax + dw:
            start_idx += 1
        if ticks[-1] > maxax - dw:
            end_idx -= 1
        # get reduction factor
        fac = int(len(ticks) / maxticks)
        ticks = ticks[start_idx:end_idx:fac]
    return ticks

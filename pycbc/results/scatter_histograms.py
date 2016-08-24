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
from pycbc.results import str_utils

def scatter_histogram(parameters, args, data, labels=None, 
                        vmin=None, vmax=None, mins=None, maxs=None):
    """Generate a figure with several scatter plots and histograms.

    Parameters
    ----------
    parameters: list
        Names of the variables to be plotted.
    args: recarray or dict
        Record array or dictionary containing the data to be plotted in the x-
        and y-axes for all variables in parameters.
    data: array
        Data to be plotted on the z-axis (color for the scatter plots).
    labels: {None, list}, optional
        A list of names for the parameters.
    vmin: {None, float}, optional
        Minimum value for the colorbar. If None, will use the minimum of data.
    vmax: {None, float}, optional
        Maximum value for the colorbar. If None, will use the maxmimum of data.
    mins: {None, dict}, optional
        Minimum value for the axis of each variable in parameters.
        If None, it will use the minimum of the corresponding variable in args.
    maxs: {None, dict}, optional
        Maximum value for the axis of each variable in parameters.
        If None, it will use the maximum of the corresponding variable in args.
    """

    if len(parameters) == 1:
        raise ValueError('You need at least two parameters.')
    if labels is None:
        labels = parameters

    # Create figure with adequate size for number of parameters.
    ndim = len(parameters)
    if ndim < 3:
        fsize = (8, 7)
    else:
        fsize = (ndim*3 - 1, ndim*3 - 2)
    fig, axes = pyplot.subplots(ndim, ndim, sharex='col',
                                figsize=fsize)

    # Select possible combinations of plots and establish rows and columns.
    combos = list(itertools.combinations(parameters,2))
    columns = [combos[0][0]]
    for ii in range(1,len(combos)):
        if combos[ii][0] != combos[ii-1][0]:
            columns.append(combos[ii][0])
    rows = []
    for ii in range(len(combos)):
        if combos[ii][0] == columns[0]:
            rows.append(combos[ii][1])

    # Sort data to get higher values on top in scatter plots
    sort_indices = data.argsort()
    data = data[sort_indices]
    args = args[sort_indices]

    # Plot histograms in diagonal
    for nrow in range(len(rows)+1):
        ax = axes[nrow, nrow]
        values = args[parameters[nrow]]
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
        fmt = '$' + str_utils.format_value(valuesM, negerror, plus_error=poserror, ndecs=2) + '$'
        ax.set_title('{} = {}'.format(labels[nrow], fmt))
        # Remove y-ticks
        pyplot.setp(ax.get_yticklabels(), visible=False)
    # Add x-label in last histogram
    ax.set_xlabel('{}'.format(labels[-1]))

    # Arguments for scatter plots
    if mins is None:
        mins = {p:args[p].min() for p in parameters}
    if maxs is None:
        maxs = {p:args[p].max() for p in parameters}
    if vmin is None:
        vmin = data.min()
    if vmax is None:
        vmax = data.max()

    # Fill figure with scatter plots
    for nrow, prow in enumerate(rows):
        for ncolumn, pcolumn in enumerate(columns):
            if (pcolumn, prow) in combos:
                ax = axes[nrow+1, ncolumn]
                plt = ax.scatter(x=args[pcolumn], y=args[prow], c=data, s=1, 
                            vmin=vmin, vmax=vmax)
                ax.set_xlim(mins[pcolumn], maxs[pcolumn])
                ax.set_ylim(mins[prow], maxs[prow])
                # Labels only on bottom and left plots
                if ncolumn == 0:
                    ax.set_ylabel('{}'.format(labels[nrow]))
                else:
                    pyplot.setp(ax.get_yticklabels(), visible=False)
                if nrow + 1 == len(rows):
                    ax.set_xlabel('{}'.format(labels[ncolumn]))
            # Make empty plots white
            else:
                ax = axes[nrow, ncolumn]
                ax.axis('off')
            # Including plots in last column
            ax = axes[nrow, len(columns)]
            ax.axis('off')

    fig.subplots_adjust(right=0.85, wspace=0.03)
    cbar_ax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    fig.colorbar(plt, cax=cbar_ax)

    return fig

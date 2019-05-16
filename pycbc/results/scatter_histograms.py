# Copyright (C) 2016 Miriam Cabero Mueller, Collin Capano
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
"""
Module to generate figures with scatter plots and histograms.
"""

import itertools
import sys

import numpy

import scipy.stats

import matplotlib

# Only if a backend is not already set ... This should really *not* be done
# here, but in the executables you should set matplotlib.use()
# This matches the check that matplotlib does internally, but this *may* be
# version dependenant. If this is a problem then remove this and control from
# the executables directly.
if 'matplotlib.backends' not in sys.modules:  # nopep8
    matplotlib.use('agg')

from matplotlib import (offsetbox, pyplot, gridspec)

from pycbc.results import str_utils
from pycbc.io import FieldArray


def create_axes_grid(parameters, labels=None, height_ratios=None,
                     width_ratios=None, no_diagonals=False):
    """Given a list of parameters, creates a figure with an axis for
    every possible combination of the parameters.

    Parameters
    ----------
    parameters : list
        Names of the variables to be plotted.
    labels : {None, dict}, optional
        A dictionary of parameters -> parameter labels.
    height_ratios : {None, list}, optional
        Set the height ratios of the axes; see `matplotlib.gridspec.GridSpec`
        for details.
    width_ratios : {None, list}, optional
        Set the width ratios of the axes; see `matplotlib.gridspec.GridSpec`
        for details.
    no_diagonals : {False, bool}, optional
        Do not produce axes for the same parameter on both axes.

    Returns
    -------
    fig : pyplot.figure
        The figure that was created.
    axis_dict : dict
        A dictionary mapping the parameter combinations to the axis and their
        location in the subplots grid; i.e., the key, values are:
        `{('param1', 'param2'): (pyplot.axes, row index, column index)}`
    """
    if labels is None:
        labels = {p: p for p in parameters}
    elif any(p not in labels for p in parameters):
        raise ValueError("labels must be provided for all parameters")
    # Create figure with adequate size for number of parameters.
    ndim = len(parameters)
    if no_diagonals:
        ndim -= 1
    if ndim < 3:
        fsize = (8, 7)
    else:
        fsize = (ndim*3 - 1, ndim*3 - 2)
    fig = pyplot.figure(figsize=fsize)
    # create the axis grid
    gs = gridspec.GridSpec(ndim, ndim, width_ratios=width_ratios,
                           height_ratios=height_ratios,
                           wspace=0.05, hspace=0.05)
    # create grid of axis numbers to easily create axes in the right locations
    axes = numpy.arange(ndim**2).reshape((ndim, ndim))

    # Select possible combinations of plots and establish rows and columns.
    combos = list(itertools.combinations(parameters, 2))
    # add the diagonals
    if not no_diagonals:
        combos += [(p, p) for p in parameters]

    # create the mapping between parameter combos and axes
    axis_dict = {}
    # cycle over all the axes, setting thing as needed
    for nrow in range(ndim):
        for ncolumn in range(ndim):
            ax = pyplot.subplot(gs[axes[nrow, ncolumn]])
            # map to a parameter index
            px = parameters[ncolumn]
            if no_diagonals:
                py = parameters[nrow+1]
            else:
                py = parameters[nrow]
            if (px, py) in combos:
                axis_dict[px, py] = (ax, nrow, ncolumn)
                # x labels only on bottom
                if nrow + 1 == ndim:
                    ax.set_xlabel('{}'.format(labels[px]), fontsize=18)
                else:
                    pyplot.setp(ax.get_xticklabels(), visible=False)
                # y labels only on left
                if ncolumn == 0:
                    ax.set_ylabel('{}'.format(labels[py]), fontsize=18)
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


def construct_kde(samples_array, use_kombine=False):
    """Constructs a KDE from the given samples.
    """
    if use_kombine:
        try:
            import kombine
        except ImportError:
            raise ImportError("kombine is not installed.")
    # construct the kde
    if use_kombine:
        kde = kombine.clustered_kde.KDE(samples_array)
    else:
        kde = scipy.stats.gaussian_kde(samples_array.T)
    return kde


def create_density_plot(xparam, yparam, samples, plot_density=True,
                        plot_contours=True, percentiles=None, cmap='viridis',
                        contour_color=None, xmin=None, xmax=None,
                        ymin=None, ymax=None, exclude_region=None,
                        fig=None, ax=None, use_kombine=False):
    """Computes and plots posterior density and confidence intervals using the
    given samples.

    Parameters
    ----------
    xparam : string
        The parameter to plot on the x-axis.
    yparam : string
        The parameter to plot on the y-axis.
    samples : dict, numpy structured array, or FieldArray
        The samples to plot.
    plot_density : {True, bool}
        Plot a color map of the density.
    plot_contours : {True, bool}
        Plot contours showing the n-th percentiles of the density.
    percentiles : {None, float or array}
        What percentile contours to draw. If None, will plot the 50th
        and 90th percentiles.
    cmap : {'viridis', string}
        The name of the colormap to use for the density plot.
    contour_color : {None, string}
        What color to make the contours. Default is white for density
        plots and black for other plots.
    xmin : {None, float}
        Minimum value to plot on x-axis.
    xmax : {None, float}
        Maximum value to plot on x-axis.
    ymin : {None, float}
        Minimum value to plot on y-axis.
    ymax : {None, float}
        Maximum value to plot on y-axis.
    exclue_region : {None, str}
        Exclude the specified region when plotting the density or contours.
        Must be a string in terms of `xparam` and `yparam` that is
        understandable by numpy's logical evaluation. For example, if
        `xparam = m_1` and `yparam = m_2`, and you want to exclude the region
        for which `m_2` is greater than `m_1`, then exclude region should be
        `'m_2 > m_1'`.
    fig : {None, pyplot.figure}
        Add the plot to the given figure. If None and ax is None, will create
        a new figure.
    ax : {None, pyplot.axes}
        Draw plot on the given axis. If None, will create a new axis from
        `fig`.
    use_kombine : {False, bool}
        Use kombine's KDE to calculate density. Otherwise, will use
        `scipy.stats.gaussian_kde.` Default is False.

    Returns
    -------
    fig : pyplot.figure
        The figure the plot was made on.
    ax : pyplot.axes
        The axes the plot was drawn on.
    """
    if percentiles is None:
        percentiles = numpy.array([50., 90.])
    percentiles = 100. - numpy.array(percentiles)
    percentiles.sort()

    if ax is None and fig is None:
        fig = pyplot.figure()
    if ax is None:
        ax = fig.add_subplot(111)

    # convert samples to array and construct kde
    xsamples = samples[xparam]
    ysamples = samples[yparam]
    arr = numpy.vstack((xsamples, ysamples)).T
    kde = construct_kde(arr, use_kombine=use_kombine)

    # construct grid to evaluate on
    if xmin is None:
        xmin = xsamples.min()
    if xmax is None:
        xmax = xsamples.max()
    if ymin is None:
        ymin = ysamples.min()
    if ymax is None:
        ymax = ysamples.max()
    npts = 100
    X, Y = numpy.mgrid[
        xmin:xmax:complex(0, npts),  # pylint:disable=invalid-slice-index
        ymin:ymax:complex(0, npts)]  # pylint:disable=invalid-slice-index
    pos = numpy.vstack([X.ravel(), Y.ravel()])
    if use_kombine:
        Z = numpy.exp(kde(pos.T).reshape(X.shape))
        draw = kde.draw
    else:
        Z = kde(pos).T.reshape(X.shape)
        draw = kde.resample

    if exclude_region is not None:
        # convert X,Y to a single FieldArray so we can use it's ability to
        # evaluate strings
        farr = FieldArray.from_kwargs(**{xparam: X, yparam: Y})
        Z[farr[exclude_region]] = 0.

    if plot_density:
        ax.imshow(numpy.rot90(Z), extent=[xmin, xmax, ymin, ymax],
                  aspect='auto', cmap=cmap, zorder=1)
        if contour_color is None:
            contour_color = 'w'

    if plot_contours:
        # compute the percentile values
        resamps = kde(draw(int(npts**2)))
        if use_kombine:
            resamps = numpy.exp(resamps)
        s = numpy.percentile(resamps, percentiles)
        if contour_color is None:
            contour_color = 'k'
        # make linewidths thicker if not plotting density for clarity
        if plot_density:
            lw = 1
        else:
            lw = 2
        ct = ax.contour(X, Y, Z, s, colors=contour_color, linewidths=lw,
                        zorder=3)
        # label contours
        lbls = ['{p}%'.format(p=int(p)) for p in (100. - percentiles)]
        fmt = dict(zip(ct.levels, lbls))
        fs = 12
        ax.clabel(ct, ct.levels, inline=True, fmt=fmt, fontsize=fs)

    return fig, ax


def create_marginalized_hist(ax, values, label, percentiles=None,
                             color='k', fillcolor='gray', linecolor='navy',
                             linestyle='-',
                             title=True, expected_value=None,
                             expected_color='red', rotated=False,
                             plot_min=None, plot_max=None):
    """Plots a 1D marginalized histogram of the given param from the given
    samples.

    Parameters
    ----------
    ax : pyplot.Axes
        The axes on which to draw the plot.
    values : array
        The parameter values to plot.
    label : str
        A label to use for the title.
    percentiles : {None, float or array}
        What percentiles to draw lines at. If None, will draw lines at
        `[5, 50, 95]` (i.e., the bounds on the upper 90th percentile and the
        median).
    color : {'k', string}
        What color to make the histogram; default is black.
    fillcolor : {'gray', string, or None}
        What color to fill the histogram with. Set to None to not fill the
        histogram. Default is 'gray'.
    linestyle : str, optional
        What line style to use for the histogram. Default is '-'.
    linecolor : {'navy', string}
        What color to use for the percentile lines. Default is 'navy'.
    title : bool, optional
        Add a title with a estimated value +/- uncertainty. The estimated value
        is the pecentile halfway between the max/min of ``percentiles``, while
        the uncertainty is given by the max/min of the ``percentiles``. If no
        percentiles are specified, defaults to quoting the median +/- 95/5
        percentiles.
    rotated : {False, bool}
        Plot the histogram on the y-axis instead of the x. Default is False.
    plot_min : {None, float}
        The minimum value to plot. If None, will default to whatever `pyplot`
        creates.
    plot_max : {None, float}
        The maximum value to plot. If None, will default to whatever `pyplot`
        creates.
    scalefac : {1., float}
        Factor to scale the default font sizes by. Default is 1 (no scaling).
    """
    if fillcolor is None:
        htype = 'step'
    else:
        htype = 'stepfilled'
    if rotated:
        orientation = 'horizontal'
    else:
        orientation = 'vertical'
    ax.hist(values, bins=50, histtype=htype, orientation=orientation,
            facecolor=fillcolor, edgecolor=color, ls=linestyle, lw=2,
            density=True)
    if percentiles is None:
        percentiles = [5., 50., 95.]
    if len(percentiles) > 0:
        plotp = numpy.percentile(values, percentiles)
    else:
        plotp = []
    for val in plotp:
        if rotated:
            ax.axhline(y=val, ls='dashed', color=linecolor, lw=2, zorder=3)
        else:
            ax.axvline(x=val, ls='dashed', color=linecolor, lw=2, zorder=3)
    # plot expected
    if expected_value is not None:
        if rotated:
            ax.axhline(expected_value, color=expected_color, lw=1.5, zorder=2)
        else:
            ax.axvline(expected_value, color=expected_color, lw=1.5, zorder=2)
    if title:
        if len(percentiles) > 0:
            minp = min(percentiles)
            maxp = max(percentiles)
            medp = (maxp + minp) / 2.
        else:
            minp = 5
            medp = 50
            maxp = 95
        values_min = numpy.percentile(values, minp)
        values_med = numpy.percentile(values, medp)
        values_max = numpy.percentile(values, maxp)
        negerror = values_med - values_min
        poserror = values_max - values_med
        fmt = '${0}$'.format(str_utils.format_value(
            values_med, negerror, plus_error=poserror))
        if rotated:
            ax.yaxis.set_label_position("right")

            # sets colored title for marginal histogram
            set_marginal_histogram_title(ax, fmt, color,
                                         label=label, rotated=rotated)

            # Remove x-ticks
            ax.set_xticks([])
            # turn off x-labels
            ax.set_xlabel('')
            # set limits
            ymin, ymax = ax.get_ylim()
            if plot_min is not None:
                ymin = plot_min
            if plot_max is not None:
                ymax = plot_max
            ax.set_ylim(ymin, ymax)

        else:

            # sets colored title for marginal histogram
            set_marginal_histogram_title(ax, fmt, color, label=label)

            # Remove y-ticks
            ax.set_yticks([])
            # turn off y-label
            ax.set_ylabel('')
            # set limits
            xmin, xmax = ax.get_xlim()
            if plot_min is not None:
                xmin = plot_min
            if plot_max is not None:
                xmax = plot_max
            ax.set_xlim(xmin, xmax)


def set_marginal_histogram_title(ax, fmt, color, label=None, rotated=False):
    """ Sets the title of the marginal histograms.

    Parameters
    ----------
    ax : Axes
        The `Axes` instance for the plot.
    fmt : str
        The string to add to the title.
    color : str
        The color of the text to add to the title.
    label : str
        If title does not exist, then include label at beginning of the string.
    rotated : bool
        If `True` then rotate the text 270 degrees for sideways title.
    """

    # get rotation angle of the title
    rotation = 270 if rotated else 0

    # get how much to displace title on axes
    xscale = 1.05 if rotated else 0.0
    if rotated:
        yscale = 1.0
    elif len(ax.get_figure().axes) > 1:
        yscale = 1.15
    else:
        yscale = 1.05

    # get class that packs text boxes vertical or horizonitally
    packer_class = offsetbox.VPacker if rotated else offsetbox.HPacker

    # if no title exists
    if not hasattr(ax, "title_boxes"):

        # create a text box
        title = "{} = {}".format(label, fmt)
        tbox1 = offsetbox.TextArea(
                   title,
                   textprops=dict(color=color, size=15, rotation=rotation,
                                  ha='left', va='bottom'))

        # save a list of text boxes as attribute for later
        ax.title_boxes = [tbox1]

        # pack text boxes
        ybox = packer_class(children=ax.title_boxes,
                            align="bottom", pad=0, sep=5)

    # else append existing title
    else:

        # delete old title
        ax.title_anchor.remove()

        # add new text box to list
        tbox1 = offsetbox.TextArea(
                   " {}".format(fmt),
                   textprops=dict(color=color, size=15, rotation=rotation,
                                  ha='left', va='bottom'))
        ax.title_boxes = ax.title_boxes + [tbox1]

        # pack text boxes
        ybox = packer_class(children=ax.title_boxes,
                            align="bottom", pad=0, sep=5)

    # add new title and keep reference to instance as an attribute
    anchored_ybox = offsetbox.AnchoredOffsetbox(
                      loc=2, child=ybox, pad=0.,
                      frameon=False, bbox_to_anchor=(xscale, yscale),
                      bbox_transform=ax.transAxes, borderpad=0.)
    ax.title_anchor = ax.add_artist(anchored_ybox)


def create_multidim_plot(parameters, samples, labels=None,
                         mins=None, maxs=None, expected_parameters=None,
                         expected_parameters_color='r',
                         plot_marginal=True, plot_scatter=True,
                         marginal_percentiles=None, contour_percentiles=None,
                         marginal_title=True, marginal_linestyle='-',
                         zvals=None, show_colorbar=True, cbar_label=None,
                         vmin=None, vmax=None, scatter_cmap='plasma',
                         plot_density=False, plot_contours=True,
                         density_cmap='viridis',
                         contour_color=None, hist_color='black',
                         line_color=None, fill_color='gray',
                         use_kombine=False, fig=None, axis_dict=None):
    """Generate a figure with several plots and histograms.

    Parameters
    ----------
    parameters: list
        Names of the variables to be plotted.
    samples : FieldArray
        A field array of the samples to plot.
    labels: dict, optional
        A dictionary mapping parameters to labels. If none provided, will just
        use the parameter strings as the labels.
    mins : {None, dict}, optional
        Minimum value for the axis of each variable in `parameters`.
        If None, it will use the minimum of the corresponding variable in
        `samples`.
    maxs : {None, dict}, optional
        Maximum value for the axis of each variable in `parameters`.
        If None, it will use the maximum of the corresponding variable in
        `samples`.
    expected_parameters : {None, dict}, optional
        Expected values of `parameters`, as a dictionary mapping parameter
        names -> values. A cross will be plotted at the location of the
        expected parameters on axes that plot any of the expected parameters.
    expected_parameters_color : {'r', string}, optional
        What color to make the expected parameters cross.
    plot_marginal : {True, bool}
        Plot the marginalized distribution on the diagonals. If False, the
        diagonal axes will be turned off.
    plot_scatter : {True, bool}
        Plot each sample point as a scatter plot.
    marginal_percentiles : {None, array}
        What percentiles to draw lines at on the 1D histograms.
        If None, will draw lines at `[5, 50, 95]` (i.e., the bounds on the
        upper 90th percentile and the median).
    marginal_title : bool, optional
        Add a title over the 1D marginal plots that gives an estimated value
        +/- uncertainty. The estimated value is the pecentile halfway between
        the max/min of ``maginal_percentiles``, while the uncertainty is given
        by the max/min of the ``marginal_percentiles. If no
        ``marginal_percentiles`` are specified, the median +/- 95/5 percentiles
        will be quoted.
    marginal_linestyle : str, optional
        What line style to use for the marginal histograms.
    contour_percentiles : {None, array}
        What percentile contours to draw on the scatter plots. If None,
        will plot the 50th and 90th percentiles.
    zvals : {None, array}
        An array to use for coloring the scatter plots. If None, scatter points
        will be the same color.
    show_colorbar : {True, bool}
        Show the colorbar of zvalues used for the scatter points. A ValueError
        will be raised if zvals is None and this is True.
    cbar_label : {None, str}
        Specify a label to add to the colorbar.
    vmin: {None, float}, optional
        Minimum value for the colorbar. If None, will use the minimum of zvals.
    vmax: {None, float}, optional
        Maximum value for the colorbar. If None, will use the maxmimum of
        zvals.
    scatter_cmap : {'plasma', string}
        The color map to use for the scatter points. Default is 'plasma'.
    plot_density : {False, bool}
        Plot the density of points as a color map.
    plot_contours : {True, bool}
        Draw contours showing the 50th and 90th percentile confidence regions.
    density_cmap : {'viridis', string}
        The color map to use for the density plot.
    contour_color : {None, string}
        The color to use for the contour lines. Defaults to white for
        density plots, navy for scatter plots without zvals, and black
        otherwise.
    use_kombine : {False, bool}
        Use kombine's KDE to calculate density. Otherwise, will use
        `scipy.stats.gaussian_kde.` Default is False.

    Returns
    -------
    fig : pyplot.figure
        The figure that was created.
    axis_dict : dict
        A dictionary mapping the parameter combinations to the axis and their
        location in the subplots grid; i.e., the key, values are:
        `{('param1', 'param2'): (pyplot.axes, row index, column index)}`
    """
    if labels is None:
        labels = {p: p for p in parameters}
    # set up the figure with a grid of axes
    # if only plotting 2 parameters, make the marginal plots smaller
    nparams = len(parameters)
    if nparams == 2:
        width_ratios = [3, 1]
        height_ratios = [1, 3]
    else:
        width_ratios = height_ratios = None

    # only plot scatter if more than one parameter
    plot_scatter = plot_scatter and nparams > 1

    # Sort zvals to get higher values on top in scatter plots
    if plot_scatter:
        if zvals is not None:
            sort_indices = zvals.argsort()
            zvals = zvals[sort_indices]
            samples = samples[sort_indices]
            if contour_color is None:
                contour_color = 'k'
        elif show_colorbar:
            raise ValueError("must provide z values to create a colorbar")
        else:
            # just make all scatter points same color
            zvals = 'gray'
            if plot_contours and contour_color is None:
                contour_color = 'navy'

    # convert samples to a dictionary to avoid re-computing derived parameters
    # every time they are needed
    samples = dict([[p, samples[p]] for p in parameters])

    # values for axis bounds
    if mins is None:
        mins = {p: samples[p].min() for p in parameters}
    else:
        # copy the dict
        mins = {p: val for p, val in mins.items()}
    if maxs is None:
        maxs = {p: samples[p].max() for p in parameters}
    else:
        # copy the dict
        maxs = {p: val for p, val in maxs.items()}

    # remove common offsets
    for pi, param in enumerate(parameters):
        values, offset = remove_common_offset(samples[param])
        if offset != 0:
            # we'll add the offset removed to the label
            labels[param] = '{} - {:d}'.format(labels[param], offset)
            samples[param] = values
            mins[param] = mins[param] - float(offset)
            maxs[param] = maxs[param] - float(offset)
        # also remove from expected parameters, if they were provided
        if expected_parameters is not None:
            try:
                expected_parameters[param] -= offset
            except KeyError:
                pass

    # create the axis grid
    if fig is None and axis_dict is None:
        fig, axis_dict = create_axes_grid(
            parameters, labels=labels,
            width_ratios=width_ratios, height_ratios=height_ratios,
            no_diagonals=not plot_marginal)

    # Diagonals...
    if plot_marginal:
        for pi, param in enumerate(parameters):
            ax, _, _ = axis_dict[param, param]
            # if only plotting 2 parameters and on the second parameter,
            # rotate the marginal plot
            rotated = nparams == 2 and pi == nparams-1
            # see if there are expected values
            if expected_parameters is not None:
                try:
                    expected_value = expected_parameters[param]
                except KeyError:
                    expected_value = None
            else:
                expected_value = None
            create_marginalized_hist(
                ax, samples[param], label=labels[param],
                color=hist_color, fillcolor=fill_color,
                linestyle=marginal_linestyle, linecolor=line_color,
                title=marginal_title, expected_value=expected_value,
                expected_color=expected_parameters_color,
                rotated=rotated, plot_min=mins[param], plot_max=maxs[param],
                percentiles=marginal_percentiles)

    # Off-diagonals...
    for px, py in axis_dict:
        if px == py:
            continue
        ax, _, _ = axis_dict[px, py]
        if plot_scatter:
            if plot_density:
                alpha = 0.3
            else:
                alpha = 1.
            plt = ax.scatter(x=samples[px], y=samples[py], c=zvals, s=5,
                             edgecolors='none', vmin=vmin, vmax=vmax,
                             cmap=scatter_cmap, alpha=alpha, zorder=2)

        if plot_contours or plot_density:
            # Exclude out-of-bound regions
            # this is a bit kludgy; should probably figure out a better
            # solution to eventually allow for more than just m_p m_s
            if (px == 'm_p' and py == 'm_s') or (py == 'm_p' and px == 'm_s'):
                exclude_region = 'm_s > m_p'
            else:
                exclude_region = None
            create_density_plot(
                px, py, samples, plot_density=plot_density,
                plot_contours=plot_contours, cmap=density_cmap,
                percentiles=contour_percentiles,
                contour_color=contour_color, xmin=mins[px], xmax=maxs[px],
                ymin=mins[py], ymax=maxs[py],
                exclude_region=exclude_region, ax=ax,
                use_kombine=use_kombine)

        if expected_parameters is not None:
            try:
                ax.axvline(expected_parameters[px], lw=1.5,
                           color=expected_parameters_color, zorder=5)
            except KeyError:
                pass
            try:
                ax.axhline(expected_parameters[py], lw=1.5,
                           color=expected_parameters_color, zorder=5)
            except KeyError:
                pass

        ax.set_xlim(mins[px], maxs[px])
        ax.set_ylim(mins[py], maxs[py])

    # adjust tick number for large number of plots
    if len(parameters) > 3:
        for px, py in axis_dict:
            ax, _, _ = axis_dict[px, py]
            ax.set_xticks(reduce_ticks(ax, 'x', maxticks=3))
            ax.set_yticks(reduce_ticks(ax, 'y', maxticks=3))

    if plot_scatter and show_colorbar:
        # compute font size based on fig size
        scale_fac = get_scale_fac(fig)
        fig.subplots_adjust(right=0.85, wspace=0.03)
        cbar_ax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
        cb = fig.colorbar(plt, cax=cbar_ax)
        if cbar_label is not None:
            cb.set_label(cbar_label, fontsize=12*scale_fac)
        cb.ax.tick_params(labelsize=8*scale_fac)

    return fig, axis_dict


def remove_common_offset(arr):
    """Given an array of data, removes a common offset > 1000, returning the
    removed value.
    """
    offset = 0
    isneg = (arr <= 0).all()
    # make sure all values have the same sign
    if isneg or (arr >= 0).all():
        # only remove offset if the minimum and maximum values are the same
        # order of magintude and > O(1000)
        minpwr = numpy.log10(abs(arr).min())
        maxpwr = numpy.log10(abs(arr).max())
        if numpy.floor(minpwr) == numpy.floor(maxpwr) and minpwr > 3:
            offset = numpy.floor(10**minpwr)
            if isneg:
                offset *= -1
            arr = arr - offset
    return arr, int(offset)


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

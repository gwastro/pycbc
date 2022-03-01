# Copyright (C) 2016 Miriam Cabero Mueller, Collin Capano
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
"""
This modules provides classes for evaluating arbitrary distributions from
a file.
"""

import h5py
import numpy
import scipy.stats
from pycbc.distributions import bounded
import pycbc.transforms

class Arbitrary(bounded.BoundedDist):
    r"""A distribution constructed from a set of parameter values using a kde.
    Bounds may be optionally provided to limit the range.

    Parameters
    ----------
    bounds : dict, optional
        Independent bounds on one or more parameters may be provided to limit
        the range of the kde.
    bandwidth : str, optional
        Set the bandwidth method for the KDE. See
        :py:func:`scipy.stats.gaussian_kde` for details. Default is "scott".
    \**params :
        The keyword arguments should provide the names of the parameters and
        a list of their parameter values. If multiple parameters are provided,
        a single kde will be produced with dimension equal to the number of
        parameters.
    """
    name = 'arbitrary'

    def __init__(self, bounds=None, bandwidth="scott", **kwargs):
        # initialize the bounds
        if bounds is None:
            bounds = {}
        bounds.update({p: None for p in kwargs if p not in bounds})
        super(Arbitrary, self).__init__(**bounds)
        # check that all parameters specified in bounds have samples
        if set(self.params) != set(kwargs.keys()):
            raise ValueError("Must provide samples for all parameters given "
                             "in the bounds dictionary")
        # if bounds are provided use logit transform to move the points
        # to +/- inifinity
        self._transforms = {}
        self._tparams = {}
        for param,bnds in self.bounds.items():
            if numpy.isfinite(bnds[1] - bnds[0]):
                tparam = 'logit'+param
                samples = kwargs[param]
                t = pycbc.transforms.Logit(param, tparam, domain=bnds)
                self._transforms[tparam] = t
                self._tparams[param] = tparam
                # remove any sample points that fall out side of the bounds
                outside = bnds.__contains__(samples)
                if outside.any():
                    samples = samples[outside]
                # transform the sample points
                kwargs[param] = t.transform({param: samples})[tparam]
            elif not (~numpy.isfinite(bnds[0]) and ~numpy.isfinite(bnds[1])):
                raise ValueError("if specifying bounds, both bounds must "
                                 "be finite")
        # build the kde
        self._kde = self.get_kde_from_arrays(*[kwargs[p] for p in self.params])
        self.set_bandwidth(bandwidth)

    @property
    def params(self):
        return self._params

    @property
    def kde(self):
        return self._kde

    def _pdf(self, **kwargs):
        """Returns the pdf at the given values. The keyword arguments must
        contain all of parameters in self's params. Unrecognized arguments are
        ignored.
        """
        for p in self._params:
            if p not in kwargs.keys():
                raise ValueError('Missing parameter {} to construct pdf.'
                                 .format(p))
        if kwargs in self:
            # transform into the kde space
            jacobian = 1.
            for param, tparam in self._tparams.items():
                t = self._transforms[tparam]
                try:
                    samples = t.transform({param: kwargs[param]})
                except ValueError as e:
                    # can get a value error if the value is exactly == to
                    # the bounds, in which case, just return 0.
                    if kwargs[param] in self.bounds[param]:
                        return 0.
                    else:
                        raise ValueError(e)
                kwargs[param] = samples[tparam]
                # update the jacobian for the transform; if p is the pdf
                # in the params frame (the one we want) and p' is the pdf
                # in the transformed frame (the one that's calculated) then:
                # p = J * p', where J is the Jacobian of going from p to p'
                jacobian *= t.jacobian(samples)
            # for scipy < 0.15.0, gaussian_kde.pdf = gaussian_kde.evaluate
            this_pdf = jacobian * self._kde.evaluate([kwargs[p]
                                                      for p in self._params])
            if len(this_pdf) == 1:
                return float(this_pdf)
            else:
                return this_pdf
        else:
            return 0.

    def _logpdf(self, **kwargs):
        """Returns the log of the pdf at the given values. The keyword
        arguments must contain all of parameters in self's params.
        Unrecognized arguments are ignored.
        """
        if kwargs not in self:
            return -numpy.inf
        else:
            return numpy.log(self._pdf(**kwargs))

    def set_bandwidth(self, set_bw="scott"):
        self._kde.set_bandwidth(set_bw)

    def rvs(self, size=1, param=None):
        """Gives a set of random values drawn from the kde.

        Parameters
        ----------
        size : {1, int}
            The number of values to generate; default is 1.
        param : {None, string}
            If provided, will just return values for the given parameter.
            Otherwise, returns random values for each parameter.

        Returns
        -------
        structured array
            The random values in a numpy structured array. If a param was
            specified, the array will only have an element corresponding to the
            given parameter. Otherwise, the array will have an element for each
            parameter in self's params.
        """
        if param is not None:
            dtype = [(param, float)]
        else:
            dtype = [(p, float) for p in self.params]
        size = int(size)
        arr = numpy.zeros(size, dtype=dtype)
        draws = self._kde.resample(size)
        draws = {param: draws[ii,:] for ii,param in enumerate(self.params)}
        for (param,_) in dtype:
            try:
                # transform back to param space
                tparam = self._tparams[param]
                tdraws = {tparam: draws[param]}
                draws[param] = self._transforms[tparam].inverse_transform(
                    tdraws)[param]
            except KeyError:
                pass
            arr[param] = draws[param]
        return arr

    @staticmethod
    def get_kde_from_arrays(*arrays):
        """Constructs a KDE from the given arrays.

        \*arrays :
            Each argument should be a 1D numpy array to construct the kde from.
            The resulting KDE will have dimension given by the number of
            parameters.
        """
        return scipy.stats.gaussian_kde(numpy.vstack(arrays))

    @classmethod
    def from_config(cls, cp, section, variable_args):
        """Raises a NotImplementedError; to load from a config file, use
        `FromFile`.
        """
        raise NotImplementedError("This class does not support loading from a "
                                  "config file. Use `FromFile` instead.")


class FromFile(Arbitrary):
    r"""A distribution that reads the values of the parameter(s) from an hdf
    file, computes the kde to construct the pdf, and draws random variables
    from it.

    Parameters
    ----------
    filename : str
        The path to an hdf file containing the values of the parameters that
        want to be used to construct the distribution. Each parameter should
        be a separate dataset in the hdf file, and all datasets should have
        the same size. For example, to give a prior for mass1 and mass2 from
        file f, f['mass1'] and f['mass2'] contain the n values for each
        parameter.
    datagroup : str, optional
        The name of the group to look in for the samples. For example, if
        ``datagroup = 'samples'``, then parameter ``param`` will be retrived
        from ``f['samples'][param]``. If none provided (the default) the data
        sets will be assumed to be in the top level directory of the file.
    \**params :
        The keyword arguments should provide the names of the parameters to be
        read from the file and (optionally) their bounds. If no parameters are
        provided, it will use all the parameters found in the file. To provide
        bounds, specify e.g. mass1=[10,100]. Otherwise, mass1=None.

    Attributes
    ----------
    norm : float
        The normalization of the multi-dimensional pdf.
    lognorm : float
        The log of the normalization.
    kde :
        The kde obtained from the values in the file.
    """
    name = 'fromfile'
    def __init__(self, filename=None, datagroup=None, **params):
        if filename is None:
            raise ValueError('A file must be specified for this distribution.')
        self._filename = filename
        self.datagroup = datagroup
        # Get the parameter names to pass to get_kde_from_file
        if len(params) == 0:
            ps = None
        else:
            ps = list(params.keys())
        param_vals, bw = self.get_arrays_from_file(filename, params=ps)
        super(FromFile, self).__init__(bounds=params, bandwidth=bw,
                                       **param_vals)

    @property
    def filename(self):
        """str: The path to the file containing values for the parameter(s).
        """
        return self._filename

    def get_arrays_from_file(self, params_file, params=None):
        """Reads the values of one or more parameters from an hdf file and
        returns as a dictionary.

        Parameters
        ----------
        params_file : str
            The hdf file that contains the values of the parameters.
        params : {None, list}
            If provided, will just retrieve the given parameter names.

        Returns
        -------
        dict
            A dictionary of the parameters mapping `param_name -> array`.
        """
        try:
            f = h5py.File(params_file, 'r')
        except:
            raise ValueError('File not found.')
        if self.datagroup is not None:
            get = f[self.datagroup]
        else:
            get = f
        if params is not None:
            if not isinstance(params, list):
                params = [params]
            for p in params:
                if p not in get.keys():
                    raise ValueError('Parameter {} is not in {}'
                                     .format(p, params_file))
        else:
            params = [str(k) for k in get.keys()]
        params_values = {p: get[p][()] for p in params}
        try:
            bandwidth = f.attrs["bandwidth"]
        except KeyError:
            bandwidth = "scott"

        f.close()
        return params_values, bandwidth

    @classmethod
    def from_config(cls, cp, section, variable_args):
        """Returns a distribution based on a configuration file.

        The parameters
        for the distribution are retrieved from the section titled
        "[`section`-`variable_args`]" in the config file.

        The file to construct the distribution from must be provided by setting
        `filename`. Boundary arguments can be provided in the same way as
        described in `get_param_bounds_from_config`.

        .. code-block:: ini

            [{section}-{tag}]
            name = fromfile
            filename = ra_prior.hdf
            min-ra = 0
            max-ra = 6.28

        Parameters
        ----------
        cp : pycbc.workflow.WorkflowConfigParser
            A parsed configuration file that contains the distribution
            options.
        section : str
            Name of the section in the configuration file.
        variable_args : str
            The names of the parameters for this distribution, separated by
            `prior.VARARGS_DELIM`. These must appear in the "tag" part
            of the section header.

        Returns
        -------
        BoundedDist
            A distribution instance from the pycbc.inference.prior module.
        """
        return bounded.bounded_from_config(cls, cp, section, variable_args,
                                           bounds_required=False)

__all__ = ['Arbitrary', 'FromFile']

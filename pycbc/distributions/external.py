# Copyright (C) 2020 Alexander Nitz, 2022 Shichao Wu
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
This modules provides classes for evaluating PDF, logPDF, CDF and inverse CDF
from external arbitrary distributions, and drawing samples from them.
"""
import logging
import importlib
import numpy as np

import scipy.integrate as scipy_integrate
import scipy.interpolate as scipy_interpolate

from pycbc import VARARGS_DELIM

logger = logging.getLogger('pycbc.distributions.external')


class External(object):
    """ Distribution defined by external cdfinv and logpdf functions

    To add to an inference configuration file:

    .. code-block:: ini

        [prior-param1+param2]
        name = external
        module = custom_mod
        logpdf = custom_function_name
        cdfinv = custom_function_name2

    Parameters
    ----------
    params : list
        list of parameter names
    custom_mod : module
        module from which logpdf and cdfinv functions can be imported
    logpdf : function
        function which returns the logpdf
    cdfinv : function
        function which applies the invcdf

    Examples
    --------
    To instantate by hand and example of function format. You must provide
    the logpdf function, and you may either provide the rvs or cdfinv function.
    If the cdfinv is provided, but not the rvs, the random values will
    be calculated using the cdfinv function.

    >>> import numpy
    >>> params = ['x', 'y']
    >>> def logpdf(x=None, y=None):
    ...     p = numpy.ones(len(x))
    ...     return p
    >>>
    >>> def cdfinv(**kwds):
    ...     return kwds
    >>> e = External(['x', 'y'], logpdf, cdfinv=cdfinv)
    >>> e.rvs(size=10)
    """
    name = "external"

    def __init__(self, params=None, logpdf=None,
                 rvs=None, cdfinv=None, **kwds):
        self.params = params
        self.logpdf = logpdf
        self.cdfinv = cdfinv
        self._rvs = rvs

        if not (rvs or cdfinv):
            raise ValueError("Must provide either rvs or cdfinv")

    def rvs(self, size=1, **kwds):
        "Draw random value"
        if self._rvs:
            return self._rvs(size=size)
        samples = {param: np.random.uniform(0, 1, size=size)
                   for param in self.params}
        return self.cdfinv(**samples)

    def apply_boundary_conditions(self, **params):
        return params

    def __call__(self, **kwds):
        return self.logpdf(**kwds)

    @classmethod
    def from_config(cls, cp, section, variable_args):
        tag = variable_args
        params = variable_args.split(VARARGS_DELIM)
        modulestr = cp.get_opt_tag(section, 'module', tag)
        mod = importlib.import_module(modulestr)

        logpdfstr = cp.get_opt_tag(section, 'logpdf', tag)
        logpdf = getattr(mod, logpdfstr)

        cdfinv = rvs = None
        if cp.has_option_tag(section, 'cdfinv', tag):
            cdfinvstr = cp.get_opt_tag(section, 'cdfinv', tag)
            cdfinv = getattr(mod, cdfinvstr)

        if cp.has_option_tag(section, 'rvs', tag):
            rvsstr = cp.get_opt_tag(section, 'rvs', tag)
            rvs = getattr(mod, rvsstr)

        return cls(params=params, logpdf=logpdf, rvs=rvs, cdfinv=cdfinv)


class DistributionFunctionFromFile(External):
    r"""Evaluating PDF, logPDF, CDF and inverse CDF from the external
        density function.

    To add to an inference configuration file:

    .. code-block:: ini

        [prior-param1]
        name = external_func_fromfile
        file_path = spin.txt
        column_index = 1

    Parameters
    ----------
    params : list
        list of parameter names
    file_path: str
        The path of the external density function's .txt file.
    column_index: int
        The column index of the density distribution. By default, the first
        should be the values of a certain parameter, such as "mass", other
        columns should be the corresponding density values (as a function of
        that parameter). If you add the name of the parameter in the first
        row, please add the '#' at the beginning.
    \**kwargs :
        All other keyword args are passed to `scipy.integrate.quad` to control
        the numerical accuracy of the inverse CDF.
        If not be provided, will use the default values in `self.__init__`.

    Notes
    -----
    This class is different from `pycbc.distributions.arbitrary.FromFile`,
    which needs samples from the hdf file to construct the PDF by using KDE.
    This class reads in any continuous functions of the parameter.
    """
    name = "external_func_fromfile"

    def __init__(self, params=None, file_path=None,
                 column_index=None, **kwargs):
        super().__init__(cdfinv=self._cdfinv, logpdf=self.logpdf)
        self.params = params
        self.data = np.loadtxt(fname=file_path, unpack=True, comments='#')
        self.column_index = int(column_index)
        self.epsabs = kwargs.get('epsabs', 1.49e-05)
        self.epsrel = kwargs.get('epsrel', 1.49e-05)
        self.x_list = np.linspace(self.data[0][0], self.data[0][-1], 1000)
        self.interp = {'pdf': callable, 'cdf': callable, 'cdfinv': callable}
        if not file_path:
            raise ValueError("Must provide the path to density function file.")

    def logpdf(self, **kwargs):
        x = kwargs.pop(self.params[0])
        return self._logpdf(x, **kwargs)

    def _pdf(self, x010, **kwargs):
        """Calculate and interpolate the PDF by using the given density
        function, then return the corresponding value at the given x."""
        if self.interp['pdf'] == callable:
            func_unnorm = scipy_interpolate.interp1d(
                self.data[0], self.data[self.column_index])
            norm_const = scipy_integrate.quad(
                func_unnorm, self.data[0][0], self.data[0][-1],
                epsabs=self.epsabs, epsrel=self.epsrel, limit=500,
                **kwargs)[0]
            self.interp['pdf'] = scipy_interpolate.interp1d(
                self.data[0], self.data[self.column_index]/norm_const,
                bounds_error=False, fill_value=0)
        pdf_val = np.float64(self.interp['pdf'](x010))
        return pdf_val

    def _logpdf(self, x010, **kwargs):
        """Calculate the logPDF by calling `pdf` function."""
        z = np.log(self._pdf(x010, **kwargs))
        return z

    def _cdf(self, x, **kwargs):
        """Calculate and interpolate the CDF, then return the corresponding
        value at the given x."""
        if self.interp['cdf'] == callable:
            cdf_list = []
            for x_val in self.x_list:
                cdf_x = scipy_integrate.quad(
                    self._pdf, self.data[0][0], x_val, epsabs=self.epsabs,
                    epsrel=self.epsrel, limit=500, **kwargs)[0]
                cdf_list.append(cdf_x)
            self.interp['cdf'] = \
                scipy_interpolate.interp1d(self.x_list, cdf_list)
        cdf_val = np.float64(self.interp['cdf'](x))
        return cdf_val

    def _cdfinv(self, **kwargs):
        """Calculate and interpolate the inverse CDF, then return the
        corresponding parameter value at the given CDF value."""
        if self.interp['cdfinv'] == callable:
            cdf_list = []
            for x_value in self.x_list:
                cdf_list.append(self._cdf(x_value))
            self.interp['cdfinv'] = \
                scipy_interpolate.interp1d(cdf_list, self.x_list)
        cdfinv_val = {self.params[0]: np.float64(
            self.interp['cdfinv'](kwargs[self.params[0]]))}
        return cdfinv_val

    @classmethod
    def from_config(cls, cp, section, variable_args):
        tag = variable_args
        params = variable_args.split(VARARGS_DELIM)
        file_path = cp.get_opt_tag(section, 'file_path', tag)
        column_index = cp.get_opt_tag(section, 'column_index', tag)
        return cls(params=params, file_path=file_path,
                   column_index=column_index)


__all__ = ['External', 'DistributionFunctionFromFile']

import numpy
import scipy.stats
import h5py
from ConfigParser import Error
import warnings
from pycbc.distributions import bounded
from pycbc.inference import boundaries

class FromFile(bounded._BoundedDist):
    """A distribution that reads the values of the parameter(s) from an hdf
    file, computes the kde to construct the pdf, and draws random variables
    from it.

    Parameters
    ----------
    file_name : str
        The path to an hdf file containing the values of the parameters that
        want to be used to construct the distribution. Each parameter should
        be a separate dataset in the hdf file, and all datasets should have
        the same size. For example, to give a prior for mass1 and mass2 from
        file f, f['mass1'] and f['mass2'] contain the n values for each
        parameter.
    \**params :
        The keyword arguments should provide the names of the parameters to be
        read from the file and (optionally) their bounds. If no parameters are
        provided, it will use all the parameters found in the file. To provide
        bounds, specify e.g. mass1=[10,100]. Otherwise, mass1=None.

    Attributes
    ----------
    name : 'fromfile'
        The name of the distribution.
    file_name : str
        The path to the file containing values for the parameter(s).
    params : list
        Parameters read from file.
    norm : float
        The normalization of the multi-dimensional pdf.
    lognorm : float
        The log of the normalization.
    kde :
        The kde obtained from the values in the file.
    """
    name = 'fromfile'
    def __init__(self, file_name=None, **params):
        if file_name is None:
            raise ValueError('A file must be specified for this distribution.')
        self._filename = file_name
        # Get the parameter names to pass to get_kde_from_file
        if len(params) == 0:
            ps = None
        else:
            ps = params.keys()
        pnames, self._kde = self.get_kde_from_file(file_name, params=ps)
        # If no parameters where given, populate with pnames
        for param in pnames:
            if param not in params:
                params[param] = None
        super(FromFile, self).__init__(**params)
        # Make sure to store parameter names in same order as given by kde function
        self._params = pnames
        # Compute the norm and save
        lower_bounds = [self.bounds[p][0] for p in pnames]
        higher_bounds = [self.bounds[p][1] for p in pnames]
        # Avoid inf because of inconsistencies in integrate_box
        RANGE_LIMIT = 2 ** 31
        for ii, bnd in enumerate(lower_bounds):
            if abs(bnd) == numpy.inf:
                lower_bounds[ii] = numpy.sign(bnd) * RANGE_LIMIT
        for ii, bnd in enumerate(higher_bounds):
            if abs(bnd) == numpy.inf:
                higher_bounds[ii] = numpy.sign(bnd) * RANGE_LIMIT
        # Array of -inf for the lower limits in integrate_box
        lower_limits = - RANGE_LIMIT * numpy.ones(shape=len(lower_bounds))
        # CDF(-inf,b) - CDF(-inf, a)
        invnorm = self._kde.integrate_box(lower_limits, higher_bounds) - \
                    self._kde.integrate_box(lower_limits, lower_bounds)
        self._norm = 1. / invnorm
        self._lognorm = numpy.log(self._norm)

    @property
    def file_name(self):
        return self._filename

    @property
    def params(self):
        return self._params

    @property
    def norm(self):
        return self._norm

    @property
    def lognorm(self):
        return self._lognorm

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
                raise ValueError('Missing parameter {} to construct pdf.'.format(p))
        if kwargs in self:
            # for scipy < 0.15.0, gaussian_kde.pdf = gaussian_kde.evaluate
            this_pdf = self._norm * self._kde.evaluate([kwargs[p]
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
        for p in self._params:
            if p not in kwargs.keys():
                raise ValueError('Missing parameter {} to construct pdf.'.format(p))
        if kwargs in self:
            # for scipy < 0.15.0,
            # gaussian_kde.logpdf = numpy.log(gaussian_kde.evaluate)
            this_logpdf = self._lognorm + \
                          numpy.log(self._kde.evaluate([kwargs[p]
                                                     for p in self._params]))
            if len(this_logpdf) == 1:
                return float(this_logpdf)
            else:
                return this_logpdf
        else:
            return -numpy.inf

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
        arr = numpy.zeros(size, dtype=dtype)
        randoms = self._kde.resample(size)
        for order, param in enumerate(dtype):
            arr[param[0]] = randoms[order]
        return arr

    @staticmethod
    def get_kde_from_file(params_file, params=None):
        """Reads the values of one or more parameters from an hdf file and
        computes the kernel density estimate (kde).

        Parameters
        ----------
        params_file : str
            The hdf file that contains the values of the parameters.
        params : {None, list}
            If provided, will just use the values for the given parameter.
            Otherwise, uses the values for each parameter in the file.
        Returns
        -------
        values
            Array with the values of the parameters.
        kde
            The kde from the parameters.
        """
        try:
            f = h5py.File(params_file, 'r')
        except:
            raise ValueError('File not found.')
        if params is not None:
            if not isinstance(params, list):
                params = [params]
            for p in params:
                if p not in f.keys():
                    raise ValueError('Parameter {} is not in {}'.format(p, params_file))
        else:
            params = [str(k) for k in f.keys()]
        params_values = {p:f[p][:] for p in params}
        f.close()
        values = numpy.vstack((params_values[p] for p in params))
        return params, scipy.stats.gaussian_kde(values)

    @classmethod
    def from_config(cls, cp, section, variable_args):
        """Returns a distribution based on a configuration file.

        The parameters
        for the distribution are retrieved from the section titled
        "[`section`-`variable_args`]" in the config file.

        The file to construct the distribution from must be provided by setting
        `file_name`. Boundary arguments can be provided in the same way as
        described in `get_param_bounds_from_config`.

        .. code-block:: ini

            [{section}-{tag}]
            name = fromfile
            file_name = ra_prior.hdf
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
        Uniform
            A distribution instance from the pycbc.inference.prior module.
        """
        return super(FromFile, cls).from_config(cp, section, variable_args,
                                                bounds_required=False)


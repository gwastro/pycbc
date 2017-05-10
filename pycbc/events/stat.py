# Copyright (C) 2016 Alex Nitz
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
""" This modules contains functions for calculating coincident ranking
statistic values
"""

import numpy
from . import events


def get_newsnr(trigs):
    """
    Calculate newsnr ('reweighted SNR') for a trigs object

    Parameters
    ----------
    trigs: dict of numpy.ndarrays
        Dictionary holding single detector trigger information.
    'chisq_dof', 'snr', and 'chisq' are required keys

    Returns
    -------
    numpy.ndarray
        Array of newsnr values
    """
    dof = 2. * trigs['chisq_dof'] - 2.
    newsnr = events.newsnr(trigs['snr'], trigs['chisq'] / dof)
    return numpy.array(newsnr, ndmin=1, dtype=numpy.float32)

def get_newsnr_sgveto(trigs):
    """
    Calculate newsnr re-weigthed by the sine-gaussian veto

    Parameters
    ----------
    trigs: dict of numpy.ndarrays
        Dictionary holding single detector trigger information.
    'chisq_dof', 'snr', and 'chisq' are required keys

    Returns
    -------
    numpy.ndarray
        Array of newsnr values
    """
    dof = 2. * trigs['chisq_dof'] - 2.
    nsnr_sg = events.newsnr_sgveto(trigs['snr'], trigs['chisq'] / dof, trigs['sg_chisq'])
    return numpy.array(nsnr_sg, ndmin=1, dtype=numpy.float32)


class Stat(object):

    """ Base class which should be extended to provide a coincident statistic"""
    def __init__(self, files):
        """Create a statistic class instance

        Parameters
        ----------
        files: list of strs
            A list containing the filenames of hdf format files used to help
        construct the coincident statistics. The files must have a 'stat'
        attribute which is used to associate them with the appropriate
        statistic class.
        """
        import h5py
        self.files = {}
        for filename in files:
            f = h5py.File(filename, 'r')
            stat = f.attrs['stat']
            self.files[stat] = f

        # Provide the dtype of the single detector method's output
        # This is used by background estimation codes that need to maintain
        # a buffer of such values.
        self.single_dtype = numpy.float32


class NewSNRStatistic(Stat):

    """ Calculate the NewSNR coincident detection statistic """

    def single(self, trigs):
        """Calculate the single detector statistic, here equal to newsnr

        Parameters
        ----------
        trigs: dict of numpy.ndarrays

        Returns
        -------
        numpy.ndarray
            The array of single detector values
        """
        return get_newsnr(trigs)

    def coinc(self, s0, s1, slide, step):
        """Calculate the coincident detection statistic.

        Parameters
        ----------
        s0: numpy.ndarray
            Single detector ranking statistic for the first detector.
        s1: numpy.ndarray
            Single detector ranking statistic for the second detector.
        slide: (unused in this statistic)
        step: (unused in this statistic)

        Returns
        -------
        numpy.ndarray
            Array of coincident ranking statistic values
        """
        return (s0**2. + s1**2.) ** 0.5


class NewSNRCutStatistic(NewSNRStatistic):

    """Same as the NewSNR statistic, but demonstrates a cut of the triggers"""

    def single(self, trigs):
        """Calculate the single detector statistic.

        Parameters
        ----------
        trigs: dict of numpy.ndarrays
            Dictionary of the single detector trigger information. 'chisq_dof',
        'snr', and 'chisq' are required keys

        Returns
        -------
        newsnr: numpy.ndarray
            Array of single detector values
        """
        newsnr = get_newsnr(trigs)
        rchisq = trigs['chisq'] / (2. * trigs['chisq_dof'] - 2.)
        newsnr[numpy.logical_and(newsnr < 10, rchisq > 2)] = -1
        return newsnr

    def coinc(self, s0, s1, slide, step):
        """Calculate the coincident detection statistic.

        Parameters
        ----------
        s0: numpy.ndarray
            Single detector ranking statistic for the first detector.
        s1: numpy.ndarray
            Single detector ranking statistic for the second detector.
        slide: (unused in this statistic)
        step: (unused in this statistic)

        Returns
        -------
        cstat: numpy.ndarray
            Array of coincident ranking statistic values
        """
        cstat = (s0**2. + s1**2.) ** 0.5
        cstat[s0==-1] = 0
        cstat[s1==-1] = 0
        return cstat


class PhaseTDStatistic(NewSNRStatistic):

    """Statistic that re-weights combined newsnr using coinc parameters.

    The weighting is based on the PDF of time delays, phase differences and
    amplitude ratios between triggers in different ifos.
    """
    def __init__(self, files):
        NewSNRStatistic.__init__(self, files)
        self.hist = self.files['phasetd_newsnr']['map'][:]

        #normalize so that peak has no effect on newsnr
        self.hist = self.hist / float(self.hist.max())
        self.hist = numpy.log(self.hist)

        # Bin boundaries are stored in the hdf file
        self.tbins = self.files['phasetd_newsnr']['tbins'][:]
        self.pbins = self.files['phasetd_newsnr']['pbins'][:]
        self.sbins = self.files['phasetd_newsnr']['sbins'][:]
        self.rbins = self.files['phasetd_newsnr']['rbins'][:]

        self.single_dtype = [('snglstat', numpy.float32),
                    ('coa_phase', numpy.float32),
                    ('end_time', numpy.float64),
                    ('sigmasq', numpy.float32),
                    ('snr', numpy.float32)]

    def single(self, trigs):
        """
        Calculate the single detector statistic and assemble other parameters

        Parameters
        ----------
        trigs: dict of numpy.ndarrays
            Dictionary holding single detector trigger information.
        'chisq_dof', 'snr', 'chisq', 'coa_phase', 'end_time', and 'sigmasq'
        are required keys.

        Returns
        -------
        numpy.ndarray
            Array of single detector parameter values
        """
        sngl_stat = get_newsnr(trigs)
        singles = numpy.zeros(len(sngl_stat), dtype=self.single_dtype)
        singles['snglstat'] = sngl_stat
        singles['coa_phase'] = trigs['coa_phase']
        singles['end_time'] = trigs['end_time']
        singles['sigmasq'] = trigs['sigmasq']
        singles['snr'] = trigs['snr']
        return numpy.array(singles, ndmin=1)

    def logsignalrate(self, s0, s1, slide, step):
        """Calculate the normalized log rate density of signals via lookup"""
        td = numpy.array(s0['end_time'] - s1['end_time'] - slide*step, ndmin=1)
        pd = numpy.array((s0['coa_phase'] - s1['coa_phase']) % \
                         (2. * numpy.pi), ndmin=1)
        rd = numpy.array((s0['sigmasq'] / s1['sigmasq']) ** 0.5, ndmin=1)
        sn0 = numpy.array(s0['snr'], ndmin=1)
        sn1 = numpy.array(s1['snr'], ndmin=1)

        snr0 = sn0 * 1
        snr1 = sn1 * 1

        snr0[rd > 1] = sn1[rd > 1]
        snr1[rd > 1] = sn0[rd > 1]
        rd[rd > 1] = 1. / rd[rd > 1]

        # Find which bin each coinc falls into
        tv = numpy.searchsorted(self.tbins, td) - 1
        pv = numpy.searchsorted(self.pbins, pd) - 1
        s0v = numpy.searchsorted(self.sbins, snr0) - 1
        s1v = numpy.searchsorted(self.sbins, snr1) - 1
        rv = numpy.searchsorted(self.rbins, rd) - 1

        # Enforce that points fits into the bin boundaries: if a point lies
        # outside the boundaries it is pushed back to the nearest bin.
        tv[tv < 0] = 0
        tv[tv >= len(self.tbins) - 1] = len(self.tbins) - 2
        pv[pv < 0] = 0
        pv[pv >= len(self.pbins) - 1] = len(self.pbins) - 2
        s0v[s0v < 0] = 0
        s0v[s0v >= len(self.sbins) - 1] = len(self.sbins) - 2
        s1v[s1v < 0] = 0
        s1v[s1v >= len(self.sbins) - 1] = len(self.sbins) - 2
        rv[rv < 0] = 0
        rv[rv >= len(self.rbins) - 1] = len(self.rbins) - 2
        
        return self.hist[tv, pv, s0v, s1v, rv]

    def coinc(self, s0, s1, slide, step):
        """
        Calculate the coincident detection statistic.

        Parameters
        ----------
        s0: numpy.ndarray
            Single detector ranking statistic for the first detector.
        s1: numpy.ndarray
            Single detector ranking statistic for the second detector.
        slide: numpy.ndarray
            Array of ints. These represent the multiple of the timeslide
        interval to bring a pair of single detector triggers into coincidence.
        step: float
            The timeslide interval in seconds.

        Returns
        -------
        coinc_stat: numpy.ndarray
            An array of the coincident ranking statistic values
        """
        rstat = s0['snglstat']**2. + s1['snglstat']**2.
        cstat = rstat + 2. * self.logsignalrate(s0, s1, slide, step)
        cstat[cstat < 0] = 0
        return cstat ** 0.5


class ExpFitStatistic(NewSNRStatistic):

    """Detection statistic using an exponential falloff noise model.

    Statistic approximates the negative log noise coinc rate density per
    template over single-ifo newsnr values.
    """

    def __init__(self, files):
        if not len(files):
            raise RuntimeError("Can't find any statistic files !")
        NewSNRStatistic.__init__(self, files)
        # the stat file attributes are hard-coded as '%{ifo}-fit_coeffs'
        parsed_attrs = [f.split('-') for f in self.files.keys()]
        self.ifos = [at[0] for at in parsed_attrs if
                     (len(at) == 2 and at[1] == 'fit_coeffs')]
        if not len(self.ifos):
            raise RuntimeError("None of the statistic files has the required "
                               "attribute called {ifo}-fit_coeffs !")
        self.fits_by_tid = {}
        self.alphamax = {}
        for i in self.ifos:
            self.fits_by_tid[i] = self.assign_fits(i)
            self.get_ref_vals(i)

        self.get_newsnr = get_newsnr

    def assign_fits(self, ifo):
        coeff_file = self.files[ifo+'-fit_coeffs']
        template_id = coeff_file['template_id'][:]
        alphas = coeff_file['fit_coeff'][:]
        lambdas = coeff_file['count_above_thresh'][:]
        # the template_ids and fit coeffs are stored in an arbitrary order
        # create new arrays in template_id order for easier recall
        tid_sort = numpy.argsort(template_id)
        return {'alpha':alphas[tid_sort], 'lambda':lambdas[tid_sort],
                'thresh':coeff_file.attrs['stat_threshold']}

    def get_ref_vals(self, ifo):
        self.alphamax[ifo] = self.fits_by_tid[ifo]['alpha'].max()

    def find_fits(self, trigs):
        """Get fit coeffs for a specific ifo and template id"""
        tnum = trigs.template_num
        # fits_by_tid is a dictionary of dictionaries of arrays
        # indexed by ifo / coefficient name / template_id
        alphai = self.fits_by_tid[trigs.ifo]['alpha'][tnum]
        lambdai = self.fits_by_tid[trigs.ifo]['lambda'][tnum]
        thresh = self.fits_by_tid[trigs.ifo]['thresh']
        return alphai, lambdai, thresh

    def lognoiserate(self, trigs):
        """
        Calculate the log noise rate density over single-ifo newsnr

        Read in single trigger information, make the newsnr statistic
        and rescale by the fitted coefficients alpha and lambda
        """
        alphai, lambdai, thresh = self.find_fits(trigs)
        newsnr = self.get_newsnr(trigs)
        # alphai is constant of proportionality between single-ifo newsnr and
        #  negative log noise likelihood in given template
        # lambdai is rate of trigs in given template compared to average
        # thresh is stat threshold used in given ifo
        lognoisel = - alphai * (newsnr - thresh) + numpy.log(alphai) + \
                      numpy.log(lambdai)
        return numpy.array(lognoisel, ndmin=1, dtype=numpy.float32)

    def single(self, trigs):
        """Single-detector statistic, here just equal to the log noise rate"""
        return self.lognoiserate(trigs)

    def coinc(self, s0, s1, slide, step):
        """Calculate the final coinc ranking statistic"""
        # Approximate log likelihood ratio by summing single-ifo negative
        # log noise likelihoods
        loglr = - s0 - s1
        # add squares of threshold stat values via idealized Gaussian formula
        threshes = [self.fits_by_tid[i]['thresh'] for i in self.ifos]
        loglr += sum([t**2. / 2. for t in threshes])
        # convert back to a coinc-SNR-like statistic
        # via log likelihood ratio \propto rho_c^2 / 2
        return (2. * loglr) ** 0.5


class ExpFitCombinedSNR(ExpFitStatistic):

    """Reworking of ExpFitStatistic designed to resemble network SNR

    Use a monotonic function of the negative log noise rate density which
    approximates combined (new)snr for coincs with similar newsnr in each ifo
    """

    def __init__(self, files):
        ExpFitStatistic.__init__(self, files)
        # for low-mass templates the exponential slope alpha \approx 6
        self.alpharef = 6.

    def use_alphamax(self):
        # take reference slope as the harmonic mean of individual ifo slopes
        inv_alphas = [1./self.alphamax[i] for i in self.ifos]
        self.alpharef = (sum(inv_alphas)/len(inv_alphas))**-1
        print self.alpharef

    def single(self, trigs):
        logr_n = self.lognoiserate(trigs)
        _, _, thresh = self.find_fits(trigs)
        # shift by log of reference slope alpha
        logr_n += -1. * numpy.log(self.alpharef)
        # add threshold and rescale by reference slope
        stat = thresh - (logr_n / self.alpharef)
        return numpy.array(stat, ndmin=1, dtype=numpy.float32)

    def coinc(self, s0, s1, slide, step):
        # scale by 1/sqrt(2) to resemble network SNR
        return (s0 + s1) / (2.**0.5)


class PhaseTDExpFitStatistic(PhaseTDStatistic, ExpFitCombinedSNR):

    """Statistic combining exponential noise model with signal histogram PDF"""

    def __init__(self, files):
        # read in both foreground PDF and background fit info
        ExpFitCombinedSNR.__init__(self, files)
        # need the self.single_dtype value from PhaseTDStatistic
        PhaseTDStatistic.__init__(self, files)

    def single(self, trigs):
        # same single-ifo stat as ExpFitCombinedSNR
        sngl_stat = ExpFitCombinedSNR.single(self, trigs)
        singles = numpy.zeros(len(sngl_stat), dtype=self.single_dtype)
        singles['snglstat'] = sngl_stat
        singles['coa_phase'] = trigs['coa_phase']
        singles['end_time'] = trigs['end_time']
        singles['sigmasq'] = trigs['sigmasq']
        singles['snr'] = trigs['snr']
        return numpy.array(singles, ndmin=1)

    def coinc(self, s0, s1, slide, step):
        # logsignalrate function inherited from PhaseTDStatistic
        logr_s = self.logsignalrate(s0, s1, slide, step)
        # rescale by ExpFitCombinedSNR reference slope as for sngl stat
        cstat = s0['snglstat'] + s1['snglstat'] + logr_s / self.alpharef
        # cut off underflowing and very small values
        cstat[cstat < 8.] = 8.
        # scale to resemble network SNR
        return cstat / (2.**0.5)


class PhaseTDExpFitSGStatistic(PhaseTDExpFitStatistic):

    """Statistic combining exponential noise model with signal histogram PDF
       and adding the sine-Gaussian veto to the single detector ranking
    """

    def __init__(self, files):
        PhaseTDExpFitStatistic.__init__(self, files)
        self.get_newsnr = get_newsnr_sgveto


class MaxContTradNewSNRStatistic(NewSNRStatistic):

    """Combination of NewSNR with the power chisq and auto chisq"""

    def single(self, trigs):
        """ Calculate the single detector statistic.

        Parameters
        ----------
        trigs: dict of numpy.ndarrays
            Dictionary of the single detector trigger information. 'chisq_dof',
        'snr', 'cont_chisq', 'cont_chisq_dof', and 'chisq' are required arrays
        for this statistic.

        Returns
        -------
        stat: numpy.ndarray
            The array of single detector values
        """
        chisq_newsnr = get_newsnr(trigs)
        rautochisq = trigs['cont_chisq'] / trigs['cont_chisq_dof']
        autochisq_newsnr = events.newsnr(trigs['snr'], rautochisq)
        return numpy.array(numpy.minimum(chisq_newsnr, autochisq_newsnr,
                           dtype=numpy.float32), ndmin=1, copy=False)


statistic_dict = {
    'newsnr': NewSNRStatistic,
    'newsnr_cut': NewSNRCutStatistic,
    'phasetd_newsnr': PhaseTDStatistic,
    'exp_fit_stat': ExpFitStatistic,
    'exp_fit_csnr': ExpFitCombinedSNR,
    'phasetd_exp_fit_stat': PhaseTDExpFitStatistic,
    'max_cont_trad_newsnr': MaxContTradNewSNRStatistic,
    'phasetd_exp_fit_stat_sgveto': PhaseTDExpFitSGStatistic
}

def get_statistic(stat):
    """
    Error-handling sugar around dict lookup

    Parameters
    ----------
    stat : string
        Name of the statistic

    Returns
    -------
    class
        Subclass of Stat base class

    Raises
    ------
    RuntimeError
        If the string is not recognized as corresponding to a Stat subclass
    """
    try:
        return statistic_dict[stat]
    except KeyError:
        raise RuntimeError('%s is not an available detection statistic' % stat)


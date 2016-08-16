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

class NewSNRStatistic(Stat):

    """ Calculate the NewSNR coincident detection statistic """

    def single(self, trigs):
        """Calculate the single detector statistic.

        Parameters
        ----------
        trigs: dict of numpy.ndarrays
            Dictionary of the single detector trigger information. 'chisq_dof',
        'snr', and 'chisq' are required arrays for this statistic.

        Returns
        -------
        stat: numpy.ndarray
            The array of single detector values
        """
        dof = 2 * trigs['chisq_dof'] - 2
        newsnr = events.newsnr(trigs['snr'], trigs['chisq'] / dof)
        return numpy.array(newsnr, ndmin=1, dtype=numpy.float32)

    def coinc(self, s1, s2, slide, step):
        """Calculate the coincident detection statistic.

        Parameters
        ----------
        s1: numpy.ndarray
            Single detector ranking statistic for the first detector.
        s2: numpy.ndarray
            Single detector ranking statistic for the second detector.
        slide: (unused in this statistic!)
        step: (unused in this statistic!)

        Returns
        -------
        coinc_stat: numpy.ndarray
            An array of the coincident ranking statitic values
        """
        return (s1**2.0 + s2**2.0) ** 0.5

class NewSNRCutStatistic(Stat):

    """Same as the NewSNR statistic, but demonstrates a cut of the triggers"""

    def single(self, trigs):
        """Calculate the single detector statistic.

        Parameters
        ----------
        trigs: dict of numpy.ndarrays
            Dictionary of the single detector trigger information. 'chisq_dof',
        'snr', and 'chisq' are requierd arrays for this statistic.

        Returns
        -------
        stat: numpy.ndarray
            The array of single detector values
        """
        dof = 2 * trigs['chisq_dof'] - 2
        rchisq = trigs['chisq'] / dof
        newsnr = events.newsnr(trigs['snr'], rchisq)
        newsnr[numpy.logical_and(newsnr < 10, rchisq > 2)] = -1
        return newsnr

    def coinc(self, s1, s2, slide, step):
        """Calculate the coincident detection statistic.

        Parameters
        ----------
        s1: numpy.ndarray
            Single detector ranking statistic for the first detector.
        s2: numpy.ndarray
            Single detector ranking statistic for the second detector.
        slide: (unused in this statistic!)
        step: (unused in this statistic!)

        Returns
        -------
        coinc_stat: numpy.ndarray
            An array of the coincident ranking statitic values
        """
        cstat = (s1**2.0 + s2**2.0) ** 0.5
        cstat[s1==-1] = 0
        cstat[s2==-1] = 0
        return cstat

class PhaseTDStatistic(NewSNRStatistic):

    """Detection statistic that re-weights the network SNR based on the
    PDF of time delays, phase difference, and amplitude ratios.
    """

    def __init__(self, files):
        NewSNRStatistic.__init__(self, files)
        self.hist = self.files['phasetd_newsnr']['map'][:]
        top = float(self.hist.max())

        #normalize so that peak has no effect on newsnr
        self.hist = self.hist / top
        self.hist = numpy.log(self.hist)

    def single(self, trigs):
        """Calculate the single detector statistic.

        Parameters
        ----------
        trigs: dict of numpy.ndarrays
            Dictionary of the single detector trigger information. 'chisq_dof',
        'snr', and 'chisq', 'coa_phase', 'end_time', and 'sigmasq'
        are required arrays for this statistic.

        Returns
        -------
        stat: numpy.ndarray
            The array of single detector values
        """
        newsnr = NewSNRStatistic.single(self, trigs)
        return numpy.array((newsnr, trigs['coa_phase'], trigs['end_time'],
                            trigs['sigmasq']**0.5, trigs['snr'])).transpose()

    def coinc(self, s1, s2, slide, step):
        """Calculate the coincident detection statistic.

        Parameters
        ----------
        s1: numpy.ndarray
            Single detector ranking statistic for the first detector.
        s2: numpy.ndarray
            Single detector ranking statistic for the second detector.
        slide: numpy.ndarray
            Array of ints. These represent the multiple of the timeslide
        interval to bring a pair of single detector triggers into coincidence.
        step: float
            The timeslide interval in seconds.

        Returns
        -------
        coinc_stat: numpy.ndarray
            An array of the coincident ranking statitic values
        """
        td = s1[:,2] - s2[:,2] - slide * step
        pd = (s1[:,1] - s2[:,1]) % (numpy.pi * 2)
        rd = s1[:, 3] / s2[:, 3]
        sn1 = s1[:,4]
        sn2 = s2[:,4]
 
        snr1 = sn1 * 1
        snr2 = sn2 * 1
        snr1[rd > 1] = sn2[rd > 1]
        snr2[rd > 1] = sn1[rd > 1]
        rd[rd > 1] = 1.0 / rd[rd > 1]

        # These are the bin boundaries stored in the hdf file
        tbins = self.files['phasetd_newsnr']['tbins'][:]
        pbins = self.files['phasetd_newsnr']['pbins'][:]
        sbins = self.files['phasetd_newsnr']['sbins'][:]
        rbins = self.files['phasetd_newsnr']['rbins'][:]

        # Find which bin each coinc falls into
        tv = numpy.searchsorted(tbins, td) - 1
        pv = numpy.searchsorted(pbins, pd) - 1
        s1v = numpy.searchsorted(sbins, snr1) - 1
        s2v = numpy.searchsorted(sbins, snr2) - 1    
        rv = numpy.searchsorted(rbins, rd) - 1  

        # The following just enforces that the point fits into
        # the bin boundaries. If a point lies outside the boundaries it is
        # pushed back to the nearest bin.
        tv[tv < 0] = 0
        tv[tv >= len(tbins) - 1] = len(tbins) - 2
        pv[pv < 0] = 0
        pv[pv >= len(pbins) - 1] = len(pbins) - 2
        s1v[s1v < 0] = 0
        s1v[s1v >= len(sbins) - 1] = len(sbins) - 2
        s2v[s2v < 0] = 0
        s2v[s2v >= len(sbins) - 1] = len(sbins) - 2
        rv[rv < 0] = 0
        rv[rv >= len(rbins) - 1] = len(rbins) - 2

        rstat = s1[:,0]**2.0 + s2[:,0]**2.0
        cstat = rstat + 2.0 * self.hist[tv, pv, s1v, s2v, rv]
        cstat[cstat < 0] = 0
        return cstat ** 0.5

class ExpFitStatistic(NewSNRStatistic):
    def __init__(self, files):
        NewSNRStatistic.__init__(self, files)
        if not len(files):
            raise RuntimeError("Can't find any statistic files !")
        # the stat file attributes are hard-coded as '%{ifo}-fit_coeffs'
        self.ifos = [f.split('-')[0] for f in self.statfiles.keys() if \
                     f.split('-')[1] == 'fit_coeffs']
        if (not len(files)) or (not len(self.ifos)):
            raise RuntimeError("None of the statistic files has the required "
                               "attribute called {ifo}-fit_coeffs !")
        self.fits_by_tid = {}
        for i in self.ifos:
           self.fits_by_tid[i] = self.assign_fits(i)

    def assign_fits(self, ifo):
        coeff_file = self.statfiles[ifo+'-fit_coeffs']
        template_id = coeff_file['template_id'][:]
        alphas = coeff_file['fit_coeff'][:]
        lambdas = coeff_file['count_above_thresh'][:]
        # the template_ids and fit coeffs are stored in an arbitrary order
        # create new arrays in template_id order for easier recall
        tid_sort = numpy.argsort(template_id)
        return {'alpha':alphas[tid_sort], 'lambda':lambdas[tid_sort],
                'thresh':coeff_file.attrs['stat_threshold']}

    def find_fits(self, trigs):
        """Get fit coeffs for a specific ifo and template id"""
        ifo = trigs.ifo
        tnum = trigs.template_num
        # fits_by_template is a dictionary of dictionaries of arrays
        # indexed by ifo / coefficient name / template_id
        alphai = self.fits_by_tid[ifo]['alpha'][tnum][0]
        lambdai = self.fits_by_tid[ifo]['lambda'][tnum][0]
        thresh = self.fits_by_tid[ifo]['thresh']
        return alphai, lambdai, thresh

    def single(self, trigs):
        """
        Calculate the parts of the coinc statistic depending on sngl parameters

        Read in single trigger information, make the newsnr statistic
        and rescale by the fitted coefficients alpha and lambda
        """
        alphai, lambdai, thresh = self.find_fits(trigs)
        dof = 2 * trigs['chisq_dof'] - 2
        newsnr = events.newsnr(trigs['snr'], trigs['chisq'] / dof)
        # alphai is constant of proportionality between single-ifo newsnr and
        #  negative log noise likelihood in given template
        # lambdai is rate of trigs in given template compared to average
        # thresh is stat threshold used in given ifo
        lognoisel = - alphai * (newsnr - thresh) + numpy.log(alphai) + \
                      numpy.log(lambdai)
        return numpy.array(lognoisel, ndmin=1, dtype=numpy.float32)

    def coinc(self, s0, s1, slide, step):
        """Calculate the final coinc ranking statistic"""
        # Approximate log likelihood ratio by summing single-ifo negative
        # log noise likelihoods
        loglr = - s0 - s1
        # add squares of threshold stat values with notional Gaussian formula
        threshes = [self.fits_by_tid[i]['thresh'] for i in self.ifos]
        loglr += sum([t**2 / 2. for t in threshes])
        # convert back to a coinc-SNR-like statistic
        # notionally, log likelihood ratio \propto rho_c^2 / 2
        return (2 * loglr) ** 0.5

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
        chisq_dof = 2 * trigs['chisq_dof'] - 2
        chisq_newsnr = events.newsnr(trigs['snr'], trigs['chisq'] / chisq_dof)
        autochisq_dof = trigs['cont_chisq_dof']
        autochisq_newsnr = events.newsnr(trigs['snr'],
                                         trigs['cont_chisq'] / autochisq_dof)
        return numpy.array(numpy.minimum(chisq_newsnr, autochisq_newsnr,
                             dtype=numpy.float32), ndmin=1, copy=False)

statistic_dict = {
    'newsnr': NewSNRStatistic,
    'newsnr_cut': NewSNRCutStatistic,
    'phasetd_newsnr': PhaseTDStatistic,
    'exp_fit_stat': ExpFitStatistic,
    'max_cont_trad_newsnr': MaxContTradNewSNRStatistic,
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


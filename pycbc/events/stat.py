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

def get_statistic(option, files):
    """ Return an appropriate coincident statistic class instance.

    Parameters
    ----------
    option: str
        Name of the statistic class to instantiate
    files: list of strs
        A list containing the filenames of hdf format files used to help
    construct the coincident statistics.

    Returns
    -------
    stat: instance of `Stat` class
        The requested class instance        
    """
    if option == 'newsnr':
        return NewSNRStatistic(files)
    elif option == 'newsnr_cut':
        return NewSNRCutStatistic(files)
    elif option == 'phasetd_newsnr':
        return PhaseTDStatistic(files)
    elif option == 'max_cont_trad_newsnr':
        return MaxContTradNewSNRStatistic(files)
    else:
        raise ValueError('%s is not an available detection statistic' % option)

class Stat(object):

    """ Base class which should be extended to provide a coincident statistic"""

    def __init__(self, files):
        """ Create a statistic class instance

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
        """ Calculate the single detector statistic.
        
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
        """ Calculate the coincident detection statistic. 

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

    """ Same as the NewSNR statistic, but demonstrates a cut of the triggers"""

    def single(self, trigs):
        """ Calculate the single detector statistic.

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
        """ Calculate the coincident detection statistic. 

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

    """ Detection statistic that re-weights the network SNR based on the
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
        """ Calculate the single detector statistic.
        
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
        """ Calculate the coincident detection statistic. 

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
        rd[rd > 1] = 1.0 / rd[rd > 1]

        # These are the bin boundaries stored in the hdf file        
        tbins = self.files['phasetd_newsnr']['tbins'][:]
        pbins = self.files['phasetd_newsnr']['pbins'][:]
        sbins = self.files['phasetd_newsnr']['sbins'][:]
        rbins = self.files['phasetd_newsnr']['rbins'][:]

        # Find which bin each coinc falls into        
        tv = numpy.searchsorted(tbins, td) - 1 
        pv = numpy.searchsorted(pbins, pd) - 1
        s1v = numpy.searchsorted(sbins, s1[:,4]) - 1
        s2v = numpy.searchsorted(sbins, s2[:,4]) - 1    
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

class MaxContTradNewSNRStatistic(NewSNRStatistic):

    """ Combination of NewSNR with the power chisq and auto chisq """
    
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

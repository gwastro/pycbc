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
""" This modules contains functions for calculating the coincident 
statistic values
"""
import h5py
import numpy
from . import events

def get_statistic(option, files):
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
    def __init__(self, files):
        self.files = {}
        for filename in files:
            f = h5py.File(filename, 'r')
            stat = f.attrs['stat']
            self.files[stat] = f

class NewSNRStatistic(Stat):
    def single(self, trigs):
        """ Read in the single detector information and make a single detector
        statistic. Results can either be a single number or a record array.
        """
        
        dof = 2 * trigs['chisq_dof'] - 2
        newsnr = events.newsnr(trigs['snr'], trigs['chisq'] / dof)
        return numpy.array(newsnr, ndmin=1, dtype=numpy.float32)

    def coinc(self, s1, s2, slide, step):
        """ Calculate the coincident statistic.
        """
        return (s1**2.0 + s2**2.0) ** 0.5

class NewSNRCutStatistic(Stat):
    def single(self, trigs):
        dof = 2 * trigs['chisq_dof'] - 2
        rchisq = trigs['chisq'] / dof
        newsnr = events.newsnr(trigs['snr'], rchisq)
        newsnr[numpy.logical_and(newsnr < 10, rchisq > 2)] = -1
        return newsnr

    def coinc(self, s1, s2, slide, step):
        cstat = (s1**2.0 + s2**2.0) ** 0.5
        cstat[s1==-1] = 0
        cstat[s2==-1] = 0
        return cstat

class PhaseTDStatistic(NewSNRStatistic):
    def __init__(self, files):
        NewSNRStatistic.__init__(self, files)
        self.hist = self.files['phasetd_newsnr']['map'][:]
        top = float(self.hist.max())     

        #normalize so that peak has no effect on newsnr
        self.hist = self.hist / top
        self.hist = numpy.log(self.hist) 

    def single(self, trigs):
        newsnr = NewSNRStatistic.single(self, trigs)
        return numpy.array((newsnr, trigs['coa_phase'], trigs['end_time'],
                            trigs['sigmasq']**0.5, trigs['snr'])).transpose()

    def coinc(self, s1, s2, slide, step):
        """ Calculate the coincident statistic.
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

        # The following just enforces that the point fits into the bin boundaries        
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
        
        m = self.hist[tv, pv, s1v, s2v, rv]
        rstat = s1[:,0]**2.0 + s2[:,0]**2.0
        cstat = rstat + 2.0 * m
        cstat[cstat < 0] = 0        
        return cstat ** 0.5

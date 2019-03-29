# Copyright (C) 2019 Alex Nitz
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

"""Implementation of a relative likelihood function
"""

import numpy, logging
import scipy.special

from pycbc import filter as pyfilter
from pycbc.waveform.spa_tmplt import spa_tmplt, spa_amplitude_factor
from pycbc.detector import Detector
from pycbc.conversions import mchirp_from_mass1_mass2

from .base_data import BaseDataModel

# In this model we only calculate terms up to a constant.
# We are primarily interested in the posterior result
def angle(v1, v2):
    p = (v1.real * v2.real + v1.imag * v2.imag)
    return numpy.arccos(p / abs(v1) / abs(v2))

def frequency_bins(hp, hp2, flow, fhigh, thr):
    ratio = hp / hp2

    kmin = int(flow / hp.delta_f)
    kmax = int(fhigh / hp.delta_f)
    krange = numpy.arange(kmin, kmax)

    edges = [kmin]
    kref = kmin

    for k in krange:
        a = angle(ratio[kref], ratio[k])
        if a > thr:
            edges.append(k-1)
            kref = k - 1

    edges.append(kmax)
    edges = numpy.array(edges)

    bins = []
    for i in range(len(edges)-1):
        bins.append((edges[i], edges[i+1]))
    return edges, bins

class Relative(BaseDataModel):
    r"""Model that assumes we know all the intrinsic parameters.

    This model assumes we know all the intrinsic parameters, and are only
    maximizing over the extrinsic ones. We also assume a dominant mode waveform
    approximant only and non-precessing.
    """
    name = 'spafocus'

    def __init__(self, data, psds,
                 mass1, mass2, spin1z, spin2z,
                 dmass1, dmass2, dphase,
                 tstart, tend,
                 low_frequency_cutoff,
                 high_frequency_cutoff,
                 sample_rate=32768,
                 **kwargs):

        super(Relative, self).__init__(data=data, **kwargs)

        low_frequency_cutoff = float(low_frequency_cutoff)
        high_frequency_cutoff = float(high_frequency_cutoff)
        mass1 = float(mass1)
        mass2 = float(mass2)
        spin1z = float(spin1z)
        spin2z = float(spin2z)
        dmass1 = float(dmass2)
        dmass2 = float(dmass2)
        dphase = float(dphase)
        self.sample_rate = int(sample_rate)
        self.tstart = float(tstart)
        self.tend = float(tend)
        self.psds = psds
        self.data = data

        # Generate reference waveform
        df = data[data.keys()[0]].delta_f
        hp = spa_tmplt(delta_f=df, distance=1,
                       mass2=mass2, mass1=mass1, spin1z=spin1z,
                       spin2z=spin2z, phase_order=-1, spin_order=-1,
                       f_lower=low_frequency_cutoff,
                       f_upper=high_frequency_cutoff + 2 * df)
        self.ampf = spa_amplitude_factor(mass1=mass1, mass2=mass2)

        # Figure out the frequency bins using another reference
        hp2 = spa_tmplt(delta_f=df, distance=1,
                       mass2=dmass1, mass1=dmass2, spin1z=spin1z,
                       spin2z=spin2z, phase_order=-1, spin_order=-1,
                       f_lower=low_frequency_cutoff,
                       f_upper=high_frequency_cutoff + 2 * df)
        self.edges, self.bins = frequency_bins(hp, hp2,
                              low_frequency_cutoff,
                              high_frequency_cutoff,
                              dphase)
        self.hreference = hp.take(self.edges)
        self.fedges = (self.edges * hp.delta_f).astype(numpy.float32)

        logging.info('Using %s bins for this model', len(self.bins))
        # Extend data and template to high sample rate
        flen = int(self.sample_rate / df) / 2 + 1
        hp.resize(flen)
        for ifo in data:
            data[ifo].resize(flen)

        # Calculate high sample rate SNR time series
        self.sh = {}
        self.shl = {}

        self.hh = {}
        self.det = {}
        for ifo in data:
            self.det[ifo] = Detector(ifo)

            # storage for the snr time series of each frequency range
            self.sh[ifo] = []
            self.shl[ifo] = []
            st = 0
            for l, h in self.bins:
                print l, h

                # Prepare frequency slice
                hpc = hp.copy()
                hpc[:l] = 0
                hpc[h:] = 0

                # Prepare linear frequency hp
                hpl = hp.copy()
                hpl[:l] = 0
                hpl[l:h] *= numpy.arange(0, h - l) / float(h - l)
                hpl[h:] = 0

                print len(hpc), len(hpl), len(hp), len(psds[ifo])

                # constant term
                snr, _, _ = pyfilter.matched_filter_core(
                    hpc.astype(numpy.complex128), data[ifo],
                    high_frequency_cutoff=high_frequency_cutoff,
                    psd=psds[ifo])
                self.sh[ifo].append(snr.time_slice(self.tstart, self.tend).numpy().copy() * 4.0 * df)

                # linear term
                snrl, _, _ = pyfilter.matched_filter_core(
                    hpl.astype(numpy.complex128), data[ifo],
                    high_frequency_cutoff=high_frequency_cutoff,
                    psd=psds[ifo])
                self.shl[ifo].append(snrl.time_slice(self.tstart, self.tend).numpy().copy() * 4.0 * df)
                st = st + snr

            self.hh[ifo] = -0.5 * pyfilter.sigmasq(
                hp.astype(numpy.complex128), psd=psds[ifo],
                low_frequency_cutoff=low_frequency_cutoff,
                high_frequency_cutoff=high_frequency_cutoff)

            self.sh[ifo] = numpy.array(self.sh[ifo]).transpose()
            self.shl[ifo] = numpy.array(self.shl[ifo]).transpose()

    def rel_sigmasq(self, m1, m2):
        return (spa_amplitude_factor(mass1=m1, mass2=m2) / self.ampf) ** 2.0

    def slow_likelihood(self, p, ifo, time):
        m1 = p['mass1']
        m2 = p['mass2']
        s1 = p['spin1z']
        s2 = p['spin2z']
        hp = spa_tmplt(f_lower=30.0, delta_f=self.psds[ifo].delta_f,
                            mass1=m1, mass2=m2,
                            distance=1,
                            spin1z=s1, spin2z=s2,
                            phase_order=-1, spin_order=-1)
        hp.resize(len(self.psds[ifo]))
        snr, _, _ = pyfilter.matched_filter_core(
                hp.astype(numpy.complex128), self.data[ifo],
                psd=self.psds[ifo],
                low_frequency_cutoff=30.0,
                high_frequency_cutoff=1000.0)
        return snr.at_time(time) * 4.0 * snr.delta_f


    def raw_likelihood(self, p, ifo, time):
        m1 = p['mass1']
        m2 = p['mass2']
        s1 = p['spin1z']
        s2 = p['spin2z']
        htarget = spa_tmplt(sample_points=self.fedges,
                            mass1=m1, mass2=m2,
                            distance=1,
                            spin1z=s1, spin2z=s2,
                            phase_order=-1, spin_order=-1)
        r = (htarget / self.hreference).conj()
        r0 = r[:-1]
        r1 = (r[1:] - r[:-1])
        idx = int((time - self.tstart) * self.sample_rate)
        return (r0*self.sh[ifo][idx] + r1*self.shl[ifo][idx]).sum()

    def _loglikelihood(self):
        return self.loglr

    def _lognl(self):
        return 0

    def _loglr(self):
        r"""Computes the log likelihood ratio

        Returns
        -------
        float
            The value of the log likelihood ratio.
        """
        # calculate <d-h|d-h> = <h|h> - 2<h|d> + <d|d> up to a constant
        p = self.current_params.copy()
        p.update(self.static_params)

        relsigsq = self.rel_sigmasq(p['mass1'], p['mass2'])
        shloglr = hhloglr = 0
        shloglr2 = 0
        for ifo in self.sh:
            fp, fc = self.det[ifo].antenna_pattern(p['ra'], p['dec'],
                                                   p['polarization'],
                                                   p['tc'])
            dt = self.det[ifo].time_delay_from_earth_center(p['ra'], p['dec'],
                                                            p['tc'])
            ip = numpy.cos(p['inclination'])
            ic = 0.5 * (1.0 + ip * ip)
            htf = (fp * ip + 1.0j * fc * ic) / p['distance']

            sh = self.raw_likelihood(p, ifo, p['tc'] + dt)
            shloglr += sh * htf
            hhloglr += self.hh[ifo] * abs(htf) ** 2.0 * relsigsq

            #sh2 = self.slow_likelihood(p, ifo, p['tc'] + dt)
            #shloglr2 += sh2 * htf


        vloglr = numpy.log(scipy.special.i0e(abs(shloglr)))
        vloglr += abs(shloglr) + hhloglr

        #vloglr2 = numpy.log(scipy.special.i0e(abs(shloglr2)))
        #loglr2 += abs(shloglr2) + hhloglr
        #if (vloglr2 > 150) or (vloglr > 150):
        #    print vloglr2, vloglr
        #print vloglr
        return float(vloglr)

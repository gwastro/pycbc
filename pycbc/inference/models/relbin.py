#!/usr/bin/env python

import numpy
import logging
from scipy.interpolate import interp1d

from pycbc.waveform.spa_tmplt import spa_tmplt
from pycbc.detector import Detector
from pycbc.types import Array

from .base_data import BaseDataModel


def setup_bins(f_full, f_lo, f_hi, chi=1.0, eps=0.5):
    """
    construct frequency binning
    f_full: full frequency grid
    [f_lo, f_hi] is the frequency range one would like to use for matched filtering
    chi, eps are tunable parameters [see Barak, Dai & Venumadhav 2018]
    return the number of bins, a list of bin-edge frequencies, and their positions in the full frequency grid
    """
    f = numpy.linspace(f_lo, f_hi, 10000)
    # f^ga power law index
    ga = numpy.array([-5./3, -2./3, 1., 5./3, 7./3])
    dalp = chi * 2.0 * numpy.pi / numpy.absolute((f_lo ** ga) - (f_hi ** ga))
    dphi = numpy.sum(numpy.array([numpy.sign(g) * d * (f ** g) for
                                  g, d in zip(ga, dalp)]), axis=0)
    dphi_diff = dphi - dphi[0]
    # now construct frequency bins
    nbin = int(dphi_diff[-1] / eps)
    dphi2f = interp1d(dphi_diff, f, kind='slinear', bounds_error=False,
                      fill_value=0.0)
    dphi_grid = numpy.linspace(dphi_diff[0], dphi_diff[-1], nbin+1)
    # frequency grid points
    fbin = dphi2f(dphi_grid)
    # indices of frequency grid points in the FFT array
    fbin_ind = numpy.array([numpy.argmin(numpy.absolute(f_full - ff)) for
                            ff in fbin])
    # make sure grid points are precise
    fbin = numpy.array([f_full[i] for i in fbin_ind])

    return nbin, fbin, fbin_ind

class Relative(BaseDataModel):
    name = "relative"
    def __init__(self, data, psds, mass1, mass2, spin1z, spin2z,
                 ra, dec, tc, low_frequency_cutoff,
                 high_frequency_cutoff, epsilon=0.5, **kwargs):
        super(Relative, self).__init__(data=data, **kwargs)
        # store data and frequencies
        f_lo = float(low_frequency_cutoff)
        f_hi = float(high_frequency_cutoff)
        self.f = numpy.array(data[data.keys()[0]].sample_frequencies)
        self.df = data[data.keys()[0]].delta_f
        self.data = data
        self.end_time = float(data[data.keys()[0]].end_time)
        self.det = {ifo: Detector(ifo) for ifo in data}
        epsilon = float(epsilon)
        # store data and psds as arrays for faster computation
        self.comp_data = {ifo: numpy.array(data[ifo].data) for ifo in data}
        self.comp_psds = {ifo: numpy.array(psds[ifo].data) for ifo in data}
        # store fiducial waveform params
        mass1 = float(mass1)
        mass2 = float(mass2)
        spin1z = float(spin1z)
        spin2z = float(spin2z)
        ra = float(ra)
        dec = float(dec)
        tc = float(tc)

        # get detector-specific arrival times relative to end of data
        dt = {ifo: self.det[ifo].time_delay_from_earth_center(ra, dec, tc) for
              ifo in data}
        self.ta = {ifo: tc + dt[ifo] - self.end_time for ifo in data}

        # generate fiducial waveform
        logging.info("Generating fiducial waveform")
        hp = spa_tmplt(f_lower=f_lo, f_upper=f_hi+self.df,
                       delta_f=self.df, mass1=mass1,
                       mass2=mass2, spin1z=spin1z,
                       spin2z=spin2z, distance=1.,
                       spin_order=-1, phase_order=-1)
        hp.resize(len(self.f))
        self.h00 = numpy.array(hp)

        # compute frequency bins
        logging.info("Computing frequency bins")
        nbin, fbin, fbin_ind = setup_bins(f_full=self.f, f_lo=f_lo, f_hi=f_hi,
                                          eps=epsilon)
        logging.info('Using %s bins for this model', nbin)
        # store bins and edges in sample and frequency space
        self.edges = fbin_ind
        self.fedges = numpy.array(fbin).astype(numpy.float32)
        self.bins = numpy.array([(self.edges[i], self.edges[i+1]) for
                                 i in range(len(self.edges) - 1)])
        self.fbins = numpy.array([(fbin[i], fbin[i+1]) for
                                  i in range(len(fbin) - 1)])
        # store low res fiducial waveform
        self.h00_sparse = self.h00.copy().take(self.edges)

        # compute summary data
        logging.info("Calculating summary data at frequency resolution %s",
                      self.df)
        self.sdat = self.summary_data()

        # store the psds and calculate the inner product weight
        logging.info("Calculating lognl")
        kmin, kmax = self.edges[0], self.edges[-1]
        self._kmin = {ifo: kmin for ifo in self.data}
        self._kmax = {ifo: kmax for ifo in self.data}
        self._data = data.copy()
        self._f_lower = {ifo: f_lo for ifo in self.data}
        self._f_upper = {ifo: f_hi for ifo in self.data}
        self._N = int(1. / (data[data.keys()[0]].delta_t * self.df))
        self._psds = {}
        self._weight = {}
        self._lognorm = {}
        self._det_lognls = {}
        self.set_psds(psds)
        logging.info("Lognl is %s", self.lognl)

    def summary_data(self):
        # timeshift the fiducial waveform for each detector
        shift = {ifo: numpy.exp(-2.0j * numpy.pi * self.f * self.ta[ifo]) for
                 ifo in self.data}
        h0 = {ifo: self.h00.copy() * shift[ifo] for ifo in self.data}
        # calculate coefficients
        sdat = {}
        for ifo in self.data:
            hd = numpy.conjugate(self.comp_data[ifo]) * h0[ifo]
            hd /= self.comp_psds[ifo]
            hh = (numpy.absolute(h0[ifo]) ** 2.0) / self.comp_psds[ifo]
            # constant terms
            a0 = numpy.array([4. * self.df * numpy.sum(hd[l:h]) for
                              l, h in self.bins])
            b0 = numpy.array([4. * self.df * numpy.sum(hh[l:h]) for
                              l, h in self.bins])
            # linear terms
            bin_centers = [0.5 * (fl + fh) for fl, fh in self.fbins]
            a1 = numpy.array([4. * self.df * \
                              numpy.sum(hd[l:h] * (self.f[l:h] - bc)) for
                              (l, h), bc in zip(self.bins, bin_centers)])
            b1 = numpy.array([4. * self.df * \
                              numpy.sum(hh[l:h] * (self.f[l:h] - bc)) for
                              (l, h), bc in zip(self.bins, bin_centers)])

            sdat[ifo] = {'a0': a0, 'a1': a1,
                         'b0': b0, 'b1': b1}
        return sdat

    def waveform_ratio(self, p, htf, dtc=0.0):
        # generate template
        hp = spa_tmplt(sample_points=self.fedges,
                       mass1=p['mass1'], mass2=p['mass2'],
                       spin1z=p['spin1z'], spin2z=p['spin2z'],
                       distance=1., spin_order=-1, phase_order=-1)
        htarget = numpy.array(hp)
        # apply antenna pattern, inclination, and distance factor
        htarget *= htf
        # compute waveform ratio and timeshift
        shift = numpy.exp(-2.0j * numpy.pi * self.fedges * dtc)
        r = htarget / self.h00_sparse * shift
        r0 = 0.5 * (r[:-1] + r[1:])
        r1 = (r[1:] - r[:-1]) / (self.fedges[1:] - self.fedges[:-1])
        return numpy.array([r0, r1], dtype=numpy.complex128)

    def binned_inner_products(self, r0, r1, sdat):
        # <h, d> is real part of the sum over bins of A0r0 + A1r1
        hd = numpy.sum(sdat['a0'] * r0 + sdat['a1'] * r1).real
        # <h, h> is real part of the sum over bins of B0|r0|^2 + 2B1Re(r1r0*)
        hh = numpy.sum(sdat['b0'] * numpy.absolute(r0) ** 2.
                       + 2. * sdat['b1'] * (r1 * numpy.conjugate(r0)).real).real

        return hd - 0.5 * hh

    def _loglikelihood(self):
        # get model params
        p = self.current_params.copy()
        p.update(self.static_params)

        # start at lognl to end at loglikelihood instead of loglr
        loglike = self.lognl
        for ifo in self.data:
            # get detector antenna pattern
            fp, fc = self.det[ifo].antenna_pattern(p['ra'], p['dec'],
                                                   p['polarization'],
                                                   p['tc'])
            ip = numpy.cos(p['inclination'])
            ic = 0.5 * (1.0 + ip * ip)
            htf = (fp * ip + 1.0j * fc * ic) / p['distance']
            # get timeshift relative to fiducial waveform
            dt = self.det[ifo].time_delay_from_earth_center(p['ra'], p['dec'],
                                                            p['tc'])
            dtc = p['tc'] + dt - self.end_time - self.ta[ifo]
            # generate template and calculate waveform ratio
            r0, r1 = self.waveform_ratio(p, htf, dtc=dtc)
            # computer inner products
            lr = self.binned_inner_products(r0, r1, self.sdat[ifo])
            # increment loglikelihood
            loglike += lr

        return float(loglike)

    def _loglr(self):
        return self.loglikelihood - self.lognl

    def _lognl(self):
        return sum([self.det_lognl(det) for det in self.data])

    def det_lognorm(self, det):
        p = self._psds[det]
        dt = self._data[det].delta_t
        kmin, kmax = self._kmin[det], self._kmax[det]
        lognorm = -float(self._N*numpy.log(numpy.pi*self._N*dt)/2. \
                         + numpy.log(p[kmin:kmax]).sum())
        self._lognorm[det] = lognorm
        return self._lognorm[det]

    def det_lognl(self, det):
        try:
            return self._det_lognls[det]
        except KeyError:
            # hasn't been calculated yet; calculate & store
            kmin = self._kmin[det]
            kmax = self._kmax[det]
            d = self._data[det]
            lognorm = self.det_lognorm(det)
            lognl = lognorm - 0.5 * d[kmin:kmax].inner(d[kmin:kmax]).real
            self._det_lognls[det] = lognl
            return self._det_lognls[det]

    def set_psds(self, psds):
        """Sets the psds, and calculates the weight and norm from them.
        The data and the low and high frequency cutoffs must be set first.
        """
        # check that the data has been set
        if self._data is None:
            raise ValueError("No data set")
        if self._f_lower is None:
            raise ValueError("low frequency cutoff not set")
        if self._f_upper is None:
            raise ValueError("high frequency cutoff not set")
        # make sure the relevant caches are cleared
        self._psds.clear()
        self._weight.clear()
        self._lognorm.clear()
        self._det_lognls.clear()
        for det, d in self._data.items():
            if psds is None:
                # No psd means assume white PSD
                p = FrequencySeries(numpy.ones(int(self._N/2+1)),
                                    delta_f=d.delta_f)
            else:
                # copy for storage
                p = psds[det].copy()
            self._psds[det] = p
            # we'll store the weight to apply to the inner product
            w = Array(numpy.zeros(len(p)))
            # only set weight in band we will analyze
            kmin = self._kmin[det]
            kmax = self._kmax[det]
            w[kmin:kmax] = numpy.sqrt(4.*p.delta_f/p[kmin:kmax])
            self._weight[det] = w
        # set the lognl and lognorm; we'll get this by just calling lognl
        _ = self.lognl

    def write_metadata(self, fp):
        """Adds writing the psds and lognl, since it's a constant.
        The lognl is written to the sample group's ``attrs``.
        Parameters
        ----------
        fp : pycbc.inference.io.BaseInferenceFile instance
            The inference file to write to.
        """
        super(Relative, self).write_metadata(fp)
        if self._psds is not None:
            fp.write_psd(self._psds)
        try:
            attrs = fp[fp.samples_group].attrs
        except KeyError:
            # group doesn't exist, create it
            fp.create_group(fp.samples_group)
            attrs = fp[fp.samples_group].attrs
        attrs['lognl'] = self.lognl
        for det in self.detectors:
            # Save lognl for each IFO as attributes in the samples group
            attrs['{}_lognl'.format(det)] = self.det_lognl(det)
            # Save each IFO's low frequency cutoff used in the likelihood
            # computation as an attribute
            fp.attrs['{}_likelihood_low_freq'.format(det)] = self._f_lower[det]
            # Save the IFO's high frequency cutoff used in the likelihood
            # computation as an attribute if one was provided the user
            if self._f_upper[det] is not None:
                fp.attrs['{}_likelihood_high_freq'.format(det)] = \
                                                        self._f_upper[det]

#!/usr/bin/env python

import numpy as np
import logging
from scipy.interpolate import interp1d

from pycbc.waveform.spa_tmplt import spa_tmplt
from pycbc.detector import Detector
from pycbc.types import Array

from .base_data import BaseDataModel


# construct frequency bins for relative binning
def setup_bins(f_full, f_lo, f_hi, chi=1.0, eps=0.5):
    """
    construct frequency binning
    f_full: full frequency grid
    [f_lo, f_hi] is the frequency range one would like to use for matched filtering
    chi, eps are tunable parameters [see Barak, Dai & Venumadhav 2018]
    return the number of bins, a list of bin-edge frequencies, and their positions in the full frequency grid
    """
    f = np.linspace(f_lo, f_hi, 10000)
    # f^ga power law index
    ga = np.array([-5.0/3.0, -2.0/3.0, 1.0, 5.0/3.0, 7.0/3.0])
    dalp = chi*2.0*np.pi/np.absolute(f_lo**ga - f_hi**ga)
    dphi = np.sum(np.array([np.sign(ga[i])*dalp[i]*f**ga[i] for i in range(len(ga))]), axis=0)
    Dphi = dphi - dphi[0]
    # now construct frequency bins
    Nbin = int(Dphi[-1]//eps)
    Dphi2f = interp1d(Dphi, f, kind='slinear', bounds_error=False, fill_value=0.0)
    Dphi_grid = np.linspace(Dphi[0], Dphi[-1], Nbin+1)
    # frequency grid points
    fbin = Dphi2f(Dphi_grid)
    # indices of frequency grid points in the FFT array
    fbin_ind = np.array([np.argmin(np.absolute(f_full - ff)) for ff in fbin])
    # make sure grid points are precise
    fbin = np.array([f_full[i] for i in fbin_ind])

    return (Nbin, fbin, fbin_ind)

class Relative(BaseDataModel):
    name = "relative"
    def __init__(self, data, psds, mass1, mass2, spin1z, spin2z,
                 ra, dec, tc, low_frequency_cutoff,
                 high_frequency_cutoff, epsilon=0.5, **kwargs):
        super(Relative, self).__init__(data=data, **kwargs)
        # define things
        f_lo = float(low_frequency_cutoff)
        f_hi = float(high_frequency_cutoff)
        self.f = np.array(data[data.keys()[0]].sample_frequencies)
        self.df = data[data.keys()[0]].delta_f
        self.fref = f_lo
        self.data = data.copy()
        # cast data and psds to arrays for faster computation
        self.comp_data = {ifo: np.array(data[ifo].data) for ifo in data}
        self.comp_psds = {ifo: np.array(psds[ifo].data) for ifo in data}
        self.det = {ifo: Detector(ifo) for ifo in data}
        mass1 = float(mass1)
        mass2 = float(mass2)
        spin1z = float(spin1z)
        spin2z = float(spin2z)
        self.ra = float(ra)
        self.dec = float(dec)
        self.tc = float(tc)
        epsilon = float(epsilon)
        self.end_time = float(data[data.keys()[0]].end_time)

        # get detector-specific arrival times relative to end of data
        dt = {ifo: self.det[ifo].time_delay_from_earth_center(
                  self.ra, self.dec, self.tc) for ifo in data}
        self.tc_loc = {ifo: self.tc + dt[ifo] - self.end_time
                       for ifo in data}
        # generate fiducial waveform
        logging.info("Generating fiducial waveform")
        hp = spa_tmplt(f_lower=f_lo, f_upper=f_hi+self.df,
                       delta_f=self.df, mass1=mass1,
                       mass2=mass2, spin1z=spin1z,
                       spin2z=spin2z, distance=1.,
                       spin_order=-1, phase_order=-1)
        hp.resize(len(self.f))
        self.h00 = np.array(hp)
        # compute frequency bins
        logging.info("Computing frequency bins")
        nbin, fbin, fbin_ind = setup_bins(f_full=self.f, f_lo=f_lo, f_hi=f_hi,
                                          eps=epsilon)
        self.edges = fbin_ind
        self.fedges = np.array(fbin).astype(np.float64)
        self.fedges32 = np.array(fbin).astype(np.float32)
        self.bins = np.array([(self.edges[i], self.edges[i+1]) for
                                  i in range(len(self.edges) - 1)])
        self.fbins = np.array([(fbin[i], fbin[i+1]) for
                                i in range(len(fbin) - 1)])
        self.h00_sparse = self.h00.copy().take(self.edges)
        logging.info('Using %s bins for this model', len(self.bins))

        # compute summary data
        logging.info("Calculating summary data")
        self.sdat = self.summary_data()

        # store the psds and calculate the inner product weight
        logging.info("Calculating lognl")
        kmin, kmax = self.edges[0], self.edges[-1]
        self._kmin = {ifo: kmin for ifo in self.data}
        self._kmax = {ifo: kmax for ifo in self.data}
        self._data = data.copy()
        self._f_lower = {ifo: f_lo for ifo in self.data}
        self._f_upper = {ifo: f_hi for ifo in self.data}
        self._N = int(1./(data[data.keys()[0]].delta_t * self.df))
        self._psds = {}
        self._weight = {}
        self._lognorm = {}
        self._det_lognls = {}
        self.set_psds(psds)
        logging.info("Lognl is {}".format(self.lognl))

    def summary_data(self):
        # timeshift the fiducial waveform for each detector
        h0 = {ifo: self.h00.copy() * np.exp(-2.0j * np.pi * self.f * self.tc_loc[ifo])
              for ifo in self.data}
        # calculate coefficients
        logging.info("Using {} seconds of data".format(1./self.df))
        sdat = {}
        for ifo in self.data:
            a0 = np.array([4.*self.df*np.sum(np.conjugate(self.comp_data[ifo][l:h]) \
                                     *h0[ifo][l:h] \
                                     /self.comp_psds[ifo][l:h]) for l, h in self.bins])
            b0 = np.array([4.*self.df*np.sum(np.absolute(h0[ifo][l:h]) ** 2.0 \
                                     /self.comp_psds[ifo][l:h]) for l, h in self.bins])
            a1 = np.array([4.*self.df*np.sum(np.conjugate(self.comp_data[ifo][l:h]) \
                                     *h0[ifo][l:h] \
                                     /self.comp_psds[ifo][l:h] \
                                     *(self.f[l:h] - 0.5*(fl+fh))) for (l, h), (fl, fh) in zip(self.bins, self.fbins)])
            b1 = np.array([4.*self.df*np.sum(np.absolute(h0[ifo][l:h]) ** 2.0 \
                                     /self.comp_psds[ifo][l:h] \
                                     *(self.f[l:h] - 0.5*(fl+fh))) for (l, h), (fl, fh) in zip(self.bins, self.fbins)])
            sdat[ifo] = {'a0': a0, 'a1': a1,
                         'b0': b0, 'b1': b1}
        return sdat

    def waveform_ratio(self, p, htf, dtc=0.0):
        # generate template
        s1z = p['spin1z']
        s2z = p['spin2z']
        hp = spa_tmplt(sample_points=self.fedges32,
                       mass1=p['mass1'],
                       mass2=p['mass2'], spin1z=s1z,
                       spin2z=s2z, distance=1.,
                       spin_order=-1, phase_order=-1)
        htarget = np.array(hp)
        # apply antenna pattern, inclination, and distance
        htarget *= htf
        # compute waveform ratio and timeshift
        r = htarget / self.h00_sparse * np.exp(-2.0j * np.pi * self.fedges * dtc)
        r0 = 0.5 * (r[:-1] + r[1:])
        r1 = (r[1:] - r[:-1]) / (self.fedges[1:] - self.fedges[:-1])
        return np.array([r0, r1], dtype=np.complex128)

    def calc_inner_products(self, r0, r1, sdata):
        # <h, d> is the sum over bins of A0r0 + A1r1
        hd = np.sum(sdata['a0']*r0 + sdata['a1']*r1)
        # <h, h> is the sum over bins of B0|r0|^2 + 2B1Re(r1r0*)
        hh = np.sum(sdata['b0']*np.absolute(r0)**2.
                    + 2.*sdata['b1']*(r1*np.conjugate(r0)).real).real

        return (hd - 0.5*hh).real

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
            ip = np.cos(p['inclination'])
            ic = 0.5 * (1.0 + ip * ip)
            htf = (fp * ip + 1.0j * fc * ic) / p['distance']
            # get timeshift relative to fiducial waveform
            dt = self.det[ifo].time_delay_from_earth_center(p['ra'], p['dec'],
                                                            p['tc'])
            dtc = p['tc'] + dt - self.end_time - self.tc_loc[ifo]
            # generate template and calculate waveform ratio
            rdata = self.waveform_ratio(p, htf, dtc=dtc)
            # computer inner products
            lr = self.calc_inner_products(rdata[0], rdata[1], self.sdat[ifo])
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
        lognorm = -float(self._N*np.log(np.pi*self._N*dt)/2. \
                         + np.log(p[kmin:kmax]).sum())
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
                p = FrequencySeries(np.ones(int(self._N/2+1)),
                                    delta_f=d.delta_f)
            else:
                # copy for storage
                p = psds[det].copy()
            self._psds[det] = p
            # we'll store the weight to apply to the inner product
            w = Array(np.zeros(len(p)))
            # only set weight in band we will analyze
            kmin = self._kmin[det]
            kmax = self._kmax[det]
            w[kmin:kmax] = np.sqrt(4.*p.delta_f/p[kmin:kmax])
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

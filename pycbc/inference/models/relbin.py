# Copyright (C) 2020  Daniel Finstad
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
"""This module provides model classes and functions for implementing
a relative binning likelihood for parameter estimation.
"""


import logging
import numpy
from scipy.interpolate import interp1d

from pycbc.waveform import get_fd_waveform_sequence
from pycbc.detector import Detector
from pycbc.types import Array

from .gaussian_noise import BaseGaussianNoise
from .relbin_cpu import likelihood_parts, likelihood_parts_v
from .tools import DistMarg


def setup_bins(f_full, f_lo, f_hi, chi=1.0, eps=0.5, gammas=None):
    """Construct frequency bins for use in a relative likelihood
    model. For details, see [Barak, Dai & Venumadhav 2018].

    Parameters
    ----------
    f_full : array
        The full resolution array of frequencies being used in the analysis.
    f_lo : float
        The starting frequency used in matched filtering. This will be the
        left edge of the first frequency bin.
    f_hi : float
        The ending frequency used in matched filtering. This will be the right
        edge of the last frequency bin.
    chi : float, optional
        Tunable parameter, see [Barak, Dai & Venumadhav 2018]
    eps : float, optional
        Tunable parameter, see [Barak, Dai & Venumadhav 2018]. Lower values
        result in larger number of bins.
    gammas : array, optional
        Frequency powerlaw indices to be used in computing bins.

    Returns
    -------
    nbin : int
        Number of bins.
    fbin : numpy.array of floats
        Bin edge frequencies.
    fbin_ind : numpy.array of ints
        Indices of bin edges in full frequency array.
    """
    f = numpy.linspace(f_lo, f_hi, 10000)
    # f^ga power law index
    ga = (
        gammas
        if gammas is not None
        else numpy.array([-5.0 / 3, -2.0 / 3, 1.0, 5.0 / 3, 7.0 / 3])
    )
    logging.info("Using powerlaw indices: %s", ga)
    dalp = chi * 2.0 * numpy.pi / numpy.absolute((f_lo ** ga) - (f_hi ** ga))
    dphi = numpy.sum(
        numpy.array([numpy.sign(g) * d * (f ** g) for g, d in zip(ga, dalp)]),
        axis=0,
    )
    dphi_diff = dphi - dphi[0]
    # now construct frequency bins
    nbin = int(dphi_diff[-1] / eps)
    dphi2f = interp1d(
        dphi_diff, f, kind="slinear", bounds_error=False, fill_value=0.0
    )
    dphi_grid = numpy.linspace(dphi_diff[0], dphi_diff[-1], nbin + 1)
    # frequency grid points
    fbin = dphi2f(dphi_grid)
    # indices of frequency grid points in the FFT array
    fbin_ind = numpy.unique(
        [numpy.argmin(numpy.absolute(f_full - ff)) for ff in fbin]
    )
    # make sure grid points are precise
    fbin = numpy.array([f_full[i] for i in fbin_ind])
    nbin = len(fbin)
    return nbin, fbin, fbin_ind


class Relative(BaseGaussianNoise, DistMarg):
    r"""Model that assumes the likelihood in a region around the peak
    is slowly varying such that a linear approximation can be made, and
    likelihoods can be calculated at a coarser frequency resolution. For
    more details on the implementation, see https://arxiv.org/abs/1806.08792.

    This model requires the use of a fiducial waveform whose parameters are
    near the peak of the likelihood. The fiducial waveform and all template
    waveforms used in likelihood calculation are currently generated using
    the SPAtmplt approximant.

    For more details on initialization parameters and definition of terms, see
    :py:class:`BaseGaussianNoise`.

    Parameters
    ----------
    variable_params : (tuple of) string(s)
        A tuple of parameter names that will be varied.
    data : dict
        A dictionary of data, in which the keys are the detector names and the
        values are the data (assumed to be unwhitened). All data must have the
        same frequency resolution.
    low_frequency_cutoff : dict
        A dictionary of starting frequencies, in which the keys are the
        detector names and the values are the starting frequencies for the
        respective detectors to be used for computing inner products.
    figucial_params : dict
        A dictionary of waveform parameters to be used for generating the
        fiducial waveform. Keys must be parameter names in the form
        'PARAM_ref' where PARAM is a recognized extrinsic parameter or
        an intrinsic parameter compatible with the chosen approximant.
    gammas : array of floats, optional
        Frequency powerlaw indices to be used in computing frequency bins.
    epsilon : float, optional
        Tuning parameter used in calculating the frequency bins. Lower values
        will result in higher resolution and more bins.
    earth_rotation: boolean, optional
        Default is False. If True, then vary the fp/fc polarization values
        as a function of frequency bin, using a predetermined PN approximation
        for the time offsets.
    \**kwargs :
        All other keyword arguments are passed to
        :py:class:`BaseGaussianNoise`.
    """
    name = "relative"

    def __init__(
        self,
        variable_params,
        data,
        low_frequency_cutoff,
        fiducial_params=None,
        gammas=None,
        epsilon=0.5,
        earth_rotation=False,
        marginalize_phase=True,
        **kwargs
    ):

        variable_params, kwargs = self.setup_distance_marginalization(
                               variable_params,
                               marginalize_phase=marginalize_phase,
                               **kwargs)

        super(Relative, self).__init__(
            variable_params, data, low_frequency_cutoff, **kwargs
        )

        self.epsilon = float(epsilon)

        # reference waveform and bin edges
        self.h00, self.h00_sparse = {}, {}
        self.f, self.df, self.end_time, self.det = {}, {}, {}, {}
        self.edges, self.fedges, self.bins, self.fbins = {}, {}, {}, {}
        self.ta = {}
        self.antenna_time = {}

        # filtered summary data for linear approximation
        self.sdat = {}

        # store data and psds as arrays for faster computation
        self.comp_data = {ifo: d.numpy() for ifo, d in self.data.items()}
        self.comp_psds = {ifo: p.numpy() for ifo, p in self.psds.items()}

        # store fiducial waveform params
        self.fid_params = self.static_params.copy()
        self.fid_params.update(fiducial_params)

        for ifo in data:
            # store data and frequencies
            d0 = self.data[ifo]
            self.f[ifo] = numpy.array(d0.sample_frequencies)
            self.df[ifo] = d0.delta_f
            self.end_time[ifo] = float(d0.end_time)
            self.det[ifo] = Detector(ifo)

            # get detector-specific arrival times relative to end of data
            dt = self.det[ifo].time_delay_from_earth_center(
                self.fid_params["ra"],
                self.fid_params["dec"],
                self.fid_params["tc"],
            )

            self.ta[ifo] = self.fid_params["tc"] + dt - self.end_time[ifo]

            # generate fiducial waveform
            f_lo = self.kmin[ifo] * self.df[ifo]
            f_hi = self.kmax[ifo] * self.df[ifo]
            logging.info(
                "%s: Generating fiducial waveform from %s to %s Hz",
                ifo,
                f_lo,
                f_hi,
            )

            # prune low frequency samples to avoid waveform errors
            fpoints = Array(self.f[ifo].astype(numpy.float64))
            fpoints = fpoints[self.kmin[ifo]:self.kmax[ifo]+1]
            fid_hp, fid_hc = get_fd_waveform_sequence(sample_points=fpoints,
                                                      **self.fid_params)

            # check for zeros at high frequencies
            # make sure only nonzero samples are included in bins
            numzeros = list(fid_hp[::-1] != 0j).index(True)
            if numzeros > 0:
                new_kmax = self.kmax[ifo] - numzeros
                f_hi = new_kmax * self.df[ifo]
                logging.info(
                    "WARNING! Fiducial waveform terminates below "
                    "high-frequency-cutoff, final bin frequency "
                    "will be %s Hz", f_hi)

            # make copy of fiducial wfs, adding back in low frequencies
            fid_hp.resize(len(self.f[ifo]))
            fid_hc.resize(len(self.f[ifo]))
            hp0 = numpy.roll(fid_hp, self.kmin[ifo])
            hc0 = numpy.roll(fid_hc, self.kmin[ifo])

            fp, fc = self.det[ifo].antenna_pattern(
                self.fid_params["ra"], self.fid_params["dec"],
                self.fid_params["polarization"], self.fid_params["tc"])

            tshift = numpy.exp(-2.0j * numpy.pi * self.f[ifo] * self.ta[ifo])
            self.h00[ifo] = numpy.array(hp0 * fp + hc0 * fc) * tshift

            # compute frequency bins
            logging.info("Computing frequency bins")
            nbin, fbin, fbin_ind = setup_bins(
                f_full=self.f[ifo],
                f_lo=f_lo,
                f_hi=f_hi,
                gammas=gammas,
                eps=self.epsilon,
            )
            logging.info("Using %s bins for this model", nbin)

            # store bins and edges in sample and frequency space
            self.edges[ifo] = fbin_ind
            self.fedges[ifo] = numpy.array(fbin).astype(numpy.float64)
            self.bins[ifo] = numpy.array(
                [
                    (self.edges[ifo][i], self.edges[ifo][i + 1])
                    for i in range(len(self.edges[ifo]) - 1)
                ]
            )
            self.fbins[ifo] = numpy.array(
                [(fbin[i], fbin[i + 1]) for i in range(len(fbin) - 1)]
            )

            # store low res copy of fiducial waveform
            self.h00_sparse[ifo] = self.h00[ifo].copy().take(self.edges[ifo])

            # compute summary data
            logging.info(
                "Calculating summary data at frequency resolution %s Hz",
                self.df[ifo],
            )
            self.sdat[ifo] = self.summary_data(ifo)

            # Calculate the times to evaluate fp/fc
            if earth_rotation is not False:
                logging.info("Enabling frequency-dependent earth rotation")
                from pycbc.waveform.spa_tmplt import spa_length_in_time

                times = spa_length_in_time(
                    phase_order=-1,
                    mass1=self.fid_params["mass1"],
                    mass2=self.fid_params["mass2"],
                    f_lower=self.fedges[ifo],
                )
                self.antenna_time[ifo] = self.fid_params["tc"] - times
                self.lik = likelihood_parts_v
            else:
                self.antenna_time[ifo] = self.fid_params["tc"]
                self.lik = likelihood_parts

        # determine the unique ifo layouts
        self.edge_unique = []
        self.ifo_map = {}
        for ifo in self.fedges:
            if len(self.edge_unique) == 0:
                self.ifo_map[ifo] = 0
                self.edge_unique.append(Array(self.fedges[ifo]))
            else:
                for i, edge in enumerate(self.edge_unique):
                    if numpy.array_equal(edge, self.fedges[ifo]):
                        self.ifo_map[ifo] = i
                        break
                else:
                    self.ifo_map[ifo] = len(self.edge_unique)
                    self.edge_unique.append(Array(self.fedges[ifo]))
        logging.info("%s unique ifo layouts", len(self.edge_unique))

    def summary_data(self, ifo):
        """Compute summary data bin coefficients encoding linear
        approximation to full resolution likelihood.

        Returns
        -------
        dict
            Dictionary containing bin coefficients a0, b0, a1, b1,
            for each frequency bin.
        """
        # calculate coefficients
        hd = numpy.conjugate(self.comp_data[ifo]) * self.h00[ifo]
        hd /= self.comp_psds[ifo]
        hh = (numpy.absolute(self.h00[ifo]) ** 2.0) / self.comp_psds[ifo]

        # constant terms
        a0 = numpy.array([
                4.0 * self.df[ifo] * numpy.sum(hd[l:h])
                for l, h in self.bins[ifo]
            ])
        b0 = numpy.array([
                4.0 * self.df[ifo] * numpy.sum(hh[l:h])
                for l, h in self.bins[ifo]
            ])

        # linear terms
        bin_lefts = [fl for fl, fh in self.fbins[ifo]]
        a1 = numpy.array([
                4.0 * self.df[ifo]
                * numpy.sum(hd[l:h] * (self.f[ifo][l:h] - bl))
                for (l, h), bl in zip(self.bins[ifo], bin_lefts)
            ])
        b1 = numpy.array([
                4.0 * self.df[ifo]
                * numpy.sum(hh[l:h] * (self.f[ifo][l:h] - bl))
                for (l, h), bl in zip(self.bins[ifo], bin_lefts)
            ])

        freqs = self.fedges[ifo]
        df = (freqs[1:] - freqs[:-1])
        a1 /= df
        b1 /= df
        return {"a0": a0, "a1": a1, "b0": b0, "b1": b1}

    def get_waveforms(self, params):
        """ Get the waveform polarizations for each ifo
        """
        wfs = []
        for edge in self.edge_unique:
            hp, hc = get_fd_waveform_sequence(sample_points=edge, **params)
            hp = hp.numpy()
            hc = hc.numpy()
            wfs.append((hp, hc))
        return {ifo: wfs[self.ifo_map[ifo]] for ifo in self.data}

    def _loglr(self):
        r"""Computes the log likelihood ratio,

        .. math::

            \log \mathcal{L}(\Theta) = \sum_i
                \left<h_i(\Theta)|d_i\right> -
                \frac{1}{2}\left<h_i(\Theta)|h_i(\Theta)\right>,

        at the current parameter values :math:`\Theta`.

        Returns
        -------
        float
            The value of the log likelihood ratio.
        """
        # get model params
        p = self.current_params.copy()
        p.update(self.static_params)
        wfs = self.get_waveforms(p)

        hh = 0.0
        hd = 0j
        for ifo in self.data:

            det = self.det[ifo]
            freqs = self.fedges[ifo]
            sdat = self.sdat[ifo]
            hp, hc = wfs[ifo]
            h00 = self.h00_sparse[ifo]
            end_time = self.end_time[ifo]
            times = self.antenna_time[ifo]

            # project waveform to detector frame
            fp, fc = det.antenna_pattern(p["ra"], p["dec"],
                                         p["polarization"], times)
            dt = det.time_delay_from_earth_center(p["ra"], p["dec"], times)
            dtc = p["tc"] + dt - end_time

            hdp, hhp = self.lik(freqs, fp, fc, dtc,
                                hp, hc, h00,
                                sdat['a0'], sdat['a1'],
                                sdat['b0'], sdat['b1'])

            hd += hdp
            hh += hhp

        return self.marginalize_loglr(hd, hh)

    def write_metadata(self, fp, group=None):
        """Adds writing the fiducial parameters and epsilon to file's attrs.

        Parameters
        ----------
        fp : pycbc.inference.io.BaseInferenceFile instance
            The inference file to write to.
        group : str, optional
            If provided, the metadata will be written to the attrs specified
            by group, i.e., to ``fp[group].attrs``. Otherwise, metadata is
            written to the top-level attrs (``fp.attrs``).
        """
        super().write_metadata(fp, group=group)
        if group is None:
            attrs = fp.attrs
        else:
            attrs = fp[group].attrs
        attrs["epsilon"] = self.epsilon
        for p, v in self.fid_params.items():
            attrs["{}_ref".format(p)] = v

    @staticmethod
    def extra_args_from_config(cp, section, skip_args=None, dtypes=None):
        """Adds reading fiducial waveform parameters from config file."""
        # add fiducial params to skip list
        skip_args += [
            option for option in cp.options(section) if option.endswith("_ref")
        ]

        # get frequency power-law indices if specified
        # NOTE these should be supplied in units of 1/3
        gammas = None
        if cp.has_option(section, "gammas"):
            skip_args.append("gammas")
            gammas = numpy.array(
                [float(g) / 3.0 for g in cp.get(section, "gammas").split()]
            )
        args = super(Relative, Relative).extra_args_from_config(
            cp, section, skip_args=skip_args, dtypes=dtypes
        )

        # get fiducial params from config
        fid_params = {
            p.replace("_ref", ""): float(cp.get("model", p))
            for p in cp.options("model")
            if p.endswith("_ref")
        }

        # add optional params with default values if not specified
        opt_params = {
            "ra": numpy.pi,
            "dec": 0.0,
            "inclination": 0.0,
            "polarization": numpy.pi,
        }
        fid_params.update(
            {p: opt_params[p] for p in opt_params if p not in fid_params}
        )
        args.update({"fiducial_params": fid_params, "gammas": gammas})
        return args

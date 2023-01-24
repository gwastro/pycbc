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
import itertools
from scipy.interpolate import interp1d

from pycbc.waveform import (get_fd_waveform_sequence,
                            get_fd_det_waveform_sequence, fd_det_sequence)
from pycbc.detector import Detector
from pycbc.types import Array

from .gaussian_noise import BaseGaussianNoise
from .relbin_cpu import (likelihood_parts, likelihood_parts_v,
                         likelihood_parts_multi, likelihood_parts_multi_v,
                         likelihood_parts_det)
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
    fbin_ind = numpy.searchsorted(f_full, fbin)
    for idx_fbin, idx_f_full in enumerate(fbin_ind):
        if idx_f_full == 0:
            curr_idx = 0
        elif idx_f_full == len(f_full):
            curr_idx = len(f_full) - 1
        else:
            abs1 = abs(f_full[idx_f_full] - fbin[idx_fbin])
            abs2 = abs(f_full[idx_f_full-1] - fbin[idx_fbin])
            if abs1 > abs2:
                curr_idx = idx_f_full - 1
            else:
                curr_idx = idx_f_full
        fbin_ind[idx_fbin] = curr_idx
    fbin_ind = numpy.unique(fbin_ind)
    return fbin_ind


class Relative(DistMarg, BaseGaussianNoise):
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

        variable_params, kwargs = self.setup_marginalization(
                               variable_params,
                               marginalize_phase=marginalize_phase,
                               **kwargs)

        super(Relative, self).__init__(
            variable_params, data, low_frequency_cutoff, **kwargs
        )

        # If the waveform handles the detector response internally, set
        # self.det_response = True
        self.no_det_response = False
        if self.static_params['approximant'] in fd_det_sequence:
            self.no_det_response = True

        # reference waveform and bin edges
        self.f, self.df, self.end_time, self.det = {}, {}, {}, {}
        self.h00, self.h00_sparse = {}, {}
        self.fedges, self.edges = {}, {}
        self.ta = {}
        self.antenna_time = {}

        # filtered summary data for linear approximation
        self.sdat = {}

        # store fiducial waveform params
        self.fid_params = self.static_params.copy()
        self.fid_params.update(fiducial_params)

        for k in self.static_params:
            if self.fid_params[k] == 'REPLACE':
               self.fid_params.pop(k)

        for ifo in data:
            # store data and frequencies
            d0 = self.data[ifo]
            self.f[ifo] = numpy.array(d0.sample_frequencies)
            self.df[ifo] = d0.delta_f
            self.end_time[ifo] = float(d0.end_time)
            # self.det[ifo] = Detector(ifo)

            # generate fiducial waveform
            f_lo = self.kmin[ifo] * self.df[ifo]
            f_hi = self.kmax[ifo] * self.df[ifo]
            logging.info(
                "%s: Generating fiducial waveform from %s to %s Hz",
                ifo, f_lo, f_hi,
            )

            # prune low frequency samples to avoid waveform errors
            fpoints = Array(self.f[ifo].astype(numpy.float64))
            fpoints = fpoints[self.kmin[ifo]:self.kmax[ifo]+1]

            if self.no_det_response:
                wave = get_fd_det_waveform_sequence(ifos=ifo,
                                                    sample_points=fpoints,
                                                    **self.fid_params)
                curr_wav = wave[ifo]
            else:
                fid_hp, fid_hc = get_fd_waveform_sequence(sample_points=fpoints,
                                                          **self.fid_params)
                curr_wav = fid_hp

            # check for zeros at low and high frequencies
            # make sure only nonzero samples are included in bins
            numzeros_lo = list(curr_wav != 0j).index(True)
            if numzeros_lo > 0:
                new_kmin = self.kmin[ifo] + numzeros_lo
                f_lo = new_kmin * self.df[ifo]
                logging.info(
                    "WARNING! Fiducial waveform starts above "
                    "low-frequency-cutoff, initial bin frequency "
                    "will be %s Hz", f_lo)
            numzeros_hi = list(curr_wav[::-1] != 0j).index(True)
            if numzeros_hi > 0:
                new_kmax = self.kmax[ifo] - numzeros_hi
                f_hi = new_kmax * self.df[ifo]
                logging.info(
                    "WARNING! Fiducial waveform terminates below "
                    "high-frequency-cutoff, final bin frequency "
                    "will be %s Hz", f_hi)

            # make copy of fiducial wfs, adding back in low frequencies
            if self.no_det_response:
                curr_wav.resize(len(self.f[ifo]))
                curr_wav = numpy.roll(curr_wav, self.kmin[ifo])
                # get detector-specific arrival times relative to end of data
                self.ta[ifo] = -self.end_time[ifo]
                tshift = numpy.exp(-2.0j * numpy.pi * self.f[ifo] * self.ta[ifo])
                h00 = numpy.array(curr_wav) * tshift
                self.h00[ifo] = h00
            else:
                fid_hp.resize(len(self.f[ifo]))
                fid_hc.resize(len(self.f[ifo]))
                hp0 = numpy.roll(fid_hp, self.kmin[ifo])
                hc0 = numpy.roll(fid_hc, self.kmin[ifo])

                self.det[ifo] = Detector(ifo)
                dt = self.det[ifo].time_delay_from_earth_center(
                    self.fid_params["ra"],
                    self.fid_params["dec"],
                    self.fid_params["tc"],
                )
                self.ta = self.fid_params["tc"] + dt - self.end_time[ifo]

                fp, fc = self.det[ifo].antenna_pattern(
                    self.fid_params["ra"], self.fid_params["dec"],
                    self.fid_params["polarization"], self.fid_params["tc"])

                tshift = numpy.exp(-2.0j * numpy.pi * self.f[ifo] * self.ta)

                h00 = (hp0 * fp + hc0 * fc) * tshift
                self.h00[ifo] = h00

            # compute frequency bins
            logging.info("Computing frequency bins")
            fbin_ind = setup_bins(
                f_full=self.f[ifo], f_lo=f_lo, f_hi=f_hi,
                gammas=gammas, eps=float(epsilon),
            )
            logging.info("Using %s bins for this model", len(fbin_ind))

            self.fedges[ifo] = self.f[ifo][fbin_ind]
            self.edges[ifo] = fbin_ind

            self.init_from_frequencies(h00, fbin_ind, ifo)
            self.antenna_time[ifo] = self.setup_antenna(earth_rotation,
                                                        self.fedges[ifo])
        self.combine_layout()

    def init_from_frequencies(self, h00, fbin_ind, ifo):
        bins = numpy.array(
            [
                (fbin_ind[i], fbin_ind[i + 1])
                for i in range(len(fbin_ind) - 1)
            ]
        )

        # store low res copy of fiducial waveform
        self.h00_sparse[ifo] = h00.copy().take(fbin_ind)

        # compute summary data
        logging.info(
            "Calculating summary data at frequency resolution %s Hz",
            self.df[ifo],
        )

        a0, a1 = self.summary_product(self.data[ifo], h00, bins, ifo)
        b0, b1 = self.summary_product(h00, h00, bins, ifo)
        self.sdat[ifo] = {"a0": a0, "a1": a1, "b0": abs(b0), "b1": abs(b1)}

    def combine_layout(self):
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

    def setup_antenna(self, earth_rotation, fedges):
        # Calculate the times to evaluate fp/fc
        if earth_rotation is not False:
            logging.info("Enabling frequency-dependent earth rotation")
            from pycbc.waveform.spa_tmplt import spa_length_in_time

            times = spa_length_in_time(
                phase_order=-1,
                mass1=self.fid_params["mass1"],
                mass2=self.fid_params["mass2"],
                f_lower=fedges,
            )
            atimes = self.fid_params["tc"] - times
            self.lik = likelihood_parts_v
            self.mlik = likelihood_parts_multi_v
        else:
            atimes = self.fid_params["tc"]
            if self.no_det_response:
                self.lik = likelihood_parts_det
            else:
                self.lik = likelihood_parts
            self.mlik = likelihood_parts_multi
        return atimes

    def summary_product(self, h1, h2, bins, ifo):
        """ Calculate the summary values for the inner product <h1|h2>
        """
        # calculate coefficients
        h12 = numpy.conjugate(h1) * h2 / self.psds[ifo]

        # constant terms
        a0 = numpy.array([
                4.0 * self.df[ifo] * h12[l:h].sum()
                for l, h in bins
            ])

        # linear terms
        a1 = numpy.array([
                4.0 / (h - l) *
                (h12[l:h] * (self.f[ifo][l:h] - self.f[ifo][l])).sum()
                for l, h in bins])

        return a0, a1

    def get_waveforms(self, params):
        """ Get the waveform polarizations for each ifo
        """
        if self.no_det_response:
            wfs = {}
            for ifo in self.data:
                wfs.update(get_fd_det_waveform_sequence(
                        ifos=ifo, sample_points=self.fedges[ifo], **params))
            return wfs

        wfs = []
        for edge in self.edge_unique:
            hp, hc = get_fd_waveform_sequence(sample_points=edge, **params)
            hp = hp.numpy()
            hc = hc.numpy()
            wfs.append((hp, hc))
        wf_ret = {ifo: wfs[self.ifo_map[ifo]] for ifo in self.data}
        return wf_ret

    @property
    def multi_signal_support(self):
        """ The list of classes that this model supports in a multi-signal
        likelihood
        """
        # Check if this model *can* be included in a multi-signal model.
        # All marginalizations must currently be disabled to work!
        if (self.marginalize_vector_params or
            self.marginalize_distance or
            self.marginalize_phase):
            logging.info("Cannot use single template model inside of"
                         "multi_signal if marginalizations are enabled")
        return [type(self)]

    def calculate_hihjs(self, models):
        """ Pre-calculate the hihj inner products on a grid
        """
        self.hihj = {}
        for m1, m2 in itertools.combinations(models, 2):
            self.hihj[(m1, m2)] = {}
            for ifo in self.data:
                h1 = m1.h00[ifo]
                h2 = m2.h00[ifo]

                # Combine the grids
                edge = numpy.unique([m1.edges[ifo], m2.edges[ifo]])

                # Remove any points where either reference is zero
                keep = numpy.where((h1[edge] != 0) | (h2[edge] != 0))[0]
                edge = edge[keep]
                fedge = m1.f[ifo][edge]

                bins = numpy.array([
                        (edge[i], edge[i + 1])
                        for i in range(len(edge) - 1)
                    ])
                a0, a1 = self.summary_product(h1, h2, bins, ifo)
                self.hihj[(m1, m2)][ifo] = a0, a1, fedge

    def multi_loglikelihood(self, models):
        """ Calculate a multi-model (signal) likelihood
        """
        models = [self] + models
        loglr = 0
        # handle sum[<d|h_i> - 0.5 <h_i|h_i>]
        for m in models:
            loglr += m.loglr

        if not hasattr(self, 'hihj'):
            self.calculate_hihjs(models)

        # finally add in the lognl term from this model
        for m1, m2 in itertools.combinations(models, 2):
            for det in self.data:
                a0, a1, fedge = self.hihj[(m1, m2)][det]

                fp, fc, dtc, hp, hc, h00 = m1._current_wf_parts[det]
                fp2, fc2, dtc2, hp2, hc2, h002 = m2._current_wf_parts[det]

                h1h2 = self.mlik(fedge,
                                 fp, fc, dtc, hp, hc, h00,
                                 fp2, fc2, dtc2, hp2, hc2, h002,
                                 a0, a1)
                loglr += - h1h2.real # This is -0.5 * re(<h1|h2> + <h2|h1>)
        return loglr + self.lognl

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

        norm = 0.0
        filt = 0j
        self._current_wf_parts = {}
        for ifo in self.data:

            freqs = self.fedges[ifo]
            sdat = self.sdat[ifo]
            h00 = self.h00_sparse[ifo]
            end_time = self.end_time[ifo]
            times = self.antenna_time[ifo]

            # project waveform to detector frame if waveform does not deal
            # with detector response. Otherwise, skip detector response.

            if self.no_det_response:
                dtc = -end_time
                channel = wfs[ifo].numpy()
                filter_i, norm_i = self.lik(freqs, dtc, channel, h00,
                                            sdat['a0'], sdat['a1'],
                                            sdat['b0'], sdat['b1'])
            else:
                hp, hc = wfs[ifo]
                det = self.det[ifo]
                fp, fc = det.antenna_pattern(p["ra"], p["dec"],
                                             p["polarization"], times)
                dt = det.time_delay_from_earth_center(p["ra"], p["dec"], times)
                dtc = p["tc"] + dt - end_time

                filter_i, norm_i = self.lik(freqs, fp, fc, dtc,
                                            hp, hc, h00,
                                            sdat['a0'], sdat['a1'],
                                            sdat['b0'], sdat['b1'])
                self._current_wf_parts[ifo] = (fp, fc, dtc, hp, hc, h00)
            filt += filter_i
            norm += norm_i
        return self.marginalize_loglr(filt, norm)

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

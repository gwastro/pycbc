# Copyright (C) 2018 Alex Nitz
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

"""This module provides model classes that assume the noise is Gaussian.
"""

import numpy

from pycbc import filter as pyfilter
from pycbc.waveform import NoWaveformError
from pycbc.types import Array

from .base_data import BaseDataModel

# In this model we only calculate terms up to a constant.
# We are primarily interested in the posterior result

class SingleTemplate(BaseDataModel):
    r"""Model that assumes we know all the intrinsic parameters.

    This model assumes we know all the intrinsic parameters, and are only
    maximizing over the extrinsic ones. We also assume a dominant mode waveform
    approximant only and non-precessing.
    """
    name = 'single_template'

    def __init__(self, variable_params, data, waveform_generator,
                 f_lower, psds=None, f_upper=None, norm=None,
                 **kwargs):
        # set up the boiler-plate attributes; note: we'll compute the
        # log evidence later
        super(GaussianNoise, self).__init__(variable_params, data,
                                            waveform_generator, **kwargs)
 
        pass

    def _lognl(self):
        """Computes the log likelihood assuming the data is noise.

        Since this is a constant for Gaussian noise, this is only computed once
        then stored.
        """
        # FIXME we won't bother calculating this for now since it is a 
        # constant in this model.
        return 0

    def _loglr(self):
        r"""Computes the log likelihood ratio,

        Returns
        -------
        float
            The value of the log likelihood ratio.
        """
        # calculate <d-h|d-h> = <h|h> - 2<h|d> + <d|d> up to a constant
        return 5.7

    def _loglikelihood(self):
        r"""Computes the log likelihood of the paramaters,

        """
        return self.loglr

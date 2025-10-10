#
#  Apapted from code in LALSimInpspiralTaylorF2.c
#
#  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
#  Copyright (C) 2012 Leo Singer, Alex Nitz
#  Adapted from code found in:
#    - LALSimInspiralTaylorF2.c
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with with program; see the file COPYING. If not, write to the
#  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA

import cupy as cp
import lal
import mako.template

taylorf2_text = mako.template.Template("""
    const float f = (i + kmin ) * delta_f;
    const float amp2 = amp * __powf(f, -7.0/6.0);
    const float v =  __powf(piM*f, 1.0/3.0);
    const float v2 = v * v;
    const float v3 = v2 * v;
    const float v4 = v2 * v2;
    const float v5 = v2 * v3;
    const float v6 = v3 * v3;
    const float v7 = v3 * v4;
    float phasing = 0.;

    float LAL_TWOPI = ${TWOPI};
    float LAL_PI_4 = ${PI_4};
    float log4 = ${LN4};
    float logv = __logf(v);

    switch (phase_order)
    {
        case -1:
        case 7:
            phasing += pfa7 * v7;
        case 6:
            phasing += (pfa6 + pfl6 * (logv + log4) ) * v6;
        case 5:
            phasing += (pfa5 + pfl5 * (logv) ) * v5;
        case 4:
            phasing += pfa4 * v4;
        case 3:
            phasing += pfa3 * v3;
        case 2:
            phasing += pfa2 * v2;
        case 0:
            phasing += 1.;
            break;
        default:
            break;
    }
    phasing *= pfaN / v5;
    phasing -=  LAL_PI_4;
    phasing -= int(phasing / (LAL_TWOPI)) * LAL_TWOPI;

    float pcos;
    float psin;
    __sincosf(phasing, &psin, &pcos);

    htilde.real(pcos * amp2);
    htilde.imag(-psin * amp2);
""").render(TWOPI=lal.TWOPI, PI_4=lal.PI_4, LN4=2*lal.LN2)


taylorf2_kernel = cp.ElementwiseKernel(
    """
        int64 kmin, int64 phase_order, float32 delta_f, float32 piM,
        float32 pfaN, float32 pfa2, float32 pfa3, float32 pfa4, float32 pfa5,
        float32 pfl5, float32 pfa6, float32 pfl6, float32 pfa7, float32 amp
    """,
    "complex64 htilde",
    taylorf2_text,
    "taylorf2_kernel",
)


def spa_tmplt_engine(htilde,  kmin,  phase_order,
                    delta_f,  piM,  pfaN,
                    pfa2,  pfa3,  pfa4,  pfa5,  pfl5,
                    pfa6,  pfl6,  pfa7, amp_factor):
    """ Calculate the spa tmplt phase
    """
    taylorf2_kernel(kmin,  phase_order,
                    delta_f,  piM,  pfaN,
                    pfa2,  pfa3,  pfa4,  pfa5,  pfl5,
                    pfa6,  pfl6,  pfa7, amp_factor, htilde.data)

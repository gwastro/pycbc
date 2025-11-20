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


# Batched version - processes multiple templates in one kernel call
taylorf2_batch_text = mako.template.Template("""
    // Determine which template and frequency bin we're computing
    int template_idx = i / freq_length;
    int freq_idx = i % freq_length;
    
    if (template_idx >= num_templates) return;
    
    // Check if this frequency index is valid for this template
    // freq_idx corresponds to the index within the template's frequency range
    // The actual frequency is (freq_idx + kmin) * delta_f
    // But we only want frequencies up to kmax
    const float f = (freq_idx + kmin[template_idx]) * delta_f;
    
    // Set output to zero if outside the template's frequency range
    if (freq_idx >= freq_spans[template_idx]) {
        htilde.real(0.0f);
        htilde.imag(0.0f);
        return;
    }
    
    const float amp2 = amp[template_idx] * __powf(f, -7.0/6.0);
    const float v = __powf(piM[template_idx] * f, 1.0/3.0);
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
    
    int po = phase_order[template_idx];

    switch (po)
    {
        case -1:
        case 7:
            phasing += pfa7[template_idx] * v7;
        case 6:
            phasing += (pfa6[template_idx] + pfl6[template_idx] * (logv + log4)) * v6;
        case 5:
            phasing += (pfa5[template_idx] + pfl5[template_idx] * logv) * v5;
        case 4:
            phasing += pfa4[template_idx] * v4;
        case 3:
            phasing += pfa3[template_idx] * v3;
        case 2:
            phasing += pfa2[template_idx] * v2;
        case 0:
            phasing += 1.;
            break;
        default:
            break;
    }
    phasing *= pfaN[template_idx] / v5;
    phasing -= LAL_PI_4;
    phasing -= int(phasing / (LAL_TWOPI)) * LAL_TWOPI;

    float pcos;
    float psin;
    __sincosf(phasing, &psin, &pcos);

    htilde.real(pcos * amp2);
    htilde.imag(-psin * amp2);
""").render(TWOPI=lal.TWOPI, PI_4=lal.PI_4, LN4=2*lal.LN2)


taylorf2_batch_kernel = cp.ElementwiseKernel(
    """
        raw int64 kmin, raw int32 freq_spans, raw int64 phase_order, float32 delta_f,
        raw float32 piM, raw float32 pfaN, raw float32 pfa2, raw float32 pfa3,
        raw float32 pfa4, raw float32 pfa5, raw float32 pfl5, raw float32 pfa6,
        raw float32 pfl6, raw float32 pfa7, raw float32 amp,
        int32 num_templates, int32 freq_length
    """,
    "complex64 htilde",
    taylorf2_batch_text,
    "taylorf2_batch_kernel",
)


def spa_tmplt_engine_batch(htilde_batch, kmin_arr, freq_spans_arr, phase_order_arr,
                          delta_f, piM_arr, pfaN_arr,
                          pfa2_arr, pfa3_arr, pfa4_arr, pfa5_arr, pfl5_arr,
                          pfa6_arr, pfl6_arr, pfa7_arr, amp_arr,
                          num_templates, freq_length):
    """
    Calculate spa tmplt phase for multiple templates in a single kernel call.
    
    Parameters
    ----------
    htilde_batch : cupy array
        Output array of shape (num_templates * freq_length,)
    kmin_arr : cupy array of int64
        Starting frequency index for each template
    freq_spans_arr : cupy array of int32
        Frequency span (kmax - kmin) for each template
    phase_order_arr : cupy array of int64
        Phase order for each template
    delta_f : float
        Frequency spacing
    piM_arr, pfaN_arr, pfa*_arr, pfl*_arr : cupy arrays of float32
        PN coefficients for each template
    amp_arr : cupy array of float32
        Amplitude factor for each template
    num_templates : int
        Number of templates to generate
    freq_length : int
        Maximum frequency span across all templates
    """
    taylorf2_batch_kernel(kmin_arr, freq_spans_arr, phase_order_arr,
                         delta_f,
                         piM_arr, pfaN_arr,
                         pfa2_arr, pfa3_arr, pfa4_arr, pfa5_arr, pfl5_arr,
                         pfa6_arr, pfl6_arr, pfa7_arr, amp_arr,
                         num_templates, freq_length,
                         htilde_batch)

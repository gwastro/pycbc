""" This module contains functions for calculating single-ifo ranking
statistic values
"""
import numpy


def effsnr(snr, reduced_x2, fac=250.):
    """Calculate the effective SNR statistic. See (S5y1 paper) for definition.
    """
    snr = numpy.array(snr, ndmin=1, dtype=numpy.float64)
    rchisq = numpy.array(reduced_x2, ndmin=1, dtype=numpy.float64)
    esnr = snr / (1 + snr ** 2 / fac) ** 0.25 / rchisq ** 0.25

    # If snr input is float, return a float. Otherwise return numpy array.
    if hasattr(snr, '__len__'):
        return esnr
    else:
        return esnr[0]


def newsnr(snr, reduced_x2, q=6., n=2.):
    """Calculate the re-weighted SNR statistic ('newSNR') from given SNR and
    reduced chi-squared values. See http://arxiv.org/abs/1208.3491 for
    definition. Previous implementation in glue/ligolw/lsctables.py
    """
    nsnr = numpy.array(snr, ndmin=1, dtype=numpy.float64)
    reduced_x2 = numpy.array(reduced_x2, ndmin=1, dtype=numpy.float64)

    # newsnr is only different from snr if reduced chisq > 1
    ind = numpy.where(reduced_x2 > 1.)[0]
    nsnr[ind] *= (0.5 * (1. + reduced_x2[ind] ** (q/n))) ** (-1./q)

    # If snr input is float, return a float. Otherwise return numpy array.
    if hasattr(snr, '__len__'):
        return nsnr
    else:
        return nsnr[0]


def newsnr_sgveto(snr, bchisq, sgchisq):
    """ Combined SNR derived from NewSNR and Sine-Gaussian Chisq"""
    nsnr = numpy.array(newsnr(snr, bchisq), ndmin=1)
    sgchisq = numpy.array(sgchisq, ndmin=1)
    t = numpy.array(sgchisq > 4, ndmin=1)
    if len(t):
        nsnr[t] = nsnr[t] / (sgchisq[t] / 4.0) ** 0.5

    # If snr input is float, return a float. Otherwise return numpy array.
    if hasattr(snr, '__len__'):
        return nsnr
    else:
        return nsnr[0]


def newsnr_sgveto_psdvar(snr, bchisq, sgchisq, psd_var_val):
    """ Combined SNR derived from NewSNR, Sine-Gaussian Chisq and PSD
    variation statistic """
    nsnr = numpy.array(newsnr_sgveto(snr, bchisq, sgchisq), ndmin=1)
    psd_var_val = numpy.array(psd_var_val, ndmin=1)
    lgc = psd_var_val >= 1.8
    nsnr[lgc] = nsnr[lgc] / numpy.sqrt(psd_var_val[lgc])

    # If snr input is float, return a float. Otherwise return numpy array.
    if hasattr(snr, '__len__'):
        return nsnr
    else:
        return nsnr[0]


def get_newsnr(trigs):
    """
    Calculate newsnr ('reweighted SNR') for a trigs object

    Parameters
    ----------
    trigs: dict of numpy.ndarrays, h5py group (or similar dict-like object)
        Dictionary-like object holding single detector trigger information.
        'chisq_dof', 'snr', and 'chisq' are required keys

    Returns
    -------
    numpy.ndarray
        Array of newsnr values
    """
    dof = 2. * trigs['chisq_dof'][:] - 2.
    nsnr = newsnr(trigs['snr'][:], trigs['chisq'][:] / dof)
    return numpy.array(nsnr, ndmin=1, dtype=numpy.float32)


def get_newsnr_sgveto(trigs):
    """
    Calculate newsnr re-weigthed by the sine-gaussian veto

    Parameters
    ----------
    trigs: dict of numpy.ndarrays, h5py group (or similar dict-like object)
        Dictionary-like object holding single detector trigger information.
        'chisq_dof', 'snr', 'sg_chisq' and 'chisq' are required keys

    Returns
    -------
    numpy.ndarray
        Array of newsnr values
    """
    dof = 2. * trigs['chisq_dof'][:] - 2.
    nsnr_sg = newsnr_sgveto(trigs['snr'][:],
                            trigs['chisq'][:] / dof,
                            trigs['sg_chisq'][:])
    return numpy.array(nsnr_sg, ndmin=1, dtype=numpy.float32)


def get_newsnr_sgveto_psdvar(trigs):
    """
    Calculate newsnr re-weighted by the sine-gaussian veto and psd variation
    statistic

    Parameters
    ----------
    trigs: dict of numpy.ndarrays
        Dictionary holding single detector trigger information.
    'chisq_dof', 'snr', 'chisq' and 'psd_var_val' are required keys

    Returns
    -------
     numpy.ndarray
        Array of newsnr values
    """
    dof = 2. * trigs['chisq_dof'][:] - 2.
    nsnr_sg_psd = \
                 newsnr_sgveto_psdvar(trigs['snr'][:], trigs['chisq'][:] / dof,
                                      trigs['sg_chisq'][:],
                                      trigs['psd_var_val'][:])
    return numpy.array(nsnr_sg_psd, ndmin=1, dtype=numpy.float32)

# Copyright (C) 2012  Alex Nitz
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
import lal
import numpy
import scipy.signal
from pycbc.types import TimeSeries

# LDAS low pass FIR filter coefficents for resampling by 2, 4 and 8      
# thse were generated using the firlp action to get the FIR coeffs       
# used by the resample action. FIR coeffs are provided for the default   
# resample action that used a Kaiser window with beta = 5 and a filter   
# order parameter n = 10. The order of the filter is 2 * n * resampRatio 
# The filter coefficents used were produced by LDAS-CIT running version 0.7.0
# of LDAS. See the LDAS dataconditioning API documentation for information.

LDAS_FIR_LP = {}

LDAS_FIR_LP[2] = numpy.array( \
    [-7.1587345178999983e-19, -1.0514587726686521e-03, 1.8542385840153159e-18,
    2.5089668121365313e-03, -3.4940298359108189e-18, -4.8948343392593557e-03,
    5.6062660030794468e-18, 8.5565590024681924e-03, -8.0907458797357440e-18,
    -1.3989973238439804e-02, 1.0780775955476128e-17, 2.2023120574050724e-02,
    -1.3459232902020823e-17, -3.4340179881765416e-02, 1.5884076076327108e-17,
    5.5288283911880384e-02, -1.7819618555760906e-17, -1.0093019739775573e-01,
    1.9068667053692095e-17, 3.1670034568032440e-01, 5.0025873529805742e-01,
    3.1670034568032440e-01, 1.9068667053692095e-17, -1.0093019739775573e-01,
    -1.7819618555760906e-17, 5.5288283911880384e-02, 1.5884076076327108e-17,
    -3.4340179881765416e-02, -1.3459232902020823e-17, 2.2023120574050724e-02,
    1.0780775955476128e-17, -1.3989973238439804e-02, -8.0907458797357440e-18,
    8.5565590024681924e-03, 5.6062660030794468e-18, -4.8948343392593557e-03,
    -3.4940298359108189e-18, 2.5089668121365313e-03, 1.8542385840153159e-18,
    -1.0514587726686521e-03, -7.1587345178999983e-19])

LDAS_FIR_LP[4] = numpy.array( \
    [-3.5797933214045194e-19, -2.8264939479410322e-04, -5.2579196542625766e-04,
    -4.7541698372288916e-04, 9.2722964970291765e-19, 7.3162852724022317e-04,
    1.2546327308624197e-03, 1.0630797667142953e-03, -1.7472228702023229e-18,
    -1.4830129609431080e-03, -2.4477084927857239e-03, -2.0065313653173720e-03,
    2.8034666665817767e-18, 2.6512725192001795e-03, 4.2787887572376141e-03,
    3.4384897676469164e-03, -4.0458544723286438e-18, -4.3947913127916887e-03,
    -6.9958192527421722e-03, -5.5549713919258352e-03, 5.3910296112353999e-18,
    6.9667290898765182e-03, 1.1012871024947616e-02, 8.6988136875756853e-03,
    -6.7304174967527627e-18, -1.0855771610536039e-02, -1.7172133746431371e-02,
    -1.3606372767203865e-02, 7.9429834019598979e-18, 1.7243084699552668e-02,
    2.7647432518244298e-02, 2.2320837133020830e-02, -8.9108698382911643e-18,
    -3.0033587538646083e-02, -5.0471105705777050e-02, -4.3494435742894542e-02,
    9.5354684261874058e-18, 7.4135704901104452e-02, 1.5836902171998732e-01,
    2.2490814275257559e-01, 2.5015914127230293e-01, 2.2490814275257559e-01,
    1.5836902171998732e-01, 7.4135704901104452e-02, 9.5354684261874058e-18,
    -4.3494435742894542e-02, -5.0471105705777050e-02, -3.0033587538646083e-02,
    -8.9108698382911643e-18, 2.2320837133020830e-02, 2.7647432518244298e-02,
    1.7243084699552668e-02, 7.9429834019598979e-18, -1.3606372767203865e-02,
    -1.7172133746431371e-02, -1.0855771610536039e-02, -6.7304174967527627e-18,
    8.6988136875756853e-03, 1.1012871024947616e-02, 6.9667290898765182e-03,
    5.3910296112353999e-18, -5.5549713919258352e-03, -6.9958192527421722e-03,
    -4.3947913127916887e-03, -4.0458544723286438e-18, 3.4384897676469164e-03,
    4.2787887572376141e-03, 2.6512725192001795e-03, 2.8034666665817767e-18,
    -2.0065313653173720e-03, -2.4477084927857239e-03, -1.4830129609431080e-03,
    -1.7472228702023229e-18, 1.0630797667142953e-03, 1.2546327308624197e-03,
    7.3162852724022317e-04, 9.2722964970291765e-19, -4.7541698372288916e-04,
    -5.2579196542625766e-04, -2.8264939479410322e-04, -3.5797933214045194e-19])

LDAS_FIR_LP[8] = numpy.array( \
    [-1.7899485045886187e-19, -6.5785565693739621e-05, -1.4132879082976897e-04,
    -2.1264395052204678e-04, -2.6290359742617416e-04, -2.7550349276906565e-04,
    -2.3771537702543715e-04, -1.4425544130102367e-04, 4.6362825333301403e-19,
    1.7887185358663154e-04, 3.6582485933411455e-04, 5.2717915649190441e-04,
    6.2733453548486414e-04, 6.3531304436406394e-04, 5.3155527927017054e-04,
    3.1369007704874736e-04, -8.7363673902677791e-19, -3.7044471629615059e-04,
    -7.4152795801187910e-04, -1.0475948732225806e-03, -1.2238896950094572e-03,
    -1.2184233937642017e-03, -1.0032947419855074e-03, -5.8331646493291016e-04,
    1.4017739341284856e-18, 6.7046275382674411e-04, 1.3256746563059471e-03,
    1.8513153796509275e-03, 2.1394563456147114e-03, 2.1081914808740686e-03,
    1.7192946812996687e-03, 9.9055884251036106e-04, -2.0229858297200544e-18,
    -1.1198042255458091e-03, -2.1974593033835155e-03, -3.0470761300287886e-03,
    -3.4980109424040981e-03, -3.4255709083983424e-03, -2.7775661451061263e-03,
    -1.5917261396652274e-03, 2.6955928804956133e-18, 1.7824568219118389e-03,
    3.4834654396768100e-03, 4.8124727581638303e-03, 5.5065950049312329e-03,
    5.3772206622082304e-03, 4.3495328232139681e-03, 2.4876883619736412e-03,
    -3.3653062207633390e-18, -2.7788423673236681e-03, -5.4280430225538204e-03,
    -7.4993135801407406e-03, -8.5863155663860897e-03, -8.3949478976005458e-03,
    -6.8033834361075291e-03, -3.9013002754523102e-03, 3.9716067341932858e-18,
    4.3911559143580796e-03, 8.6217920704845935e-03, 1.1985707616989357e-02,
    1.3824116659430464e-02, 1.3633243414968320e-02, 1.1160741825101029e-02,
    6.4757499732487761e-03, -4.4555639696470439e-18, -7.5079556703908984e-03,
    -1.5017228726807877e-02, -2.1332840162384223e-02, -2.5236283794044537e-02,
    -2.5641711624359045e-02, -2.1747847773897343e-02, -1.3165154002435106e-02,
    4.7678723092621354e-18, 1.7117561087684679e-02, 3.7068926111156336e-02,
    5.8343464299929801e-02, 7.9186804418539564e-02, 9.7784206656366959e-02,
    1.1245732857891018e-01, 1.2185045954177413e-01, 1.2508319353304170e-01,
    1.2185045954177413e-01, 1.1245732857891018e-01, 9.7784206656366959e-02,
    7.9186804418539564e-02, 5.8343464299929801e-02, 3.7068926111156336e-02,
    1.7117561087684679e-02, 4.7678723092621354e-18, -1.3165154002435106e-02,
    -2.1747847773897343e-02, -2.5641711624359045e-02, -2.5236283794044537e-02,
    -2.1332840162384223e-02, -1.5017228726807877e-02, -7.5079556703908984e-03,
    -4.4555639696470439e-18, 6.4757499732487761e-03, 1.1160741825101029e-02,
    1.3633243414968320e-02, 1.3824116659430464e-02, 1.1985707616989357e-02,
    8.6217920704845935e-03, 4.3911559143580796e-03, 3.9716067341932858e-18,
    -3.9013002754523102e-03, -6.8033834361075291e-03, -8.3949478976005458e-03,
    -8.5863155663860897e-03, -7.4993135801407406e-03, -5.4280430225538204e-03,
    -2.7788423673236681e-03, -3.3653062207633390e-18, 2.4876883619736412e-03,
    4.3495328232139681e-03, 5.3772206622082304e-03, 5.5065950049312329e-03,
    4.8124727581638303e-03, 3.4834654396768100e-03, 1.7824568219118389e-03,
    2.6955928804956133e-18, -1.5917261396652274e-03, -2.7775661451061263e-03,
    -3.4255709083983424e-03, -3.4980109424040981e-03, -3.0470761300287886e-03,
    -2.1974593033835155e-03, -1.1198042255458091e-03, -2.0229858297200544e-18,
    9.9055884251036106e-04, 1.7192946812996687e-03, 2.1081914808740686e-03,
    2.1394563456147114e-03, 1.8513153796509275e-03, 1.3256746563059471e-03,
    6.7046275382674411e-04, 1.4017739341284856e-18, -5.8331646493291016e-04,
    -1.0032947419855074e-03, -1.2184233937642017e-03, -1.2238896950094572e-03,
    -1.0475948732225806e-03, -7.4152795801187910e-04, -3.7044471629615059e-04,
    -8.7363673902677791e-19, 3.1369007704874736e-04, 5.3155527927017054e-04,
    6.3531304436406394e-04, 6.2733453548486414e-04, 5.2717915649190441e-04,
    3.6582485933411455e-04, 1.7887185358663154e-04, 4.6362825333301403e-19,
    -1.4425544130102367e-04, -2.3771537702543715e-04, -2.7550349276906565e-04,
    -2.6290359742617416e-04, -2.1264395052204678e-04, -1.4132879082976897e-04,
    -6.5785565693739621e-05, -1.7899485045886187e-19])


_resample_func = {numpy.dtype('float32'): lal.ResampleREAL4TimeSeries,
                 numpy.dtype('float64'): lal.ResampleREAL8TimeSeries}

def resample_to_delta_t(timeseries, delta_t, method='butterworth'):
    """Resmple the time_series to delta_t

    Resamples the TimeSeries instance time_series to the given time step, 
    delta_t. Only powers of two and real valued time series are supported 
    at this time. Additional restrictions may apply to particular filter
    methods.

    Parameters
    ----------
    time_series: TimeSeries
        The time series to be resampled
    delta_t: float
        The desired time step 

    Returns
    -------
    Time Series: TimeSeries
        A TimeSeries that has been resampled to delta_t.

    Raises
    ------
    TypeError: 
        time_series is not an instance of TimeSeries.
    TypeError: 
        time_series is not real valued

    Examples
    --------

    >>> h_plus_sampled = resample_to_delta_t(h_plus, 1.0/2048)
    """

    if not isinstance(timeseries,TimeSeries):
        raise TypeError("Can only resample time series")

    if timeseries.kind is not 'real':
        raise TypeError("Time series must be real")

    if timeseries.delta_t == delta_t:
        return timeseries * 1

    if method == 'butterworth':
        lal_data = timeseries.lal()
        _resample_func[timeseries.dtype](lal_data, delta_t)
        data = lal_data.data.data 
        
    elif method == 'ldas':  
        factor = int(delta_t / timeseries.delta_t)
        
        if factor == 8:
            timeseries = resample_to_delta_t(timeseries, timeseries.delta_t * 4.0, method='ldas')
            factor = 2
        elif factor == 16:
            timeseries = resample_to_delta_t(timeseries, timeseries.delta_t * 4.0, method='ldas')
            factor = 4 
        elif factor == 32:
            timeseries = resample_to_delta_t(timeseries, timeseries.delta_t * 8.0, method='ldas')
            factor = 4 
        elif factor == 64:
            timeseries = resample_to_delta_t(timeseries, timeseries.delta_t * 16.0, method='ldas')
            factor = 4 

        try:
            filter_coefficients = LDAS_FIR_LP[factor]
        except:
            raise ValueError('Unsupported resample factor, %s, given' %factor)
            
        # apply the filter
        series = scipy.signal.lfilter(filter_coefficients, 1.0, 
                                      timeseries.numpy())
        
        # reverse the time shift caused by the filter
        corruption_length = len(filter_coefficients)
        data = numpy.zeros(len(timeseries))
        data[:len(data)-corruption_length/2] = series[corruption_length/2:]
        
        # zero out corrupted region
        data[0:corruption_length/2] = 0
        data[len(data)-corruption_length/2:] = 0       

        # Decimate the time series
        data = data[::factor] * 1
        
    else:
        raise ValueError('Invalid resampling method: %s' % method)
        
    return TimeSeries(data, delta_t = delta_t,
                      dtype=timeseries.dtype, 
                      epoch=timeseries._epoch)
       

_highpass_func = {numpy.dtype('float32'): lal.HighPassREAL4TimeSeries,
                 numpy.dtype('float64'): lal.HighPassREAL8TimeSeries}

def highpass(timeseries, frequency, filter_order=8, attenuation=0.1):
    """Return a new timeseries that is highpassed.

    Return a new time series that is highpassed above the `frequency`. 

    Parameters
    ----------
    Time Series: TimeSeries
        The time series to be high-passed.
    frequency: float
        The frequency below which is suppressed. 
    filter_order: {8, int}, optional
        The order of the filter to use when high-passing the time series.
    attenuation: {0.1, float}, optional
        The attenuation of the filter. 

    Returns
    -------
    Time Series: TimeSeries
        A  new TimeSeries that has been high-passed. 

    Raises
    ------
    TypeError: 
        time_series is not an instance of TimeSeries.
    TypeError: 
        time_series is not real valued

    """

    if not isinstance(timeseries, TimeSeries):
        raise TypeError("Can only resample time series")

    if timeseries.kind is not 'real':
        raise TypeError("Time series must be real")

    lal_data = timeseries.lal()
    _highpass_func[timeseries.dtype](lal_data, frequency, 
                                     1-attenuation, filter_order)

    return TimeSeries(lal_data.data.data, delta_t = lal_data.deltaT,
                      dtype=timeseries.dtype, epoch=timeseries._epoch)

    

__all__ = ['resample_to_delta_t', 'highpass']


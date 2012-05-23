from pycbc.types import zeros,TimeSeries,FrequencySeries,float32,complex64
from pycbc.fft import fft,ifft
from math import log,ceil,sqrt

def get_padded_frequencyseries(vec):
    if not isinstance(vec,TimeSeries):
        raise TypeError("Can only return padded frequency series from a time series")
    else:
        power = ceil(log(len(vec),2))+1
        N = 2 ** power
        n = N/2+1
        
        vec_pad = TimeSeries(zeros(N),delta_t=vec.delta_t,dtype=float32)
        vec_pad[0:len(vec)] = vec
        
        vectilde = FrequencySeries(zeros(n),delta_f=1, dtype=complex64)
        
        fft(vec_pad,vectilde)
        
        return vectilde

def get_frequencyseries(vec):
    if isinstance(vec,FrequencySeries):
        return vec
    if isinstance(vec,TimeSeries):
        N = len(vec)
        n = N/2+1    
        vectilde = FrequencySeries(zeros(n),delta_f=1, dtype=complex64)
        fft(vec,vectilde)
        
        return vectilde
    else:
        raise TypeError("Can only convert a TimeSeries to a FrequencySeries")
        

def sigmasq(htilde,psd = None,low_frequency_cutoff=None,high_frequency_cutoff=None):
    N = (len(htilde)-1) * 2 
    norm = 4.0 / (N * N * htilde.delta_f)
    moment = htilde.conj()*htilde
    kmin,kmax = get_cutoff_indices(low_frequency_cutoff,high_frequency_cutoff,htilde.delta_f,N)
    if psd is not None:
        moment /= psd
    return moment[kmin:kmax].sum() * norm
    
def get_cutoff_indices(flow,fhigh,df,N):
    if flow:
        kmin = int(flow / df)
    else:
        kmin = 1
    if fhigh:
        kmax = int(fhigh / df )
    else:
        kmax = N/2 + 1
        
    return kmin,kmax
    
def matchedfilter(template,data,psd=None,low_frequency_cutoff=None,high_frequency_cutoff=None):

    # Get the Inputs in terms of htilde and stilde
    htilde = get_frequencyseries(template)
    stilde = get_frequencyseries(data)

    # Calculate the length we need for the temporary memory 
    # Check that this is a power of two?
    N = (len(htilde)-1) * 2 
    
    kmin,kmax = get_cutoff_indices(low_frequency_cutoff,high_frequency_cutoff,stilde.delta_f,N)
   
    # Create workspace memory
    q = zeros(N,dtype=complex64)
    qtilde = zeros(N,dtype=complex64)
   
    #Weighted Correlation
    qtilde[kmin:kmax] = htilde.conj()[kmin:kmax] * stilde[kmin:kmax]
    
    if psd is not None:
        qtilde[kmin:kmax] /= psd[kmin:kmax]

    #Inverse FFT
    ifft(qtilde,q) 

    #Calculate the Normalization
    norm = sqrt(((4.0 / (N * N * stilde.delta_f)) **2) / sigmasq(htilde,psd,low_frequency_cutoff,high_frequency_cutoff) )

    #return complex snr
    return q,norm
    
    
def match(vec1,vec2,psd=None,low_frequency_cutoff=None,high_frequency_cutoff=None):
    htilde = get_frequencyseries(vec1)
    stilde = get_frequencyseries(vec2)
    snr,norm = matchedfilter(htilde,stilde,psd,low_frequency_cutoff,high_frequency_cutoff)
    maxsnrsq = (snr.conj()*snr).max()
    return sqrt(maxsnrsq/sigmasq(stilde,psd,low_frequency_cutoff,high_frequency_cutoff))*norm
    

from numpy import *


# The following functions all convert between physical parameters and f-stat values
# In particular, they do not need anything about a detector or network.


def params_to_a(d, cosi, psi, phi=0, d0=1.):
    """
    Calculate the F-stat A params given the physical parameters and a choice of
    d0 to set the overall scaling
    :param d: distance to source
    :param cosi: cos(inclination) of source
    :param psi: polarization of source
    :param phi: coalescence phase of source
    :param d0: overall scaling of A's
    """
    a_plus = d0 / d * (1. + cosi ** 2) / 2
    a_cross = d0 / d * cosi
    a = zeros(5)
    a[0] = d0
    a[1] = a_plus * cos(2 * phi) * cos(2 * psi) - a_cross * sin(2 * phi) * sin(2 * psi)
    a[2] = a_plus * cos(2 * phi) * sin(2 * psi) + a_cross * sin(2 * phi) * cos(2 * psi)
    a[3] = - a_plus * sin(2 * phi) * cos(2 * psi) - a_cross * cos(2 * phi) * sin(2 * psi)
    a[4] = - a_plus * sin(2 * phi) * sin(2 * psi) + a_cross * cos(2 * phi) * cos(2 * psi)
    return a


def a_to_params(a):
    """
    Calculate the physical parameters based upon the F-stat A parameters.
    :param a: array of f-stat params, entry zero assumed to be d0
    """
    # these variables are what they say [ (a_plus +/- a_cross)^2 ]
    ap_plus_ac_2 = (a[1] + a[4]) ** 2 + (a[2] - a[3]) ** 2
    ap_minus_ac_2 = (a[1] - a[4]) ** 2 + (a[2] + a[3]) ** 2
    a_plus = 0.5 * (sqrt(ap_plus_ac_2) + sqrt(ap_minus_ac_2))
    a_cross = 0.5 * (sqrt(ap_plus_ac_2) - sqrt(ap_minus_ac_2))
    amp = a_plus + sqrt(a_plus ** 2 - a_cross ** 2)
    cosi = a_cross / amp
    d = a[0] / amp
    psi = 0.5 * arctan2(a_plus * a[2] + a_cross * a[3], a_plus * a[1] - a_cross * a[4])
    phi = 0.5 * arctan2(-a_plus * a[3] - a_cross * a[2], a_plus * a[1] - a_cross * a[4])
    return d, cosi, psi, phi


def a_to_circular(a):
    """
    Calculate the circular F-stat A parameters given in Whelan et al 2013
    :param a: array of f-stat params, entry zero assumed to be d0
    """
    a_circ = zeros(5)
    a_circ[0] = a[0]
    a_circ[1] = 0.5 * (a[1] + a[4])
    a_circ[2] = 0.5 * (a[2] - a[3])
    a_circ[3] = 0.5 * (a[1] - a[4])
    a_circ[4] = - 0.5 * (a[2] + a[3])
    return a_circ


def a_to_circ_amp(a):
    """
    Calculate the amplitudes of left/right circularly polarized waveforms
    from the F-stat A parameters
    :param a: array of f-stat params, entry zero assumed to be d0
    """
    a_circ = a_to_circular(a)
    ar_hat = sqrt(a_circ[1] ** 2 + a_circ[2] ** 2)
    al_hat = sqrt(a_circ[3] ** 2 + a_circ[4] ** 2)
    return ar_hat, al_hat


def phase_diff(a):
    """
    Calculate the phase difference (not sure what of)
    :param a: array of f-stat params, entry zero assumed to be d0
    """
    return arctan2(a[1] * a[4] - a[2] * a[3], a[1] * a[2] + a[3] * a[4])


def amp_ratio(a):
    """
    Calculate the amplitude ratio
    :param a: array of f-stat params, entry zero assumed to be d0
    """
    return sqrt((a[1] ** 2 + a[3] ** 2) / (a[2] ** 2 + a[4] ** 2))


# The following functions calculate SNRs, likelihoods, etc for a signal, given a network.
# They all work in the dominant polarization (i.e. assuming that the network is described
# by F+, Fx and they're orthogonal)

def expected_snr(a, f_plus, f_cross):
    """
    Calculate the SNR for a given set of A parameters and network sensitivity.
    :param a: the F-stat A parameters
    :param f_plus: F_plus sensitivity
    :param f_cross: F_cross sensitivity
    """
    f = array([0, f_plus, f_cross, f_plus, f_cross])
    snrsq = sum(f ** 2 * a ** 2)
    return sqrt(snrsq)


def set_snr(a, f_plus, f_cross, snr):
    """
    rescale distance to give desired SNR, return rescaled as and distance
    :param a: the F-stat A parameters
    :param f_plus: F_plus sensitivity
    :param f_cross: F_cross sensitivity
    :param snr: the desired SNR
    """
    s = expected_snr(a, f_plus, f_cross)
    a_scale = a * snr / s
    a_scale[0] = a[0]
    d_scale = a_to_params(a_scale)[0]
    return a_scale, d_scale


def lost_snrsq(a_hat, a, f_plus, f_cross):
    """
    Calculate the difference in SNRSQ between the true parameters a_hat
    and the template a, network sensitivity f_plus, f_cross
    :param a_hat: the observed F-stat A parameters
    :param a: the "template" F-stat A parameters
    :param f_plus: sensitivity to plus polarization
    :param f_cross: sensitivity to cross polarization
    """
    f = array([0, f_plus, f_cross, f_plus, f_cross])
    snrsq = sum(f ** 2 * (a_hat - a) ** 2)
    return snrsq


def circ_project(a, f_plus, f_cross, hand):
    """
    Project the f-stat A parameters to those that would be observed by restricting
    to left or right circular polarization
    :param f_plus: sensitivity to plus polarization
    :param f_cross: sensitivity to cross polarization
    :param hand: one of "left", "right"
    """
    if hand == "right":
        x = 1
    elif hand == "left":
        x = -1
    else:
        raise ValueError("hand must be either left or right")

    f = array([f_plus, f_cross, f_plus, f_cross])

    fa = (a[1:] * f)
    # matrix that projects FA onto circular configuration
    p = array([[f_plus ** 2, 0, 0, x * f_plus * f_cross],
               [0, f_cross ** 2, -x * f_plus * f_cross, 0],
               [0, -x * f_plus * f_cross, f_plus ** 2, 0],
               [x * f_plus * f_cross, 0, 0, f_cross ** 2]])
    p /= (f_plus ** 2 + f_cross ** 2)
    fa_proj = inner(p, fa)
    a_proj = zeros_like(a)
    a_proj[0] = a[0]
    a_proj[1:] = fa_proj / f
    a_proj[isnan(a_proj)] = 0
    return a_proj

# These functions allow us to go from SNRs to F-stat params and back

def snr_f_to_a(z, f_sig):
    """
    Given the complex SNR and the detector sensitivities, calculate the f-stat A parameters
    :param z: array containing complex snrs for the detectors
    :param f_sig: sensitivity of detectors sigma * (F+, Fx)
    """
    m = zeros((2, 2))
    for f in f_sig:
        m += outer(f, f)
    s_h = inner(z, f_sig.transpose())
    a_max = inner(s_h, linalg.inv(m))
    a = array([1.0, a_max[0].real, a_max[1].real, a_max[0].imag, a_max[1].imag])
    return a

def a_f_to_snr(a, f_plus, f_cross):
    """
    Given the F-stat a parameters and f_plus, f_cross for a detector, calculate the SNR.
    From the Harry-Fairhurst paper, the waveform is $h = A^{\mu} h_{\mu}$ where 
    $h_1 = F_{+} h_{0}$, $h_2 = F_{\times} h_{0}$, $h_{3} = F_{+} h_{\pi/2}$, 
    $h_{4} = F_{\times} h_{\pi/2}$ and $z = (s | h_{0}) + i (s | h_{\pi/2})$.  
    Given the $A^{\mu}$, the expected SNR is:
    $z = F_{+} A^{1} + F_{\times} A^{2} + i( F_{+} A^{3} + F_{\times} A^{4})$
    :param a: F-stat parameters
    :param f_plus: Sensitivity to plus polarization
    :param f_cross: Sensitivity to cross polarization
    """
    z = f_plus * (a[1] + 1j * a[3]) + f_cross * (a[2] + 1j * a[4]) 
    return z

# Calculate the dominant polarization F+, Fx and the angle between this and the original frame.

def dominant_polarization(f_sig):
    """
    Given the detector responses, compute the dominant polarization F+, Fx and the
    angle that we have to rotate through to get to D.P.
    :param f_sig: sensitivity of detectors: sigma * (F+, Fx)
    """
    m = zeros((2, 2))
    for f in f_sig:
        m += outer(f, f)
    f_cross, f_plus = sort(sqrt(linalg.eig(m)[0])) #This returns eigenvalues
 #   f_cross, f_plus = linalg.eig(m)[1] #This returns eigenvectors
    chi = 1./4 * arctan2(2 * m[0,1], m[0,0] - m[1,1])
    return f_plus, f_cross, chi 


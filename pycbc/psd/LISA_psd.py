L = 0.03342293561

def LISA_psd_X(S_acc, facc, S_IMS=1e-40):
    """Returns the nosie PSD X, PSD XY
    Parameters
    ----------
    S_acc : PSD of the acceleration noise
    facc : frequencies
    S_IMS : PSD of the Interferometric Measurement System
    Returns
    -------
    numpy.ndarray
        PSD X, PSD XY, in order
    """
    omega = 2*np.pi*facc
    Sx = 16*(sin(omega*L)**2)*(2*(1 + cos(omega*L)**2)*S_acc + S_IMS)
    Sxy = -8*(sin(omega*L)**2)*(cos(omega*L)*(S_IMS + 4*S_acc))
    return np.array([Sx, Sxy])

def LISA_psd_AET(S_acc, facc, S_IMS=1e-40):
    """Returns the nosie PSD A, E, T
    Parameters
    ----------
    S_acc : PSD of the acceleration noise
    facc : frequencies
    S_IMS : PSD of the Interferometric Measurement System
    Returns
    -------
    numpy.ndarray
        PSD A, E, T in order
    """
    omega = 2*np.pi*facc
    Sa = 8*sin(omega*L)**2*4*((1 + cos(omega*L) + \
            cos(omega*L)**2)*S_acc + (2 + cos(omega*L))*S_IMS)
    Se = np.copy(Sa)
    St = 32*(sin(omega*L)**2)*(sin(omega*L/2)**2)* \
            (4*sin(omega*L/2)**2*S_acc + S_IMS)
    return np.array([Sa, Se, St])

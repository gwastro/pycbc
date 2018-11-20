import numpy as np
from numpy import log
import copy, h5py
from scipy.interpolate import interp1d
from scipy.integrate import quad
from astropy.cosmology import WMAP9 as cosmo

from pycbc.conversions import mchirp_from_mass1_mass2 as m1m2tomch

_mch_BNS = 1.4/2**.2
_redshifts, _d_lum, _I = np.arange(0., 5., 0.01), [], []
_save_params = ['mass1', 'mass2', 'spin1z', 'spin2z', 'spin1y', 'spin2y',
                                'spin1x', 'spin2x', 'distance', 'end_time']

for zz in _redshifts:
    _d_lum.append(cosmo.luminosity_distance(zz).value)
_dlum_interp = interp1d(_d_lum, _redshifts)

def read_injections(sim_files, m_dist, s_dist, d_dist):
    ''' Read all the injections from the files in the provided folder.
        The files must belong to individual set i.e. no files that combine
        all the injections in a run.
        Identify injection strategies and finds parameter boundaries.
        Collect injection according to GPS.

        Parameters
        ----------
        sim_files: list
           List containign names of the simulation files
        m_dist: list
           The mass distribution used in the simulation runs
        s_dist: list
           The spin distribution used in the simulation runs
        d_dist: list
           The distance distribution used in the simulation runs

        Returns
        -------
        injections: dictionary
           Contains the organized information about the injections
    '''

    injections = {}
    min_d, max_d = 1e12, 0
    nf = len(sim_files)
    for i in range(nf):

        key = str(i)
        injections[key] = process_injections(sim_files[i])
        injections[key]['file_name'] = sim_files[i]
        injections[key]['m_dist'] = m_dist[i]
        injections[key]['s_dist'] = s_dist[i]
        injections[key]['d_dist'] = d_dist[i]

        mass1, mass2 = injections[key]['mass1'], injections[key]['mass2']
        distance = injections[key]['distance']

        mchirp = m1m2tomch(mass1, mass2)
        injections[key]['chirp_mass'] = mchirp
        injections[key]['total_mass'] = mass1 + mass2

        injections[key]['mtot_range'] = [min(mass1 + mass2), max(mass1 + mass2)]
        injections[key]['m1_range'] = [min(mass1), max(mass1)]
        injections[key]['m2_range'] = [min(mass2), max(mass2)]
        injections[key]['d_range'] = [min(distance), max(distance)]

        min_d, max_d = min(min_d, min(distance)), max(max_d, max(distance))

    injections['z_range'] = [dlum_to_z(min_d), dlum_to_z(max_d)]

    return injections

def estimate_vt(injections, mchirp_sampler, model_pdf, **kwargs):
    #Try including ifar threshold
    '''Based on injection strategy and the desired astro model estimate the injected volume.
       Scale injections and estimate sensitive volume.

       Parameters
       ----------
       injections: dictionary
            Dictionary obtained after reading injections from read_injections
       mchirp_sampler: function
            Sampler for producing chirp mass samples for the astro model.
       model_pdf: function
            The PDF for astro model in mass1-mass2-spin1z-spin2z space.
            This is easily extendible to include precession
       kwargs: key words
            Inputs for thresholds and astrophysical models

       Returns
       -------
       injection_chunks: dictionary
        The input dictionary with VT and VT error included with the injections
    '''

    thr_var = kwargs.get('thr_var')
    thr_val = kwargs.get('thr_val')

    nsamples = 1000000 #Used to calculate injected astro volume
    injections = copy.deepcopy(injections)
    min_z, max_z = injections['z_range']
    V = quad(contracted_dVdc, 0., max_z)[0]

    z_astro = astro_redshifts(min_z, max_z, nsamples)
    astro_lum_dist = cosmo.luminosity_distance(z_astro).value

    mch_astro = np.array(mchirp_sampler(nsamples = nsamples, **kwargs))
    mch_astro_det = mch_astro * (1. + z_astro)
    idx_within = np.zeros(nsamples)

    for key in injections.keys():

        if key == 'z_range':
            # This is repeated down again and is so
            continue

        mchirp = injections[key]['chirp_mass']
        min_mchirp, max_mchirp = min(mchirp),  max(mchirp)
        distance = injections[key]['distance']

        if injections[key]['d_dist'] == 'uniform':
            d_min, d_max = min(distance), max(distance)
        elif injections[key]['d_dist'] == 'dchirp':
            d_fid_min = min(distance / (mchirp/_mch_BNS)**(5/6.))
            d_fid_max = max(distance / (mchirp/_mch_BNS)**(5/6.))

            d_min = d_fid_min * (mch_astro_det/_mch_BNS)**(5/6.)
            d_max = d_fid_max * (mch_astro_det/_mch_BNS)**(5/6.)

        bound = np.sign((max_mchirp-mch_astro_det)*(mch_astro_det-min_mchirp))
        bound += np.sign((d_max - astro_lum_dist)*(astro_lum_dist - d_min))

        idx = np.where(bound == 2)
        idx_within[idx] = 1

    inj_V0 = 4*np.pi*V*len(idx_within[idx_within == 1])/float(nsamples)
    injections['inj_astro_vol'] = inj_V0

    # Estimate the sensitive volume
    z_range = injections['z_range']
    V_min = quad(contracted_dVdc, 0., z_range[0])[0]
    V_max = quad(contracted_dVdc, 0., z_range[1])[0]

    thr_falloff, i_inj, i_det, i_det_sq = [], 0, 0, 0
    gps_min, gps_max = 1e15, 0
    keys = injections.keys()
    for key in keys:

        if key == 'z_range' or key == 'inj_astro_vol':
            continue

        data = injections[key]
        distance = data['distance']
        mass1, mass2 = data['mass1'], data['mass2']
        spin1z, spin2z = data['spin1z'], data['spin2z']
        mchirp = data['chirp_mass']
        gps_min = min(gps_min, min(data['end_time']))
        gps_max = max(gps_max, max(data['end_time']))

        z_inj = dlum_to_z(distance)
        m1_sc, m2_sc = mass1/(1 + z_inj), mass2/(1 + z_inj)
        p_out = model_pdf(m1_sc, m2_sc, spin1z, spin2z)
        p_out *= pdf_z_astro(z_inj, V_min, V_max)

        p_in = 0
        J = cosmo.luminosity_distance(z_inj + 0.0005).value
        J -= cosmo.luminosity_distance(z_inj - 0.0005).value
        J = abs(J)/0.001 # A quick way to get dD_l/dz

        # Sum probability of injections from j-th set for all the strategies
        for key2 in keys:

            if key2 == 'z_range' or key2 == 'inj_astro_vol':
                continue

            dt_j = injections[key2]
            dist_j = dt_j['distance']
            m1_j, m2_j = dt_j['mass1'], dt_j['mass2']
            s1x_2, s2x_2 = dt_j['spin1x'], dt_j['spin2x']
            s1y_2, s2y_2 = dt_j['spin1y'], dt_j['spin2y']
            s1z_2, s2z_2 = dt_j['spin1z'], dt_j['spin2z']
            s1 = np.sqrt(s1x_2**2 + s1y_2**2 + s1z_2**2)
            s2 = np.sqrt(s2x_2**2 + s2y_2**2 + s2z_2**2)
            mch_j = dt_j['chirp_mass']

        #Get probability density for injections in mass-distance space
            if dt_j['m_dist'] == 'totalMass':
                lomass, himass = min(min(m1_j), min(m2_j), max(max(m1_j), max(m2_j)))
                lomass_2, himass_2 = lomass, himass
            elif dt_j['m_dist'] == 'componentMass' or dt_j['m_dist'] == 'log':
                lomass, himass = min(m1_j), max(m1_j)
                lomass_2, himass_2 = min(m2_j), max(m2_j)

            if dt_j['d_dist'] == 'dchirp':
                l_dist = min(dist_j / (mch_j/_mch_BNS)**(5/6.))
                h_dist = max(dist_j / (mch_j/_mch_BNS)**(5/6.))
            elif dt_j['d_dist'] == 'uniform':
                l_dist, h_dist = min(dist_j), max(dist_j)

            mdist = dt_j['m_dist']
            prob_mass = inj_mass_pdf(mdist, mass1, mass2,
                                          lomass, himass, lomass_2, himass_2)

            ddist = dt_j['d_dist']
            prob_dist = inj_distance_pdf(ddist, distance, l_dist,
                                                              h_dist, mchirp)

            hspin1, hspin2 = max(s1), max(s2)
            prob_spin = inj_spin_pdf(dt_j['s_dist'], hspin1, spin1z)
            prob_spin *= inj_spin_pdf(dt_j['s_dist'], hspin2, spin2z)

            p_in += prob_mass * prob_dist * prob_spin * J * (1 + z_inj)**2

        p_in[p_in == 0] = 1e12
        p_out_in = p_out/p_in

        i_inj += np.sum(p_out_in)
        i_det += np.sum((p_out_in)[data[thr_var] > thr_val])
        i_det_sq += np.sum((p_out_in)[data[thr_var] > thr_val]**2)

        idx_thr = np.where(data[thr_var] > thr_val)
        thrs = data[thr_var][idx_thr]
        ratios = p_out_in[idx_thr]/max(p_out_in[idx_thr])
        rndn = np.random.uniform(0, 1, len(ratios))
        idx_ratio = np.where(ratios > rndn)
        thr_falloff.append(thrs[idx_ratio])

    inj_V0 = injections['inj_astro_vol']
    injections['ninj'] = i_inj
    injections['ndet'] = i_det
    injections['ndetsq'] = i_det_sq
    injections['VT'] = ((inj_V0*i_det/i_inj) * (gps_max - gps_min)/31557600)
    injections['VT_err'] = injections['VT'] * np.sqrt(i_det_sq)/i_det
    injections['thr_falloff'] = np.hstack(np.array(thr_falloff).flat)

    return injections

def process_injections(hdffile):
    """Function to read in the injection file and
       extract the found injections and all injections

       Parameters
       ----------
       hdffile: hdf file
           File for which injections are to be processed

       Returns
       -------
       data: dictionary
           Dictionary containing injection read from the input file
    """
    data = {}

    with h5py.File(hdffile, 'r') as inp:
        found_index = inp['found_after_vetoes/injection_index'][:]

        for param in _save_params:
            data[param] = inp['injections/'+param][:]

        ifar = np.zeros_like(data[_save_params[0]])
        ifar[found_index] = inp['found_after_vetoes/ifar'][:]

        data['ifar'] = ifar

        stat = np.zeros_like(data[_save_params[0]])
        stat[found_index] = inp['found_after_vetoes/stat'][:]

        data['stat'] = stat

    return data

def dlum_to_z(dl):
    ''' Get the redshift for a luminosity distance

        Parameters
        ----------
        dl: array
           The array of luminosity distances

        Returns
        -------
        array
           The redshift values corresponding to the luminosity distances
    '''

    return _dlum_interp(dl)

def astro_redshifts(min_z, max_z, nsamples):
    '''Sample the redshifts for sources, with redshift
            independent rate, using standard cosmology

       Parameters
       ----------
       min_z: float
            Minimum redshift
       max_z: float
            Maximum redshift
       nsamples: int
            Number of samples

       Returns
       -------
       z_astro: array
            nsamples of redshift, between min_z, max_z, by standard cosmology
    '''

    dz, fac = 0.001, 3.0
    # use interpolation instead of directly estimating all the pdfz for rndz
    V = quad(contracted_dVdc, 0., max_z)[0]
    zbins = np.arange(min_z, max_z + dz/2., dz)
    zcenter = (zbins[:-1] + zbins[1:]) / 2
    pdfz = cosmo.differential_comoving_volume(zcenter).value/(1+zcenter)/V

    int_pdf = interp1d(zcenter, pdfz, bounds_error=False, fill_value=0)

    rndz = np.random.uniform(min_z, max_z, int(fac*nsamples))
    pdf_zs = int_pdf(rndz)
    maxpdf = max(pdf_zs)
    rndn = np.random.uniform(0, 1, int(fac*nsamples)) * maxpdf
    diff = pdf_zs - rndn
    idx = np.where(diff > 0)
    z_astro = rndz[idx]

    np.random.shuffle(z_astro)
    z_astro.resize(nsamples)

    return z_astro

def pdf_z_astro(z, V_min, V_max):
    ''' Get the probability density for the rate of events
        at a redshift assuming standard cosmology
    '''
    return contracted_dVdc(z)/(V_max - V_min)

def contracted_dVdc(z):
    #Return the time-dilated differential comoving volume
    return cosmo.differential_comoving_volume(z).value/(1+z)

##### Defining current standard strategies used for making injections #####

def inj_mass_pdf(key, mass1, mass2, lomass, himass, lomass_2 = 0, himass_2 = 0):

    '''Estimate the probability density based on the injection strategy

       Parameters
       ----------
       key: string
          Injection strategy
       mass1: array
          First mass of the injections
       mass2: array
          Second mass of the injections
       lomass: float
          Lower value of the mass distributions
       himass: float
          higher value of the mass distribution

       Returns
       -------
       pdf: array
          Probability density of the injections
    '''

    mass1, mass2 = np.array(mass1), np.array(mass2)

    if key == 'totalMass':
        # Returns the PDF of mass when total mass is uniformly distributed.
        # Both the component masses have the same distribution for this case.

        # Parameters
        # ----------
        # lomass: lower component mass
        # himass: higher component mass

        bound = np.sign((lomass + himass) - (mass1 + mass2))
        bound += np.sign((himass - mass1)*(mass1 - lomass))
        bound += np.sign((himass - mass2)*(mass2 - lomass))
        idx = np.where(bound != 3)
        pdf = 1./(himass - lomass)/(mass1 + mass2 - 2 * lomass)
        pdf[idx] = 0

        return pdf

    if key == 'componentMass':
        # Returns the PDF of mass when component mass is uniformly
        # distributed. Component masses are independent for this case.

        # Parameters
        # ----------
        # lomass: lower component mass
        # himass: higher component mass

        bound = np.sign((himass - mass1)*(mass1 - lomass))
        bound += np.sign((himass_2 - mass2)*(mass2 - lomass_2))
        idx = np.where(bound != 2)
        pdf = np.ones_like(mass1) / (himass - lomass) / (himass_2 - lomass_2)
        pdf[idx] = 0

        return pdf

    if key == 'log':
        # Returns the PDF of mass when component mass is uniform in log.
        # Component masses are independent for this case.

        # Parameters
        # ----------
        # lomass: lower component mass
        # himass: higher component mass

        bound = np.sign((himass - mass1)*(mass1 - lomass))
        bound += np.sign((himass_2 - mass2)*(mass2 - lomass_2))
        idx = np.where(bound != 2)
        pdf = 1 / (log(himass) - log(lomass)) / (log(himass_2) - log(lomass_2))
        pdf /= (mass1 * mass2)
        pdf[idx] = 0

        return pdf

def inj_spin_pdf(key, high_spin, spinz):
    ''' Estimate the probability density of the
           injections for the spin distribution.

        Parameters
        ----------
        key: string
          Injections strategy
        high_spin: float
          Maximum spin used in the strategy
        spinz: array
          Spin of the injections (for one component)
    '''

    # If the data comes from disable_spin simulation
    if spinz[0] == 0:
        return np.ones_like(spinz)

    spinz = np.array(spinz)

    bound = np.sign(np.absolute(high_spin) - np.absolute(spinz))
    bound += np.sign(1 - np.absolute(spinz))

    if key == 'precessing':
        # Returns the PDF of spins when total spin is
        # isotropically distributed. Both the component
        # masses have the same distribution for this case.

        pdf = (np.log(high_spin - np.log(abs(spinz)))/high_spin/2)
        idx = np.where(bound != 2)
        pdf[idx] = 0

        return pdf

    if key == 'aligned':
        # Returns the PDF of mass when spins are aligned and uniformly
        # distributed. Component spins are independent for this case.

        pdf = (np.ones_like(spinz) / 2 / high_spin)
        idx = np.where(bound != 2)
        pdf[idx] = 0

        return pdf

    if key == 'disable_spin':
        # Returns unit array

        pdf = np.ones_like(spinz)

        return pdf

def inj_distance_pdf(key, distance, low_dist, high_dist, mchirp = 1):
    ''' Estimate the probability density of the
        injections for the distance distribution.

        Parameters
        ----------
        key: string
          Injections strategy
        distance: array
          Array of distances
        low_dist: float
          Lower value of distance used in the injection strategy
        high_dist: float
          Higher value of distance used in the injection strategy
    '''

    distance = np.array(distance)

    if key == 'uniform':
        # Returns the PDF at a distance when
        # distance is uniformly distributed.

        pdf = np.ones_like(distance)/(high_dist - low_dist)
        bound = np.sign((high_dist - distance)*(distance - low_dist))
        idx = np.where(bound != 1)
        pdf[idx] = 0
        return pdf

    if key == 'dchirp':
        # Returns the PDF at a distance when distance is uniformly
        # distributed but scaled by the chirp mass

        weight = (mchirp/_mch_BNS)**(5./6)
        pdf = np.ones_like(distance) / weight / (high_dist - low_dist)
        bound = np.sign((weight*high_dist - distance)*(distance - weight*low_dist))
        idx = np.where(bound != 1)
        pdf[idx] = 0
        return pdf

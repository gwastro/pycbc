import numpy as np
import glob, sys, copy, h5py
from scipy import interpolate
from scipy.integrate import quad
from scipy.stats import kstest, ks_2samp
from astropy.cosmology import WMAP9 as cosmo

from pycbc.conversions import mchirp_from_mass1_mass2

_mch_BNS = 1.4/2**.2
_redshifts, _d_lum, _I = np.arange(0., 5., 0.01), [], []
_save_params = ['mass1', 'mass2', 'spin1z', 'spin2z', 'spin1y', 'spin2y', 'spin1x', 'spin2x', 'distance', 'end_time']

for zz in _redshifts:
    _d_lum.append(cosmo.luminosity_distance(zz).value)
_dlum_interp = interpolate.interp1d(_d_lum, _redshifts)

def read_injections(folder_name):
    ''' Read all the injections from the files in the provided folder.
        The files must belong to individual set i.e. no files that combine 
        all the injections in a run.
        Identify injection strategies and finds parameter boundaries.
        Collect injection according to GPS.

        Parameters
        ----------
        folder_name: string
           Location of folder containing simulation files

        Returns
        -------
        chunk_data: dictionary
           Contains the organized information about the injections
    '''
    
    injections = {}
    all_files = glob.glob(folder_name)
    
    nf, min_d, max_d = len(all_files), 1e12, 0
    for i in range(nf):
        injections[str(i)] = process_injections(all_files[i])
        injections[str(i)]['file_name'] = all_files[i]
        
        mass1, mass2 = injections[str(i)]['mass1'], injections[str(i)]['mass2']   
        mchirp = mchirp_from_mass1_mass2(mass1, mass2)
        injections[str(i)]['chirp_mass'] = mchirp
        injections[str(i)]['total_mass'] = mass1 + mass2
    
        min_m1, max_m1 = min(mass1), max(mass1)
        min_m2, max_m2 = min(mass2), max(mass2)
        min_totM, max_totM = min(mass1+mass2), max(mass1+mass2)        
     
        # Identify the mass distribution
    
        def cdf_componentMass(xx): #Assumes both the masses have the same distribution
            return (xx - min_m2)/(max_m2 - min_m2)

        def cdf_totalMass(xx):
            return (xx - min_totM)/(max_totM - min_totM)
 
        def cdf_log(xx): #Assumes both the masses have the same distribution
            return (np.log(xx) - np.log(min_m2))/(np.log(max_m2) - np.log(min_m2))        
    
        p_thr = 0.0001
        stat, pvalue = kstest(mass1, cdf_componentMass)
        if p_thr < pvalue:
            injections[str(i)]['m_dist'] = 'componentMass'        
        
        stat, pvalue = kstest(mass1 + mass2, cdf_totalMass)    
        if p_thr < pvalue:
            injections[str(i)]['m_dist'] = 'totalMass'  
            
        stat, pvalue = kstest(mass1, cdf_log)
        if p_thr < pvalue:
            injections[str(i)]['m_dist'] = 'log'             
        
        spin1z, spin2z = injections[str(i)]['spin1z'], injections[str(i)]['spin2z']   
        spin1y, spin2y = injections[str(i)]['spin1y'], injections[str(i)]['spin2y']
        spin1x, spin2x = injections[str(i)]['spin1x'], injections[str(i)]['spin2x']   
        spin1, spin2 = np.sqrt(spin1z**2 + spin1y**2 + spin1x**2), np.sqrt(spin2z**2 + spin2y**2 + spin2x**2) 
        injections[str(i)]['spin1'] = spin1
        injections[str(i)]['spin2'] = spin2
    
        min_s1z, max_s1z = min(spin1z), max(spin1z)
        min_s2z, max_s2z = min(spin2z), max(spin2z)
        min_totS, max_totS = min(spin1), max(spin2)        
        
        # Identify the spin distribution  
        if spin1[0] == 0 and spin1z[0] == 0:
            injections[str(i)]['s_dist'] = 'disable_spin'    
    
        if spin1y[0] == 0 and spin2y[0] == 0 and spin1x[0] == 0 and spin2x[0] == 0 and max_totS > 0.: 
            injections[str(i)]['s_dist'] = 'aligned'
        
        if spin1y[0] != 0 or spin2y[0] != 0 or spin1x[0] != 0 and spin2x[0] != 0:
            injections[str(i)]['s_dist'] = 'precessing'   
        
        # Identify the distance distribution
        distance = injections[str(i)]['distance']
        fiducial_distance = distance / (mchirp/_mch_BNS)**(5/6.)
        injections[str(i)]['fiducial_distance'] = fiducial_distance
        min_dist, max_dist = min(distance), max(distance)
        min_fid_dist, max_fid_dist = min(fiducial_distance), max(fiducial_distance)
    
        def cdf_distance(xx):
            return (xx - min_dist)/(max_dist - min_dist)

        def cdf_fid_distance(xx):
            return (xx - min_fid_dist)/(max_fid_dist - min_fid_dist)    
    
        stat, pvalue = kstest(distance, cdf_distance)
        if p_thr < pvalue:
            injections[str(i)]['d_dist'] = 'uniform'    
        
        stat, pvalue = kstest(fiducial_distance, cdf_fid_distance)
        if p_thr < pvalue:
            injections[str(i)]['d_dist'] = 'dchirp'    
        
        if not(set(['m_dist', 's_dist', 'd_dist']) < set(injections[str(i)].keys())):
            sys.exit("Couldn't indentify the Mass, Spin or Distance distribution at the given threshold!")  
            
            
        d1, d2 = min(distance), max(distance)
        min_d, max_d = min(min_d, d1), max(max_d, d2)            
        
    # Get the injections ranges
    variables = ['total_mass', 'mass1', 'mass2',
                 'spin1z', 'spin2z', 'spin1', 'spin2', 
                 'distance', 'fiducial_distance', 'end_time']
    keys = ['mtot_range', 'm1_range', 'm2_range', 
            's1z_range', 's2z_range', 's1_range', 's2_range', 
            'd_range', 'fid_range', 'gps']
    threshold = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.8]
    for var, key, thr in zip(variables, keys, threshold):
        for i in range(nf):        
            value = injections[str(i)][var]
            bins, min_value, max_value = [], 1.e12, 0.
            for j in range(i, nf):
        
                if key not in injections[str(j)].keys():
                    value_2 = injections[str(j)][var]
                    stat, pvalue = ks_2samp(value, value_2)
        
                    if pvalue > thr:  
                        bins.append(j)
                        min_value, max_value = min(min_value, min(value_2)), max(max_value, max(value_2))    
            for j in bins:
                if key not in ['s1z_range', 's2z_range']:
                    #low, high = 10**(2 - round(np.log10(min_value))), 10**(2 - round(np.log10(max_value)))
                    injections[str(j)][key] = [round(min_value - 0.5), round(max_value + 0.5)]
                else:
                    injections[str(j)][key] = [min_value, max_value]

    for i in range(nf): 
        if injections[str(i)]['d_dist'] == 'dchirp':
            injections[str(i)]['d_range'] = injections[str(i)]['fid_range']
        del injections[str(i)]['fiducial_distance']
        del injections[str(i)]['fid_range']
        del injections[str(i)]['total_mass']
    
    nchunks, chunk_data, chunk_gps = -1, {}, []
    for key_files_1 in injections.keys():
    
        if injections[key_files_1]['gps'] in chunk_gps:
            continue
        chunk_gps.append(injections[key_files_1]['gps'])
        nchunks += 1
    
        key = str(nchunks)
        chunk_data[key] = {}
        chunk_data[key]['gps'] = injections[key_files_1]['gps']   
        for data_key in injections[key_files_1].keys():
            if data_key == 'gps':
                continue        
            chunk_data[key][data_key] = []    
        for key_files_2 in injections.keys():
            if injections[key_files_1]['gps'] == injections[key_files_2]['gps']:
            
                for data_key in injections[key_files_1].keys():
                    if data_key == 'gps':
                        continue
                    chunk_data[key][data_key].append(injections[key_files_2][data_key])

    chunk_data['z_range'] = [round(10*dlum_to_z(min_d) - 0.5)/10., round(10*dlum_to_z(max_d) + 0.5)/10.]
                    
    return chunk_data

def estimate_vt(inj_chunks, mchirp_sampler, model_pdf, **kwargs): #Need to make this flexible to include ifar threshold
    '''Based on injection strategy and the desired astro model estimate the injected volume. 
       Scale injections and estimate sensitive volume.

       Parameters
       ----------       
       inj_chunks: dictionary 
            Dictionary obtained after reading injections using function read_injections 
       mchirp_sampler: function
            Sampler for producing chirp mass samples for the astro model.
       model_pdf: function
            The PDF for astro model in mass1-mass2-spin1z-spin2z space. This is easily extendible to include precession
       kwargs: key words
            Inputs for thresholds and astrophysical models

       Returns
       -------
       injection_chunks: dictionary
            The input dictionary with VT and VT error included for each chunk
    '''
    
    min_mass = kwargs.get('min_mass', 5.)
    max_mass = kwargs.get('max_mass', 95.)
    max_mtotal = kwargs.get('max_mtotal', min_mass + max_mass)
    thr_var = kwargs.get('thr_var')
    thr_val = kwargs.get('thr_val')

    injection_chunks = copy.deepcopy(inj_chunks)
    nsamples, min_z, max_z = 1000000, injection_chunks['z_range'][0], injection_chunks['z_range'][1]
    for key in injection_chunks.keys():
        
        if key == 'z_range': # This is repeated down again and is sort of hard-coded
            continue
            
        data = injection_chunks[key]
        inj_distr = zip(data['mtot_range'], data['m1_range'], data['m2_range'], #Hard coded
                             data['d_range'], zip(data['m_dist'], data['d_dist']))
        
        unique_distr = np.unique(inj_distr, axis = 0)
        injection_chunks[key]['inj_distr'] = inj_distr
        injection_chunks[key]['unique_distr'] = unique_distr
        
        if 'inj_astro_vol' in injection_chunks[key].keys():
            continue        
        
        z_astro = astro_redshifts(min_z, max_z, nsamples)
        astro_lum_dist = cosmo.luminosity_distance(z_astro).value
        V = quad(lambda z: cosmo.differential_comoving_volume(z).value/(1+z), 0., max_z)[0]
        
        bound, mch_astro_det = [], np.array(mchirp_sampler(nsamples = nsamples, **kwargs)) * (1. + z_astro)
        
        for uu in unique_distr:
            if uu[4][0] == 'totalMass':
                min_mchirp = mchirp_from_mass1_mass2(float(uu[1][0]), float(uu[1][0]))
                max_mchirp = mchirp_from_mass1_mass2(float(uu[0][1])/2., float(uu[0][1])/2.)
            if uu[4][0] == 'componentMass' or uu[4][0] == 'log':
                min_mchirp = mchirp_from_mass1_mass2(float(uu[1][0]), float(uu[2][0]))
                max_mchirp = mchirp_from_mass1_mass2(float(uu[1][1]), float(uu[2][1]))         
            if uu[4][1] == 'uniform':
                i_dmin, i_dmax = float(uu[3][0]) * np.ones_like(mch_astro_det), float(uu[3][1]) * np.ones_like(mch_astro_det)
            if uu[4][1] == 'dchirp':
                factor = (mch_astro_det/_mch_BNS)**(5./6)
                i_dmin, i_dmax = float(uu[3][0]) * factor, float(uu[3][1]) * factor
            
            bb = np.sign((max_mchirp - mch_astro_det) * (mch_astro_det - min_mchirp))
            bb += np.sign((i_dmin - astro_lum_dist) * (astro_lum_dist - i_dmax))
    
            bound.append(bb)        
        
        i_within = sum(2 in set(xx) for xx in np.array(bound).T)
        inj_V0 = 4*np.pi*V*i_within/float(nsamples)
        
        injection_chunks[key]['inj_astro_vol'] = inj_V0
        
        for kee in injection_chunks.keys():
            
            if kee == 'z_range':
                continue            
            
            if 'inj_astro_vol' in injection_chunks[kee].keys():
                continue                 
               
            inj_distr_2 = zip(injection_chunks[kee]['mtot_range'], injection_chunks[kee]['m1_range'], 
                              injection_chunks[kee]['m2_range'], injection_chunks[kee]['d_range'], 
                                  zip(injection_chunks[kee]['m_dist'], injection_chunks[kee]['d_dist']))
        
            unique_distr_2 = np.unique(inj_distr_2, axis = 0)   
            
            if np.array_equal(unique_distr, unique_distr_2):
                injection_chunks[kee]['inj_astro_vol'] = inj_V0               
                
    # Estimate the sensitive volume
    z_range = injection_chunks['z_range']
    V_min = quad(lambda z: cosmo.differential_comoving_volume(z).value/(1+z), 0., z_range[0])[0]
    V_max = quad(lambda z: cosmo.differential_comoving_volume(z).value/(1+z), 0., z_range[1])[0]
    def pdf_z_astro(z):
        ''' Get the probability density for the rate of events at a redshift assuming standard cosmology'''
        return cosmo.differential_comoving_volume(z).value/(1+z)/(V_max - V_min)    
    
    for key in injection_chunks.keys():    
        
        if key == 'z_range':
            continue
            
        data = injection_chunks[key]
        i_det, i_inj, i_det_sq, thr_falloff = 0, 0, 0, []
        narrays = len(data['distance'])
        for i in range(narrays):
            mchirp, mass1, mass2 = data['chirp_mass'][i], data['mass1'][i], data['mass2'][i]
            distance, spin1z, spin2z = data['distance'][i], data['spin1z'][i], data['spin2z'][i]
            
            inj_distr = data['inj_distr']
            unique_distr = data['unique_distr']            
            prob_dist, prob_mass = [], []
            
            for uu in unique_distr:
        
                'Get probability density for injections in mass-distance space'
                if uu[4][0] == 'totalMass':
                    low_mass = .5*float(uu[0][0])
                    high_mass = float(uu[0][1]) - low_mass
                    low_mass_2, high_mass_2 = low_mass, high_mass
                if uu[4][0] == 'componentMass' or uu[4][0] == 'log':
                    low_mass, high_mass = float(uu[1][0]), float(uu[1][1])
                    low_mass_2, high_mass_2 = float(uu[2][0]), float(uu[2][1])
            
                mm = inj_mass_pdf(uu[4][0], mass1, mass2, low_mass, high_mass, low_mass_2, high_mass_2)
        
                prob_mass.append(mm)
                dd = inj_distance_pdf(uu[4][1], distance, float(uu[3][0]), float(uu[3][1]), mchirp)
                prob_dist.append(dd)
        
            prob_mass = np.array(prob_mass)
            prob_dist = np.array(prob_dist)         
    
            z_inj = dlum_to_z(distance)
            m1_sc, m2_sc = mass1/(1 + z_inj), mass2/(1 + z_inj)
            p_out = model_pdf(m1_sc, m2_sc, spin1z, spin2z) * pdf_z_astro(z_inj)
        
            p_in = 0
            J = abs(cosmo.luminosity_distance(z_inj + 0.0005).value - cosmo.luminosity_distance(z_inj - 0.0005).value)/0.001
            for ii in range(narrays):
                inj_d = inj_distr[ii]
                hspin1, hspin2 = data['s1_range'][ii][1], data['s2_range'][ii][1]
                prob_spin = inj_spin_pdf(data['s_dist'][ii], hspin1, spin1z) * inj_spin_pdf(data['s_dist'][ii], hspin2, spin2z)
                
                all_idx = np.array([np.array_equal(inj_d, uu) for uu in unique_distr]).astype(float)
                idx = np.where(all_idx == 1)[0][0]
                p_in += prob_mass[idx] * prob_dist[idx] * prob_spin * J * (1 + z_inj)**2 
    
            p_out_in = p_out/p_in
            i_inj += np.sum(p_out_in)
            i_det += np.sum((p_out_in)[data[thr_var][i] > thr_val])
            i_det_sq += np.sum((p_out_in)[data[thr_var][i] > thr_val]**2)
            
            idx_thr = np.where(data[thr_var][i] > thr_val)
            thrs = data[thr_var][i][idx_thr]
            ratios = p_out_in[idx_thr]/max(p_out_in[idx_thr])  
            rndn = np.random.uniform(0, 1, len(ratios))
            idx_ratio = np.where(ratios > rndn)                              
            thr_falloff.append(thrs[idx_ratio])
            
        inj_V0 = injection_chunks[key]['inj_astro_vol']
        injection_chunks[key]['ninj'] = i_inj
        injection_chunks[key]['ndet'] = i_det
        injection_chunks[key]['ndetsq'] = i_det_sq
        injection_chunks[key]['VT'] = ((inj_V0*i_det/i_inj) * (data['gps'][1] - data['gps'][0])/31557600)
        injection_chunks[key]['VT_err'] = injection_chunks[key]['VT'] * np.sqrt(i_det_sq)/i_det
        injection_chunks[key]['thr_falloff'] = np.hstack(np.array(thr_falloff).flat)
    
    return injection_chunks

def process_injections(hdffile, inj_type = 'pipeline'):
    """Function to read in the injection file and extract the found injections and all injections

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

def merge_injections(all_inj):
    """Merge injections across chunks.

       Parameters
       ----------
       all_inj: dictionary
           dictionary containing injections for various chunks

       Returns
       -------
       data: dictionary
           Dictionary containing the merged injections
    """
    injs = {}
    for data in np.append(_save_params, 'found'):
        injs[data] = np.concatenate([inj[data] for inj in tuple(all_inj.values())])
   
    return injs
        
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
    '''Sample the redshifts for sources, with redshift independent rate, using standard cosmology

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
            nsamples of redshift, between min_z, max_z, confirming standard cosmology
    '''
    
    dz, fac = 0.001, 3.0 # use interpolation instead of directly estimating all the pdfz for rndz
    V = quad(lambda z: cosmo.differential_comoving_volume(z).value/(1+z), 0., max_z)[0]
    zs = np.arange(min_z, max_z + dz/2., dz)
    zbins = np.arange(min_z, max_z + dz/2., dz)
    zcenter = (zbins[:-1] + zbins[1:]) / 2
    pdfz = cosmo.differential_comoving_volume(zcenter).value/(1+zcenter)/V

    int_pdf = interpolate.interp1d(zcenter, pdfz, bounds_error=False, fill_value=0)

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

def get_summed_vt(dictionary):
    '''Read the VT results produced from estimate_vt and combine them.

       Parameters
       ----------
       dictionary: dictionary
             Dictionary obtained from estimate_vt function

       Returns
       -------
       sum_vt: float
          Sum of VTs for all the chunks in the dictionary
       float
          Error on summed VT estimated by propogating error
    '''
    sum_vt, sum_vt_err = 0, 0
    for key in dictionary.keys():
    
        if key == 'z_range':
            continue
    
        sum_vt += dictionary[key]['VT']
        sum_vt_err += dictionary[key]['VT_err']**2
    
    return sum_vt, np.sqrt(sum_vt_err) 

def get_accumulated_falloff(dictionary):
    '''collect SNRs from all chunks to get collected fall off.

       Parameters
       ----------
       dictionary: dictionary
             Dictionary obtained from estimate_vt function

       Returns
       -------
       float
          array contaning all the SNRs over the chunks
    '''
    
    thrs = []
    for key in dictionary.keys():
        
        if key == 'z_range':
            continue
        thrs.append(dictionary[key]['thr_falloff'])
        
    return np.hstack(np.array(thrs).flat)

########## Defining current standard strategies used for making injections ##########

def inj_mass_pdf(key, mass1, mass2, low_mass, high_mass, low_mass_2 = 0, high_mass_2 = 0):
    
    '''Estimate the probability density based on the injection strategy

       Parameters
       ----------
       key: string
          Injection strategy
       mass1: array
          First mass of the injections
       mass2: array
          Second mass of the injections
       low_mass: float
          Lower value of the mass distributions 
       high_mass: float
          higher value of the mass distribution

       Returns
       -------
       pdf: array
          Probability density of the injections 
    '''

    mass1, mass2 = np.array(mass1), np.array(mass2)
    
    if key == 'totalMass':
        ''' Returns the PDF of mass when total mass is uniformly distributed. 
            Both the component masses have the same distribution for this case.

            Parameters
            ----------
            low_mass: lower component mass
            high_mass: higher component mass
        '''
        
        bound = np.sign((low_mass + high_mass) - (mass1 + mass2)) 
        bound += np.sign((high_mass - mass1)*(mass1 - low_mass)) + np.sign((high_mass - mass2)*(mass2 - low_mass))
        idx = np.where(bound != 3)
        pdf = 1./(high_mass - low_mass)/(mass1 + mass2 - 2 * low_mass)
        pdf[idx] = 0
        
        return pdf
    
    if key == 'componentMass':
        ''' Returns the PDF of mass when component mass is uniformly distributed. 
            Component masses are independent for this case.

            Parameters
            ----------
            low_mass: lower component mass
            high_mass: higher component mass
         '''
        
        bound = np.sign((high_mass - mass1)*(mass1 - low_mass)) + np.sign((high_mass_2 - mass2)*(mass2 - low_mass_2))
        idx = np.where(bound != 2)
        pdf = np.ones_like(mass1) / (high_mass - low_mass) / (high_mass_2 - low_mass_2)
        pdf[idx] = 0
        
        return pdf
    
    if key == 'log':
        ''' Returns the PDF of mass when component mass is uniform in log. 
            Component masses are independent for this case.

            Parameters
            ----------
            low_mass: lower component mass
            high_mass: higher component mass
         '''  
        
        bound = np.sign((high_mass - mass1)*(mass1 - low_mass)) + np.sign((high_mass_2 - mass2)*(mass2 - low_mass_2))
        idx = np.where(bound != 2)
        pdf = 1 / (np.log(high_mass) - np.log(low_mass)) / (np.log(high_mass_2) - np.log(low_mass_2)) / mass1 / mass2
        pdf[idx] = 0
        
        return pdf   
    
def inj_spin_pdf(key, high_spin, spinz):
    ''' Estimate the probability density of the injections for the spin distribution.

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
        ''' Returns the PDF of spins when total spin is isotropically distributed. 
            Both the component masses have the same distribution for this case.
         '''

        pdf = (np.log(high_spin - np.log(abs(spinz)))/high_spin/2)
        idx = np.where(bound != 2)
        pdf[idx] = 0
        
        return pdf
    
    if key == 'aligned':
        ''' Returns the PDF of mass when spins are aligned and uniformly distributed. 
            Component spins are independent for this case.
        '''

        pdf = (np.ones_like(spinz) / 2 / high_spin)
        idx = np.where(bound != 2)
        pdf[idx] = 0
        
        return pdf   
    
    if key == 'disable_spin':
        ''' Returns unit array '''    
        pdf = np.ones_like(spinz)
        
        return pdf    
    
def inj_distance_pdf(key, distance, low_dist, high_dist, mchirp = 1):
    ''' Estimate the probability density of the injections for the distance distribution.

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
        ''' Returns the PDF at a distance when distance is uniformly distributed.
        '''

        pdf = np.ones_like(distance)/(high_dist - low_dist)
        bound = np.sign((high_dist - distance)*(distance - low_dist)) 
        idx = np.where(bound != 1)
        pdf[idx] = 0
        return pdf

    if key == 'dchirp':
        ''' Returns the PDF at a distance when distance is uniformly distributed but scaled by the chirp mass
        '''
      
        weight = (mchirp/_mch_BNS)**(5./6)
        pdf = np.ones_like(distance) / weight / (high_dist - low_dist)
        bound = np.sign((weight*high_dist - distance)*(distance - weight*low_dist)) 
        idx = np.where(bound != 1)
        pdf[idx] = 0        
        return pdf

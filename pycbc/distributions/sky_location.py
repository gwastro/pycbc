# Copyright (C) 2016  Collin Capano
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

"""This modules provides classes for evaluating sky distributions in
right ascension and declination.
"""

import logging
import numpy # changer les np en numpy
import time #pas besoin enfaite?


from scipy.spatial.transform import Rotation

from pycbc.distributions import angular
from pycbc import VARARGS_DELIM
from pycbc.io import FieldArray

logger = logging.getLogger('pycbc.distributions.sky_location')


class UniformSky(angular.UniformSolidAngle):
    """A distribution that is uniform on the sky. This is the same as
    UniformSolidAngle, except that the polar angle varies from pi/2 (the north
    pole) to -pi/2 (the south pole) instead of 0 to pi. Also, the default
    names are "dec" (declination) for the polar angle and "ra" (right
    ascension) for the azimuthal angle, instead of "theta" and "phi".
    """
    name = 'uniform_sky'
    _polardistcls = angular.CosAngle
    _default_polar_angle = 'dec'
    _default_azimuthal_angle = 'ra'


class FisherSky:
    """A distribution that returns a random angle drawn from an approximate
    `Von_Mises-Fisher distribution`_. Assumes that the Fisher concentration
    parameter is large, so that we can draw the samples from a simple
    rotationally-invariant distribution centered at the North Pole (which
    factors as a uniform distribution for the right ascension, and a Rayleigh
    distribution for the declination, as described in
    `Fabrycky and Winn 2009 ApJ 696 1230`) and then rotate the samples to be
    centered around the specified mean position. As in UniformSky, the
    declination varies from π/2 to -π/2 and the right ascension varies from
    0 to 2π.

    .. _Von_Mises-Fisher distribution:
        http://en.wikipedia.org/wiki/Von_Mises-Fisher_distribution

    .. _Fabrycky and Winn 2009 ApJ 696 1230:
        https://doi.org/10.1088/0004-637X/696/2/1230

    .. _Briggs et al 1999 ApJS 122 503:
        https://doi.org/10.1086/313221

    Parameters
    ----------
    mean_ra: float
        RA of the center of the distribution.
    mean_dec: float
        Declination of the center of the distribution.
    sigma: float
        Spread of the distribution. For the precise interpretation, see Eq 8
        of `Briggs et al 1999 ApJS 122 503`_. This should be smaller than
        about 20 deg for the approximation to be valid.
    angle_unit: str
        Unit for the angle parameters: either "deg" or "rad".
    """
    name = 'fisher_sky'
    _params = ['ra', 'dec']

    def __init__(self, **params):
        if params['angle_unit'] not in ['deg', 'rad']:
            raise ValueError("Only deg or rad is allowed as angle unit")
        mean_ra = params['mean_ra']
        mean_dec = params['mean_dec']
        sigma = params['sigma']
        if params['angle_unit'] == 'deg':
            mean_ra = numpy.deg2rad(mean_ra)
            mean_dec = numpy.deg2rad(mean_dec)
            sigma = numpy.deg2rad(sigma)
        if mean_ra < 0 or mean_ra > 2 * numpy.pi:
            raise ValueError(
                f'The mean RA must be between 0 and 2π, {mean_ra} rad given'
            )
        if mean_dec < -numpy.pi/2 or mean_dec > numpy.pi/2:
            raise ValueError(
                'The mean declination must be between '
                f'-π/2 and π/2, {mean_dec} rad given'
            )
        if sigma <= 0 or sigma > 2 * numpy.pi:
            raise ValueError(
                'Sigma must be positive and smaller than 2π '
                '(preferably much smaller)'
            )
        if sigma > 0.35:
            logger.warning(
                'Warning: sigma = %s rad is probably too large for the '
                'Fisher approximation to be valid', sigma
            )
        self.rayleigh_scale = 0.66 * sigma
        # Prepare a rotation that puts the North Pole at the mean position
        self.rotation = Rotation.from_euler(
            'yz',
            [numpy.pi / 2 - mean_dec, mean_ra]
        )

    @property
    def params(self):
        return self._params

    @classmethod
    def from_config(cls, cp, section, variable_args):
        tag = variable_args
        variable_args = variable_args.split(VARARGS_DELIM)
        if set(variable_args) != set(cls._params):
            raise ValueError("Not all parameters used by this distribution "
                             "included in tag portion of section name")
        mean_ra = float(cp.get_opt_tag(section, 'mean_ra', tag))
        mean_dec = float(cp.get_opt_tag(section, 'mean_dec', tag))
        sigma = float(cp.get_opt_tag(section, 'sigma', tag))
        angle_unit = cp.get_opt_tag(section, 'angle_unit', tag)
        return cls(
            mean_ra=mean_ra,
            mean_dec=mean_dec,
            sigma=sigma,
            angle_unit=angle_unit
        )

    def rvs(self, size):
        # Draw samples from a distribution centered on the North pole
        np_ra = numpy.random.uniform(
            low=0,
            high=(2*numpy.pi),
            size=size
        )
        np_dec = numpy.random.rayleigh(
            scale=self.rayleigh_scale,
            size=size
        )

        # Convert the samples to intermediate cartesian representation
        np_cart = numpy.empty(shape=(size, 3))
        np_cart[:, 0] = numpy.cos(np_ra) * numpy.sin(np_dec)
        np_cart[:, 1] = numpy.sin(np_ra) * numpy.sin(np_dec)
        np_cart[:, 2] = numpy.cos(np_dec)

        # Rotate the samples according to our pre-built rotation
        rot_cart = self.rotation.apply(np_cart)

        # Convert the samples back to spherical coordinates.
        # Some unpleasant conditional operations are needed
        # to get the correct angle convention.
        rot_radec = FieldArray(
            size,
            dtype=[
                ('ra', '<f8'),
                ('dec', '<f8')
            ]
        )
        rot_radec['ra'] = numpy.arctan2(rot_cart[:, 1], rot_cart[:, 0])
        neg_mask = rot_radec['ra'] < 0
        rot_radec['ra'][neg_mask] += 2 * numpy.pi
        rot_radec['dec'] = numpy.arcsin(rot_cart[:, 2])
        return rot_radec




class HealpixSky:
    """Extract the distribution of a healpix map by using the rejection 
    sampling method on the smallest acceptable piece of the celestial sphere.
    Assume the map is not an empty map.
    As in UniformSky and FisherSky, the declination varies from π/2 to -π/2 
    and the right ascension varies from 0 to 2π.

    Parameters
    ----------
    healpix_file : str
        path to a fits file containing probability distribution encoded in a 
        HEALPix scheme
    coverage : {float, 0.9999}
        percentage of the map covered by the method
        
    rasterization_nside : {int, 64}
        nside of the rasterized map used to determine the 
        boundaries of the input map.
    """
    name = 'healpix_sky'
    _params = ['ra', 'dec']
    
    def __init__(self, **params): 
        import healpy 
        import mhealpy
        
        #give the boundaries of a map m in radian, and following the 
        #radec convention delta in [-π/2, π/2] alpha in [0,2π]
        def boundaries(file_name,nside,coverage): 
            healpix_map = mhealpy.HealpixMap.read_map(file_name)
            nside = min(nside,healpix_map.nside) #min si jamais nside < m.nside pour pas ralentir inutilement
            
            delta_max= -numpy.pi/2
            delta_min= numpy.pi/2
            alpha_max= 0
            alpha_min = 2*numpy.pi
            
            rasterized_map = healpix_map.rasterize(scheme = 'NESTED',nside = nside) # scheme pas important ? RING par default #LE MIN ETAIT SEULEMENT ICI
            data = rasterized_map.data
            non_zero_data = data[data != 0]
            renormalization_constant = non_zero_data.sum() 
            
            normalized_data = data/renormalization_constant
            sort_normalized_data = - numpy.sort(-non_zero_data/renormalization_constant) # TRI DECROISSANT DE non_zero_data
            
            N = len(non_zero_data)
            map_coverage = 0
            pix =  numpy.zeros(N,dtype=int)
            delta = numpy.zeros(N)
            alpha =numpy.zeros(N)
           
            for i in range(N):
                j = numpy.argmax(normalized_data == sort_normalized_data[i]) 
                pix[i] = rasterized_map.uniq[j]-4*nside*nside
                delta[i],alpha[i] =  healpy.pix2ang(
                    nside = nside,                                                         #  IL Y AVAIT PEUT ETRE UN PROBLEME AVEC NSIDE QUI ETAIT DIFF SI LE MIN ETAIT UTILISE PLUS HAUT
                    ipix = pix[i],
                    nest = rasterized_map.is_nested,
                    lonlat = False
                )
                # converts colatitude to latitude
                delta[i] = numpy.pi/2 - delta[i] 
                
                map_coverage += sort_normalized_data[i]
                delta_max,delta_min = max(delta_max,delta[i]),min(delta_min,delta[i])
                alpha_max,alpha_min = max(alpha_max,alpha[i]),min(alpha_min,alpha[i])
                    
                if map_coverage > coverage : 
                    break
            #ajout de la securité : taille_pix ~ np.pi/4nside-1
            secu = 2*numpy.pi/(4*nside -1)
            
            delta_max = min(delta_max + secu, numpy.pi/2)
            delta_min = max(delta_min - secu, -numpy.pi/2)
            alpha_max = min(alpha_max + secu, 2*numpy.pi)
            alpha_min = max(alpha_min - secu, 0)
     
            return delta_min,delta_max,alpha_min,alpha_max
       
        
        file_name = params['healpix_file']
        if 'coverage' in params:
            coverage = params['coverage']
        else :
            coverage = 0.9999    
        if coverage > 1 or coverage < 0 : 
            raise ValueError(
                'Coverage must be between 0 and 1'
                )
        if 'rasterization_nside' in params:
            rasterization_nside = params['rasterization_nside']
        else :
            rasterization_nside = 64 # or 128
        if bin(rasterization_nside).count('1') != 1 :
            raise ValueError(
                'Nside must be a power of 2'
                '(preferably between 16 and 256 for better efficiency)'
                )
        self.healpix_map = mhealpy.HealpixMap.read_map(file_name)
        self.boundaries = boundaries(file_name,rasterization_nside,coverage) # self.m à la place de file_name ?
        
    
        

    @property
    def params(self):
        return self._params

    @classmethod
    def from_config(cls, cp, section, variable_args):
        tag = variable_args
        variable_args = variable_args.split(VARARGS_DELIM)
        if set(variable_args) != set(cls._params):
            raise ValueError("Not all parameters used by this distribution "
                             "included in tag portion of section name")
        healpix_file = str(cp.get_opt_tag(section, 'healpix_file', tag))
        coverage = 0.9999
        if cp.has_option_tag(section,'coverage',tag):
            coverage = float(cp.get_opt_tag(section, 'coverage', tag))
        
        rasterization_nside = 64
        if cp.has_option_tag(section,'rasterization_nside',tag):
            rasterization_nside = int(cp.get_opt_tag(section, 'rasterization_nside', tag))
        return cls(
            healpix_file=healpix_file,
            coverage=coverage,
            rasterization_nside=rasterization_nside
        )

    def rvs(self, size):
        
       
        
        # gives arrays of random coordinates in radians within the boundaries
        # alpha being the right ascention and delta the declination
        def uniform_distribution_sphere_partial(size,boundaries):
            
            delta_min,delta_max,alpha_min,alpha_max = boundaries
            
            # angles must be in radians and following the radec convention
            u = numpy.random.uniform(0, 1, size)
            v = numpy.random.uniform(0, 1, size)
            
            delta = numpy.arcsin( 
                ( numpy.sin(delta_max) - numpy.sin(delta_min) ) * u  
                + numpy.sin(delta_min) 
                )  
            alpha = (alpha_max - alpha_min) * v + alpha_min
            

            return alpha, delta  
        #delta in [-pi/2, pi/2] alpha in [0,2pi]
        
        #gives the points accepted by the rejection sampling method
        def simple_rejection_sampling(healpix_map,size,boundaries):                                                  

            data = healpix_map.data
            random_data = numpy.random.uniform(0, data.max(), size)

            alpha,delta = uniform_distribution_sphere_partial(size,boundaries)                                
            theta,phi = numpy.pi/2 -delta,alpha                                           
           
            d_data =  numpy.array( healpix_map.get_interp_val(theta,phi,lonlat = False) )  #  NOUVELLE VERSION DE MHEALPY NECESSAIRE # ARRAY NECESSAIRE POUR NOUVELLE MAPS
            
            dist_theta = theta[ d_data > random_data ] 
            dist_phi = phi[ d_data > random_data ]

            return (dist_theta,dist_phi)
       
        #sampling method
        theta,phi = numpy.array([]),  numpy.array([])
        while len(theta) < size:
            new_theta,new_phi = simple_rejection_sampling(self.healpix_map, #pas de maj
                                                  size,
                                                  self.boundaries
                                                  )
            theta = numpy.concatenate((theta,new_theta), axis = 0) # concatene theta et THETA
            phi   = numpy.concatenate((phi, new_phi), axis = 0)    # concatene phi et PHI

        if len(theta) > size:
            theta = theta[:size]
            phi = phi[:size]
        #end of sampling method
        radec = FieldArray(
            size,
            dtype=[
                ('ra', '<f8'),
                ('dec', '<f8')
            ]
        )
        radec['ra'] = phi
        radec['dec'] = numpy.pi/2 - theta
        return radec



__all__ = ['UniformSky', 'FisherSky', 'HealpixSky']

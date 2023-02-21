class Fisher(UniformSky):
    """A distribution that returns a random (ra, dec) angle drawn from the
    Fisher distribution. Assume that the concentration parameter (kappa)
    is large so that we can use a Rayleigh distribution about the north
    pole and rotate it to be centered at the (ra, dec) coordinate mu.
    Assume kappa = 1 / sigma**2 (kappa should be in units of steradians)
    As in UniformSky, the declination (dec) varies from pi/2 to-pi/2
    and right ascension (ra) varies from 0 to 2pi. And the angles
    should be provided in (ra,dec) format in radians (mu_radians=True),
    rather than factors of pi, or in degrees (mu_radians=False).
    References:
      * http://en.wikipedia.org/wiki/Von_Mises-Fisher_distribution
      * http://arxiv.org/pdf/0902.0737v1 (states the Rayleigh limit)
    """
    name = 'fisher'
    _polardistcls = angular.CosAngle
    _default_polar_angle = 'dec'
    _default_azimuthal_angle = 'ra'

    def __init__(self, ra, dec, kappa, mu_radians=True):
        self.kappa = kappa
        if kappa >= 500:
            if mu_radians:
                mu_values = numpy.array(self.ra, self.dec)
            else:
                mu_values = numpy.array(numpy.deg2rad([self.ra,
                                                            self.dec]))
            mu_values = decra2polaz(self.dec, self.ra)
        else:
            raise ValueError("Kappa too low, minimum should be 500")
            
    @property 
    def rvs_polaz(self, size):
        """
        Randomly draw multiple samples from the Fisher distribution
        and returns (polar, azimuthal) angle values.
        """
        arr = numpy.array([
            numpy.random.rayleigh(scale=1./numpy.sqrt(self.kappa),
                                  size=size),
            numpy.random.uniform(low=0,
                                 high=2*numpy.pi,
                                 size=size)]).reshape((2, size)).T
        alpha, beta = new_z_to_euler(mu_values)
        return rotate_euler(arr, alpha, beta, 0)
    @property
    def rvs_radec(self, size):
        """
        Randomly draw multiple samples from the Fisher distribution
        and returns (ra, dec) values
        """
        rot_eu = self.rvs_polaz(size)
        ra_a = rot_eu[:, 1]
        dec_p = rot_eu[:, 0]
        right_ascension, declination = polaz2radec(dec_p, ra_a)
        return right_ascension, declination

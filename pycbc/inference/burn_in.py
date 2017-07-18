

def max_posterior(sampler, fp):
    """Burn in based on samples being within dim/2 of maximum posterior.

    Parameters
    ----------
    chain_posteriors : array
        (extra dim) x nwalkers x niterations array of log posterior values.
    dim : int
        Number of dimensions of the run (= the number of variable args).

    Returns
    -------
    burnin_idx :
        Array of indices giving the burn-in index for each chain.
    """
    # get the posteriors
    chain_posteriors = sampler.read_samples(fp, ['loglr', 'prior'],
        samples_group=fp.stats_group, thin_interval=1, thin_start=0,
        thin_end=None, flatten=False)
    dim = len(sampler.variable_args)
    # find the posterior to compare against
    maxP = chain_posteriors.max()
    criteria = maxP - dim/2
    nwalkers = chain_posteriors.shape[-2]
    niterations = chain_posteriors.shape[-1]
    burnin_idx = numpy.repeat(niterations, nwalkers).astype(int)
    for ii in range(nwalkers):
        chain = chain_posteriors[...,ii,:]
        idx = numpy.where(chain >= criteria)[-1]
        if idx.size != 0:
            burnin_idx[ii] = idx[0] 
    return burnin_idx


def half_chain(sampler, fp):
    """Takes the second half of the iterations as post-burn in.
    """ 
    nwalkers = sampler.nwalkers
    niterations = fp.niterations
    return numpy.repeat(niterations/2, nwalkers).astype(int)


def kombine(sampler, fp):
    """Uses kombine's burn_in function.
    """
    if sampler.name != 'kombine':
        raise ValueError("can only use kombine burn in with kombine")
    sampler.burn_in()
    return numpy.repeat(sampler.burn_in_iterations,
                        sampler.nwalkers).astype(int)


burnin_functions = {
    'max_posterior': posterior_based,
    'half_chain': half_chain,
    }

class BurnIn(object):
    """Class to estimate the number of burn in iterations.
    """

    def __init__(self, burnin_functions):
        self.burnin_functions = {fname: burnin_funcs[fname]
                                 for fname in burnin_functions}

    def evaluate(self, sampler, fp):
        """Evaluates sampler's chains to find burn in.
        """
        return numpy.vstack([func(sampler, fp)
                for func in self.burnin_functions.values()]).max(axis=0)

    def update_burnin(self, sampler, fp):
        """Evaluates burn in saves updated indices to the given file.
        """
        burnidx = self.evaluate(sampler, fp)
        sampler.write_burnin_iterations(fp, burninidx)


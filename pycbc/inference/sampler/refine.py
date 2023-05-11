""" Sampler that uses kde refinement of an existing posterior estimate.
"""

import logging
import numpy
import numpy.random

from scipy.special import logsumexp
from scipy.stats import gaussian_kde
from scipy.stats import entropy as sentropy

from pycbc.inference import models
from pycbc.pool import choose_pool
from pycbc.inference.io import loadfile

from .base import setup_output, initial_dist_from_config
from .dummy import DummySampler


def call_model(params):
    models._global_instance.update(**params)
    return (models._global_instance.logposterior,
            models._global_instance.loglikelihood)


def resample_equal(samples, logwt, seed=0):
    weights = numpy.exp(logwt - logsumexp(logwt))
    N = len(weights)
    positions = (numpy.random.random() + numpy.arange(N)) / N
    idx = numpy.zeros(N, dtype=int)
    cumulative_sum = numpy.cumsum(weights)
    cumulative_sum /= cumulative_sum[-1]
    i, j = 0, 0
    while i < N:
        if positions[i] < cumulative_sum[j]:
            idx[i] = j
            i += 1
        else:
            j += 1
    try:
        rng = numpy.random.default_rng(seed)
    except AttributeError:
        # numpy pre-1.17 uses RandomState
        # Py27: delete this after we drop python 2.7 support
        rng = numpy.random.RandomState(seed)
    rng.shuffle(idx)
    return {p: samples[p][idx] for p in samples}


class RefineSampler(DummySampler):
    """Sampler for kde drawn refinement of existing posterior estimate

    Parameters
    ----------
    model : Model
        An instance of a model from ``pycbc.inference.models``.
    num_samples: int
        The number of samples to draw from the kde at the conclusion
    iterative_kde_samples: int
        The number of samples to add to the kde during each iterations
    min_refinement_steps: int
        The minimum number of iterations to take.
    max_refinement_steps: The maximum number of refinment steps to take.
    entropy: float
        The target entropy between iterative kdes
    dlogz: float
        The target evidence difference between iterative kde updates
    kde: kde
        The inital kde to use.
    """
    name = 'refine'

    def __init__(self, model, *args, nprocesses=1, use_mpi=False,
                 num_samples=int(1e5),
                 iterative_kde_samples=int(1e3),
                 min_refinement_steps=5,
                 max_refinement_steps=40,
                 entropy=0.001,
                 dlogz=0.01,
                 kde=None,
                 **kwargs):
        super().__init__(model, *args)

        self.model = model
        self.kde = kde
        self.vparam = model.variable_params
        models._global_instance = model
        self.num_samples = int(num_samples)
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        self._samples = {}

        self.num_samples = int(num_samples)
        self.iterative_kde_samples = int(iterative_kde_samples)
        self.min_refinement_steps = int(min_refinement_steps)
        self.max_refinement_steps = int(max_refinement_steps)
        self.entropy = float(entropy)
        self.dlogz_target = float(dlogz)

    def draw_samples(self, size):
        """Draw new samples within the model priors"""
        logging.info('getting from kde')
        ksamples = self.kde.resample(size=size)
        params = {k: ksamples[i, :] for i, k in enumerate(self.vparam)}
        logging.info('checking prior')
        keep = self.model.prior_distribution.contains(params)
        return ksamples[:, keep]

    @staticmethod
    def compare_kde(kde1, kde2, size=int(1e4)):
        """ Calculate information difference between two kde distributions
        """
        s = kde1.resample(size=size)
        return sentropy(kde1.pdf(s), kde2.pdf(s))

    def converged(self, step, kde_new, factor):
        """ Check that kde is converged by comparing to previous iteration
        """
        if not hasattr(self, 'old_logz'):
            self.old_logz = numpy.inf

        entropy_diff = self.compare_kde(self.kde, kde_new)

        # Compare how the logz changes when adding new samples
        # this is guaranteed to decrease as old samples included
        logz = logsumexp(factor) - numpy.log(len(factor))
        dlogz = logz - self.old_logz
        self.old_logz = logz

        # compare evidence subsets agree
        choice2 = numpy.random.choice(factor, len(factor) // 2)
        logz2 = logsumexp(choice2) - numpy.log(len(choice2))
        choice3 = numpy.random.choice(factor, len(factor) // 2)
        logz3 = logsumexp(choice3) - numpy.log(len(choice3))
        dlogz2 = logz3 - logz2

        logging.info('%s: Checking convergence: dlogz_iter=%.4f,'
                     'dlogz_half=%.4f, entropy=%.4f',
                     step,  dlogz, dlogz2, entropy_diff)
        if (entropy_diff < self.entropy and step >= self.min_refinement_steps
            and max(abs(dlogz), abs(dlogz2)) < self.dlogz_target):
            return True
        else:
            return False

    @classmethod
    def from_config(cls, cp, model, output_file=None, nprocesses=1,
                    use_mpi=False):
        """This should initialize the sampler given a config file.
        """
        kwargs = {k: cp.get('sampler', k) for k in cp.options('sampler')}
        obj = cls(model, nprocesses=nprocesses, use_mpi=use_mpi, **kwargs)
        obj.set_start_from_config(cp)
        setup_output(obj, output_file, check_nsamples=False, validate=False)
        return obj

    def set_start_from_config(self, cp):
        """Sets the initial state of the sampler from config file
        """
        if cp.has_option('sampler', 'start-file'):
            start_file = cp.get('sampler', 'start-file')
            logging.info("Using file %s for initial positions", start_file)
            samples = loadfile(start_file, 'r').read_samples(self.vparam)
        else:
            init_prior = initial_dist_from_config(
                cp, self.model.variable_params, self.model.static_params)
            if init_prior is not None:
                samples = init_prior.rvs(size=self.iterative_kde_samples)
            else:
                p = self.model.prior_distribution
                samples = p.rvs(size=self.iterative_kde_samples)

        ksamples = numpy.array([samples[v] for v in self.vparam])
        self.kde = gaussian_kde(ksamples)

    def run_samples(self, ksamples):
        """ Calculate the likelihoods and weights for a set of samples
        """
        # Calculate likelihood for each sample
        args = []
        for i in range(len(ksamples[0])):
            param = {k: ksamples[j][i] for j, k in enumerate(self.vparam)}
            args.append(param)

        result = self.pool.map(call_model, args)
        logp = numpy.array([r[0] for r in result])
        logl = numpy.array([r[1] for r in result])
        logw = logp - numpy.log(self.kde.pdf(ksamples))

        k = logp != - numpy.inf
        ksamples = ksamples[:, k]
        logp, logl, logw = logp[k], logl[k], logw[k]
        return ksamples, logp, logl, logw

    def run(self):
        """ Iterative sample from kde and update based on likelihood values
        """
        total_samples = None
        total_logp = None
        total_logw = None
        total_logl = None

        for r in range(self.max_refinement_steps):
            logging.info('calculating likelihoods...')
            ksamples = self.draw_samples(self.iterative_kde_samples)
            ksamples, logp, logl, logw = self.run_samples(ksamples)

            logging.info('..done')

            if total_samples is not None:
                total_samples = numpy.concatenate([total_samples,
                                                   ksamples], axis=1)
                total_logp = numpy.concatenate([total_logp, logp])
                total_logw = numpy.concatenate([total_logw, logw])
                total_logl = numpy.concatenate([total_logl, logl])
            else:
                total_samples = ksamples
                total_logp = logp
                total_logw = logw
                total_logl = logl

            logging.info('setting up next kde iteration..')
            ntotal_logw = total_logw - logsumexp(total_logw)
            kde_new = gaussian_kde(total_samples,
                                   weights=numpy.exp(ntotal_logw))
            logging.info('done')
            if self.converged(r, kde_new, total_logl + total_logw):
                break

            self.kde = kde_new

        logging.info('Drawing final samples')
        ksamples = self.draw_samples(self.num_samples)
        logging.info('Calculating final likelihoods')
        ksamples, logp, logl, logw = self.run_samples(ksamples)
        self._samples = {k: ksamples[j,:] for j, k in enumerate(self.vparam)}
        self._samples['loglikelihood'] = logl
        logging.info("Reweighting to equal samples")
        self._samples = resample_equal(self._samples, logw)

""" Direct monte carlo sampling using pregenerated mapping files that
encode the intrinsic parameter space.
"""
import numpy, h5py, logging, tqdm
from numpy.random import choice
from scipy.special import logsumexp
from pycbc.io import FieldArray

from pycbc.inference.io import PosteriorFile
from pycbc.inference import models
from pycbc.pool import choose_pool
from .dummy import DummySampler

from .base import (BaseSampler, setup_output)

def call_likelihood(params):
    """ Accessor to update the global model and call its reconstruction
    routine.
    """
    models._global_instance.update(**params)
    return models._global_instance.loglikelihood

def call_rlikelihood(params):
    """ Accessor to update the global model and call its reconstruction
    routine.
    """
    models._global_instance.update(**params)
    logl = models._global_instance.loglikelihood
    models._global_instance.update(**params)
    rec = models._global_instance.reconstruct()
    return logl, rec

class NetBank(DummySampler):
    """Direct monte-carlo sampler using a preconstructed parameter space
    mapping file.

    Parameters
    ----------
    model : Model
        An instance of a model from ``pycbc.inference.models``.
    """
    name = 'net_bank'

    def __init__(self, model, *args, nprocesses=1, use_mpi=False,
                 mapfile=None,
                 loglr_region=25,
                 target_likelihood_calls=1e5,
                 resample_draw_size=2e5,
                 reconstruct = False,
                 **kwargs):
        super().__init__(model, *args)

        self.mapfile = mapfile

        models._global_instance = model
        self.model = model
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        self._samples = {}

        self.target_likelihood_calls = int(target_likelihood_calls)
        self.reconstruct = bool(reconstruct)
        self.loglr_region = float(loglr_region)
        self.resample_draw_size = int(float(resample_draw_size))

    def run(self):
        logging.info('Retrieving params of parameter space nodes')
        with h5py.File(self.mapfile, 'r') as f:
            bparams = {p: f['bank'][p][:] for p in self.variable_params}
            dmap = {}
            num_nodes = len(bparams[list(bparams.keys())[0]])
            lengths = numpy.array([len(f['map'][str(x)]) for x in range(num_nodes)])
            dtype = f['map']['0'].dtype

        logging.info('Calculating likelihood at nodes')
        args = []
        for i in range(num_nodes):
            pset = {p: bparams[p][i] for p in self.model.variable_params}
            args.append(pset)

        node_loglrs = list(tqdm.tqdm(self.pool.imap(call_likelihood, args),
                                     total=len(args)))
        node_loglrs = numpy.array(node_loglrs)
        loglr_bound = node_loglrs.max() - self.loglr_region

        logging.info('Drawing proposal samples from node regions')
        logw = node_loglrs + numpy.log(lengths)
        passed = numpy.where(node_loglrs > loglr_bound)[0]
        logw2 = logw[passed]
        logw2 -= logsumexp(logw2)
        weight = numpy.exp(logw2)

        logging.info("...reading template bins")
        with h5py.File(self.mapfile, 'r') as f:
            for i in passed:
                dmap[i] = f['map'][str(i)][:]

        logging.info("...draw template bins")
        drawcount = (weight * self.target_likelihood_calls).astype(int)

        dorder = node_loglrs[passed].argsort()[::-1]
        remainder = 0
        for i in dorder:
            c = drawcount[i]
            l = lengths[passed][i]
            if c > l:
                drawcount[i] = l
                remainder += c - l
            elif c < l:
                asize = min(l - c, remainder)
                drawcount[i] += asize
                remainder -= asize
        drawweight = weight / drawcount
        total_draw = drawcount.sum()

        logging.info('...drawn random points within bins')
        psamp = FieldArray(total_draw, dtype=dtype)
        pweight = numpy.zeros(total_draw, dtype=float)
        j = 0
        for i, c, w in zip(passed, drawcount, drawweight):
            bdraw = choice(dmap[i], size=c, replace=False)
            psamp[j:j+len(bdraw)] = FieldArray.from_records(bdraw, dtype=dtype)
            pweight[j:j+len(bdraw)] = - node_loglrs[i] + numpy.log(w)
            j += len(bdraw)

        logging.info("Possible unique values %s", lengths[passed].sum())
        logging.info("Templates drawn from %s", len(passed))
        logging.info("Unique values first draw %s", len(psamp))

        # Calculate the likelihood values for the unique parameter space
        # points
        args = []
        for i, s in enumerate(psamp):
            pset = {p: psamp[p][i] for p in self.model.variable_params}
            args.append(pset)

        if self.reconstruct:
            logging.info('Reconstructing any marginalized params...')
            res = list(tqdm.tqdm(self.pool.imap(call_rlikelihood, args),
                       total=len(args)))
            loglr_samp = numpy.array([r[0] for r in res])
            rec = [r[1] for r in res]
        else:
            loglr_samp = numpy.array(list(tqdm.tqdm(self.pool.imap(call_likelihood, args),
                                     total=len(args))))

        # Draw samples based on the actual likelihood relative to the
        # initial weights
        logw3 = loglr_samp + pweight
        logw3 -= logsumexp(logw3)
        weight2 = numpy.exp(logw3)

        ess = 1.0 / (weight2 ** 2.0).sum()
        logging.info("ESS = %s", ess)

        draw2 = choice(len(psamp), size=int(ess * 5),
                       replace=True, p=weight2)
        logging.info("Unique values second draw %s",
                     len(numpy.unique(psamp[draw2])))

        # Prepare the equally weighted output samples
        fsamp = FieldArray(len(draw2), dtype=dtype)
        for i in range(len(draw2)):
            fsamp[i] = psamp[draw2[i]]

        self._samples = {p: fsamp[p] for p in self.model.variable_params}
        self._samples['loglikelihood'] = loglr_samp[draw2]

        if self.reconstruct:
            for k in rec[0]:
                values = [r[k] for r in rec]
                self._samples[k] = numpy.array(values)[draw2]

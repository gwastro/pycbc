""" Direct monte carlo sampling using pregenerated mapping files that
encode the intrinsic parameter space.
"""
import numpy, h5py, logging
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
                 node_draw_size=1e6,
                 resample_draw_size=2e5,
                 **kwargs):
        super().__init__(model, *args)

        self.mapfile = mapfile

        models._global_instance = model
        self.model = model
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        self._samples = {}
          
        self.loglr_region = float(loglr_region)
        self.node_draw_size = int(float(node_draw_size))
        self.resample_draw_size = int(float(resample_draw_size))

    def run(self):
        logging.info('Retrieving params of parameter space nodes')
        with h5py.File(self.mapfile, 'r') as f:
            bparams = {p: f['bank'][p][:] for p in self.variable_params}
            dmap = []
            num_nodes = len(bparams[list(bparams.keys())[0]])
            for i in range(num_nodes):
                dmap.append(f['map'][str(i)][:])
        dtype = dmap[0].dtype     
        lengths = numpy.array([len(x) for x in dmap])  
        
        logging.info('Calculating likelihood at nodes')
        args = []
        for i in range(num_nodes):
            pset = {p: bparams[p][i] for p in self.model.variable_params}
            args.append(pset)
            
        node_loglrs = numpy.array(self.pool.map(call_likelihood, args))
        loglr_bound = node_loglrs.max() - self.loglr_region

        logging.info('Drawing proposal samples from node regions')
        logw = node_loglrs + numpy.log(lengths)
        passed = numpy.where(node_loglrs > loglr_bound)[0]
        logw2 = logw[passed]
        logw2 -= logsumexp(logw2)
        weight = numpy.exp(logw2)

        logging.info("...draw template bins")
        draw = choice(passed, size=self.node_draw_size, replace=True, p=weight)

        logging.info('...drawn random points within bins')
        psamp = FieldArray(self.node_draw_size, dtype=dtype)
        for i in range(len(psamp)):
            psamp[i] = choice(dmap[draw[i]])

        upsamp, expand = numpy.unique(psamp, return_inverse=True)
        logging.info("Possible unique values %s", lengths[passed].sum())
        logging.info("Templates drawn from %s", len(numpy.unique(draw)))
        logging.info("Unique values first draw %s", len(upsamp))

        # Calculate the likelihood values for the unique parameter space
        # points
        args = []
        for i, s in enumerate(upsamp):
            pset = {p: upsamp[p][i] for p in self.model.variable_params}
            args.append(pset)
        loglr_samp = numpy.array(self.pool.map(call_likelihood, args))

        # Draw samples based on the actual likelihood relative to the
        # initial weights
        logw3 = loglr_samp[expand] - numpy.array(node_loglrs)[draw]
        logw3 -= logsumexp(logw3)
        weight2 = numpy.exp(logw3)

        draw2 = choice(len(psamp), size=self.resample_draw_size,
                       replace=True, p=weight2)
        logging.info("Unique values second draw %s",
                     len(numpy.unique(psamp[draw2])))
        logging.info("ESS = %s", 1.0 / (weight2 ** 2.0).sum())

        # Prepare the equally weighted output samples
        fsamp = FieldArray(len(draw2), dtype=dtype)
        for i in range(len(draw2)):
            fsamp[i] = psamp[draw2[i]]

        self._samples = {p: fsamp[p] for p in self.model.variable_params}
        self._samples['loglikelihood'] = loglr_samp[expand][draw2]

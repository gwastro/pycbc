""" Dummy class when no actual sampling is needed, but we may want to do
some reconstruction supported by the likelihood model.
"""

import numpy, tqdm, h5py, logging
from numpy.random import choice
from scipy.special import logsumexp
from pycbc.io import FieldArray

from pycbc.inference.io import PosteriorFile
from pycbc.inference import models
from pycbc.pool import choose_pool
from .dummy import DummySampler

from .base import (BaseSampler, setup_output)

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
                 num_samples=1000, 
                 bank=None,
                 mapfile=None,
                 **kwargs):
        super().__init__(model, *args)

        self.bank = bank
        self.mapfile = mapfile

        models._global_instance = model
        self.model = model
        
        self.num_samples = int(num_samples)
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        self._samples = {}
          

    def run(self):
        logging.info('Retrieving params of parameter space nodes')
        bparams = {}
        with h5py.File(self.bank, 'r') as f:
            for p in self.model.variable_params:
                bparams[p] = f[p][:]
 
        logging.info('formatting input')
        f = h5py.File(self.mapfile, 'r')

        idx = f['idx'][:]
        match = f['match'][:]
        params = f['params'][:]

        dmap = []
        for i in tqdm.tqdm(range(len(bparams[p]))):
            dmap.append([])
            select =  params[numpy.where(idx == i)[0]]
            dmap[i] = select  

        self.dmap = dmap
        self.lengths = numpy.array([len(x) for x in dmap])
        f.close()   
        
        logging.info('Calculating likelihood at nodes')
        node_loglrs = []
        for i in tqdm.tqdm(range(len(bparams[p]))):
            pset = {p: bparams[p][i] for p in self.model.variable_params}
            self.model.update(**pset)
            node_loglrs.append(self.model.loglr)
        node_loglrs = numpy.array(node_loglrs)
        loglr_bound = node_loglrs.max() - 25

        logging.info('Drawing proposal samples from node regions')
        logw = node_loglrs + numpy.log(self.lengths)
        passed = numpy.where(node_loglrs > loglr_bound)[0]
        logw2 = logw[passed]
        logw2 -= logsumexp(logw2)
        weight = numpy.exp(logw2)

        size = int(1e6)
        logging.info("...draw template bins")
        draw = choice(passed, size=size, replace=True, p=weight)

        logging.info('...drawn random points within bins')
        psamp = FieldArray(size, dtype=params.dtype)
        for i in range(len(psamp)):
            psamp[i] = choice(dmap[draw[i]])

        logw3 = numpy.zeros(size)
        upsamp, expand = numpy.unique(psamp, return_inverse=True)

        logging.info("Possible unique values %s", self.lengths[passed].sum())
        logging.info("Templates drawn from %s", len(numpy.unique(draw)))
        logging.info("Unique values first draw %s", len(upsamp))

        for i, s in tqdm.tqdm(enumerate(upsamp), total=len(upsamp), position=0, leave=True):
            self.model.update(**{p: upsamp[p][i] for p in self.model.variable_params}) 
            ll =  self.model.logl
            logw3[i] += ll

        logw3 = logw3[expand] - numpy.array(node_loglrs)[draw]

        logw3 -= logsumexp(logw3)
        weight2 = numpy.exp(logw3)

        draw2 = choice(len(psamp), size=200000, replace=True, p=weight2)
        logging.info("Unique values second draw %s", len(numpy.unique(psamp[draw2])))
        logging.info("ESS = %s", 1.0 / (weight2 ** 2.0).sum())

        fsamp = FieldArray(len(draw2), dtype=params.dtype)
        for i in range(len(draw2)):
            fsamp[i] = psamp[draw2[i]]

        self._samples = {p: fsamp[p] for p in self.model.variable_params}

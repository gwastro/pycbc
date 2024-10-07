""" Direct monte carlo sampling using pregenerated mapping files that
encode the intrinsic parameter space.
"""
import logging
import tqdm
import h5py
import numpy
import numpy.random
from scipy.special import logsumexp
from pycbc.io import FieldArray

from pycbc.inference import models
from pycbc.pool import choose_pool
from .dummy import DummySampler


def call_likelihood(params):
    """ Accessor to update the global model
    """
    models._global_instance.update(**params)
    return models._global_instance.loglikelihood


class OutOfSamples(Exception):
    """Exception if we ran out of samples"""


class GameSampler(DummySampler):
    """Direct importance sampling using a preconstructed parameter space
    mapping file.

    Parameters
    ----------
    model : Model
        An instance of a model from ``pycbc.inference.models``.
    mapfile : str
        Path to the pre-generated file containing the pre-mapped prior volume
    loglr_region: int
        Only use regions from the prior volume tiling that are within
        this loglr difference of the maximum tile.
    target_likelihood_calls: int
        Try to use this many likelihood calls in each round of the analysis.
    rounds: int
        The number of iterations to use before terminated.
    """
    name = 'games'

    def __init__(self, model, *args, nprocesses=1, use_mpi=False,
                 mapfile=None,
                 loglr_region=25,
                 target_likelihood_calls=1e5,
                 rounds=1,
                 **kwargs):
        super().__init__(model, *args)

        self.meta = {}
        self.mapfile = mapfile
        self.rounds = int(rounds)
        self.dmap = {}
        self.draw = {}

        models._global_instance = model
        self.model = model
        self.pool = choose_pool(mpi=use_mpi, processes=nprocesses)
        self._samples = {}

        self.target_likelihood_calls = int(target_likelihood_calls)
        self.loglr_region = float(loglr_region)

    def draw_samples_from_bin(self, i, size):
        """ Get samples from the binned prior space """
        if i not in self.draw:
            self.draw[i] = numpy.arange(0, len(self.dmap[i]))

        if size > len(self.draw[i]):
            raise OutOfSamples

        numpy.random.shuffle(self.draw[i])
        selected = self.draw[i][:size]
        self.draw[i] = self.draw[i][size:]

        if size > 0:
            remain = len(self.draw[i])
            logging.info('Drew %i, %i remains in bin %i', size, remain, i)

        return self.dmap[i][selected]

    def sample_round(self, bin_weight, node_idx, lengths):
        """ Sample from the posterior using pre-binned sets of points and
        the weighting factor of each bin.

        bin_weight: Array
            The weighting importance factor of each bin of the prior space
        node_idx: Array
            The set of ids into the prebinned prior volume to use. This should
            map to the given weights.
        lengths: Array
            The size of each bin, used to self-normalize
        """
        logging.info("...draw template bins")
        drawcount = (bin_weight * self.target_likelihood_calls).astype(int)

        dorder = bin_weight.argsort()[::-1]
        remainder = 0
        for i in dorder:
            bincount = drawcount[i]
            binlen = lengths[i]
            if bincount > binlen:
                drawcount[i] = binlen
                remainder += bincount - binlen
            elif bincount < binlen:
                asize = min(binlen - bincount, remainder)
                drawcount[i] += asize
                remainder -= asize

        drawweight = bin_weight / drawcount
        total_draw = drawcount.sum()

        logging.info('...drawn random points within bins')
        psamp = FieldArray(total_draw, dtype=self.dtype)
        pweight = numpy.zeros(total_draw, dtype=float)
        bin_id = numpy.zeros(total_draw, dtype=int)
        j = 0
        for i, (c, w) in enumerate(zip(drawcount, drawweight)):
            bdraw = self.draw_samples_from_bin(node_idx[i], c)
            psamp[j:j+len(bdraw)] = FieldArray.from_records(bdraw,
                                                            dtype=self.dtype)
            pweight[j:j+len(bdraw)] = numpy.log(bin_weight[i]) - numpy.log(w)
            bin_id[j:j+len(bdraw)] = i
            j += len(bdraw)

        logging.info("Possible unique values %s", lengths.sum())
        logging.info("Templates drawn from %s", len(lengths))
        logging.info("Unique values first draw %s", len(psamp))

        # Calculate the likelihood values for the unique parameter space
        # points
        args = []
        for i in range(len(psamp)):
            pset = {p: psamp[p][i] for p in self.model.variable_params}
            args.append(pset)

        loglr_samp = list(tqdm.tqdm(self.pool.imap(call_likelihood, args),
                                    total=len(args)))
        loglr_samp = numpy.array(loglr_samp)

        # Calculate the weights from the actual likelihood relative to the
        # initial weights
        logw3 = loglr_samp + numpy.log(lengths[bin_id]) - pweight
        logw3 -= logsumexp(logw3)
        weight2 = numpy.exp(logw3)
        return psamp, loglr_samp, weight2, bin_id

    def run(self):
        """ Produce posterior samples """
        logging.info('Retrieving params of parameter space nodes')
        with h5py.File(self.mapfile, 'r') as mapfile:
            bparams = {p: mapfile['bank'][p][:] for p in self.variable_params}
            num_nodes = len(bparams[list(bparams.keys())[0]])
            lengths = numpy.array([len(mapfile['map'][str(x)])
                                   for x in range(num_nodes)])
            self.dtype = mapfile['map']['0'].dtype

        logging.info('Calculating likelihood at nodes')
        args = []
        for i in range(num_nodes):
            pset = {p: bparams[p][i] for p in self.model.variable_params}
            args.append(pset)

        node_loglrs = list(tqdm.tqdm(self.pool.imap(call_likelihood, args),
                                     total=len(args)))
        node_loglrs = numpy.array(node_loglrs)
        loglr_bound = node_loglrs[~numpy.isnan(node_loglrs)].max()
        loglr_bound -= self.loglr_region

        logging.info('Drawing proposal samples from node regions')
        logw = node_loglrs + numpy.log(lengths)
        passed = (node_loglrs > loglr_bound) & ~numpy.isnan(node_loglrs)
        passed = numpy.where(passed)[0]
        logw2 = logw[passed]
        logw2 -= logsumexp(logw2)
        weight = numpy.exp(logw2)

        logging.info("...reading template bins")
        with h5py.File(self.mapfile, 'r') as mapfile:
            for i in passed:
                self.dmap[i] = mapfile['map'][str(i)][:]

        # Sample from posterior
        psamp = None
        loglr_samp = None
        weight2 = None
        bin_ids = None

        weight = lengths[passed] / lengths[passed].sum()

        for i in range(self.rounds):
            try:
                psamp_v, loglr_samp_v, weight2_v, bin_id = \
                    self.sample_round(weight / weight.sum(),
                                      passed, lengths[passed])
            except OutOfSamples:
                logging.info("No more samples to draw from")
                break

            for j, v in zip(bin_id, weight2_v):
                weight[j] += v

            if psamp is None:
                psamp = psamp_v
                loglr_samp = loglr_samp_v
                weight2 = weight2_v
                bin_ids = bin_id
            else:
                psamp = numpy.concatenate([psamp_v, psamp])
                loglr_samp = numpy.concatenate([loglr_samp_v, loglr_samp])
                weight2 = numpy.concatenate([weight2_v, weight2])
                bin_ids = numpy.concatenate([bin_id, bin_ids])

            ess = 1.0 / ((weight2/weight2.sum()) ** 2.0).sum()
            logging.info("ESS = %s", ess)

        # Prepare the equally weighted output samples
        self.meta['ncalls'] = len(weight2)
        self.meta['ess'] = ess

        weight2 /= weight2.sum()
        draw2 = numpy.random.choice(len(psamp), size=int(ess * 1),
                                    replace=True, p=weight2)
        logging.info("Unique values second draw %s",
                     len(numpy.unique(psamp[draw2])))

        fsamp = FieldArray(len(draw2), dtype=self.dtype)
        for i, v in enumerate(draw2):
            fsamp[i] = psamp[v]

        self._samples = {p: fsamp[p] for p in self.model.variable_params}
        self._samples['loglikelihood'] = loglr_samp[draw2]

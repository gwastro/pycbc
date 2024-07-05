""" Direct monte carlo sampling using pregenerated mapping files that
encode the intrinsic parameter space.
"""
import numpy, h5py, logging, tqdm
import numpy.random
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

class OutOfSamples(Exception):
    pass

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
                 reconstruct = False,
                 rounds = 1,
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
        self.reconstruct = bool(reconstruct)
        self.loglr_region = float(loglr_region)

    def draw_samples_from_bin(self, i, size):
        if i not in self.draw:
            self.draw[i] = numpy.arange(0, len(self.dmap[i]))

        if size > len(self.draw[i]):
            raise OutOfSamples

        numpy.random.shuffle(self.draw[i])
        selected = self.draw[i][:size]
        self.draw[i] = self.draw[i][size:]

        if size > 0:
            remain = len(self.draw[i])
            logging.info(f'Drew {size}, {remain} remains in bin {i}')

        return self.dmap[i][selected]

    def sample_round(self, bin_weight, node_idx, lengths):
        logging.info("...draw template bins")
        drawcount = (bin_weight * self.target_likelihood_calls).astype(int)

        dorder = bin_weight.argsort()[::-1]
        remainder = 0
        for i in dorder:
            c = drawcount[i]
            l = lengths[i]
            if c > l:
                drawcount[i] = l
                remainder += c - l
            elif c < l:
                asize = min(l - c, remainder)
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
            psamp[j:j+len(bdraw)] = FieldArray.from_records(bdraw, dtype=self.dtype)
            pweight[j:j+len(bdraw)] = numpy.log(bin_weight[i]) - numpy.log(w)
            bin_id[j:j+len(bdraw)] = i
            j += len(bdraw)

        logging.info("Possible unique values %s", lengths.sum())
        logging.info("Templates drawn from %s", len(lengths))
        logging.info("Unique values first draw %s", len(psamp))

        # Calculate the likelihood values for the unique parameter space
        # points
        args = []
        for i, s in enumerate(psamp):
            pset = {p: psamp[p][i] for p in self.model.variable_params}
            args.append(pset)

        loglr_samp = numpy.array(list(tqdm.tqdm(self.pool.imap(call_likelihood, args),
                                     total=len(args))))

        # Calculate the weights from the actual likelihood relative to the
        # initial weights
        logw3 = loglr_samp + numpy.log(lengths[bin_id]) - pweight
        logw3 -= logsumexp(logw3)
        weight2 = numpy.exp(logw3)
        return psamp, loglr_samp, weight2, bin_id

    def calculate_evidence(self, node_loglrs, node_lengths, loglrsamp, bin_id):
        logz = node_loglrs * 1

        #s = bin_id.argsort()
        #bin_id = bin_id[s]
        #loglrsamp = loglrsamp[s]

        bid = numpy.unique(bin_id)
        for i in bid:
            x = numpy.where(i == bin_id)[0]
            logz[i] = logsumexp(loglrsamp[x]) - numpy.log(len(x))

        logz = (logz + numpy.log(node_lengths) - numpy.log(node_lengths.sum()))
        return logsumexp(logz)

    def run(self):
        logging.info('Retrieving params of parameter space nodes')
        with h5py.File(self.mapfile, 'r') as f:
            bparams = {p: f['bank'][p][:] for p in self.variable_params}
            num_nodes = len(bparams[list(bparams.keys())[0]])
            lengths = numpy.array([len(f['map'][str(x)]) for x in range(num_nodes)])
            self.dtype = f['map']['0'].dtype

        logging.info('Calculating likelihood at nodes')
        args = []
        for i in range(num_nodes):
            pset = {p: bparams[p][i] for p in self.model.variable_params}
            args.append(pset)

        node_loglrs = list(tqdm.tqdm(self.pool.imap(call_likelihood, args),
                                     total=len(args)))
        node_loglrs = numpy.array(node_loglrs)
        loglr_bound = node_loglrs[~numpy.isnan(node_loglrs)].max() - self.loglr_region

        logging.info('Drawing proposal samples from node regions')
        logw = node_loglrs + numpy.log(lengths)
        passed = numpy.where((node_loglrs > loglr_bound)  & ~numpy.isnan(node_loglrs))[0]
        logw2 = logw[passed]
        logw2 -= logsumexp(logw2)
        weight = numpy.exp(logw2)
        node_loglrsp = node_loglrs[passed]

        logging.info("...reading template bins")
        with h5py.File(self.mapfile, 'r') as f:
            for i in passed:
                self.dmap[i] = f['map'][str(i)][:]

        # Sample from posterior

        psamp = None
        loglr_samp = None
        weight2 = None
        bin_ids = None

        weight = lengths[passed] / lengths[passed].sum()

        for i in range(self.rounds):
            try:
                psamp_v, loglr_samp_v, weight2_v, bin_id = self.sample_round(weight/weight.sum(), passed, lengths[passed])
            except OutOfSamples:
                logging.info("No more samples to draw from")
                break

            for i, v in zip(bin_id, weight2_v):
                weight[i] += v

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

            uniq = numpy.unique(psamp)

            ess = 1.0 / ((weight2/weight2.sum()) ** 2.0).sum()
            logging.info("ESS = %s", ess)

            logz = self.calculate_evidence(node_loglrs, lengths, loglr_samp, bin_ids)
            dlogz = (loglr_samp * weight2).sum() / weight2.sum() + numpy.log(len(bin_id)) - numpy.log(lengths.sum())
            logging.info("dlogz = %s, %s", dlogz, logz)

        # Prepare the equally weighted output samples
        self.meta['ncalls'] = len(weight2)
        self.meta['ess'] = ess

        weight2 /= weight2.sum()
        draw2 = numpy.random.choice(len(psamp), size=int(ess * 1),
                       replace=True, p=weight2)
        logging.info("Unique values second draw %s",
                     len(numpy.unique(psamp[draw2])))

        fsamp = FieldArray(len(draw2), dtype=self.dtype)
        for i in range(len(draw2)):
            fsamp[i] = psamp[draw2[i]]

        self._samples = {p: fsamp[p] for p in self.model.variable_params}
        self._samples['loglikelihood'] = loglr_samp[draw2]

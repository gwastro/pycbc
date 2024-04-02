# Copyright (C) 2022  Collin Capano
#               2023  Alex Nitz & Shichao Wu
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

"""Hierarchical model definitions."""

import shlex
import logging
import numpy
from pycbc import transforms
from pycbc.workflow import WorkflowConfigParser
from .base import BaseModel

#
# =============================================================================
#
#                       Hierarhical model definition
#
# =============================================================================
#


class HierarchicalModel(BaseModel):
    r"""Model that is a combination of other models.

    Sub-models are treated as being independent of each other, although
    they can share parameters. In other words, the hierarchical likelihood is:

    .. math::

        p(\mathbf{D}|\mathbf{\vartheta}, \mathbf{H}) =
            \prod_{I}^{K} p(\mathbf{d}_I|\mathbf{\vartheta}, H_{I})

    Submodels are provided as a dictionary upon initialization with a unique
    label assigned to each model, e.g., ``{'event1' -> model1, 'event2' ->
    model2}``. Variable and static parameters that are specific to each
    submodel should be prepended with ``{label}__``, where ``{label}__`` is the
    label associated with the given submodel. Shared parameters across multiple
    models have no labels prepended. To specify shared models over a subset of
    models, separate models with an underscore.  For example,
    ``event1_event2__foo`` will result in ``foo`` being common between models
    ``event1`` and ``event2``. For more details on parameter naming see
    :py:class:`HierarchicalParam
    <pycbc.inference.models.hierarchical.HierarchicalParam>`.

    All waveform and sampling transforms, as well as prior evaluation, are
    handled by this model, not the sub-models. Parameters created by waveform
    transforms should therefore also have sub-model names prepended to them,
    to indicate which models they should be provided to for likelihood
    evaluation.

    Parameters
    ----------
    variable_params: (tuple of) string(s)
        A tuple of parameter names that will be varied.
    submodels: dict
        Dictionary of model labels -> model instances of all the submodels.
    \**kwargs :
        All other keyword arguments are passed to
        :py:class:`BaseModel <pycbc.inference.models.base.BaseModel>`.
    """
    name = 'hierarchical'

    def __init__(self, variable_params, submodels, **kwargs):
        # sub models is assumed to be a dict of model labels -> model instances
        self.submodels = submodels
        # initialize standard attributes
        super().__init__(variable_params, **kwargs)
        # store a map of model labels -> parameters for quick look up later
        self.param_map = map_params(self.hvariable_params)
        # add any parameters created by waveform transforms
        if self.waveform_transforms is not None:
            derived_params = set()
            derived_params.update(*[t.outputs
                                    for t in self.waveform_transforms])
            # convert to hierarchical params
            derived_params = map_params(hpiter(derived_params,
                                               list(self.submodels.keys())))
            for lbl, pset in derived_params.items():
                self.param_map[lbl].update(pset)
        # make sure the static parameters of all submodels are set correctly
        self.static_param_map = map_params(self.hstatic_params.keys())
        # also create a map of model label -> extra stats created by each model
        # stats are prepended with the model label. We'll include the
        # loglikelihood returned by each submodel in the extra stats.
        self.extra_stats_map = {}
        self.__extra_stats = []
        for lbl, model in self.submodels.items():
            model.static_params = {p.subname: self.static_params[p.fullname]
                                   for p in self.static_param_map[lbl]}
            self.extra_stats_map.update(map_params([
                HierarchicalParam.from_subname(lbl, p)
                for p in model._extra_stats+['loglikelihood']]))
            self.__extra_stats += self.extra_stats_map[lbl]
            # also make sure the model's sampling transforms and waveform
            # transforms are not set, as these are handled by the hierarchical
            # model, except for `joint_primary_marginalized` model, because
            # this specific model needs to allow its submodels to handle
            # transform with prefix on the submodel's level
            if self.name != "joint_primary_marginalized":
                if model.sampling_transforms is not None:
                    raise ValueError("Model {} has sampling transforms "
                                     "set; in a hierarchical analysis, "
                                     "these are handled by the "
                                     "hierarchical model".format(lbl))
                if model.waveform_transforms is not None:
                    raise ValueError("Model {} has waveform transforms "
                                     "set; in a hierarchical analysis, "
                                     "these are handled by the "
                                     "hierarchical model".format(lbl))

    @property
    def hvariable_params(self):
        """The variable params as a tuple of :py:class:`HierarchicalParam`
        instances.
        """
        return self._variable_params

    @property
    def variable_params(self):
        # converts variable params back to a set of strings before returning
        return tuple(p.fullname for p in self._variable_params)

    @variable_params.setter
    def variable_params(self, variable_params):
        # overrides BaseModel's variable params to store the variable params
        # as HierarchicalParam instances
        if isinstance(variable_params, str):
            variable_params = [variable_params]
        self._variable_params = tuple(HierarchicalParam(p, self.submodels)
                                      for p in variable_params)

    @property
    def hstatic_params(self):
        """The static params with :py:class:`HierarchicalParam` instances used
        as dictionary keys.
        """
        return self._static_params

    @property
    def static_params(self):
        # converts the static param keys back to strings
        return {p.fullname: val for p, val in self._static_params.items()}

    @static_params.setter
    def static_params(self, static_params):
        if static_params is None:
            static_params = {}
        self._static_params = {HierarchicalParam(p, self.submodels): val
                               for p, val in static_params.items()}

    @property
    def _extra_stats(self):
        return [p.fullname for p in self.__extra_stats]

    @property
    def _hextra_stats(self):
        """The extra stats as :py:class:`HierarchicalParam` instances."""
        return self.__extra_stats

    def _loglikelihood(self):
        # takes the sum of the constitutent models' loglikelihoods
        logl = 0.
        for lbl, model in self.submodels.items():
            # update the model with the current params. This is done here
            # instead of in `update` because waveform transforms are not
            # applied until the loglikelihood function is called
            model.update(**{p.subname: self.current_params[p.fullname]
                            for p in self.param_map[lbl]})
            # now get the loglikelihood from the model
            sublogl = model.loglikelihood
            # store the extra stats
            mstats = model.current_stats
            for stat in self.extra_stats_map[lbl]:
                setattr(self._current_stats, stat, mstats[stat.subname])
            # add to the total loglikelihood
            logl += sublogl
        return logl

    def write_metadata(self, fp, group=None):
        """Adds data to the metadata that's written.

        Parameters
        ----------
        fp : pycbc.inference.io.BaseInferenceFile instance
            The inference file to write to.
        group : str, optional
            If provided, the metadata will be written to the attrs specified
            by group, i.e., to ``fp[group].attrs``. Otherwise, metadata is
            written to the top-level attrs (``fp.attrs``).

        """
        # write information about self
        super().write_metadata(fp, group=group)
        # write information about each submodel into a different group for
        # each one
        if group is None or group == '/':
            prefix = ''
        else:
            prefix = group+'/'
        for lbl, model in self.submodels.items():
            model.write_metadata(fp, group=prefix+lbl)

        # if all submodels support it, write a combined lognl parameter
        try:
            sampattrs = fp.getattrs(group=fp.samples_group)
            lognl = [self.submodels[k].lognl for k in self.submodels]
            sampattrs['{}lognl'.format(prefix)] = sum(lognl)
        except AttributeError:
            pass

    @classmethod
    def from_config(cls, cp, **kwargs):
        r"""Initializes an instance of this class from the given config file.

        Sub-models are initialized before initializing this class. The model
        section must have a ``submodels`` argument that lists the names of all
        the submodels to generate as a space-separated list. Each sub-model
        should have its own ``[{label}__model]`` section that sets up the
        model for that sub-model. For example:

        .. code-block:: ini

            [model]
            name = hierarchical
            submodels = event1 event2

            [event1__model]
            <event1 model options>

            [event2__model]
            <event2 model options>

        Similarly, all other sections that are specific to a model should start
        with the model's label. All sections starting with a model's label will
        be passed to that model's ``from_config`` method with the label removed
        from the section name. For example, if a sub-model requires a data
        section to be specified, it should be titled ``[{label}__data]``. Upon
        initialization, the ``{label}__`` will be stripped from the section
        header and passed to the model.

        No model labels should preceed the ``variable_params``,
        ``static_params``, ``waveform_transforms``, or ``sampling_transforms``
        sections.  Instead, the parameters specified in these sections should
        follow the naming conventions described in :py:class:`HierachicalParam`
        to determine which sub-model(s) they belong to. (Sampling parameters
        can follow any naming convention, as they are only handled by the
        hierarchical model.) This is because the hierarchical model handles
        all transforms, communication with the sampler, file IO, and prior
        calculation. Only sub-model's loglikelihood functions are called.

        Metadata for each sub-model is written to the output hdf file under
        groups given by the sub-model label. For example, if we have two
        submodels labelled ``event1`` and ``event2``, there will be groups
        with the same names in the top level of the output that contain that
        model's subdata. For instance, if event1 used the ``gaussian_noise``
        model, the GW data and PSDs will be found in ``event1/data`` and the
        low frequency cutoff used for that model will be in the ``attrs`` of
        the ``event1`` group.

        Parameters
        ----------
        cp : WorkflowConfigParser
            Config file parser to read.
        \**kwargs :
            All additional keyword arguments are passed to the class. Any
            provided keyword will override what is in the config file.
        """
        # we need the read from config function from the init; to prevent
        # circular imports, we import it here
        from pycbc.inference.models import read_from_config
        # get the submodels
        submodel_lbls = shlex.split(cp.get('model', 'submodels'))
        # sort parameters by model
        vparam_map = map_params(hpiter(cp.options('variable_params'),
                                       submodel_lbls))
        sparam_map = map_params(hpiter(cp.options('static_params'),
                                       submodel_lbls))

        # we'll need any waveform transforms for the initializing sub-models,
        # as the underlying models will receive the output of those transforms
        if any(cp.get_subsections('waveform_transforms')):
            waveform_transforms = transforms.read_transforms_from_config(
                cp, 'waveform_transforms')
            wfoutputs = set.union(*[t.outputs
                                    for t in waveform_transforms])
            wfparam_map = map_params(hpiter(wfoutputs, submodel_lbls))
        else:
            wfparam_map = {lbl: [] for lbl in submodel_lbls}
        # initialize the models
        submodels = {}
        logging.info("Loading submodels")
        for lbl in submodel_lbls:
            logging.info("============= %s =============", lbl)
            # create a config parser to pass to the model
            subcp = WorkflowConfigParser()
            # copy sections over that start with the model label (this should
            # include the [model] section for that model)
            copy_sections = [
                HierarchicalParam(sec, submodel_lbls)
                for sec in cp.sections() if lbl in
                sec.split('-')[0].split(HierarchicalParam.delim, 1)[0]]
            for sec in copy_sections:
                # check that the user isn't trying to set variable or static
                # params for the model (we won't worry about waveform or
                # sampling transforms here, since that is checked for in the
                # __init__)
                if sec.subname in ['variable_params', 'static_params']:
                    raise ValueError("Section {} found in the config file; "
                                     "[variable_params] and [static_params] "
                                     "sections should not include model "
                                     "labels. To specify parameters unique to "
                                     "one or more sub-models, prepend the "
                                     "individual parameter names with the "
                                     "model label. See HierarchicalParam for "
                                     "details.".format(sec))
                subcp.add_section(sec.subname)
                for opt, val in cp.items(sec):
                    subcp.set(sec.subname, opt, val)
            # set the static params
            subcp.add_section('static_params')
            for param in sparam_map[lbl]:
                subcp.set('static_params', param.subname,
                          cp.get('static_params', param.fullname))
            # set the variable params: for now we'll just set all the
            # variable params as static params
            # so that the model doesn't raise an error looking for
            # prior sections. We'll then manually set the variable
            # params after the model is initialized

            subcp.add_section('variable_params')
            for param in vparam_map[lbl]:
                subcp.set('static_params', param.subname, 'REPLACE')
            # add the outputs from the waveform transforms
            for param in wfparam_map[lbl]:
                subcp.set('static_params', param.subname, 'REPLACE')

            # initialize
            submodel = read_from_config(subcp)
            # move the static params back to variable
            for p in vparam_map[lbl]:
                submodel.static_params.pop(p.subname)
            submodel.variable_params = tuple(p.subname
                                             for p in vparam_map[lbl])
            # remove the waveform transform parameters
            for p in wfparam_map[lbl]:
                submodel.static_params.pop(p.subname)
            # store
            submodels[lbl] = submodel
            logging.info("")
        # now load the model
        logging.info("Loading hierarchical model")
        return super().from_config(cp, submodels=submodels)


class HierarchicalParam(str):
    """Sub-class of str for hierarchical parameter names.

    This adds attributes that keep track of the model label(s) the parameter
    is associated with, along with the name that is passed to the models.

    The following conventions are used for parsing parameter names:

      * Model labels and parameter names are separated by the ``delim`` class
        attribute, which by default is ``__``, e.g., ``event1__mass``.
      * Multiple model labels can be provided by separating the model labels
        with the ``model_delim`` class attribute, which by default is ``_``,
        e.g., ``event1_event2__mass``. Note that this means that individual
        model labels cannot contain ``_``, else they'll be parsed as separate
        models.
      * Parameters that have no model labels prepended to them (i.e., there
        is no ``__`` in the name) are common to all models.

    These parsing rules are applied by the :py:meth:`HierarchicalParam.parse`
    method.

    Parameters
    ----------
    fullname : str
        Name of the hierarchical parameter. Should have format
        ``{model1}[_{model2}[_{...}]]__{param}``.
    possible_models : set of str
        The possible sub-models a parameter can belong to. Should a set of
        model labels.

    Attributes
    ----------
    fullname : str
        The full name of the parameter, including model labels. For example,
        ``e1_e2__foo``.
    models : set
        The model labels the parameter is associated with. For example,
        ``e1_e2__foo`` yields models ``e1, e2``.
    subname : str
        The name of the parameter without the model labels prepended to it.
        For example, ``e1_e2__foo`` yields ``foo``.
    """
    delim = '__'
    model_delim = '_'

    def __new__(cls, fullname, possible_models):
        fullname = str(fullname)
        obj = str.__new__(cls, fullname)
        obj.fullname = fullname
        models, subp = HierarchicalParam.parse(fullname, possible_models)
        obj.models = models
        obj.subname = subp
        return obj

    @classmethod
    def from_subname(cls, model_label, subname):
        """Creates a HierarchicalParam from the given subname and model label.
        """
        return cls(cls.delim.join([model_label, subname]), set([model_label]))

    @classmethod
    def parse(cls, fullname, possible_models):
        """Parses the full parameter name into the models the parameter is
        associated with and the parameter name that is passed to the models.

        Parameters
        ----------
        fullname : str
            The full name of the parameter, which includes both the model
            label(s) and the parameter name.
        possible_models : set
            Set of model labels the parameter can be associated with.

        Returns
        -------
        models : list
            List of the model labels the parameter is associated with.
        subp : str
            Parameter name that is passed to the models. This is the parameter
            name with the model label(s) stripped from it.
        """
        # make sure possible models is a set
        possible_models = set(possible_models)
        p = fullname.split(cls.delim, 1)
        if len(p) == 1:
            # is a global fullname, associate with all
            subp = fullname
            models = possible_models.copy()
        else:
            models, subp = p
            # convert into set of model label(s)
            models = set(models.split(cls.model_delim))
            # make sure the given labels are in the list of possible models
            unknown = models - possible_models
            if any(unknown):
                raise ValueError('unrecognized model label(s) {} present in '
                                 'parameter {}'.format(', '.join(unknown),
                                                       fullname))
        return models, subp


def hpiter(params, possible_models):
    """Turns a list of parameter strings into a list of HierarchicalParams.

    Parameters
    ----------
    params : list of str
        List of parameter names.
    possible_models : set
        Set of model labels the parameters can be associated with.

    Returns
    -------
    iterator :
        Iterator of :py:class:`HierarchicalParam` instances.
    """
    return map(lambda x: HierarchicalParam(x, possible_models), params)


def map_params(params):
    """Creates a map of models -> parameters.

    Parameters
    ----------
    params : list of HierarchicalParam instances
        The list of hierarchical parameter names to parse.

    Returns
    -------
    dict :
        Dictionary of model labels -> associated parameters.
    """
    param_map = {}
    for p in params:
        for lbl in p.models:
            try:
                param_map[lbl].update([p])
            except KeyError:
                param_map[lbl] = set([p])
    return param_map


class MultiSignalModel(HierarchicalModel):
    """ Model for multiple signals which share data

    Sub models are treated as if the signals overlap in data. This requires
    constituent models to implement a specific method to handle this case.
    All models must be of the same type or the specific model is responsible
    for implement cross-compatibility with another model. Each model h_i is
    responsible for calculating its own loglikelihood ratio for itself, and
    must also implement a method to calculate crossterms of the form
    <h_i | h_j> which arise from the full calculation of <d - h|d - h>.
    This model inherits from the HierarchicalModel so the syntax for
    configuration files is the same. The primary model is used to determine
    the noise terms <d | d>, which by default will be the first model used.
    """
    name = 'multi_signal'

    def __init__(self, variable_params, submodels, **kwargs):
        super().__init__(variable_params, submodels, **kwargs)

        # Check what models each model supports
        support = {}
        ctypes = set()  # The set of models we need to completely support
        for lbl in self.submodels:
            model = self.submodels[lbl]

            ctypes.add(type(model))
            if hasattr(model, 'multi_signal_support'):
                support[lbl] = set(model.multi_signal_support)

        # pick the primary model if it supports the set of constituent models
        for lbl in support:
            if ctypes <= support[lbl]:
                self.primary_model = lbl
                logging.info('MultiSignalModel: PrimaryModel == %s', lbl)
                break
        else:
            # Oh, no, we don't support this combo!
            raise RuntimeError("It looks like the combination of models, {},"
                               "for the MultiSignal model isn't supported by"
                               "any of the constituent models.".format(ctypes))

        self.other_models = self.submodels.copy()
        self.other_models.pop(self.primary_model)
        self.other_models = list(self.other_models.values())

    def write_metadata(self, fp, group=None):
        """Adds metadata to the output files

        Parameters
        ----------
        fp : pycbc.inference.io.BaseInferenceFile instance
            The inference file to write to.
        group : str, optional
            If provided, the metadata will be written to the attrs specified
            by group, i.e., to ``fp[group].attrs``. Otherwise, metadata is
            written to the top-level attrs (``fp.attrs``).
        """
        super().write_metadata(fp, group=group)
        sampattrs = fp.getattrs(group=fp.samples_group)
        # if a group is specified, prepend the lognl names with it
        if group is None or group == '/':
            prefix = ''
        else:
            prefix = group.replace('/', '__')
            if not prefix.endswith('__'):
                prefix += '__'
        try:
            model = self.submodels[self.primary_model]
            sampattrs['{}lognl'.format(prefix)] = model.lognl
        except AttributeError:
            pass

    def _loglikelihood(self):
        for lbl, model in self.submodels.items():
            # Update the parameters of each
            model.update(**{p.subname: self.current_params[p.fullname]
                            for p in self.param_map[lbl]})

        # Calculate the combined loglikelihood
        p = self.primary_model
        logl = self.submodels[p].multi_loglikelihood(self.other_models)

        # store any extra stats from the submodels
        for lbl, model in self.submodels.items():
            mstats = model.current_stats
            for stat in self.extra_stats_map[lbl]:
                setattr(self._current_stats, stat, mstats[stat.subname])
        return logl


class JointPrimaryMarginalizedModel(HierarchicalModel):
    """ Hierarchical heterodyne likelihood for coherent multiband
    parameter estimation which combines data from space-borne and
    ground-based GW detectors coherently. Currently, this only
    supports LISA as the space-borne GW detector.

    Sub models are treated as if the same GW source (such as a GW
    from stellar-mass BBH) is observed in different frequency bands by
    space-borne and ground-based GW detectors, then transform all
    the parameters into the same frame in the sub model level, use
    `HierarchicalModel` to get the joint likelihood, and marginalize
    over all the extrinsic parameters supported by `RelativeTimeDom`
    or its variants. Note that LISA submodel only supports the `Relative`
    for now, for ground-based detectors, please use `RelativeTimeDom`
    or its variants.

    Although this likelihood model is used for multiband parameter
    estimation, users can still use it for other purposes, such as
    GW + EM parameter estimation, in this case, please use `RelativeTimeDom`
    or its variants for the GW data, for the likelihood of EM data,
    there is no restrictions.
    """
    name = 'joint_primary_marginalized'

    def __init__(self, variable_params, submodels, **kwargs):
        super().__init__(variable_params, submodels, **kwargs)

        # assume the ground-based submodel as the primary model
        self.primary_model = self.submodels[kwargs['primary_lbl'][0]]
        self.primary_lbl = kwargs['primary_lbl'][0]
        self.other_models = self.submodels.copy()
        self.other_models.pop(kwargs['primary_lbl'][0])
        self.other_models = list(self.other_models.values())

    def write_metadata(self, fp, group=None):
        """Adds metadata to the output files

        Parameters
        ----------
        fp : pycbc.inference.io.BaseInferenceFile instance
            The inference file to write to.
        group : str, optional
            If provided, the metadata will be written to the attrs specified
            by group, i.e., to ``fp[group].attrs``. Otherwise, metadata is
            written to the top-level attrs (``fp.attrs``).
        """
        super().write_metadata(fp, group=group)
        sampattrs = fp.getattrs(group=fp.samples_group)
        # if a group is specified, prepend the lognl names with it
        if group is None or group == '/':
            prefix = ''
        else:
            prefix = group.replace('/', '__')
            if not prefix.endswith('__'):
                prefix += '__'
        try:
            for lbl, model in self.submodels.items():
                sampattrs['{}lognl'.format(prefix + '%s__' % lbl)
                          ] = model.lognl
        except AttributeError:
            pass

    def total_loglr(self):
        r"""Computes the total log likelihood ratio,

        .. math::

            \log \mathcal{L}(\Theta) = \sum_i
                \left<h_i(\Theta)|d_i\right> -
                \frac{1}{2}\left<h_i(\Theta)|h_i(\Theta)\right>,

        at the current parameter values :math:`\Theta`.

        Returns
        -------
        float
            The value of the log likelihood ratio.
        """
        # calculate <d-h|d-h> = <h|h> - 2<h|d> + <d|d> up to a constant

        # note that for SOBHB signals, ground-based detectors dominant SNR
        # and accuracy of (tc, ra, dec)
        self.primary_model.return_sh_hh = True
        sh_primary, hh_primary = self.primary_model.loglr
        self.primary_model.return_sh_hh = False

        margin_names_vector = list(
            self.primary_model.marginalize_vector_params.keys())
        if 'logw_partial' in margin_names_vector:
            margin_names_vector.remove('logw_partial')

        margin_params = {}
        nums = 1
        for key, value in self.primary_model.current_params.items():
            # add marginalize_vector_params
            if key in margin_names_vector:
                margin_params[key] = value
                if isinstance(value, numpy.ndarray):
                    nums = len(value)
        # add distance if it has been marginalized,
        # use numpy array for it is just let it has the same
        # shape as marginalize_vector_params, here we assume
        # self.primary_model.current_params['distance'] is a number
        if self.primary_model.distance_marginalization:
            margin_params['distance'] = numpy.full(
                nums, self.primary_model.current_params['distance'])

        # add likelihood contribution from space-borne detectors, we
        # calculate sh/hh for each marginalized parameter point
        sh_others = numpy.full(nums, 0 + 0.0j)
        hh_others = numpy.zeros(nums)

        # update parameters in other_models
        for _, other_model in enumerate(self.other_models):
            # not using self.primary_model.current_params, because others_model
            # may have its own static parameters
            current_params_other = other_model.current_params.copy()
            for i in range(nums):
                current_params_other.update(
                    {key: value[i] if isinstance(value, numpy.ndarray) else
                        value for key, value in margin_params.items()})
                other_model.update(**current_params_other)
                other_model.return_sh_hh = True
                sh_others[i], hh_others[i] = other_model.loglr
                other_model.return_sh_hh = False

        if nums == 1:
            sh_others = sh_others[0]
        sh_total = sh_primary + sh_others
        hh_total = hh_primary + hh_others

        # calculate marginalize_vector_weights
        self.primary_model.marginalize_vector_weights = \
            - numpy.log(self.primary_model.vsamples)
        loglr = self.primary_model.marginalize_loglr(sh_total, hh_total)
        return loglr

    def others_lognl(self):
        """Calculate the combined lognl from all others sub-models."""
        total_others_lognl = 0
        for model in self.other_models:
            total_others_lognl += model.lognl
        return total_others_lognl

    def update_all_models(self, **params):
        """This update method is also useful for loglr checking,
        the original update method in base module can't update
        parameters in submodels correctly in loglr checking."""
        for lbl, model in self.submodels.items():
            if self.param_map != {}:
                p = {params.subname: self.current_params[params.fullname]
                     for params in self.param_map[lbl]}
            else:
                # dummy sampler doesn't have real variables,
                # which means self.param_map is {}
                p = {}
            p.update(params)
            model.update(**p)

    def _loglikelihood(self):
        self.update_all_models()

        # calculate the combined loglikelihood
        logl = self.total_loglr() + self.primary_model.lognl + \
            self.others_lognl()

        # store any extra stats from the submodels
        for lbl, model in self.submodels.items():
            mstats = model.current_stats
            for stat in self.extra_stats_map[lbl]:
                setattr(self._current_stats, stat, mstats[stat.subname])
        return logl

    @classmethod
    def from_config(cls, cp, **kwargs):
        r"""Initializes an instance of this class from the given config file.
        For more details, see `from_config` in `HierarchicalModel`.

        Parameters
        ----------
        cp : WorkflowConfigParser
            Config file parser to read.
        \**kwargs :
            All additional keyword arguments are passed to the class. Any
            provided keyword will override what is in the config file.
        """
        # we need the read from config function from the init; to prevent
        # circular imports, we import it here
        from pycbc.inference.models import read_from_config
        # get the submodels
        kwargs['primary_lbl'] = shlex.split(cp.get('model', 'primary_model'))
        kwargs['others_lbls'] = shlex.split(cp.get('model', 'other_models'))
        submodel_lbls = kwargs['primary_lbl'] + kwargs['others_lbls']
        # sort parameters by model
        vparam_map = map_params(hpiter(cp.options('variable_params'),
                                       submodel_lbls))
        sparam_map = map_params(hpiter(cp.options('static_params'),
                                       submodel_lbls))

        # we'll need any waveform transforms for the initializing sub-models,
        # as the underlying models will receive the output of those transforms

        # if `waveform_transforms` section doesn't have the prefix of
        # sub-model's name, then add this `waveform_transforms` section
        # into top level, if not, add it into sub-models' config
        if any(cp.get_subsections('waveform_transforms')):
            waveform_transforms = transforms.read_transforms_from_config(
                cp, 'waveform_transforms')
            wfoutputs = set.union(*[t.outputs
                                    for t in waveform_transforms])
            wfparam_map = map_params(hpiter(wfoutputs, submodel_lbls))
        else:
            wfparam_map = {lbl: [] for lbl in submodel_lbls}
        # initialize the models
        submodels = {}
        logging.info("Loading submodels")
        for lbl in submodel_lbls:
            logging.info("============= %s =============", lbl)
            # create a config parser to pass to the model
            subcp = WorkflowConfigParser()
            # copy sections over that start with the model label (this should
            # include the [model] section for that model)
            copy_sections = [
                HierarchicalParam(sec, submodel_lbls)
                for sec in cp.sections() if lbl in
                sec.split('-')[0].split(HierarchicalParam.delim, 1)[0]]
            for sec in copy_sections:
                # check that the user isn't trying to set variable or static
                # params for the model (we won't worry about waveform or
                # sampling transforms here, since that is checked for in the
                # __init__)
                if sec.subname in ['variable_params', 'static_params']:
                    raise ValueError("Section {} found in the config file; "
                                     "[variable_params] and [static_params] "
                                     "sections should not include model "
                                     "labels. To specify parameters unique to "
                                     "one or more sub-models, prepend the "
                                     "individual parameter names with the "
                                     "model label. See HierarchicalParam for "
                                     "details.".format(sec))
                subcp.add_section(sec.subname)
                for opt, val in cp.items(sec):
                    subcp.set(sec.subname, opt, val)
            # set the static params
            subcp.add_section('static_params')
            for param in sparam_map[lbl]:
                subcp.set('static_params', param.subname,
                          cp.get('static_params', param.fullname))

            # set the variable params: different from the standard
            # hierarchical model, in this multiband model, all sub-models
            # has the same variable parameters, so we don't need to worry
            # about the unique variable issue. Besides, the primary model
            # needs to do marginalization, so we must set variable_params
            # and prior section before initializing it.

            subcp.add_section('variable_params')
            for param in vparam_map[lbl]:
                if lbl in kwargs['primary_lbl']:
                    subcp.set('variable_params', param.subname,
                              cp.get('variable_params', param.fullname))
                else:
                    subcp.set('static_params', param.subname, 'REPLACE')

            for section in cp.sections():
                # the primary model needs prior of marginlized parameters
                if 'prior-' in section and lbl in kwargs['primary_lbl']:
                    prior_section = '%s' % section
                    subcp[prior_section] = cp[prior_section]
                # if `waveform_transforms` has a prefix,
                # add it into sub-models' config
                elif '%s_waveform_transforms' % lbl in section:
                    transforms_section = '%s' % section
                    subcp[transforms_section] = cp[transforms_section]
                else:
                    pass

            # similar to the standard hierarchical model,
            # add the outputs from the waveform transforms if sub-model
            # doesn't need marginalization
            if lbl not in kwargs['primary_lbl']:
                for param in wfparam_map[lbl]:
                    subcp.set('static_params', param.subname, 'REPLACE')

            # save the vitual config file to disk for later check
            with open('%s.ini' % lbl, 'w', encoding='utf-8') as file:
                subcp.write(file)

            # initialize
            submodel = read_from_config(subcp)

            if lbl not in kwargs['primary_lbl']:
                # similar to the standard hierarchical model,
                # move the static params back to variable if sub-model
                # doesn't need marginalization
                for p in vparam_map[lbl]:
                    submodel.static_params.pop(p.subname)
                submodel.variable_params = tuple(p.subname
                                                 for p in vparam_map[lbl])
                # similar to the standard hierarchical model,
                # remove the waveform transform parameters if sub-model
                # doesn't need marginalization
                for p in wfparam_map[lbl]:
                    submodel.static_params.pop(p.subname)
            submodels[lbl] = submodel
            logging.info("")

        # remove all marginalized parameters from the top-level model's
        # `variable_params` and `prior` sections
        # here we ignore `coa_phase`, because if it's been marginalized,
        # it will not be listed in `variable_params` and `prior` sections
        primary_model = submodels[kwargs['primary_lbl'][0]]
        marginalized_params = primary_model.marginalize_vector_params.copy()
        if 'logw_partial' in marginalized_params:
            marginalized_params.pop('logw_partial')
            marginalized_params = list(marginalized_params.keys())
        else:
            marginalized_params = []
        # this may also include 'f_ref', 'f_lower', 'approximant',
        # but doesn't matter
        marginalized_params += list(primary_model.static_params.keys())
        for p in primary_model.static_params.keys():
            p_full = '%s__%s' % (kwargs['primary_lbl'][0], p)
            if p_full not in cp['static_params']:
                cp['static_params'][p_full] = "%s" % \
                    primary_model.static_params[p]

        for section in cp.sections():
            if 'prior-' in section:
                p = section.split('-')[-1]
                if p in marginalized_params:
                    cp['variable_params'].pop(p)
                    cp.pop(section)

        # now load the model
        logging.info("Loading joint_primary_marginalized model")
        return super(HierarchicalModel, cls).from_config(
                cp, submodels=submodels, **kwargs)

    def reconstruct(self, rec=None, seed=None):
        """ Reconstruct marginalized parameters by using the primary
        model's reconstruct method, total_loglr, and others_lognl.
        """
        if seed:
            numpy.random.seed(seed)

        if rec is None:
            rec = {}

        def get_loglr():
            self.update_all_models(**rec)
            return self.total_loglr()

        rec = self.primary_model.reconstruct(
                rec=rec, seed=seed, set_loglr=get_loglr)
        # the primary model's reconstruct doesn't know lognl in other models
        rec['loglikelihood'] += self.others_lognl()
        return rec

.. _hierachical_model:

-----------------------------------
Using the hierarchical model
-----------------------------------

.. |HierarchicalModel| replace:: :py:class:`HierarchicalModel <pycbc.inference.models.hierarchical.HierarchicalModel>`


The |HierarchicalModel| is for performing a joint Bayesian analysis over two or
more events that share one or more common parameters. For example, you may wish
to analyze two different binary neutron star mergers assuming the neutron stars
share the same equation of state, or you want to analyze two signals that
arrived at different times under the assumption that they are lensed copies of
the same event.

The events are assumed to not share the same data, and may not even be
observations from the same type of detector. What type of data is read for each
event is determined by the model used for that event. Upon initialization,
the hierarchical model passes the relevant parameters and sections in the
config file to each sub-model's ``from_config`` method. During the analysis, the
sub-models' loglikelihood functions are called and summed over. In that regard,
the hierarchical model treats the sub-models as black boxes.

To set up the hierarchical model, you provide a list of sub-model labels in
the ``[model]`` section of your config file that represent each event to
analyze. For each sub-model label provided, you provide a ``[{label}__model]``
section that in turn specifies how to initialize that event's model.

In the ``[variable_params]`` and ``[static_params]`` sections you specify
which parameters belong to which event by prepending the parameter name with
the event's label, followed by a ``__``; i.e., ``{label}__{param}``. To specify
that a parameter is common to a sub-set of events, you prepend each event's
label (separated by a ``_``). Parameters that have no event labels prepended to
them are treated as common to all events.

For example, say we have three events ``a``, ``b``, and ``c``, who's models
take two parameters, ``foo`` and ``bar``. Say we want ``foo`` to be common
between ``a`` and ``b``, but unique to ``c``, and ``bar`` to common all events.
Our analysis would therefore have three parameters, ``a_b__foo``, ``c__foo``,
and ``bar``.

Additional sections that are required by each model should have the model label
prepended to them in the same manner. In the above example, if the models for
events ``a``, ``b``, and ``c`` require a ``data`` section, then the config file
should have sections ``[a__data]``, ``[b__data]``, and ``[c__data]``.

When sub-models are initialized, the event labels are stripped from the
parameter names (and section headers), then passed to the sub-model. Sub-models
do not need any special features to be run in a hierarchical analysis; any
model that inherits from :py:class:`BaseModel
<pycbc.inference.models.base.BaseModel>` can be used as a sub-model.

^^^^^^^^^^^^^^^
Lensing example
^^^^^^^^^^^^^^^

To illustrate how to use the |HierarchicalModel| here we perform a simple
lensing analysis of two binary black hole simulations. The two simulations
are injected into fake noise ~one week apart, with slightly different apparent
sky locations (the difference may not be physically possible; this is only
for illustration purposes). We will analyze the two events allowing them to
have different sky locations and coalescence times, but sharing the same
masses.

First, we need to create the simulations using ``pycbc_create_injections``.
To create the two injections we'll use the configuration files:

.. literalinclude:: ../../../examples/inference/hierarchical/event1_inj.ini
    :language: ini

:download:`Download <../../../examples/inference/hierarchical/event1_inj.ini>`

.. literalinclude:: ../../../examples/inference/hierarchical/event2_inj.ini
    :language: ini

:download:`Download <../../../examples/inference/hierarchical/event2_inj.ini>`

Create the injection hdf files by running:

.. literalinclude:: ../../../examples/inference/hierarchical/make_injections.sh
    :language: bash

:download:`Download <../../../examples/inference/hierarchical/make_injections.sh>`

Now we'll setup the configuration files to run the hierarchical analysis on
these two events. First we'll set the ``[model]`` section to use the
|HierarchicalModel|, and tell it to analyze two events called ``event1`` and
``event2``:

.. literalinclude:: ../../../examples/inference/hierarchical/model.ini
    :language: ini

:download:`Download <../../../examples/inference/hierarchical/model.ini>`

We now need to specify the model to use for each event. For this analysis
we'll use the :py:class:`Relative <pycbc.inference.models.relbin.Relative>`
model for both events:

.. literalinclude::  ../../../examples/inference/hierarchical/model-event1_relbin.ini
    :language: ini

:download:`Download <../../../examples/inference/hierarchical/model-event1_relbin.ini>`

.. literalinclude::  ../../../examples/inference/hierarchical/model-event2_relbin.ini
    :language: ini

:download:`Download <../../../examples/inference/hierarchical/model-event2_relbin.ini>`

We also need to provide data sections for each event:

.. literalinclude:: ../../../examples/inference/hierarchical/data.ini
    :language: ini

:download:`Download <../../../examples/inference/hierarchical/data.ini>`

Now the prior:

.. literalinclude:: ../../../examples/inference/hierarchical/prior.ini
    :language: ini

:download:`Download <../../../examples/inference/hierarchical/prior.ini>`

And finally the sampler configuration:

.. literalinclude:: ../../../examples/inference/hierarchical/dynesty.ini
    :language: ini

:download:`Download <../../../examples/inference/hierarchical/dynesty.ini>`

Now run:

.. literalinclude:: ../../../examples/inference/hierarchical/run.sh
   :language: bash

:download:`Download <../../../examples/inference/hierarchical/run.sh>`

When it is done, you will have a file called ``hierarchical.hdf`` which
contains the results.

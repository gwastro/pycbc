===========
Sampler API
===========

The :py:mod:`pycbc.inference.sampler` module is the interface between
``pycbc_inference`` and the sampling engines, such as ``emcee``. Below, we
provide an overview of the general structure of the sampler classes, how it
interacts with ``pycbc_inference``, and how to add support for new samplers. We
also provide inheritance diagrams for all of the currently supported samplers.

---------------------
Overview & Guidelines
---------------------

The following guidelines apply to the sampler classes:

1. All sampler classes must inherit from the 
   :py:class:`BaseSampler <pycbc.inference.sampler.BaseSampler>` class. This is
   an :py:mod:`abstract base class <abc>` that defines methods that all
   samplers must implement, as they will be used by ``pycbc_inference``. (See
   `this tutorial <https://www.python-course.eu/python3_abstract_classes.php>`_
   for a primer on abstract base classes.)
2. Duplicate code should be avoided. If multiple samplers have common methods,
   those methods should be added to one or more support classes that the
   samplers can inherit from, in addition to
   :py:class:`BaseSampler <pycbc.inference.sampler.base.BaseSampler>`. For
   example, all MCMC samplers need to be able to compute an autocorrelation
   length. That functionality is provided for single-temperature MCMCs in
   :py:class:`MCMCAutocorrSupport <pycbc.inference.sampler.MCMCAutocoorSupport>`
   class. These support classes may themselves be abstract base classes which
   add more required methods to samplers that inherit from them.
3. Inheritance is kept to one level. For example, if we have sampler class
   ``Foo(Bar, BaseSampler)``, both ``Bar`` and ``BaseSampler`` do not inherit
   from any parent classes, only ``object``.
4. To avoid confusion, only inherited abstract methods should be overridden.
6. All sampler classes need a corresponding class in the
   :py:mod:`pycbc.inference.io` module for handling reading and writing. See
   :doc:`io` for more details on IO classes.

As mentioned above, the
:py:class:`BaseSampler <pycbc.inference.sampler.base.BaseSampler>` class is an
abstract base class. It defines a collection of abstract methods and properties
that all samplers must override in order to function properly. These are (click
on the names to see their documentation):

 * :py:meth:`from_config <pycbc.inference.sampler.base.BaseSampler.from_config>`
 * :py:attr:`io <pycbc.inference.sampler.base.BaseSampler.io>`
 * :py:meth:`set_initial_conditions <pycbc.inference.sampler.base.BaseSampler.set_initial_conditions>`
 * :py:meth:`run <pycbc.inference.sampler.base.BaseSampler.run>`
 * :py:meth:`checkpoint <pycbc.inference.sampler.base.BaseSampler.checkpoint>`
 * :py:attr:`samples <pycbc.inference.sampler.base.BaseSampler.samples>`
 * :py:attr:`model_stats <pycbc.inference.sampler.base.BaseSampler.model_stats>`
 * :py:meth:`finalize <pycbc.inference.sampler.base.BaseSampler.finalize>`

---------------------
Detailed example
---------------------

Let's examine the
:py:class:`EmceeEnsembleSampler <pycbc.inference.sampler.emcee.EmceeEnsembleSampler>`
class to see how these guidelines apply in practice. Here is its inheritance
structure (click on the names of the classes to see their documentation):

.. inheritance-diagram:: pycbc.inference.sampler.emcee
    :parts: 2

In addition to :py:class:`BaseSampler <pycbc.inference.sampler.base.BaseSampler>`,
:py:class:`EmceeEnsembleSampler <pycbc.inference.sampler.emcee.EmceeEnsembleSampler>`
inherits from :py:class:`BaseMCMC <pycbc.inference.sampler.base_mcmc.BaseMCMC>` and
:py:class:`MCMCAutocorrSupport <pycbc.inference.sampler.base_mcmc.MCMCAutocorrSupport>`.
Inspecting :py:class:`BaseMCMC <pycbc.inference.sampler.base_mcmc.BaseMCMC>`,
we see that it implements several of the methods that
:py:class:`BaseSampler <pycbc.inference.sampler.base.BaseSampler>` requires: namely,
:py:meth:`set_initial_conditions <pycbc.inference.sampler.bsae_mcmc.BaseMCMC.set_initial_conditions>`,
:py:meth:`run <pycbc.inference.sampler.base_mcmc.BaseMCMC.run>`, and
:py:meth:`checkpoint <pycbc.inference.sampler.base_mcmc.BaseMCMC.checkpoint>`.
This is because the steps taken in these functions are common across MCMC
samplers. For example, in :py:meth:`run <pycbc.inference.sampler.base_mcmc.BaseMCMC.run>`,
the sampler is run for blocks of iterations (specified by
:py:attr:`checkpoint_interval <pycbc.inference.sampler.base_mcmc.checkpoint_interval>`)
until the convergence criteria has been met (which is determined by
:py:meth:`set_target <pycbc.inference.sampler.base_mcmc.set_target>`). This,
generally, is what all MCMC samplers do.

*How* an MCMC sampler is run for some number of iterations is unique to each
sampling engine. Thus, in :py:meth:`run <pycbc.inference.sampler.base_mcmc.BaseMCMC.run>`,
a call to :py:meth:`run_mcmc <pycbc.inference.sampler.base_mcmc.BaseMCMC.run_mcmc>`
is made. This is an abstract method: i.e., :py:class:`BaseMCMC <pycbc.inference.sampler.base_mcmc.BaseMCMC>`
is itself an abstract base class. Since
:py:class:`EmceeEnsembleSampler <pycbc.inference.sampler.emcee.EmceeEnsembleSampler>`
inheritance is ``EmceeEnsembleSampler(MCMCAutocorrSupport, BaseMCMC, BaseSampler)``,
:py:class:`BaseMCMC <pycbc.inference.sampler.base_mcmc.BaseMCMC>` fulfils
:py:class:`BaseSampler's <pycbc.inference.sampler.base.BaseSampler>` requirement
that a ``run`` method be implemented, but replaces it with the requirement that
a ``run_mcmc`` method be implemented. Thus, looking at its source code, we see
that :py:class:`EmceeEmsembleSampler <pycbc.inference.sampler.emcee.EmceeEnsembleSampler>`
implements a :py:meth:`run_mcmc <pycbc.inference.sampler.emcee.EmceeEnsembleSampler.run_mcmc>`
method. Likewise,

(note
that BaseSampler is furthest to the right: when multiple classes are specified


iterations are looped over, with results being dumped to 
dumped to file looking at the source code for
:py:meth:`run <pycbc.inference.sampler.base_mcmc.BaseMCMC.run>`, we see that
ther

---------------------
Inheritance diagrams
---------------------

Here are inheritance diagrams for the rest of the samplers currently supported:

* ``emcee_pt``:

.. inheritance-diagram:: pycbc.inference.sampler.emcee_pt
    :parts: 2

===========
Sampler API
===========

.. Following are useful substituions for classes and modules
.. BaseSampler:
.. |BaseSampler| replace:: :py:class:`BaseSampler <pycbc.inference.sampler.base.BaseSampler>`
.. BaseMCMC:
.. |BaseMCMC| replace:: :py:class:`BaseMCMC <pycbc.inference.sampler.base_mcmc.BaseMCMC>`
.. |BaseMCMC.set_initial_conditions| replace:: :py:meth:`set_initial_conditions <pycbc.inference.sampler.base_mcmc.BaseMCMC.set_initial_conditions>`
.. |BaseMCMC.run| replace:: :py:meth:`run <pycbc.inference.sampler.base_mcmc.BaseMCMC.run>`
.. |BaseMCMC.run_mcmc| replace:: :py:meth:`BaseMCMC.run_mcmc <pycbc.inference.sampler.base_mcmc.BaseMCMC.run_mcmc>`
.. |BaseMCMC.checkpoint| replace:: :py:meth:`checkpoint <pycbc.inference.sampler.base_mcmc.BaseMCMC.checkpoint>`
.. |BaseMCMC.checkpoint_interval| replace:: :py:attr:`checkpoint_interval <pycbc.inference.sampler.base_mcmc.BaseMCMC.checkpoint_interval>`
.. |BaseMCMC.set_target| replace:: :py:meth:`set_target <pycbc.inference.sampler.base_mcmc.BaseMCMC.set_target>`
.. |BaseMCMC.compute_acf| replace:: :py:meth:`compute_acf <pycbc.inference.sampler.base_mcmc.BaseMCMC.compute_acf>`
.. |BaseMCMC.compute_acl| replace:: :py:meth:`compute_acl <pycbc.inference.sampler.base_mcmc.BaseMCMC.compute_acl>`
.. MCMCAutocorrSupport
.. |MCMCAutocorrSupport| replace:: :py:class:`MCMCAutocorrSupport <pycbc.inference.sampler.base_mcmc.MCMCAutocorrSupport>`
.. EmceEnsembleSampler
.. |EmceeEnsembleSampler| replace:: :py:class:`EmceeEnsembleSampler <pycbc.inference.sampler.emcee.EmceeEnsembleSampler>`
.. |EmceeEnsembleSampler.run_mcmc| replace:: :py:meth:`run_mcmc <pycbc.inference.sampler.emcee.EmceeEnsembleSampler.run_mcmc>`
.. MultiTemperedAutocorrSupport
.. |MultiTemperedAutocorrSupport| replace:: :py:class:`MultiTemperedAutocorrSupport <pycbc.sampler.base_multitemper.MultiTemperedAutocorrSuppport>`
.. EmceePTSampler
.. |EmceePTSampler| replace:: :py:class:`EmceePTSampler <pycbc.sampler.emcee_pt.EmceePTSampler>`

The :py:mod:`pycbc.inference.sampler` module is the interface between
``pycbc_inference`` and the sampling engines, such as ``emcee``. Below, we
provide an overview of the general structure of the sampler classes, how it
interacts with ``pycbc_inference``, and how to add support for new samplers. We
also provide inheritance diagrams for all of the currently supported samplers.

.. _sampler_api-overview:

---------------------
Overview & Guidelines
---------------------

The following guidelines apply to the sampler classes:

1. All sampler classes must inherit from the |BaseSampler| class. This is
   an :py:mod:`abstract base class <abc>` that defines methods that all
   samplers must implement, as they will be used by ``pycbc_inference``. (See
   `this tutorial <https://www.python-course.eu/python3_abstract_classes.php>`_
   for a primer on abstract base classes.)
2. All sampler classes must have a ``name`` attribute that is unique across
   all sampler classes in :py:mod:`pycbc.inference.sampler`. This name is used
   to reference the sampler throughout the code, and is how the user specifies
   which sampler to use in their config file. For example,
   |EmceeEnsembleSampler|'s name is ``'emcee'``.
3. Duplicate code should be avoided. If multiple samplers have common methods,
   those methods should be added to one or more support classes that the
   samplers can inherit from, in addition to |BaseSampler|. For example, all
   MCMC samplers need to be able to compute an autocorrelation length. That
   functionality is provided for single-temperature MCMCs in the
   |MCMCAutocorrSupport| class. These support classes may themselves be
   abstract base classes which add more required methods to samplers that
   inherit from them.
4. Inheritance is kept to one level. For example, if we have sampler class
   ``Foo(Bar, BaseSampler)``, both ``Bar`` and ``BaseSampler`` do not inherit
   from any parent classes, only ``object``.
5. To avoid confusion, only inherited abstract methods should be overridden.
6. All sampler classes need a corresponding class in the
   :py:mod:`pycbc.inference.io` module for handling reading and writing. See
   :doc:`io` for more details on IO classes.

As mentioned above, the |BaseSampler| class is an abstract base class. It
defines a collection of abstract methods and properties that all samplers must
override in order to function properly. These are (click on the names to see
their documentation):

 * :py:meth:`from_config <pycbc.inference.sampler.base.BaseSampler.from_config>`
 * :py:attr:`io <pycbc.inference.sampler.base.BaseSampler.io>`
 * :py:meth:`set_initial_conditions <pycbc.inference.sampler.base.BaseSampler.set_initial_conditions>`
 * :py:meth:`run <pycbc.inference.sampler.base.BaseSampler.run>`
 * :py:meth:`checkpoint <pycbc.inference.sampler.base.BaseSampler.checkpoint>`
 * :py:attr:`samples <pycbc.inference.sampler.base.BaseSampler.samples>`
 * :py:attr:`model_stats <pycbc.inference.sampler.base.BaseSampler.model_stats>`
 * :py:meth:`finalize <pycbc.inference.sampler.base.BaseSampler.finalize>`

----------------
Detailed example
----------------

Let's examine the |EmceeEnsembleSampler| class to see how these guidelines
apply in practice. Here is its inheritance structure (click on the names of the
classes to see their documentation):

.. inheritance-diagram:: pycbc.inference.sampler.emcee
    :parts: 2

|

In addition to |BaseSampler|, |EmceeEnsembleSampler| inherits from |BaseMCMC|
and |MCMCAutocorrSupport|.  Inspecting |BaseMCMC|, we see that it implements
several of the methods that |BaseSampler| requires: namely,
|BaseMCMC.set_initial_conditions|, |BaseMCMC.run|, and |BaseMCMC.checkpoint|.
This is because the steps taken in these functions are common across MCMC
samplers. For example, in |BaseMCMC.run|, the sampler is run for blocks of
iterations (specified by |BaseMCMC.checkpoint_interval|) until the convergence
criteria has been met (which is determined by |BaseMCMC.set_target|). This,
generally, is what all MCMC samplers do.

*How* an MCMC sampler is run for some number of iterations is unique to each
sampling engine. To accommodate this, |BaseMCMC| adds |BaseMCMC.run_mcmc|,
which it calls from within its |BaseMCMC.run|.  |BaseMCMC.run_mcmc| is an
abstract method -- |BaseMCMC| is itself an abstract base class. Since
|EmceeEnsembleSampler| inherits from |BaseSampler|, followed by |BaseMCMC| (see
:ref:`note <python-inheritance-note>`), |BaseMCMC| fulfils |BaseSampler|'s
requirement that a ``run`` method be implemented, but replaces it with the
requirement that the class define a ``run_mcmc`` method.  As a result,
|EmceeEnsembleSampler| has its own |EmceeEnsembleSampler.run_mcmc|; this is
where the call to the underlying sampling engine (external to pycbc) is made.

.. _python-inheritance-note: 
.. note::
   In python, the order of inheritance is determined by the order the parents
   are given in the class definition, from right to left. For example,
   |EmceeEnsembleSampler| is defined as:

   .. code-block:: python

      class EmceeEnsembleSampler(MCMCAutocorrSupport, BaseMCMC, BaseSampler):

   This means that methods introduced by |BaseSampler| will be overridden by
   |BaseMCMC|, which in turn will be overridden by |MCMCAutocorrSupport|.  For
   this reason, all sampler class definitions must have |BaseSampler| listed
   last.

All MCMC samplers need to be able to compute an autocorrelation function (ACF)
and length (ACL). This is used to determine how to thin the chains to obtain
independent samples. Consequently, |BaseMCMC| also adds abstract base methods
|BaseMCMC.compute_acf| and |BaseMCMC.compute_acl|; these are called by its
|BaseMCMC.checkpoint| method.  The |MCMCAutocorrSupport| class provides these
functions. These functions are provided in a class separate from |BaseMCMC|
because not all MCMC samplers estimate ACF/Ls in the same way. For example,
multi-tempered samplers need to compute ACF/Ls separately for each temperature
chain. Consequently, there is an equivalent class
|MultiTemperedAutocorrSupport| that offers the same functions for
multi-tempered MCMCs. This class is used by, e.g., |EmceePTSampler| (see its
:ref:`inheritance diagram <inheritance-emcee_pt>`, below). By making the
compute ACF/L functions abstract base methods in |BaseMCMC|, both single and
multi-tempered MCMC samplers can inherit from |BaseMCMC|.


We see that by separating functionality out into support classes
and using multiple inheritance, we are able to provide support for all of the
unique features of different samplers, while keeping the base API that
``pycbc_inference`` interacts with simple.

---------------------
Inheritance diagrams
---------------------

Here are inheritance diagrams for all of the currently supported samplers:

.. include:: ../_include/sampler_inheritance_diagrams.rst

.. _how-to-add-a-sampler:

--------------------
How to add a sampler
--------------------

To add support for a new sampler, do the following:

1. Create a file in ``pycbc/inference/sampler`` for the new sampler's class.
2. Add the new class definition. The class must inherit from at least
   |BaseSampler|.
3. Give a name attribute to the class that is unique across the supported
   sampler classes.
4. :ref:`Add an IO class <how-to-add-a-sampler-io>` for the sampler to the
   :py:mod:`inference.io <pycbc.inference.io>` modules. Set your new class's
   ``io`` attribute to point to this new class.
5. Add any other methods you need to satisfy the |BaseSampler|'s required
   methods. When doing so, try to follow the guidelines above: do not duplicate
   code, and try to use support classes that offer functionality that you need.
   If you think some of the methods will be useful for more than just your
   sampler, create a new support class and add those methods to it.  However,
   if you're unsure what is available or you would have to make changes to the
   support classes that may break other samplers, just add the methods you need
   to your new class definition. Fixing code duplication or rearranging support
   classes can be done through the review process when you wish to add your new
   sampler to the main gwastro repository.
6. Add the sampler to the ``samplers`` dictionary in
   ``pycbc/inference/sampler/__init__.py`` so that ``pycbc_inference`` is aware
   of it to use.
7. When you're satisfied that your new sampler works,
   `file a pull request <https://github.com/gwastro/pycbc/blob/master/CONTRIBUTING.md#open-a-pull-request>`_ to get it into the main gwastro repostiory. Thank you for
   your contributions!

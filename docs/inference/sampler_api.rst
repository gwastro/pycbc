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
   :py:class:`BaseSampler <pycbc.inference.sampler.base.BaseSampler>` class. This is
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
5. All sampler classes need a corresponding class in the
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

----------------
Detailed example
----------------

Let's examine the
:py:class:`EmceeEnsembleSampler <pycbc.inference.sampler.emcee.EmceeEnsembleSampler>`
class to see how these guidelines apply in practice. Here is its inheritance
structure (click on the names of the classes to see their documentation):

.. _inheritance-emcee:
.. inheritance-diagram:: pycbc.inference.sampler.emcee
    :parts: 2
|

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
inherits from :py:class:`BaseSampler <pycbc.inference.sampler.base.BaseSampler>`
followed by :py:class:`BaseMCMC <pycbc.inference.sampler.base_mcmc.BaseMCMC>`
(see :ref:`note <python-inheritance-note>`),
:py:class:`BaseMCMC <pycbc.inference.sampler.base_mcmc.BaseMCMC>` fulfils
:py:class:`BaseSampler's <pycbc.inference.sampler.base.BaseSampler>` requirement
that a ``run`` method be implemented, but replaces it with the requirement that
a ``run_mcmc`` method be implemented. This is why
:py:class:`EmceeEmsembleSampler <pycbc.inference.sampler.emcee.EmceeEnsembleSampler>`
implements a :py:meth:`run_mcmc <pycbc.inference.sampler.emcee.EmceeEnsembleSampler.run_mcmc>`
method.

.. _python-inheritance-note: 
.. note::
   In python, the order of inheritance when a class inherits from multiple
   parents is determined by the order the parents are given in the class
   definition, from right to left. For example,
   :py:class:`EmceeEnsembleSampler <pycbc.inference.sampler.emcee.EmceeEnsembleSampler>`
   is defined as:

   .. code-block:: python

      class EmceeEnsembleSampler(MCMCAutocorrSupport, BaseMCMC, BaseSampler):

   This means that methods introduced by
   :py:class:`BaseSampler <pycbc.inference.sampler.base.BaseSampler>`
   will be overridden by
   :py:class:`BaseMCMC <pycbc.inference.sampler.base_mcmc.BaseMCMC>`, which in
   turn will be overridden by
   :py:class:`MCMCAutocorrSupport <pycbc.inference.sampler.base_mcmc.MCMCAutocorrSupport>`.
   For this reason, all sampler class definitions must have
   :py:class:`BaseSampler <pycbc.inference.sampler.base.BaseSampler>`
   listed last.

All MCMC samplers need to be able to compute an autocorrelation
function (ACF) and length (ACL). This is used to determine how to thin the chains
to obtain independent samples. Consequently, :py:class:`BaseMCMC <pycbc.inference.sampler.base_mcmc.BaseMCMC>`
also adds abstract base methods :py:meth:`compute_acf <pycbc.inference.sampler.base_mcmc.BaseMCMC.compute_acf>`
and :py:meth:`compute_acl <pycbc.inference.sampler.base_mcmc.BaseMCMC.compute_acl>`; these
are called by its :py:class:`checkpoint <pycbc.inference.sampler.base_mcmc.BaseMCMC.checkpoint>` method.
The :py:class:`MCMCAutocorrSupport <pycbc.inference.sampler.base_mcmc.MCMCAutocorrSupport>`
provides these functions. These
functions are provided in a class separate from :py:class:`BaseMCMC <pycbc.inference.sampler.base_mcmc.BaseMCMC>`
because not all MCMC samplers estimate ACF/Ls in the same way. For example,
multi-tempered samplers need to compute ACF/Ls separately for each temperature
chain. Consequently, there is an equivalent class,
:py:class:`MultiTemperedAutocorrSupport <pycbc.sampler.base_multitemper.MultiTemperedAutocorrSuppport>`
which offers the same functions for multi-tempered MCMCs. This class is used by,
e.g., :py:class:`EmceePTSampler <pycbc.sampler.emcee_pt.EmceePTSampler>` (see its
:ref:`inheritance diagram <inheritance-emcee_pt>`, below). By making the
compute ACF/L functions abstract base methods in
:py:class:`BaseMCMC <pycbc.inference.sampler.base_mcmc.BaseMCMC>`, both single and
multi-tempered MCMC samplers can inherit from
:py:class:`BaseMCMC <pycbc.inference.sampler.base_mcmc.BaseMCMC>`.


We see that by separating functionality out into support classes
and using multiple inheritance, we are able to provide support for all of the
unique features of different samplers, while keeping the base API that
``pycbc_inference`` interacts with simple.

---------------------
Inheritance diagrams
---------------------

Here are inheritance diagrams for the rest of the samplers currently supported:

.. _inheritance-emcee_pt:
* ``emcee_pt``:

.. inheritance-diagram:: pycbc.inference.sampler.emcee_pt
    :parts: 2
|

--------------------
How to add a sampler
--------------------

To add support for a new sampler, do the following:

1. Create a file in ``pycbc/inference/sampler`` for the new sampler's class.
2. Add the new class definition. The class must inherit from at least
   :py:class:`BaseSampler <pycbc.inference.sampler.base.BaseSampler>`.
3. Give a name attribute to the class that is unique across the supported
   sampler classes.
4. Add an IO class for the sampler to the :py:mod:`inference.io <pycbc.inference.io>`
   modules. Set your new class's ``io`` attribute to point to this new class.

5. Add any other methods you need to satisfy the
   :py:class:`BaseSampler's <pycbc.inference.sampler.base.BaseSampler>`
   required methods. When doing so, try to follow the guidelines above: do not
   duplicate code, and try to use support classes that offer functionality that
   you need. If you think some of the methods will be useful for more than just
   your sampler, create a new support class and add those methods to it.
   However, if you're unsure what is available or you would have to make
   changes to the support classes that may break other samplers, just add the
   methods you need to your new class definition. Fixing code duplication or
   rearranging support classes to accommodate can be done through the review
   process when you wish to add your new sampler to main gwastro repository.

6. Add the sampler to the ``samplers`` dictionary in
   ``pycbc/inference/sampler/__init__.py`` so that ``pycbc_inference`` is aware
   of it to use.
7. When you're satisfied that your new sampler works,
   `file a pull request <https://github.com/gwastro/pycbc/blob/master/CONTRIBUTING.md#open-a-pull-request>`_ to get it into the main gwastro repostiory. Thank you for
   your contributions!

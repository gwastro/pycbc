.. _waveform_plugin:

------------------------------------------------------
Making new waveform approximants available to PyCBC
------------------------------------------------------

=================================================
Adding a custom waveform model within a script
=================================================

By example, the following script shows how to write a waveform model
in the form required for PyCBC. We can also make this new waveform directly
accessible by using the :py:func:`~pycbc.waveform.plugin.add_custom_waveform` function.
If you are developing in a notebook or self-contained script, this may be
what you want to do. However, if you want to make your waveform available
to pycbc-based executables such as PyCBC Inference, also read the next
section.

There are two kinds of models you can make. In the example below, we
make a time-domain model. You can also make a freuqency-domain model. The only
difference is that your function should return an instance of :py:class:`~pycbc.types.frequencyseries.FrequencySeries` and
the required sample step option is `delta_f` instead of `delta_t`.

Each waveform generation function must take only keyword arguments, and
should be able to take an arbitrary number of them. You may add new parameters
as you like. These will be automatically useable by PyCBC Inference and
other pycbc codes.

Each waveform model must have an associate `approximant` name, which identifies
the model and distinguishes it from any other. If the name has already been
used, you should select a different name. By default, an error will be raised
unless overridden.

.. plot:: ../examples/waveform/add_waveform.py
   :include-source:

=================================================
Creating a plugin for PyCBC
=================================================

To make a waveform model universally available to PyCBC so it can be called
from PyCBC Inference, or the pycbc-based searched codes, you can create
a plugin package which advertises your model. PyCBC will automatically
detect your package and make your waveform model available for use.

The steps are:

 * Create a waveform model just like as in the above example
 * Create a python package for your module
 * In your packages setup.py advertise that it contains a PyCBC compatible
   waveform model in it's `entry_points` option.

Your `setup.py` should look like the following, the key addition being the `entry_points`
parameter passed to setup.py.

.. code-block:: python

    setup (
        name = 'pycbc-revchirp',
        version = VERSION,
        description = 'An example waveform plugin for PyCBC',
        long_description = open('descr.rst').read(),
        author = 'The PyCBC team',
        author_email = 'alex.nitz@gmail.org',
        url = 'http://www.pycbc.org/',
        download_url = 'https://github.com/gwastro/revchirp/tarball/v%s' % VERSION,
        keywords = ['pycbc', 'signal processing', 'gravitational waves'],
        install_requires = ['pycbc'],
        py_modules = ['revchirp'],
        entry_points = {"pycbc.waveform.td":"revchirp = revchirp:reverse_chirp_td",
                        "pycbc.waveform.fd":"revchirp = revchirp:reverse_chirp_fd"},
        classifiers=[
            'Programming Language :: Python',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.6',
            'Intended Audience :: Science/Research',
            'Natural Language :: English',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Astronomy',
            'Topic :: Scientific/Engineering :: Physics',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        ],
    )

The format for the `entry_points` is `"capability":"approximant_name = module_path:function_name"`.
The module path may include dots if the module is within a package or sub-package. The
valid `capbility` is `pycbc.waveform.td` and `pycbc.waveform.fd` for time and frequency
domain waveform models,  respectively.

For a complete working minimal example of a PyCBC waveform plugin, see the
example package on github to
`make a reversed-chirp waveform <https://github.com/gwastro/example-waveform-plugin>`_ .

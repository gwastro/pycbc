""" Utilities for handling waveform plugins
"""


def add_custom_waveform(approximant, function, domain, force=False):
    """ Make custom waveform available to pycbc

    Parameters
    ----------
    approximant : str
        The name of the waveform
    function : function
        The function to generate the waveform
     domain : str
        Either 'frequency' or 'time' to indicate the domain of the waveform.
    """
    from pycbc.waveform.waveform import cpu_fd, cpu_td

    if domain == 'time':
        if not force and (approximant in cpu_td):
            raise RuntimeError("Can't load plugin waveform {}, the name is"
                               " already in use.".format(approximant))
        cpu_td[approximant] = function
    elif domain == 'frequency':
        if not force and (approximant in cpu_fd):
            raise RuntimeError("Can't load plugin waveform {}, the name is"
                               " already in use.".format(approximant))
        cpu_fd[approximant] = function
    else:
        raise ValueError("Invalid domain ({}), should be "
                         "'time' or 'frequency'".format(domain))


def add_length_estimator(approximant, function):
    """ Add length estimator for an approximant

    Parameters
    ----------
    approximant : str
        Name of approximant
    function : function
        A function which takes kwargs and returns the waveform length
    """
    from pycbc.waveform.waveform import _filter_time_lengths
    if approximant in _filter_time_lengths:
        raise RuntimeError("Can't load length estimator {}, the name is"
                           " already in use.".format(approximant))
    _filter_time_lengths[approximant] = function


def retrieve_waveform_plugins():
    """ Process external waveform plugins
    """
    import pkg_resources

    # Check for fd waveforms
    for plugin in pkg_resources.iter_entry_points('pycbc.waveform.fd'):
        add_custom_waveform(plugin.name, plugin.resolve(), 'frequency')

    # Check for td waveforms
    for plugin in pkg_resources.iter_entry_points('pycbc.waveform.td'):
        add_custom_waveform(plugin.name, plugin.resolve(), 'time')

    # Check for wavveform length estimates
    for plugin in pkg_resources.iter_entry_points('pycbc.waveform.length'):
        add_length_estimator(plugin.name, plugin.resolve())

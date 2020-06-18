""" Utilities for handling waveform plugins
"""

def add_custom_waveform(approximant, function, domain):
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
        cpu_td[approximant] = function
    elif domain == 'frequency':
        cpu_fd[approximant] = function
    else:
        raise ValueError("Invalid domain ({}), should be "
                         "'time' or 'frequency'".format(domain))

def add_length_estimator(approximant, function):
    from pycbc.waveform.waveform import _filter_time_lengths
    _filter_time_lengths[approximant] = function

def retrieve_waveform_plugins():
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

# Copyright (C) 2016  Collin Capano
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
"""
This modules provides functions for formatting values into strings for display.
"""

import numpy

mjax_header = """
<script type="text/x-mathjax-config">
  MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$']]}});
</script>
<script type="text/javascript"
    src="//cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>
"""


def mathjax_html_header():
    """Standard header to use for html pages to display latex math.

    Returns
    -------
    header: str
        The necessary html head needed to use latex on an html page.
    """
    return mjax_header

def drop_trailing_zeros(num):
    """
    Drops the trailing zeros in a float that is printed.
    """
    txt = '%f' %(num)
    txt = txt.rstrip('0')
    if txt.endswith('.'):
        txt = txt[:-1]
    return txt

def get_signum(val, err, max_sig=numpy.inf):
    """
    Given an error, returns a string for val formated to the appropriate
    number of significant figures.
    """
    coeff, pwr = ('%e' % err).split('e')
    if pwr.startswith('-'):
        pwr = int(pwr[1:])
        if round(float(coeff)) == 10.:
            pwr -= 1
        pwr = min(pwr, max_sig)
        tmplt = '%.' + str(pwr+1) + 'f'
        return tmplt % val
    else:
        pwr = int(pwr[1:])
        if round(float(coeff)) == 10.:
            pwr += 1
        # if the error is large, we can sometimes get 0;
        # adjust the round until we don't get 0 (assuming the actual
        # value isn't 0)
        return_val = round(val, -pwr+1)
        if val != 0.:
            loop_count = 0
            max_recursion = 100
            while return_val == 0.:
                pwr -= 1
                return_val = round(val, -pwr+1)
                loop_count += 1
                if loop_count > max_recursion:
                    raise ValueError("Maximum recursion depth hit! Input " +\
                        "values are: val = %f, err = %f" %(val, err))
        return drop_trailing_zeros(return_val)

def format_value(value, error, plus_error=None, use_scientific_notation=3,
        include_error=True, use_relative_error=False, ndecs=None):
    """Given a numerical value and some bound on it, formats the number into a
    string such that the value is rounded to the nearest significant figure,
    which is determined by the error = abs(value-bound).

    Note: if either use_scientific_notation or include_error are True, the
    returned string will include LaTeX characters.

    Parameters
    ----------
    value : float
        The value to format.
    error : float
        The uncertainty in the value. This is used to determine the
        number of significant figures to print. If the value has no
        uncertainty, you can just do value*1e-k, where k+1 is the number
        of significant figures you want.
    plus_error : {None, float}
        The upper uncertainty on the value; i.e., what you need to add to the
        value to get its upper bound. If provided, ``error`` is assumed to be
        the negative; i.e., value +plus_error -error. The number of
        significant figures printed is determined from min(error,
        plus_error).
    use_scientific_notation : int, optional
        If ``abs(log10(value))`` is greater than the given, the return string
        will be formated to "\%.1f \\times 10^{p}", where p is the powers of 10
        needed for the leading number in the value to be in the singles spot.
        Otherwise will return "\%.(p+1)f". Default is 3.  To turn off, set to
        ``numpy.inf``. Note: using scientific notation assumes that the
        returned value will be enclosed in LaTeX math mode.
    include_error : {True, bool}
        Include the error in the return string; the output will be formated
        val \\pm err, where err is the error rounded to the same
        power of 10 as val. Otherwise, just the formatted value will
        be returned. If plus_error is provided then the return text will be
        formatted as ``val^{+plus_error}_{-error}``.
    use_relative_error : {False, bool}
        If include_error, the error will be formatted as a percentage of the
        the value.
    ndecs: {None, int}
        Number of values after the decimal point. If not provided,
        it will default to the number of values in the error.

    Returns
    -------
    string
        The value (and error, if include_error is True) formatted as a string.


    Examples
    --------
    Given a value and its uncertainty:

    >>> val, err
    (3.9278372067613837e-22, 2.2351435286500487e-23)

    Format with error quoted:

    >>> format_value(val, err)
    '3.93 \\pm 0.22\\times 10^{-22}'

    Quote error as a relative error:

    >>> format_value(val, err, use_relative_error=True)
    '3.93 \\times 10^{-22} \\pm5.6\\%'

    Format without the error and without scientific notation:

    >>> format_value(val, err, use_scientific_notation=float('inf'),
                     include_error=False)
    '0.000000000000000000000393'

    Given an plus error:

    >>> err_plus
    8.2700310560051804e-24

    Format with both bounds quoted:

    >>> format_value(val, err, plus_error=err_plus)
    '3.928^{+0.083}_{-0.224}\\times 10^{-22}'

    Format with both bounds quoted as a relative error:

    >>> format_value(val, err, plus_error=err_plus, use_relative_error=True)
    '3.928\\times 10^{-22}\\,^{+2.1\\%}_{-5.7\\%}'

    """
    minus_sign = '-' if value < 0. else ''
    value = abs(value)
    minus_err = abs(error)
    if plus_error is None:
        plus_err = minus_err
    else:
        plus_err = abs(plus_error)
    error = min(minus_err, plus_err)
    if value == 0. or abs(numpy.log10(value)) < use_scientific_notation:
        conversion_factor = 0.
    else:
        conversion_factor = numpy.floor(numpy.log10(value))
    value = value * 10**(-conversion_factor)
    error = error * 10**(-conversion_factor)
    if conversion_factor == 0.:
        powfactor = ''
    elif conversion_factor == 1.:
        powfactor = r'\times 10'
    else:
        powfactor = r'\times 10^{%i}' %(int(conversion_factor))

    if ndecs is not None:
        decs = value * 10**(-ndecs)
    else:
        decs = error
    # now round the the appropriate number of sig figs
    valtxt = get_signum(value, decs)
    valtxt = '{}{}'.format(minus_sign, valtxt)

    if include_error:
        if plus_error is None:
            errtxt = get_signum(error, error)
            if use_relative_error and float(valtxt) != 0.:
                relative_err = 100.*float(errtxt)/float(valtxt)
                # we round the relative error to the nearest 1% using
                # get_signum; Note that if the relative error is < 1%,
                # get_signum will automatically increase the number of values
                # after the decimal until it gets to the first non-zero value
                relative_err = get_signum(relative_err, 1.)
                txt = r'%s %s \pm%s\%%' %(valtxt, powfactor, relative_err)
            else:
                txt = r'%s \pm %s%s' %(valtxt, errtxt, powfactor)
        else:
            plus_err = plus_err * 10**(-conversion_factor)
            minus_err = minus_err * 10**(-conversion_factor)
            minus_err_txt = get_signum(minus_err, decs)
            plus_err_txt = get_signum(plus_err, decs)
            if use_relative_error and float(valtxt) != 0.:
                # same as above, but with plus and minus
                rel_plus_err = get_signum(
                    100.*float(plus_err_txt)/float(valtxt), 1.)
                rel_minus_err = get_signum(
                    100.*float(minus_err_txt)/float(valtxt), 1.)
                txt = r'%s%s\,^{+%s\%%}_{-%s\%%}' %(valtxt, powfactor,
                    rel_plus_err, rel_minus_err)
            else:
                txt = r'%s^{+%s}_{-%s}%s' %(valtxt, plus_err_txt,
                    minus_err_txt, powfactor)
    else:
        txt = r'%s%s' %(valtxt, powfactor)
    return txt

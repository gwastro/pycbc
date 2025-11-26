# Copyright (C) 2012 Alex Nitz, Tito Dal Canton
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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
"""Utilities to read PSDs from files.
"""

import logging
import numpy
import xml.etree.ElementTree as ET
import io

import scipy.interpolate
from igwn_ligolw import utils as ligolw_utils

from pycbc.types import FrequencySeries

logger = logging.getLogger('pycbc.psd.read')

def from_numpy_arrays(freq_data, noise_data, length, delta_f, low_freq_cutoff):
    """Interpolate n PSD (as two 1-dimensional arrays of frequency and data)
    to the desired length, delta_f and low frequency cutoff.

    Parameters
    ----------
    freq_data : array
        Array of frequencies in hertz.
    noise_data : array
        PSD values corresponding to frequencies in freq_arr.
    length : int
        Length of the frequency series in samples.
    delta_f : float
        Frequency resolution of the frequency series in hertz.
    low_freq_cutoff : float
        Frequencies below this value (in hertz) are set to zero.

    Returns
    -------
    psd : FrequencySeries
        The generated frequency series.
    """
    # Only include points above the low frequency cutoff
    if freq_data[0] > low_freq_cutoff:
        raise ValueError(
            f'Lowest frequency in input PSD data ({freq_data[0]} Hz) is '
            f'higher than requested low-frequency cutoff ({low_freq_cutoff} Hz)'
        )

    kmin = int(low_freq_cutoff / delta_f)
    flow = kmin * delta_f

    data_start = (0 if freq_data[0]==low_freq_cutoff else numpy.searchsorted(freq_data, flow) - 1)
    data_start = max(0, data_start)
    # If the cutoff is exactly in the file, start there
    if freq_data[data_start+1] == low_freq_cutoff:
        data_start += 1

    freq_data = freq_data[data_start:]
    noise_data = noise_data[data_start:]

    if (length - 1) * delta_f > freq_data[-1]:
        logger.warning('Requested number of samples exceeds the highest '
                       'available frequency in the input data, '
                       'will use max available frequency instead. '
                       '(requested %f Hz, available %f Hz)',
                       (length - 1) * delta_f, freq_data[-1])
        length = int(freq_data[-1]/delta_f + 1)

    flog = numpy.log(freq_data)
    slog = numpy.log(noise_data)

    psd_interp = scipy.interpolate.interp1d(
        flog, slog, fill_value=(slog[0], slog[-1]), bounds_error=False)
    psd = numpy.zeros(length, dtype=numpy.float64)

    vals = numpy.log(numpy.arange(kmin, length) * delta_f)
    psd[kmin:] =  numpy.exp(psd_interp(vals))

    return FrequencySeries(psd, delta_f=delta_f)


def from_txt(filename, length, delta_f, low_freq_cutoff, is_asd_file=True):
    """Read an ASCII file containing one-sided ASD or PSD data and generate
    a frequency series with the corresponding PSD. The ASD or PSD data is
    interpolated in order to match the desired frequency resolution.

    Parameters
    ----------
    filename : string
        Path to a two-column ASCII file. The first column must contain
        the frequency in hertz (positive frequencies only) and the second column
        must contain the amplitude density OR power spectral density.
    length : int
        Length of the frequency series in samples.
    delta_f : float
        Frequency resolution of the frequency series in hertz.
    low_freq_cutoff : float
        Frequencies below this value (in hertz) are set to zero.
    is_asd_file : Boolean
        If false assume that the second column holds power spectral density.
        If true assume that the second column holds amplitude spectral density.
        Default: True

    Returns
    -------
    psd : FrequencySeries
        The generated frequency series.

    Raises
    ------
    ValueError
        If the ASCII file contains negative, infinite or NaN frequencies
        or amplitude densities.
    """
    file_data = numpy.loadtxt(filename)
    if (file_data < 0).any() or \
                            numpy.logical_not(numpy.isfinite(file_data)).any():
        raise ValueError('Invalid data in ' + filename)

    freq_data = file_data[:, 0]
    noise_data = file_data[:, 1]
    if is_asd_file:
        noise_data = noise_data ** 2

    return from_numpy_arrays(freq_data, noise_data, length, delta_f,
                             low_freq_cutoff)


def _strip_ns(tag):
    # remove namespace if present
    return tag.split('}', 1)[-1] if '}' in tag else tag

def _extract_array_text(elem):
    # find a child element that likely contains the numerical array
    # common candidates: 'Array', 'data', 'ArrayData'
    for child in elem.iter():
        tag = _strip_ns(child.tag)
        if tag.lower() in ('array', 'arraydata', 'data', 'values') and child.text:
            return child.text
    # fallback: if element itself contains whitespace-separated numbers
    if elem.text and any(ch.isdigit() for ch in elem.text):
        return elem.text
    return None

def _parse_psd_xmldoc(xml_bytes, root_name):
    """
    Parse bytes into ElementTree and extract PSD entries.

    Parameters
    ----------
    xml_bytes: byte string
        The xml document as a bytestring
    root_name: string
        If given use this as the root name for the PSD XML file.

    Returns
    -------
    dict: The PSDs

    """

    # Avoid importing lal.series here by parsing the LIGOLW XML directly.
    # Use igwn_ligolw's load_fileobj to transparently handle compression,
    # then parse the resulting bytes with ElementTree. This keeps the
    # runtime dependency on lal.series optional.
    try:
        tree = ET.parse(io.BytesIO(xml_bytes))
        root = tree.getroot()
    except Exception:
        # Try parsing as text
        root = ET.fromstring(xml_bytes.decode('utf-8'))

    psd_entries = {}
    # Find all elements whose tag matches root_name (ignoring namespace)
    for elem in root.iter():
        if _strip_ns(elem.tag) == root_name:
            # identify name / ifo string
            name = None
            # common places: attribute 'Name', child <Name>, child <ifo>
            if 'Name' in elem.attrib:
                name = elem.attrib.get('Name')
            if not name:
                name_tag = elem.find('.//{*}Name')
                if name_tag is not None and name_tag.text:
                    name = name_tag.text.strip()
            if not name:
                ifo_tag = elem.find('.//{*}ifo')
                if ifo_tag is not None and ifo_tag.text:
                    name = ifo_tag.text.strip()
            # deltaF may appear as child <deltaF> or attribute
            delta_f = None
            if 'deltaF' in elem.attrib:
                try:
                    delta_f = float(elem.attrib['deltaF'])
                except Exception:
                    delta_f = None
            if delta_f is None:
                df_tag = elem.find('.//{*}deltaF')
                if df_tag is not None and df_tag.text:
                    try:
                        delta_f = float(df_tag.text.strip())
                    except Exception:
                        delta_f = None

            # extract numerical array text and convert to numpy array
            arr_text = _extract_array_text(elem)
            noise_data = None
            if arr_text:
                # the array text may contain newlines, commas or extra spaces
                # parse floats by splitting on whitespace and commas
                import re

                tokens = re.split('[,\s]+', arr_text.strip())
                # filter out any empty tokens
                tokens = [t for t in tokens if t]
                try:
                    noise_data = numpy.array([float(t) for t in tokens], dtype=numpy.float64)
                except Exception:
                    noise_data = None

            if noise_data is None or delta_f is None:
                # skip entries we couldn't parse
                continue

            # if name still None, try to fall back to an index-based name
            if name is None:
                name = f'psd_{len(psd_entries)}'

            psd_entries[name] = {'data': noise_data, 'delta_f': delta_f}

    return psd_entries

def from_xml(filename, length, delta_f, low_freq_cutoff, ifo_string=None,
             root_name='psd'):
    """Read a LIGOLW XML file containing one-sided PSD data and generate
    a frequency series with the corresponding PSD. The data is interpolated in
    order to match the desired frequency resolution.

    Parameters
    ----------
    filename : string
        Path to a LIGOLW XML file.
    length : int
        Length of the frequency series in samples.
    delta_f : float
        Frequency resolution of the frequency series in hertz.
    low_freq_cutoff : float
        Frequencies below this value (in hertz) are set to zero.
    ifo_string : string
        Use the PSD in the file's PSD dictionary with this ifo string.
        If not given, and only one PSD present in the file, return that; if not
        given, and multiple (or zero) PSDs present, then an exception is raised.
    root_name : string (default='psd')
        If given use this as the root name for the PSD XML file. If this means
        nothing to you, then it is probably safe to ignore this option.

    Returns
    -------
    psd : FrequencySeries
        The generated frequency series.
    """
    

    with open(filename, 'rb') as fp:
        # load_fileobj handles compression transparently and returns raw bytes
        xml_bytes = ligolw_utils.load_fileobj(fp, compress='auto')
        psd_dict = _parse_psd_xmldoc(xml_bytes, root_name=root_name)

    if ifo_string is not None:
        try:
            psd_freq_series = psd_dict[ifo_string]
        except KeyError:
            raise KeyError(f"PSD with ifo string '{ifo_string}' not found in XML file")
    elif len(psd_dict.keys()) == 1:
        psd_freq_series = psd_dict[next(iter(psd_dict.keys()))]
    else:
        raise ValueError(
            "No ifo string given and input XML file does not contain exactly one PSD. Specify which PSD you want to use."
        )

    noise_data = psd_freq_series['data']
    freq_data = numpy.arange(len(noise_data)) * psd_freq_series['deltaF']

    return from_numpy_arrays(freq_data, noise_data, length, delta_f,
                             low_freq_cutoff)

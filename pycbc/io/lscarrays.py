# Copyright (C) 2015  Collin Capano
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
#                           Preamble
#
# =============================================================================
#
"""
This modules provides definitions of, and helper functions for, LSCArrays.
LSCArrays are wrappers of numpy recarrays with additional functionality useful for storing and retrieving data created by a search for gravitational waves.
"""

import os, sys, types
import re
import copy
import hashlib
import inspect
import numpy
import numpy.lib.recfunctions as recfunctions

import lal
import lalsimulation as lalsim
from glue.ligolw import types as ligolw_types

#
# =============================================================================
#
#                           Data type mappings
#
# =============================================================================
#
# add ligolw_types to numpy typeDict
numpy.typeDict.update(ligolw_types.ToNumPyType)

# Annoyingly, numpy has no way to store NaNs in an integer field to indicate
# the equivalent of None. This can be problematic for fields that store ids:
# if an array has an id field with value 0, it isn't clear if this is because
# the id is the first element, or if no id was set. To clear up the ambiguity,
# we define here an integer to indicate 'id not set'. 
ID_NOT_SET = -1
EMPTY_OBJECT = None

def set_default_empty(array):
    if array.dtype.names is None:
        # scalar dtype, just set
        if array.dtype.str[1] == 'i':
            # integer, set to ID_NOT_SET
            array[:] = ID_NOT_SET
        elif array.dtype.str[1] == 'O':
            # object, set to EMPTY_OBJECT
            array[:] = EMPTY_OBJECT
    else:
        for name in array.dtype.names:
            set_default_empty(array[name])

def default_empty(shape, dtype):
    """Numpy's empty array can have random values in it. To prevent that, we
    define here a default emtpy array. This default empty is a numpy.zeros
    array, except that objects are set to None, and all ints to ID_NOT_SET.
    """
    default = numpy.zeros(shape, dtype=dtype)
    set_default_empty(default)
    return default

# set default data types
_default_types_status = {
    'default_strlen': 50,
    'ilwd_as_int': True,
    'lstring_as_obj': False
}

def lstring_as_obj(true_or_false=None):
    """Toggles whether lstrings should be treated as strings or as objects.
    When lscarrays is first loaded, the default is True.

    Parameters
    ----------
    true_or_false : {None|bool}
        Pass True to map lstrings to objects; False otherwise. If None
        provided, just returns the current state.

    Return
    ------
    current_stat : bool
        The current state of lstring_as_obj.

    Examples
    --------
    >>> from pycbc.io import lscarrays
    >>> lscarrays.lstring_as_obj()
        True
    >>> lscarrays.LSCArray.from_arrays([numpy.zeros(10)], dtype=[('foo', 'lstring')])
    LSCArray([(0.0,), (0.0,), (0.0,), (0.0,), (0.0,), (0.0,), (0.0,), (0.0,),
           (0.0,), (0.0,)], 
          dtype=[('foo', 'O')])
    >>> lscarrays.lstring_as_obj(False)
        False
    >>> lscarrays.LSCArray.from_arrays([numpy.zeros(10)], dtype=[('foo', 'lstring')])
    LSCArray([('0.0',), ('0.0',), ('0.0',), ('0.0',), ('0.0',), ('0.0',),
           ('0.0',), ('0.0',), ('0.0',), ('0.0',)], 
          dtype=[('foo', 'S50')])
    """
    if true_or_false is not None:
        _default_types_status['lstring_as_obj'] = true_or_false
        # update the typeDict
        numpy.typeDict[u'lstring'] = numpy.object_ \
            if _default_types_status['lstring_as_obj'] \
            else 'S%i' % _default_types_status['default_strlen']
    return _default_types_status['lstring_as_obj']

def ilwd_as_int(true_or_false=None):
    """Similar to lstring_as_obj, sets whether or not ilwd:chars should be
    treated as strings or as ints. Default is True.
    """
    if true_or_false is not None:
        _default_types_status['ilwd_as_int'] = true_or_false
        numpy.typeDict[u'ilwd:char'] = int \
            if _default_types_status['ilwd_as_int'] \
            else 'S%i' % default_strlen
    return _default_types_status['ilwd_as_int']


def default_strlen(strlen=None):
    """Sets the default string length for lstring and ilwd:char, if they are
    treated as strings. Default is 50.
    """
    if strlen is not None:
        _default_types_status['default_strlen'] = strlen
        # update the typeDicts as needed
        lstring_as_obj(_default_types_status['lstring_as_obj'])
        set_ilwd_as_int(_default_types_status['ilwd_as_int'])
    return _default_types_status['default_strlen']

# set the defaults
lstring_as_obj(True)
ilwd_as_int(True)


#
# =============================================================================
#
#                           Helper functions 
#
# =============================================================================
#
def get_dtype_descr(dtype):
    """Numpy's ``dtype.descr`` will return empty void fields if a dtype has
    offsets specified. This function tries to fix that by not including
    fields that have no names and are void types.
    """
    return [dt for dt in dtype.descr if not (dt[0] == '' and dt[1][1] == 'V')]


def combine_fields(dtypes):
    """Combines the fields in the list of given dtypes into a single dtype.

    Parameters
    ----------
    dtypes : (list of) numpy.dtype(s)
        Either a numpy.dtype, or a list of numpy.dtypes.

    Returns
    -------
    numpy.dtype
        A new dtype combining the fields in the list of dtypes.
    """
    if not isinstance(dtypes, list):
        dtypes = [dtypes]
    # Note: incase any of the dtypes have offsets, we won't include any fields
    # that have no names and are void
    new_dt = numpy.dtype([dt for dtype in dtypes \
        for dt in get_dtype_descr(dtype)])
    return new_dt


def merge_arrays(merge_list, names=None, flatten=True, outtype=None):
    """Merges the given arrays into a single array. The arrays must all have
    the same shape. If one or more of the given arrays has multiple fields,
    all of the fields will be included as separate fields in the new array.

    Parameters
    ----------
    merge_list : list of arrays
        The list of arrays to merge.
    names : {None | sequence of strings}
        Optional, the names of the fields in the output array. If flatten is
        True, must be the same length as the total number of fields in
        merge_list.  Otherise, must be the same length as the number of
        arrays in merge_list.  If None provided, and flatten is True, names
        used will be the same as the name of the fields in the given arrays.
        If the datatype has no name, or flatten is False, the new field will
        be `fi` where i is the index of the array in arrays.
    flatten : bool
        Make all of the fields in the given arrays separate fields in the
        new array. Otherwise, each array will be added as a field. If an
        array has fields, they will be subfields in the output array. Default
        is True.
    outtype : {None | class}
        Cast the new array to the given type. Default is to return a
        numpy structured array.

    Returns
    -------
    new array : {numpy.ndarray | outtype}
        A new array with all of the fields in all of the arrays merged into
        a single array.
    """
    if not all(merge_list[0].shape == arr.shape for arr in merge_list):
        raise ValueError("all of the arrays in merge_list must have the " +
            "same shape")
    if flatten:
        new_dt = combine_fields([arr.dtype for arr in merge_list])
    else:
        new_dt = numpy.dtype([('f%i' %ii, arr.dtype.descr) \
            for ii,arr in enumerate(merge_list)])
    new_arr = numpy.empty(merge_list[0].shape, dtype=new_dt)
    # ii is a counter to keep track of which fields from the new array
    # go with which arrays in merge list
    ii = 0
    for arr in merge_list:
        if arr.dtype.names is None:
            new_arr[new_dt.names[ii]] = arr
            ii += 1
        else:
            for field in arr.dtype.names:
                new_arr[field] = arr[field]
                ii += 1
    # set the names if desired
    if names is not None:
        new_arr.dtype.names = names
    # ditto the outtype
    if outtype is not None:
        new_arr = new_arr.view(type=outtype)
    return new_arr


def add_fields(input_array, arrays, names=None, assubarray=False):
    """Adds the given array(s) as new field(s) to the given input array.
    Returns a new instance of the input_array with the new fields added.

    Parameters
    ----------
    input_array : instance of a numpy.ndarray or numpy recarray
        The array to to add the fields to.
    arrays : (list of) numpy array(s)
        The arrays to add. If adding multiple arrays, must be a list;
        if adding a single array, can just be that array.
    names : (list of) strings
        Optional, the name(s) of the new fields in the output array. If
        adding multiple fields, must be a list of strings with the same
        length as the list of arrays. If None provided, names used will
        be the same as the name of the datatype in the given arrays.
        If the datatype has no name, the new field will be ``'fi'`` where
        i is the index of the array in arrays.
    assubarray : bool
        Add the list of arrays as a single subarray field. If True, and names
        provided, names should be a string or a length-1 sequence. Default is
        False, in which case each array will be added as a separate field.

    Returns
    -------
    new_array : new instance of `input_array`
        A copy of the `input_array` with the desired fields added.
    """
    if not isinstance(arrays, list):
        arrays = [arrays]
    # set the names
    if names is not None:
        if isinstance(names, str) or isinstance(names, unicode):
            names = [names]
        # check if any names are subarray names; if so, we have to add them
        # separately
        subarray_names = [name for name in names if len(name.split('.')) > 1]
    else:
        subarray_names = []
    if any(subarray_names):
        subarrays = [arrays[ii] for ii,name in enumerate(names) \
            if name in subarray_names]
        # group together by subarray
        groups = {}
        for name,arr in zip(subarray_names, subarrays):
            key = name.split('.')[0]
            subkey = '.'.join(name.split('.')[1:])
            try:
                groups[key].append((subkey, arr))
            except KeyError:
                groups[key] = [(subkey, arr)]
        # now cycle over the groups, adding all of the fields in each group
        # as a subarray
        for group_name in groups:
            # we'll create a dictionary out of the subarray field names ->
            # subarrays
            thisdict = dict(groups[group_name])
            # check if the input array has this field; if so, remove it, then
            # add it back with the other new arrays
            if group_name in input_array.fieldnames:
                # get the data
                new_subarray = input_array[group_name]
                # add the new fields to the subarray
                new_subarray = add_fields(new_subarray, thisdict.values(),
                    thisdict.keys())
                # remove the original from the input array
                input_array = input_array.without_fields(group_name)
            else:
                new_subarray = thisdict.values()
            # add the new subarray to input_array as a subarray
            input_array = add_fields(input_array, new_subarray,
                names=group_name, assubarray=True)
            # set the subarray names
            input_array[group_name].dtype.names = thisdict.keys()
        # remove the subarray names from names 
        keep_idx = [ii for ii,name in enumerate(names) \
            if name not in subarray_names]
        names = [names[ii] for ii in keep_idx]
        # if there's nothing left, just return
        if names == []:
            return input_array
        # also remove the subarray arrays
        arrays = [arrays[ii] for ii in keep_idx] 
    if assubarray:
        # merge all of the arrays into a single array
        if len(arrays) > 1:
            arrays = [merge_arrays(arrays, flatten=True)]
        # now merge all the fields as a single subarray
        merged_arr = numpy.empty(len(arrays[0]),
            dtype=[('f0', arrays[0].dtype.descr)])
        merged_arr['f0'] = arrays[0]
        arrays = [merged_arr]
    merge_list = [input_array] + arrays
    if names is not None:
        names = list(input_array.dtype.names) + names
    # merge into a single array
    return merge_arrays(merge_list, names=names, flatten=True,
        outtype=type(input_array))


def get_fields(input_array, names, copy=False, outtype=None):
    """Given an array with named fields, creates a new LSCArray with the given
    fields. Only fields in input_array that are in the list of names will be
    extracted. All the fields listed in names must be present in the
    input_array. The method for creating a view is from:
    <http://stackoverflow.com/a/21819324/1366472>

    Parameters
    ----------
    input_array : array-like
        Anything that can be cast as a numpy array. The resulting array
        must have a dtype with at least one field in names.
    names : {strings|list of strings}
        List of field names for the returned array; may be either a single
        name or a list of names. All of the names must be a field in the input
        array.
    copy : bool, optional
        If True, will force a copy of the input array rather than a view.
    outtype : optional
        Make the output array be the given type. If None, the type of the
        output will be the same as the input.

    Returns
    -------
    output : {type(input_array)|outtype}
        An view or copy of the input array with the given names.
    """
    if outtype is None:
        outtype = type(input_array)
    if isinstance(names, str) or isinstance(names, unicode):
        names = [names]
    if copy:
        drop_fields = [name for name in input_array.dtype.names \
            if name not in names]
        return recfunctions.rec_drop_fields(input_array, drop_fields).view(
            type=outtype)
    else:
        new_dtype = numpy.dtype({name: input_array.dtype.fields[name] \
            for name in names})
        return numpy.ndarray(input_array.shape, new_dtype,
            input_array, 0, input_array.strides).view(type=outtype)
    

def build_lookup_table(input_array):
    """Given an array, builds a dictionary of the values of the array.
    """
    unique_vals, unique_idx, map_idx = numpy.unique(input_array,
        return_index=True, return_inverse=True)
    # if the number of unique values is the same as the number
    # of elements in self, this is a one-to-one mapping
    if unique_vals.shape == input_array.shape:
        try:
            return dict(zip(unique_vals, unique_idx))
        # if unique vals is a list (as in the case when looking up multiple
        # fields), we'll get a TypeError; in that case, cast to tuple
        except TypeError:
            return dict(zip(tuple(unique_vals.tolist()), unique_idx))
    # else, we need to be set each key value pointing to the
    # appropriate indices
    else:
        try:
            return dict([[unique_vals[ii], numpy.where(map_idx == ii)] \
                for ii in range(len(unique_vals))
                ])
        except TypeError:
            return dict([[tuple(unique_vals[ii]), numpy.where(map_idx == ii)]\
                for ii in range(len(unique_vals))
                ])

def copy_attributes(this_array, other_array):
    """Copies attributes of `other_array` that are not in `this_array` to
    `this_array`. Checks if each attribute is a property, a method, or
    a plain attribute, adding each attribute accordingly. Returns a view of
    `this_array` with all of the attributes added.
    """
    # create a view of this array so as not to add the attributes to the 
    # original
    new_array = this_array.view(type=type(this_array))
    # get attributes in other that are not in this array, and figure out what
    # are properties, what are methods, and what are just plain attributes
    new_attrs = set(dir(other_array)) - set(dir(this_array))
    properties = {}
    methods = {}
    for attrname in new_attrs:
        # check if the attribute is a property
        try:
            attr = getattr(type(other_array), attrname)
            isproperty = isinstance(attr, property)
        except AttributeError:
            # we'll get an attribute error if attrname is just an attribute
            # of this instance; in that case, it cannot be a property
            attr = getattr(other_array, attrname)
            isproperty = False
        if isproperty:
            properties[attrname] = attr.fget
        elif inspect.ismethod(attr):
            # is an instance method
            methods[attrname] = attr
        else:
            # just an ordinary attribute, just add it
            setattr(new_array, attrname, attr)
    # now add each thing accordingly
    if properties != {}:
        new_array = new_array.add_properties(properties.keys(),
            properties.values())
    if methods != {}:
        new_array.add_methods(methods.keys(), methods.values())
    return new_array

def join_arrays(this_array, other_array, map_field, expand_field_name=None,
        other_map_field=None, get_fields=None, map_indices=None,
        copy_methods=True):
    """Joins `other_array` to `this_array` using the provided map fields.
    The information from `other_array` is added to `this_array` as a
    subarray. For a given element in `this_array`, all elements in
    `other_array` with
    ``this_array[map_field] == other_array[(other_)map_field]`` are
    retrieved. If multiple elements in `other_array` map to a single
    element in `this_array`, the expanded sub-array will have a shape equal
    to the maximum number of elements in `other_array` that map to a single
    element in `this_array`. Instance methods and properties of
    `other_array` that are not in `this_array` will also be copied.

    Parameters
    ----------
    this_array : any subclass of numpy array
        The array to add the fields to.
    other_array : LSCArray or similar
        The array that information is retrieved from. Must have `lookup` and
        `with_fields` methods.
    map_field : string
        The name of the field in `this_array` to use for mapping.
    expand_field_name : {None|string}
        If provided, all of the fields from `other_array` will be added
        as a subfield with name `expand_field_name` to the output array.
        Otherwise, all requested fields (see `get_fields`) are added as
        fields to `this_array`.
    other_map_field : {None | string}
        The name of the field in `other_array` to use for mapping. If None,
        `map_field` will be used.
    get_fields : {None | (list of) strings}
        Optionally specify what fields to retrieve from `other_array`. If
        None provided, will get all the fields in `other_array`.
    map_indices : {None | array of ints}
        If provided, will only map rows in `this_array` that have indices in
        the given array of indices. Any rows that are skipped will have a
        zeroed element in the new fields of the returned array. If None (the
        default), all rows in `this_array` are mapped.
    copy_methods : bool
        Copy instance methods and properties of `other_array` that are not
        in `this_array` to `this_array`. Default is True.

    Returns
    -------
    new_array : type(this_array)
        A copy of `this_array` with the mapped fields added.
    """
    if other_map_field is None:
        other_map_field = map_field
    # strip off the expand field in this_array, if is in the array
    if expand_field_name is not None:
        new_array = this_array.without_fields(expand_field_name)
    else:
        new_array = this_array
    # get a view of the other_array with just the desired fields;
    # note: this will also give a clean lookup table
    if isinstance(get_fields, str) or isinstance(get_fields, unicode):
        get_fields = [get_fields]
    elif get_fields is None:
        get_fields = list(other_array.dtype.names)
    # ensure the map field is included
    if other_map_field not in get_fields:
        get_fields.append(other_map_field)
    # XXX: I've found that copying other array is necessary here, otherwise
    # I get corrupted values in expanded_info, below.
    other_array = other_array.with_fields(get_fields, copy=True)
    # set an empty default in case one or map values in this array is not
    # found in other array
    # Note: for some strange reason, running dtype.descr will yield void fields
    # if with_fields is not a copy, so we'll strip those out
    other_dtdescr = get_dtype_descr(other_array.dtype)
    default = default_empty(1, dtype=other_dtdescr).view(
        type=type(other_array))
    if map_indices is not None:
        # only map rows whose indices are listed in map inidices
        mask = numpy.zeros(this_array.size, dtype=bool)
        mask[map_indices] = True
        expanded_info = [
            other_array.lookup(
                other_map_field, this_array[map_field][ii], default
            ) if applymap else default \
            for ii,applymap in enumerate(mask)]
    else:
        expanded_info = [other_array.lookup(other_map_field, mapval, default) \
            for mapval in this_array[map_field]]
    # need to know the maximum size of the subarray
    maxlen = numpy.array([x.size for x in expanded_info]).max()
    # convert to LSCArray
    if expand_field_name is not None:
        expanded_info = LSCArray.from_records(expanded_info,
            dtype=[(expand_field_name, other_dtdescr, maxlen)]).flatten()
    else:
        expanded_info = LSCArray.from_records(expanded_info,
            dtype=[(name, dt, maxlen) for (name,dt) in other_dtdescr])
        # if the map fields are the same name, remove from the expanded info
        # to avoid name collisions
        if map_field == other_map_field:
            expanded_info = expanded_info.without_fields(map_field)
    # add to this_array
    new_array = new_array.add_fields(expanded_info)
    # copy over any new methods
    if copy_methods:
        new_array = copy_attributes(new_array, other_array)
    return new_array

def file_checksum(filename, hasher=hashlib.sha256, buffersize=65536):
    """Provides a checksum of the given file. Modified from:
    <http://stackoverflow.com/a/3431835/1366472>.

    Parameters
    ----------
    filename : string
        Name of the file to checksum.
    hasher : {hashlib.sha256 | hashlib algorithm}
        The algorithm used for computing the hash.
    buffersize : {65536 | int}
        The number of bytes to read in from the file at a time.

    Returns
    -------
    checksum : hex
        The checksum of the file.
    """
    hasher = hasher()
    f = open(filename, 'rb')
    buf = f.read(buffersize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = f.read(buffersize)
    f.close()
    return hasher.hexdigest()

def get_all_field_names(dtype):
    """Given a numpy dtype, returns a list of all the field and subfield names
    present. Subfield names use `.` notation; e.g., if field `foo` has subfield
    `bar`, will return `foo.bar`.

    Parameters
    ----------
    dtype : numpy.dtype
        The dtype to get the field names from.

    Returns
    -------
    names : list
        The list of all of the field names.
    """
    names = []
    for name in dtype.names:
        if len(dtype[name]) > 0:
            subnames = ['%s.%s' %(name, subname) \
                for subname in get_all_field_names(dtype[name])]
            names.extend(subnames)
        else:
            names.append(name)
    return names

#
#   Argument syntax parsing
#
# this parser will pull out sufields as separate identifiers from their parent
# field; e.g., foo.bar --> ['foo', 'bar']
_pyparser = re.compile(r'(?P<identifier>[\w_][\w\d_]*)')
# this parser treats subfields as one identifier with their parent field;
# e.g., foo.bar --> ['foo.bar']
_fieldparser = re.compile(r'(?P<identifier>[\w_][.\w\d_]*)')
def get_vars_from_arg(arg):
    """Given a python string, gets the names of any identifiers use in it.
    For example, if ``arg = '3*narf/foo.bar'``, this will return
    ``set(['narf', 'foo', 'bar'])``.
    """
    return set(_pyparser.findall(arg))

def get_fields_from_arg(arg):
    """Given a python string, gets LSCArray field names used in it. This
    differs from get_vars_from_arg in that any identifier with a '.' in it
    will be treated as one identifier. For example, if
    ``arg = '3*narf/foo.bar'``, this will return ``set(['narf', 'foo.bar'])``.
    """
    return set(_fieldparser.findall(arg))

#
# =============================================================================
#
#                           Base LSCArray definitions
#
# =============================================================================
#
class LSCArray(numpy.recarray):
    """
    Subclass of numpy.recarray that adds additional functionality.

    Initialization is done the same way as numpy.recarray, with the addition
    that a "name" attribute can be passed to name the output array. When you
    initialize an array it creates a new zeroed array. This is similar to
    numpy.recarray, except that ``numpy.recarray(shape)`` will create an empty
    array, whereas here the default is to zero all of the elements (see
    ``default_zero`` for definition of zero for different data types). If you
    prefer an empty array, set ``zero=False`` when initializing.
    
    You cannot pass an array or sequence as input as you do with numpy.array.
    To initialize an LSCArray from an already existing arrays, use the
    ``LSCArray.from_arrays`` class method. To initialize from a list of
    tuples, use ``LSCArray.from_records``. See the docstring for those methods
    for details. For more information on initalizing an empty array, see
    ``numpy.recarray`` help.

    Some additional features:

    * **Arbitrary functions**:
    You can retrive functions on fields in the same manner that you access
    individual fields. For example, if you have an LSCArray ``x`` with fields
    ``a`` and ``b``, you can access each field with ``x['a'], x['b']``.
    You can also do ``x['a*b/(a+b)**2.']``, ``x[cos(a)*sin(b)]``, etc. Boolean
    operations are also possible, e.g., ``x['(a < 3) & (b < 2)']``. Syntax
    for functions is python, and any numpy ufunc can be used to operate
    on the fields. Note that while fields may be accessed as
    attributes (e.g, field ``a`` can be accessed via ``x['a']`` or ``x.a``),
    functions on multiple fields may not (``x.a+b`` does not work, for obvious
    reasons).

    * **Subfields and '.' indexing**:
    Structured arrays, which are the base class for recarrays and, by
    inheritance, LSCArrays, allows for fields to themselves have fields. For
    example, an array ``x`` may have fields ``a`` and ``b``, with ``b`` having
    subfields ``c`` and ``d``. You can access subfields using other index
    notation or attribute notation. So, the subfields ``d`` may be retrieved
    via ``x['b']['d']``, ``x.b.d``, ``x['b'].d`` or ``x['b.d']``. Likewise,
    functions can be carried out on the subfields, as they can on fields. If
    ``d`` is a float field, we could get the log of it via ``x['log(b.d)']``.
    There is no limit to the number of subfields. So, ``c`` could also have
    subfield ``c0``, which would be accessed via ``x.c.c0``, or any of the
    other methods.

    .. warning::
        Record arrays also allow you to set values of a field using attribute
        notation. However, this can lead to unexpected results if you
        accidently misspell the attribute. For example, if ``x`` has field
        ``foo``, and you misspell this when setting, e.g., you try to do
        ``x.fooo = numpy.arange(x.size)``, ``foo`` will not be set, nor will
        you get an error. Instead, the attribute ``fooo`` will be added to
        ``x``. If you tried to do this using index notation, however ---
        ``x['fooo'] = numpy.arange(x.size)`` --- you will
        get an ``AttributeError`` as you might expect. For this reason, it is
        recommended that you always use index notation when *setting* values;
        you can use either index or attribute notation when *retrieving*
        values.

    * **Properties and methods as fields**:
    If a propety or instance method is defined for a class that inherits from
    LSCArray, those can be accessed in the same way as fields are. For example,
    define ``Foo`` as:

    .. code-block:: python

        class Foo(LSCArray):
            @property
            def bar(self):
                return self['a']**2.

            def narf(self, y):
                return self['a'] + y

    Then if we have an instance: ``foo = Foo(100, dtype=[('a', float)])``.
    The ``bar`` and ``narf`` attributes may be accessed via field notation:
    ``foo.bar``, ``foo['bar']``, ``foo.narf(10)`` and ``foo['narf(10)']``.

    * **Add/drop fields**:
    Fields may be added or dropped to an already intialized array using
    ``add_fields`` and ``with[out]_fields``; see the documentation for
    those functions for details.


    * **Appending arrays**:
    Two instances of LSCArrays may be appended together if they have the
    same fields using the append instance method. If the array has an id
    field that you want incremented with the append to prevent collisions,
    you can indicate this in the append. For example, say you have array ``x``
    and ``y`` both with fields ``event_id`` that go from ``0-1000``. To append
    ``y`` to ``x``, while keeping ``event_id`` unique:
    
    .. code-block:: python

        z = x.append(y, remap_ids='event_id')

    This will result in the elements from ``y`` having event ids 1001-2001
    in ``z``. A dictionary mapping the new ids to original is stored in the
    ``id_maps`` attribute. See ``append`` for more details.

    * **Joining arrays**
    You can join two arrays using a common field with the ``join`` method.
    One-to-one, one-to-many, many-to-one, and many-to-many joins are supported.
    For example, if ``x`` has fields ``event_id, a`` and ``y`` has fields
    ``event_id, b``, then ``z = x.join(y, 'event_id')`` will have fields
    ``event_id, a, b``, with ``b`` set such that ``x.event_id == y.event_id``.
    In addition, all properties and methods of ``y`` that are not in ``x``
    will be copied to ``z``. For example, if ``y`` had property ``foo`` that
    operated on field ``b``, ``z`` will inherit that property. See ``join``
    for details.


    * **Lookup tables**:
    A lookup function is provided that allows you to quickly get all rows in
    the array for which a paricular field matches a particular value, e.g.,
    ``x.lookup('a', 10.)`` will return all rows in ``x`` for which ``x['a'] ==
    10.``.  This is done by building an internal dictionary using the
    requested field as a key the first time it is requested. Since this relies
    on the order of the results, the internal lookup table is not passed to
    views or copies of the array, and it cleared is whenever an in-place sort
    is carried out.  The lookup table does increase memory overhead. Also, if
    you change a value in the field that is used as key after the lookup table
    is created, you will get spurious results. For these reasons, a
    clear_lookup method is also provided. See ``lookup`` and ``clear_lookup``
    for details.

    Parameters
    ----------
    shape : {int | tuple}
        The shape of the new array.
    name : {None | str}
        Optional, what to name the new array. The array's ``name`` attribute
        is set to this.

    For details on other keyword arguments, see ``numpy.recarray`` help.

    Attributes
    ----------
    name : str
        Instance attribute. The name of the array.

    Examples
    --------
    .. note:: For some predefined arrays with default fields, see the other
        array classes defined below. For utilities that make loading arrays
        from data sources easier, see ``lscarray_utils.py`` in the various
        sub-directories of ``io``.

    Create an empty array with four rows and two fields called `foo` and
    `bar`, both of which are floats:

    >>> x = LSCArray(4, dtype=[('foo', float), ('bar', float)])

    Set/retrieve a fields using index or attribute syntax:

    >>> x['foo'] = [1.,2.,3.,4.]
    >>> x['bar'] = [5.,6.,7.,8.]
    >>> x
    LSCArray([(1.0, 5.0), (2.0, 6.0), (3.0, 7.0), (4.0, 8.0)], 
          dtype=[('foo', '<f8'), ('bar', '<f8')])
    >>> x.foo
        array([ 1.,  2.,  3.,  4.])
    >>> x['bar']
        array([ 5.,  6.,  7.,  8.])
    
    Get the names of the fields:

    >>> x.fieldnames
        ('foo', 'bar')

    Rename the fields to `a` and `b`:

    >>> x.dtype.names = ['a', 'b']
    >>> x.fieldnames
        ('a', 'b')

    Retrieve a function of the fields as if it were a field:

    >>> x['sin(a/b)']
    array([ 0.19866933,  0.3271947 ,  0.41557185,  0.47942554])

    Create an array with subfields:

    >>> x = LSCArray(4, dtype=[('foo', [('cat', float), ('hat', int)]), ('bar', float)])
    >>> x.all_fieldnames
        ['foo.cat', 'foo.hat', 'bar']

    Load from a list of arrays (in this case, from an hdf5 file):

    >>> bankhdf = h5py.File('bank/H1L1-BANK2HDF-1117400416-928800.hdf')
    >>> bankhdf.keys()
        [u'mass1', u'mass2', u'spin1z', u'spin2z', u'template_hash']
    >>> templates = LSCArray.from_arrays(bankhdf.values(), names=bankhdf.keys())
    >>> templates.fieldnames
        ('mass1', 'mass2', 'spin1z', 'spin2z', 'template_hash')
    >>> templates.mass1
    array([ 1.71731389,  1.10231435,  2.99999857, ...,  1.67488706,
            1.00531888,  2.11106491], dtype=float32)

    Sort by a field without having to worry about also sorting the other
    fields:

    >>> templates[['mass1', 'mass2']]
    array([(1.7173138856887817, 1.2124452590942383),
           (1.1023143529891968, 1.0074082612991333),
           (2.9999985694885254, 1.0578444004058838), ...,
           (1.6748870611190796, 1.1758257150650024),
           (1.0053188800811768, 1.0020891427993774),
           (2.111064910888672, 1.0143394470214844)], 
          dtype=[('mass1', '<f4'), ('mass2', '<f4')])
    >>> templates.sort(order='mass1')
    >>> templates[['mass1', 'mass2']]
    array([(1.000025987625122, 1.0000133514404297),
           (1.0002814531326294, 1.0002814531326294),
           (1.0005437135696411, 1.0005437135696411), ...,
           (2.999999523162842, 1.371169090270996),
           (2.999999523162842, 1.4072519540786743), (3.0, 1.4617927074432373)], 
          dtype=[('mass1', '<f4'), ('mass2', '<f4')])

    Convert a LIGOLW xml table:

    >>> type(sim_table)
        glue.ligolw.lsctables.SimInspiralTable
    >>> sim_array = LSCArray.from_ligolw_table(sim_table)
    >>> sim_array.mass1
    array([ 2.27440691,  1.85058105,  1.61507106, ...,  2.0504961 ,
            2.33554196,  2.02732205], dtype=float32)
    >>> sim_array.waveform
    array([u'SpinTaylorT2', u'SpinTaylorT2', u'SpinTaylorT2', ...,
           u'SpinTaylorT2', u'SpinTaylorT2', u'SpinTaylorT2'], dtype=object)
    
    Only view a few of the fields:

    >>> sim_array.with_fields(['simulation_id', 'mass1', 'mass2'])
    LSCArray([(0, 2.274406909942627, 2.6340370178222656),
           (1, 1.8505810499191284, 2.8336880207061768),
           (2, 1.6150710582733154, 2.2336490154266357), ...,
           (11607, 2.0504961013793945, 2.6019821166992188),
           (11608, 2.3355419635772705, 1.2164380550384521),
           (11609, 2.0273220539093018, 2.2453839778900146)], 
          dtype={'names':['simulation_id','mass1','mass2'], 'formats':['<i8','<f4','<f4'], 'offsets':[200,236,240], 'itemsize':244})

    ...or just retrieve a few of the fields to begin with:

    >>> sim_array = LSCArray.from_ligolw_table(sim_table, columns=['simulation_id', 'mass1', 'mass2'])
    >>> sim_array
    LSCArray([(0, 2.274406909942627, 2.6340370178222656),
           (1, 1.8505810499191284, 2.8336880207061768),
           (2, 1.6150710582733154, 2.2336490154266357), ...,
           (11607, 2.0504961013793945, 2.6019821166992188),
           (11608, 2.3355419635772705, 1.2164380550384521),
           (11609, 2.0273220539093018, 2.2453839778900146)], 
          dtype=[('simulation_id', '<i8'), ('mass1', '<f4'), ('mass2', '<f4')])

    Add a field to the array:

    >>> optimal_snrs = numpy.random.uniform(4.,40., size=len(sim_array))
    >>> sim_array = sim_array.add_fields(optimal_snrs, 'optimal_snrs')
    >>> sim_array.fieldnames
        ('simulation_id', 'mass1', 'mass2', 'optimal_snrs')

    Notes
    -----
    Input arrays with variable-length strings in one or more fields can be
    tricky to deal with. Numpy arrays are designed to use fixed-length
    datasets, so that quick memory access can be achieved. To deal with
    variable-length strings, there are two options: 1. set the data type to
    object, or 2. set the data type to a string with a fixed length larger
    than the longest string in the input array.
    
    The first option, using objects, essentially causes the array to store a
    pointer to the string.  This is the most flexible option, as it allows
    strings in the array to be updated to any length. However, operations on
    object fields are slower, as numpy cannot take advantage of its fast
    memory striding abilities (see `this question/answer on stackoverflow
    <http://stackoverflow.com/a/14639568/1366472>`_ for details). Also,
    numpy's support of object arrays is more limited.  In particular, prior
    to version 1.9.2, you cannot create a view of an array that changes the 
    dtype if the array has any fields that are objects, even if the view does
    not affect the object fields. (This has since been relaxed.)

    The second option, using strings of a fixed length, solves the issues
    with object fields. However, if you try to change one of the strings
    after the array is created, the string will be truncated at whatever
    string length is used. Additionally, if you choose too large of a string
    length, you can substantially increase the memory overhead for large
    arrays.

    This class offers the option to set strings designated as 'lstring' as
    objects, or as fixed-length strings. To toggle what it does use
    ``lscarrays.set_lstring_as_obj`` (see the docstring for that function for
    more details).
    """
    __persistent_attributes__ = ['name', 'source_files', 'id_maps']

    def __new__(cls, shape, name=None, zero=True, **kwargs):
        """Initializes a new empty array.
        """
        obj = super(LSCArray, cls).__new__(cls, shape, **kwargs).view(
            type=cls)
        obj.name = name
        obj.source_files = None
        obj.id_maps = None
        obj.__persistent_attributes__ = cls.__persistent_attributes__
        # zero out the array if desired
        if zero:
            default = default_empty(1, dtype=obj.dtype)
            obj[:] = default
        return obj


    def __array_finalize__(self, obj):
        """Default values are set here.
        """
        if obj is None:
            return
        # copy persisitent attributes
        [setattr(self, attr, getattr(obj, attr, None)) for attr in \
            self.__persistent_attributes__]
        # numpy has some issues with dtype field names that are unicode,
        # so we'll force them to strings here
        if obj.dtype.names is not None and \
                any([isinstance(name, unicode) for name in obj.dtype.names]):
            obj.dtype.names = map(str, obj.dtype.names)
        self.__lookuptable__ = {}


    def __copy_attributes__(self, other, default=None):
        """Copies the values of all of the attributes listed in
        `self.__persistent_attributes__` to other.
        """
        [setattr(other, attr, copy.deepcopy(getattr(self, attr, default))) \
            for attr in self.__persistent_attributes__]


    def __setitem__(self, item, values):
        """Wrap's recarray's setitem to allow attribute-like indexing when
        setting values.
        """
        try:
            return super(LSCArray, self).__setitem__(item, values)
        except ValueError:
            # we'll get a ValueError if a subarray is being referenced using
            # '.'; so we'll try to parse it out here
            fields = item.split('.')
            if len(fields) > 1:
                for field in fields[:-1]:
                    self = self[field]
                item = fields[-1]
            # now try again
            return super(LSCArray, self).__setitem__(item, values)


    def __getsubitem__(self, item):
        """Gets a subfield using `field.subfield` notation.
        """
        try:
            return super(LSCArray, self).__getitem__(item)
        except ValueError as err:
            subitems = item.split('.')
            if len(subitems) > 1:
                return super(LSCArray, self).__getitem__(subitems[0]
                    ).__getsubitem__('.'.join(subitems[1:]))
            else:
                raise ValueError(err.message)


    def __getitem__(self, item):
        """Wraps recarray's  `__getitem__` so that math functions on fields and
        attributes can be retrieved. Any function in numpy's library may be
        used.
        """
        try:
            return self.__getsubitem__(item)
        except ValueError:
            # arg isn't a simple argument of row, so we'll have to eval it
            # get the set of fields & attributes we will need
            #itemvars = get_fields_from_arg(item)
            itemvars = get_vars_from_arg(item)
            # pull out the fields: note, by getting the parent fields, we
            # also get the sub fields name
            item_dict = dict([ [fn,
                super(LSCArray, self).__getitem__(fn)] \
                for fn in set(self.fieldnames).intersection(itemvars)])
            # add any aliases
            item_dict.update({alias: item_dict[name] \
                for alias,name in self.aliases.items()
                    if name in item_dict})
            # pull out any attributes needed
            itemvars = get_fields_from_arg(item)
            item_dict.update({attr: getattr(self, attr) for attr in \
                set(dir(self)).intersection(itemvars)})
            # add numpy functions
            item_dict.update(numpy.__dict__)
            return eval(item, {"__builtins__": None}, item_dict)

    def addattr(self, attrname, value=None, persistent=True):
        """Adds an attribute to self. If persistent is True, the attribute will
        be made a persistent attribute. Persistent attributes are copied
        whenever a view or copy of this array is created. Otherwise, new views
        or copies of this will not have the attribute.
        """
        setattr(self, attrname, value)
        # add as persistent
        if persistent and attrname not in self.__persistent_attributes__:
            self.__persistent_attributes__.append(attrname)

    def add_methods(self, names, methods):
        """Adds the given method(s) as instance method(s) of self. The
        method(s) must take `self` as a first argument.
        """
        if isinstance(names, str) or isinstance(names, unicode):
            names = [names]
            methods = [methods]
        for name,method in zip(names, methods):
            setattr(self, name, types.MethodType(method, self))
        
    def add_properties(self, names, methods):
        """Returns a view of self with the given methods added as properties.

        From: <http://stackoverflow.com/a/2954373/1366472>.
        """
        cls = type(self)
        if not hasattr(cls, '__perinstance'):
            cls = type(cls.__name__, (cls,), {})
            cls.__perinstance = True
        if isinstance(names, str) or isinstance(names, unicode):
            names = [names]
            methods = [methods]
        for name,method in zip(names, methods):
            setattr(cls, name, property(method))
        return self.view(type=cls)

    @classmethod
    def from_arrays(cls, arrays, name=None, **kwargs):
        """Creates a new instance of self from the given (list of) array(s).
        This is done by calling numpy.rec.fromarrays on the given arrays with
        the given kwargs. The type of the returned array is cast to this
        class, and the name (if provided) is set.

        Parameters
        ----------
        arrays : (list of) numpy array(s)
            A list of the arrays to create the LSCArray from.
        name : {None|str}
            What the output array should be named.

        For other keyword parameters, see the numpy.rec.fromarrays help.

        Returns
        -------
        array : instance of this class
            An array that is an instance of this class in which the field
            data is from the given array(s).
        """
        obj = numpy.rec.fromarrays(arrays, **kwargs).view(type=cls)
        obj.name = name
        return obj

    @classmethod
    def from_records(cls, records, name=None, **kwargs):
        """
        Creates a new instance of self from the given (list of) record(s). A
        "record" is a tuple in which each element is the value of one field
        in the resulting record array. This is done by calling
        `numpy.rec.fromrecords` on the given records with the given kwargs.
        The type of the returned array is cast to this class, and the name
        (if provided) is set.

        Parameters
        ----------
        records : (list of) tuple(s)
            A list of the tuples to create the LSCArray from.
        name : {None|str}
            What the output array should be named.

        For other keyword parameters, see the `numpy.rec.fromrecords` help.

        Returns
        -------
        array : instance of this class
            An array that is an instance of this class in which the field
            data is from the given record(s).
        """
        obj = numpy.rec.fromrecords(records, **kwargs).view(
            type=cls)
        obj.name = name
        return obj

    @classmethod
    def from_ligolw_table(cls, table, columns=None, cast_to_dtypes=None):
        """Converts the given ligolw table into an LSCArray. The `tableName`
        attribute is copied to the array's `name`.

        Parameters
        ----------
        table : LIGOLw table instance
            The table to convert.
        columns : {None|list}
            Optionally specify a list of columns to retrieve. All of the
            columns must be in the table's validcolumns attribute. If None
            provided, all the columns in the table will be converted.
        dtype : {None | dict}
            Override the columns' dtypes using the given dictionary. The
            dictionary should be keyed by the column names, with the values
            a tuple that can be understood by numpy.dtype. For example, to
            cast a ligolw column called "foo" to a field called "bar" with
            type float, cast_to_dtypes would be: ``{"foo": ("bar", float)}``.

        Returns
        -------
        array : LSCArray
            The input table as an LSCArray.
        """
        name = table.tableName.split(':')[0]
        if columns is None:
            # get all the columns
            columns = table.validcolumns
        else:
            # note: this will raise a KeyError if one or more columns is
            # not in the table's validcolumns
            columns = {col: table.validcolumns[col] \
                for col in columns}
        if cast_to_dtypes is not None:
            dtype = [cast_to_dtypes[col] for col in columns]
        else:
            dtype = columns.items()
        # get the values
        if _default_types_status['ilwd_as_int']:
            input_array = [tuple(
                        getattr(row, col) if dt != 'ilwd:char' \
                        else int(getattr(row, col)) \
                    for col,dt in columns.items()) \
                for row in table]
        else:
            input_array = [tuple(getattr(row, col) for col in columns) \
                for row in table]
        # return the values as an instance of cls
        return cls.from_records(input_array, dtype=dtype,
            name=name)

    @property
    def fieldnames(self):
        """Returns a tuple listing the field names in self. Equivalent to
        `array.dtype.names`, where `array` is self.
        """
        return self.dtype.names

    @property
    def all_fieldnames(self):
        """Returns a list of all of the field names in self, including
        subfields. Subfields are named `field.subfield`.
        """
        return get_all_field_names(self.dtype)

    @property
    def aliases(self):
        """Returns a dictionary of the aliases, or "titles", of the field names
        in self. An alias can be specified by passing a tuple in the name
        part of the dtype. For example, if an array is created with
        ``dtype=[(('foo', 'bar'), float)]``, the array will have a field
        called `bar` that has alias `foo` that can be accessed using
        either `arr['foo']` or `arr['bar']`. Note that the first string
        in the dtype is the alias, the second the name. This function returns
        a dictionary in which the aliases are the keys and the names are the
        values. Only fields that have aliases are returned.
        """
        return dict(c[0] for c in self.dtype.descr if isinstance(c[0], tuple))

    def add_source_file(self, filename, id_field=None, id_ranges=None,
            include_checksum=False, force=False):
        """Adds a source file to self.source_files. If id_field is provided,
        will also include information about the ids that the source file
        is valid for. This information is stored in a dictionary as:
        ``self.source_files[filename][id_field] = id_ranges``

        If include_checksum is True, a checksum of the file will be computed
        and added to:
        ``self.source_files[filename]['_checksum'] = checksum``.

        If id_field is None and include_checksum is False,
        ``self.source_files[filename]`` will just be an empty dictionary.

        Parameters
        ----------
        filename : string
            The name of the file.
        id_field : {None | string}
            The name of an id field in self for which this file is valid. To
            refer to a subarray field, use "``.``", e.g.,
            ``recovered.event_id``.
        id_ranges : {None | array}
            The id values for which the file is valid. Can either be an array
            of the individual id values, or just the minimum and maximum
            values, to specify a range. If provided, must also provide
            id_field. If id_field provided and id_ranges is None, the ranges
            will be ``[self[id_field].min(), self[id_field].max()]``.
        include_checksum : bool
            If True, will calculate a checksum for the file and add it to
            the checksum fields. Default is False.
        force : bool
            If the filename is already in self.source_files, overwrite it.
            Otherwise, a ValueError will be raised. Default is False.
        """
        if self.source_files is None:
            self.source_files = {}
        if id_field is None and id_ranges is not None:
            raise ValueError("if providing an id_ranges must " + \
                "also provide id_field")
        if filename in self.source_files and not force:
            raise ValueError("filename %s already in source_files" %(filename))
        self.source_files[filename] = {}
        if id_field is not None:
            if id_ranges is None:
                id_ranges = numpy.array([self[id_field].min(),
                    self[id_field].max()])
            self.source_files[filename][id_field] = id_ranges
        if include_checksum:
            self.source_files[filename]['_checksum'] = file_checksum(filename)


    def add_fields(self, arrays, names=None, assubarray=False):
        """
        Adds the given arrays as new fields to self.  Returns a new instance
        with the new fields added. Note: this array does not change; the
        returned array is a new copy. The internal lookup table of the new
        array will also be empty.

        Parameters
        ----------
        arrays : (list of) numpy array(s)
            The arrays to add. If adding multiple arrays, must be a list;
            if adding a single array, can just be that array.
        names : (list of) strings
            Optional, the name(s) of the new fields in the output array. If
            adding multiple fields, must be a list of strings with the same
            length as the list of arrays. If None provided, names used will
            be the same as the name of the datatype in the given arrays.
            If the datatype has no name, the new field will be ``'fi'`` where
            i is the index of the array in arrays.
        assubarray : bool
            Add the list of arrays as a single subarray field. If True, and
            names provided, names should be a string or a length-1 sequence.
            Default is False, in which case each array will be added as a
            separate field.

        Returns
        -------
        new_array : new instance of this array
            A copy of this array with the desired fields added.
        """
        newself = add_fields(self, arrays, names=names, assubarray=assubarray)
        self.__copy_attributes__(newself)
        return newself


    def with_fields(self, names, copy=False):
        """
        Get a view/copy of this array with only the specified fields.

        Parameters
        ----------
        names : {strings|list of strings}
            List of field names for the returned array; may be either a
            single name or a list of names.
        copy : bool, optional
            If True, will return a copy of the array rather than a view.
            Default is False.

        Returns
        -------
        lscarray : new instance of this array
            The view or copy of the array with only the specified fields.
        """
        newself = get_fields(self, names, copy=copy)  
        # copy relevant attributes
        self.__copy_attributes__(newself)
        return newself


    def without_fields(self, names, copy=False):
        """
        Get a view/copy of this array without the specified fields.

        Parameters
        ----------
        names : {strings|list of strings}
            List of field names for the returned array; may be either a
            single name or a list of names.
        copy : bool, optional
            If True, will return a copy of the array rather than a view.
            Default is False.

        Returns
        -------
        lscarray : new instance of this array
            The view or copy of the array without the specified fields.
        """
        # check if only given a single name
        if isinstance(names, str) or isinstance(names, unicode):
            names = [names]
        # cycle over self's fields, excluding any names in names
        keep_names = [name for name in self.dtype.names if name not in names]
        newself = get_fields(self, keep_names, copy=copy)  
        # copy relevant attributes
        self.__copy_attributes__(newself)
        return newself


    def lookup(self, field, value, default=KeyError):
        """
        Returns the elements in self for which the given field(s) matches the
        given value(s). If this is the first time that the field(s) has been
        requested, an internal dictionary is built first, then the elements
        returned. If the value(s) cannot be found, a KeyError is raised, or
        a default is returned if provided. To lookup multiple fields, pass
        the field names and the corresponding values as tuples.
        
        .. note::
            Every time a lookup on a new field is done, an array of indices
            is created that maps every unique value in the field to the
            indices in self that match that value. Since this mapping will
            change if the array is sorted, or if a slice of self is created,
            this look up table is not passed to new views of self. To reduce
            memory overhead, run clear_lookup on a field if you will no
            longer need to look up self using that field. See clear_look up
            for details.

        .. warning::
            If you change the value of an item that is used as a look
            up item, the internal lookup dictionary will not give correct
            values. Always run clear_lookup if you change the value of
            a field that you intend to use as a lookup key.

        Parameters
        ----------
        field : (tuple of) string(s)
            The name(s) of the field to look up.
        value : (tuple of) value(s)
            The value(s) in the fields to get.
        default : {KeyError | value}
            Optionally specify a value to return is the value is not found
            in self's field. Otherwise, a KeyError is raised.

        Returns
        -------
        matching : type(self)
            The rows in self that match the requested value.
        """
        try:
            return self[self.__lookuptable__[field][value]]
        except KeyError:
            # build the look up for this field
            if field not in self.__lookuptable__:
                # if field is a joint field, convert to list
                if isinstance(field, tuple):
                    colvals = self[list(field)]
                else:
                    colvals = self[field]
                self.__lookuptable__[field] = build_lookup_table(colvals)
            # try again
            try:
                return self[self.__lookuptable__[field][value]]
            except KeyError as err:
                if inspect.isclass(default) and issubclass(default, Exception):
                    raise default(err.message)
                else:
                    return default


    def clear_lookup(self, field=None):
        """
        Clears the internal lookup table used to provide fast lookup of
        the given field. If no field is specified, the entire table will
        be deleted. Run this if you will no longer need to look up self
        using a particular field (or, if field is None, any look ups at
        all). Note: if you clear the lookup table of a field, then try to
        lookup that field again, a new lookup table will be created.
        """
        if field is None:
            self.__lookuptable__.clear()
        else:
            self.__lookuptable__[field].clear()
            del self.__lookuptable__[field]


    def sort(self, *args, **kwargs):
        """
        Clears self's lookup table before sorting. See numpy.ndarray.sort for
        help.
        """
        self.clear_lookup()
        super(LSCArray, self).sort(*args, **kwargs)


    def join(self, other, map_field, expand_field_name=None,
            other_map_field=None, get_fields=None, map_indices=None):
        """
        Join another array to this array such that:
        ``self[map_field]`` == ``other[other_map_field]``. The fields from
        ``other`` are added as a sub-array to self with field name
        ``expand_field_name``. Any attributes in ``other`` that are not
        attributes of this array are copied to the output array. Therefore, if
        ``other`` array has properties/methods that use the fields
        that are added to this array, the properties/methods will continue to
        work on the output array.

        Parameters
        ----------
        other : LSCArray or similar
            The array that information is retrieved from. Must have
            ``lookup`` and ``with_fields`` methods.
        map_field : string
            The name of the field in self to use for mapping.
        expand_field_name : string
            The name of the field that will be added to self. The information
            from ``other`` will be contained as a subarray under this field.
        other_map_field : {None|string}
            The name of the field in ``other`` to use for mapping. If None
            provided, ``map_field`` will be used.
        get_fields : {None | (list of) strings}
            Optionally specify what fields to retrieve from ``other_array``.
            If None provided, will get all the fields in ``other_array``.
        map_indices : {None | array of ints}
            If provided, will only map rows in this array that have indices
            in the given array of indices. Any rows that are skipped will
            have a zeroed element in the expand field of the returned array.
            If None (the default), all rows in this array are mapped.

        Returns
        -------
        new_array : type(self)
            A copy of this array with the mapped information added to
            ``expand_field_name``.
        """
        return join_arrays(self, other, map_field,
            expand_field_name=expand_field_name,
            other_map_field=other_map_field, get_fields=get_fields,
            map_indices=map_indices)

    def append(self, other, remap_ids=None):
        """
        Appends another array to this array, returning a copy of the combined
        array. If a list of remap_ids is provided, those ids will be updated
        in the combined array to prevent collisions. The persistent attributes
        of the new array (e.g., name) will be retrieved from this array.

        Parameters
        ----------
        other : array
            The array to append. Must have the same dtype as this array.
        remap_ids : {None | (list of) string(s)}
            A (list of) field name(s) to update. For every field name provided,
            the values in other array that are added to self will be increased
            by the maximum value of the field in self. The ids that are updated
            and the padding by which they were updated is stored to the new
            array's ``id_maps`` attribute.

        Returns
        -------
        combined : array
            A copy of self with other appended.

        Examples
        --------
        Append two sim_inspiral arrays together:

        >>> sims1.name, sims1.size, sims1.simulation_id, sims1.process_id, sims1.source_files, sims1.id_maps
        ('sim_inspiral',
         11610,
         array([    0,     1,     2, ..., 11607, 11608, 11609]),
         array([0, 0, 0, ..., 0, 0, 0]),
         {'inj_files/HL-INJECTIONS_BNS0INJ-1117400416-928800.xml': {'_checksum': '2a7a0f88f88476cd2ae85d03f9772a8da0d49ee734352cee79cfd64943d99bbe',
           'simulation_id': array([    0, 11609])}},
         None)
        >>> sims2.name, sims2.size, sims2.simulation_id, sims2.process_id, sims2.source_files, sims2.id_maps
        ('other_sim_inspiral',
         11610,
         array([    0,     1,     2, ..., 11607, 11608, 11609]),
         array([0, 0, 0, ..., 0, 0, 0]),
         {'inj_files/HL-INJECTIONS_BNS1INJ-1117400416-928800.xml': {'_checksum': 'cab57e614de403fd58affa3994f5920c56280c2c283871e8b730717338629fa9',
           'simulation_id': array([    0, 11609])}},
         None)
        >>> sims = sims1.append(sims2, ['simulation_id', 'process_id'])
        >>> sims.name, sims.size, sims.simulation_id, sims.process_id, sims.source_files, sims.id_maps
        ('sim_inspiral',
         23220,
         array([    0,     1,     2, ..., 23217, 23218, 23219]),
         array([0, 0, 0, ..., 1, 1, 1]),
         {'inj_files/HL-INJECTIONS_BNS0INJ-1117400416-928800.xml': {'_checksum': '2a7a0f88f88476cd2ae85d03f9772a8da0d49ee734352cee79cfd64943d99bbe',
           'simulation_id': array([    0, 11609])},
          'inj_files/HL-INJECTIONS_BNS1INJ-1117400416-928800.xml': {'_checksum': 'cab57e614de403fd58affa3994f5920c56280c2c283871e8b730717338629fa9',
           'simulation_id': array([11610, 23219])}},
         {'process_id': [(1, 1, 1)], 'simulation_id': [(11610, 23219, 11610)]})
        """
        newarr = numpy.append(self, other).view(type=type(self))
        # copy the persistent attributes from self
        self.__copy_attributes__(newarr)
        # remap ids if desired
        if remap_ids is not None:
            if isinstance(remap_ids, str) or isinstance(remap_ids, unicode):
                remap_ids = [remap_ids]
            if newarr.id_maps is None:
                newarr.id_maps = {}
            for idfield in remap_ids:
                # get the maximum value from self and use that as a padding
                pad = self[idfield].max() + 1
                newarr[idfield][self.size:] += pad
                # store in id map: we store the range of values of the new
                # id values, and the pad used; this can be used to get the
                # value of the original ids
                new_ids = newarr[idfield][self.size:]
                idmap = (new_ids.min(), new_ids.max(), pad)
                try:
                    newarr.id_maps[idfield].append(idmap)
                except KeyError:
                    newarr.id_maps[idfield] = []
                    newarr.id_maps[idfield].append(idmap)
        if hasattr(other, 'source_files') and other.source_files is not None:
            for filename in other.source_files:
                source_info = copy.deepcopy(other.source_files[filename])
                newarr.source_files[filename] = source_info 
                # update any ids that were remapped
                if remap_ids is not None:
                    update_fields = [id_field for id_field in source_info \
                        if id_field in remap_ids]
                    for id_field in update_fields:
                        newarr.source_files[filename][id_field] = \
                              source_info[id_field] + \
                              newarr.id_maps[id_field][-1][-1]
        return newarr


def aliases_from_fields(fields):
    """
    Given a dictionary of fields, will return a dictionary mapping the aliases
    to the names.
    """
    return dict(c for c in fields if isinstance(c, tuple))


def fields_from_names(fields, names=None):
    """
    Given a dictionary of fields and a list of names, will return a dictionary
    consisting of the fields specified by names. Names can be either the names
    of fields, or their aliases.
    """
    if names is None:
        return fields
    if isinstance(names, str) or isinstance(names, unicode):
        names = [names]
    aliases_to_names = aliases_from_fields(fields)
    names_to_aliases = dict(zip(aliases_to_names.values(),
        aliases_to_names.keys()))
    outfields = {}
    for name in names:
        try:
            outfields[name] = fields[name]
        except KeyError:
            if name in aliases_to_names:
                key = (name, aliases_to_names[name])
            elif name in names_to_aliases:
                key = (names_to_aliases[name], name)
            else:
                raise KeyError('default fields has no field %s' % name)
            outfields[key] = fields[key]
    return outfields

class _LSCArrayWithDefaults(LSCArray):
    """
    Subclasses LSCArray, adding class method ``default_fields`` and class
    attribute ``default_name``. If no name is provided when initialized, the
    ``default_name`` will be used. Likewise, if no dtype is provided when
    initalized, the default fields will be used. If names are provided, they
    must be names or aliases that are in default fields, else a KeyError is
    raised. Non-default fields can be created by specifying dtype directly. 

    The default ``default_name`` is None and ``default_fields`` returns an
    empty dictionary. This class is mostly meant to be subclassed by other
    classes, so they can add their own defaults.
    """
    default_name = None

    @classmethod
    def default_fields(cls, **kwargs):
        """
        The default fields. This function should be overridden by subclasses
        to return a dictionary of the desired default fields. Allows for
        key word arguments to be passed to it, for classes that need to be
        able to alter properties of some of the default fields.
        """
        return {}

    def __new__(cls, shape, name=None, field_args={}, **kwargs):
        """
        Makes use of cls.default_fields and cls.default_name.
        """
        if 'names' in kwargs and 'dtype' in kwargs:
            raise ValueError("Please provide names or dtype, not both")
        fields = cls.default_fields(**field_args)
        if 'names' in kwargs:
            names = kwargs.pop('names')
            if isinstance(names, str) or isinstance(names, unicode):
                names = [names]
            kwargs['dtype'] = fields_from_names(fields, names).items()
        if 'dtype' not in kwargs:
            kwargs['dtype'] = fields.items()
        return super(_LSCArrayWithDefaults, cls).__new__(cls, shape,
            name=None, **kwargs)

    def add_default_fields(self, names, **kwargs):
        """
        Adds one or more empty default fields to self.

        Parameters
        ----------
        names : (list of) string(s)
            The names of the fields to add. Must be a field in self's default
            fields.
        
        Other keyword args are any arguments passed to self's default fields.

        Returns
        -------
        new array : instance of this array
            A copy of this array with the field added.
        """
        if isinstance(names, str) or isinstance(names, unicode):
            names = [names]
        fields = fields_from_names(self.default_fields(**kwargs), names)
        arrays = []
        names = []
        for name,dt in fields.items():
            arrays.append(default_empty(self.size, dtype=[(name, dt)])) 
            names.append(name)
        return self.add_fields(arrays, names)


#
# =============================================================================
#
#                           Default LSCArrays
#
# =============================================================================
#

class Waveform(_LSCArrayWithDefaults):
    """
    Subclasses LSCArrayWithDefaults, with default name `waveform`. The
    `_static_fields` define fields needed to describe a CBC waveform in the
    radiation frame. These are returned by `default_fields`.

    Also adds various common functions decorated as properties, such as
    ``mtotal = mass1+mass2``.
    """
    default_name = 'waveform'

    # _static_fields are fields that are automatically inherited by
    # subclasses of this class
    _static_fields = {
        # intrinsic parameters
        ('m1', 'mass1'): float,
        ('m2', 'mass2'): float,
        ('s1x', 'spin1x'): float,
        ('s1y', 'spin1y'): float,
        ('s1z', 'spin1z'): float,
        ('s2x', 'spin2x'): float,
        ('s2y', 'spin2y'): float,
        ('s2z', 'spin2z'): float,
        'lambda1': float,
        'lambda2': float,
        'quadparam1': float,
        'quadparam2': float,
        'eccentricity': float,
        'argument_periapsis': float,
        # extrinsic parameters
        'phi_ref': float,
        ('inc', 'inclination'): float,
        # waveform parameters
        'sample_rate': int,
        'segment_length': int,
        'f_min': float,
        'f_ref': float,
        'f_max': float,
        'duration': float,
        'amp_order': int,
        'phase_order': int,
        'spin_order': int,
        'tidal_order': int,
        ('apprx', 'approximant'): 'lstring',
        'taper': 'lstring',
        'frame_axis': 'lstring',
        'modes_choice': 'lstring',
        }

    @classmethod
    def default_fields(cls):
        return cls._static_fields 

    # some other derived parameters
    @property
    def mtotal(self):
        return self.mass1 + self.mass2

    @property
    def mtotal_s(self):
        return lal.MTSUN_SI*self.mtotal

    @property
    def q(self):
        return self.mass1 / self.mass2

    @property
    def eta(self):
        return self.mass1*self.mass2 / self.mtotal**2.

    @property
    def mchirp(self):
        return self.eta**(3./5)*self.mtotal

    def tau0(self, f0=None):
        """
        Returns tau0. If f0 is not specified, uses self.f_min.
        """
        if f0 is None:
            f0 = self.f_min
        return (5./(256 * numpy.pi * f0 * self.eta)) * \
            (numpy.pi * self.mtotal_s * f0)**(-5./3.)
   
    def v0(self, f0=None):
        """
        Returns the velocity at f0, as a fraction of c. If f0 is not
        specified, uses self.f_min.
        """
        if f0 is None:
            f0 = self.f_min
        return (2*numpy.pi* f0 * self.mtotal_s)**(1./3)

    @property
    def s1(self):
        return numpy.array([self.spin1x, self.spin1y, self.spin1z]).T

    @property
    def s2(self):
        return numpy.array([self.spin2x, self.spin2y, self.spin2z]).T

    @property
    def s1mag(self):
        return numpy.sqrt((self.s1**2).sum(axis=1))

    @property
    def s2mag(self):
        return numpy.sqrt((self.s2**2).sum(axis=1))


class TmpltInspiral(Waveform):
    """
    Subclasses Waveform, with default name `tmplt_inspiral`. Adds fields
    `template_id` and `process_id` as ids and an `ifo` field for listing what
    ifo(s) the template is filtered in. The `ifo` field may have more than one
    ifo listed. The maxium length of this field is specified at creating using
    the `nifos` key word.
   
    **Default fields:**
    %s

    Examples
    --------
    Create a TmpltInspiral array from an hdf bank file:

    >>> bankhdf = h5py.File('H1L1-BANK2HDF-1117400416-928800.hdf', 'r')
    >>> templates = TmpltInspiral.from_arrays(bankhdf.values(), names=bankhdf.keys())
    >>> templates = templates.add_fields(numpy.arange(len(templates)), names='template_id')
    >>> templates.fieldnames
        ('mass1', 'mass2', 'spin1z', 'spin2z', 'template_hash', 'template_id')
    >>> templates.mass1
    array([ 1.71731389,  1.10231435,  2.99999857, ...,  1.67488706,
            1.00531888,  2.11106491], dtype=float32)
    >>> templates[['mass1', 'mass2']]
    array([(1.7173138856887817, 1.2124452590942383),
           (1.1023143529891968, 1.0074082612991333),
           (2.9999985694885254, 1.0578444004058838), ...,
           (1.6748870611190796, 1.1758257150650024),
           (1.0053188800811768, 1.0020891427993774),
           (2.111064910888672, 1.0143394470214844)], 
          dtype=[('mass1', '<f4'), ('mass2', '<f4')])
    """
    default_name = 'tmplt_inspiral'

    @classmethod
    def default_fields(cls, nifos=1):
        """
        Admits argument nifos, which is used to set the size of the ifo
        sub-array.
        """
        fields = {
            'template_id': int,
            'process_id': int,
            # Note: ifo is a subarray with length set by nifos
            'ifo': ('S2', nifos),
            }
        # Note: cls._static_fields is Waveform's static fields
        return dict(cls._static_fields.items() + fields.items())

    def __new__(cls, shape, name=None, nifos=1, **kwargs):
        """
        Adds nifos to initialization.
        """
        return super(TmpltInspiral, cls).__new__(cls, shape, name=name, 
            field_args={'nifos': nifos}, **kwargs)


def end_time_from_float(end_time):
    """Given an end time as a float, returns the seconds and nanoseconds.
    """
    try:
        end_time_s = end_time.astype(int)
        end_time_ns = (end_time % 1 * 1e9).astype(int)
    except AttributeError:
        end_time_s = int(end_time)
        end_time_ns = (end_time % 1 * 1e9).astype(int)
    return end_time_s, end_time_ns


class SnglEvent(_LSCArrayWithDefaults):
    """
    Has default fields for storing information about an event found in a single
    detector. One of these fields is `ranking_stat`. This can be aliased to
    any other string at initialization by setting the keyword argument
    `ranking_stat_alias`; the default is `new_snr`.

    By default, the array is initialized with a `template_id` field. If given
    a `TmpltInspiral` array (which has a `template_id` field), the waveform
    parameters of the events will be added to the array. See `expand_templates`
    for details.

    End times are stored in the `end_time_s` field (for seconds) and the
    `end_time_ns` field (for nanoseconds). The `end_time` property converts
    these into a floats.

    Examples
    --------
    Create a new empty SnglEvent array of length 10 and default fields:

    >>> sngls = SnglEvent(10)
    >>> sngls.all_fieldnames
    ['end_time_s', 'detector', 'cont_chisq_dof', 'event_id', 'ranking_stat',
     'bank_chisq', 'chisq', 'chisq_dof', 'cont_chisq', 'process_id', 'snr',
     'bank_chisq_dof', 'end_time_ns', 'sigma', 'template_id']
    >>> sngls.aliases
        {'ifo': 'detector', 'new_snr': 'ranking_stat'}

    Create a SngEvent array from an hdfcoinc merged file:

    .. code-block:: python
        hdf = h5py.File('BNS1INJ/H1-HDF_TRIGGER_MERGE_BNS1INJ-1117400416-928800.hdf', 'r')
        hsngls = SnglEvent(len(hdf['end_time']))
        hsngls['event_id'] = numpy.arange(hsngls.size)
        hsngls['snr'] = hdf['snr']
        hsngls['chisq'] = hdf['chisq']
        hsngls['chisq_dof'] = 2.*hdf['chisq_dof'].value - 2
        hsngls['end_time_s'] = hdf['end_time'].value.astype(int)
        hsngls['end_time_ns'] = (hdf['end_time'].value % 1 * 1e9).astype(int)
        hsngls['template_id'] = hdf['template_id']
        hsngls['ifo'] = 'H1'
        hsngls['sigma'] = numpy.sqrt(hdf['sigmasq'])
        hsngls['ranking_stat'] = hsngls.snr
        reweight_idx = numpy.where(hsngls['chisq/chisq_dof > 1'])
        hsngls.ranking_stat[reweight_idx] = hsngls['snr'][reweight_idx] / ((1. + hsngls['chisq/chisq_dof'][reweight_idx]**3.)/2.)**(1./6)

    Get the ranking stat by its alias:

    >>> hsngls.ranking_stat
    array([ 9.59688939,  5.25923089,  5.67155045, ...,  5.08446111,
            7.24139452,  7.55083953])
    >>> hsngls.ranking_stat_alias
        'new_snr'
    >>> hsngls.new_snr
    array([ 9.59688939,  5.25923089,  5.67155045, ...,  5.08446111,
            7.24139452,  7.55083953])

    Expand the template field (see TmpltInspiral help for how to create
    `templates` from an hdf bank file) (note that properties such as `mchirp`,
    are also copied from the `templates` array to the `hsngls` array):

    >>> hsngls.fieldnames
    ('end_time_s', 'detector', 'cont_chisq_dof', 'event_id', 'ranking_stat',
     'bank_chisq', 'chisq', 'chisq_dof', 'cont_chisq', 'process_id', 'snr',
     'bank_chisq_dof', 'end_time_ns', 'sigma', 'template_id')
    >>> templates.fieldnames
        ('mass1', 'mass2', 'spin1z', 'spin2z', 'template_hash', 'template_id')
    >>> 'mchirp' in dir(hsngls)
        False
    >>> hsngls = hsngls.expand_templates(templates)
    >>> hsngls.fieldnames
    ('end_time_s', 'detector', 'cont_chisq_dof', 'event_id', 'ranking_stat',
     'bank_chisq', 'chisq', 'chisq_dof', 'cont_chisq', 'process_id', 'snr',
     'bank_chisq_dof', 'end_time_ns', 'sigma', 'template_id', 'mass1', 'mass2',
     'spin1z', 'spin2z', 'template_hash')
    >>> hsngls.ranking_stat, hsngls.template_id, hsngls.mass1
    (array([ 9.59688939,  5.25923089,  5.67155045, ...,  5.08446111,
             7.24139452,  7.55083953]),
     array([   0,    0,    0, ..., 4030, 4030, 4030]),
     array([ 1.71731389,  1.71731389,  1.71731389, ...,  2.11106491,
             2.11106491,  2.11106491], dtype=float32))
    >>> 'mchirp' in dir(hsngls)
        True
    >>> hsngls.mchirp
    array([ 1.25239313,  1.25239313,  1.25239313, ...,  1.25727332,
            1.25727332,  1.25727332], dtype=float32)

    Sort by the ranking stat (note that all other fields are sorted too):

    >>> hsngls.sort(order='ranking_stat')
    >>> hsngls.event_id, hsngls.ranking_stat, hsngls.template_id, hsngls.mass1
    (array([253515,  39853, 155364, ..., 166129, 211314,   5742]),
     array([   5.00000016,    5.00000954,    5.00000962, ...,  116.15787474,
             124.83552448,  155.72376072]),
     array([3785,  590, 2314, ..., 2473, 3159,   86]),
     array([ 2.99986959,  1.85319662,  1.41802227, ...,  2.99999189,
             1.87226772,  2.88923597], dtype=float32))
    """
    default_name = 'sngl_event'
    ranking_stat_alias = 'new_snr'

    # we define the following as static parameters because they are inherited
    # by CoincEvent
    _static_fields = {
        # ids
        'process_id': int,
        'event_id': int,
        'template_id': int,
        # end times
        'end_time_s': "int_4s",
        'end_time_ns': "int_4s",
        }

    @classmethod
    def default_fields(cls, ranking_stat_alias='new_snr'):
        """
        The ranking stat alias can be set; default is ``'new_snr'``.
        """
        fields = {
            ('ifo', 'detector'): 'S2',
            'sigma': float,
            (ranking_stat_alias, 'ranking_stat'): float,
            'snr': float,
            'chisq': float,
            'chisq_dof': float,
            'bank_chisq': float,
            'bank_chisq_dof': float,
            'cont_chisq': float,
            'cont_chisq_dof': float,
            }
        return dict(cls._static_fields.items() + fields.items())

    def __new__(cls, shape, name=None, ranking_stat_alias='new_snr', **kwargs):
        """
        Adds ranking_stat_alias to initialization.
        """
        field_args = {'ranking_stat_alias': ranking_stat_alias}
        return super(SnglEvent, cls).__new__(cls, shape, name=name,
            field_args=field_args, **kwargs)

    @property
    def end_time(self):
        return self.end_time_s + 1e-9*self.end_time_ns

    def expand_templates(self, tmplt_inspiral_array, get_fields=None,
            selfs_map_field='template_id', tmplts_map_field='template_id'):
        """
        Given an array of templates, adds the fields from the templates to this
        array. This is done by getting all rows in the `tmplt_inspiral_array`
        such that:
        ``self[selfs_map_field] == tmplt_inspiral_array[tmplts_map_field]``.

        Parameters
        ----------
        tmplt_inspiral_array : (any subclass of) LSCArray
            The array of templates with additional fields to add.
        get_fields : {None | (list of) strings}
            The names of the fields to get from the tmplt_inspiral_array.
            If ``None``, all fields will be retrieved.
        selfs_map_field : {'template.template_id' | string}
            The name of the field in self's current template sub-array to use
            to match to templates in the tmplt_inspiral_array. Default is
            `template_id`.
        tmplts_map_field : {'template_id' | string}
            The name of the field in the tmplt_inspiral_array to match.
            Default is `template_id`.

        Returns
        -------
        new_array : new instances of this array
            A copy of this array with the template sub-array containing all
            of the fields specified by `get_fields`. All methods and properties
            of the `tmpl_inspiral_array` that are not methods/properties of
            this array are added to `new_array`.
        """
        return self.join(tmplt_inspiral_array, selfs_map_field,
            other_map_field=tmplts_map_field, get_fields=get_fields)


class CoincEvent(SnglEvent):
    """
    Has default fields for storing information about an event found in multiple
    detectors. Like `SnglEvent`, one of these fields is `ranking_stat`, which
    can be aliased to another string at initialization (default is `new_snr`). 
    Also stores end times in a `_s` and `_ns` fields, but these are added
    together using the `end_time` property. This also has a `template_id`
    column, and it inherits the `expand_templates` method for adding paramteer
    information (see `SnglEvent` for details).

    Unlike `SnglEvent`, `CoincEvent` initialization has an additional keyword
    `detectors`. This is a list of the names of possible detectors that could
    form events. The provided names are added as fields to CoincEvent. By
    default, each detector field has sub-fields `event_id`, `end_time_s`, and
    `end_time_ns`. These can be used to store information about the single
    events that formed the coincidence. For example, if
    ``detectors=['H1', 'L1']``, the array will have fields `H1` and `L1`, each 
    with their own `event_id` and `end_time_(n)s` fields. If a SnglEvent array
    is provided, the detector fields can be expanded to include all of the
    single event parameters, similar to `expand_templates`; see `expand_sngls`
    for details.  Note that the `end_time` property will also work on the
    detector end times. For example if `CoincEvent` array `foo` has detector
    field `H1`, ``foo.H1.end_time`` will return the H1 `end_time_s` and
    `end_time_ns` as a float.
    
    If one or more events in the array were not detected in all of the
    detectors provided, the `event_id` for the detector fields they were not
    detected is is set to `lscarrays.ID_NOT_SET`. The property `detected_in`
    will return an array of strings listing the detectors that each event was
    found in.

    Other differences from SnglEvent are the default fields for statistics:
    this has `ifar(_exc)` and `ifap(_exc)` fields for storing inverse false
    alarm rates and inverse false alarm probabilities, respectively. The
    properties `far` and `fap` return the inverse of these. Inverse fars
    (faps) are stored because a new empty array will initialize these to 0.
    Thus the initialized fars and faps = infinity.


    Examples
    --------
    Create a new empty CoincEvent array of length 10 and default fields (note
    that the default detectors are `detector1` and `detector2`):

    >>> coincs = CoincEvent(10)
    >>> coincs.all_fieldnames
    ['end_time_s', 'event_id', 'ranking_stat', 'ifar_exc', 'ifap_exc',
     'detector1.event_id', 'detector1.end_time_s', 'detector1.end_time_ns',
     'detector2.event_id', 'detector2.end_time_s', 'detector2.end_time_ns',
     'process_id', 'ifar', 'ifap', 'snr', 'end_time_ns', 'template_id']

    Create an array from a statmap hdfcoinc file:

    .. code-block:: python
        statmap = h5py.File('BNS1INJ_coinc/H1L1-HDFINJFIND_BNS1INJ_INJ_INJ-1117400416-928800.hdf', 'r')
        fg = statmap['found']
        coincs = CoincEvent(len(fg['fap']), detectors=[statmap.attrs['detector_1'], statmap.attrs['detector_2']], names=['ifap', 'ifar', 'ranking_stat', 'template_id', 'event_id'])
        coincs['ifap'] = 1./fg['fap'].value
        coincs['ifar'] = fg['ifar']
        coincs['ranking_stat'] = fg['stat']
        coincs['template_id'] = fg['template_id']
        coincs['event_id'] = numpy.arange(len(coincs))
        coincs[statmap.attrs['detector_1']]['event_id'] = fg['trigger_id1']
        coincs[statmap.attrs['detector_1']]['end_time_s'] = fg['time1'].value.astype(int)
        coincs[statmap.attrs['detector_1']]['end_time_ns'] = (fg['time1'].value % 1 * 1e9).astype(int)
        coincs[statmap.attrs['detector_2']]['event_id'] = fg['trigger_id2']
        coincs[statmap.attrs['detector_2']]['end_time_s'] = fg['time2'].value.astype(int)
        coincs[statmap.attrs['detector_2']]['end_time_ns'] = (fg['time2'].value % 1 * 1e9).astype(int)

    Add information about the single events (see SnglEvent help for how to
    create `hsngls` and `lsngls` in this example):

    >>> coincs.detectors
        ('H1', 'L1')
    >>> coincs.H1.fieldnames
        ('event_id', 'end_time_s', 'end_time_ns')
    >>> coincs = coincs.expand_sngls(hsngls)
    >>> coincs.H1.fieldnames
    ('end_time_s', 'cont_chisq_dof', 'event_id', 'ranking_stat', 'bank_chisq',
     'chisq', 'chisq_dof', 'cont_chisq', 'process_id', 'snr', 'bank_chisq_dof',
     'end_time_ns', 'sigma', 'template_id', 'mass1', 'mass2', 'spin1z',
     'spin2z', 'template_hash')
>>> coincs.L1.fieldnames
    ('event_id', 'end_time_s', 'end_time_ns')
>>> coincs = coincs.expand_sngls(lsngls)
>>> coincs.L1.fieldnames
    ('end_time_s', 'cont_chisq_dof', 'event_id', 'ranking_stat', 'bank_chisq',
     'chisq', 'chisq_dof', 'cont_chisq', 'process_id', 'snr', 'bank_chisq_dof',
     'end_time_ns', 'sigma', 'template_id')

    >>> coincs.ranking_stat, coincs['H1']['ranking_stat'], coincs.L1.ranking_stat
    (array([ 76.14434814, 8.63374901, ..., 139.97875977, 18.55389023]),
     array([ 44.2959137,  6.21183914, ..., 124.83552448, 14.78488204]),
     array([ 61.93411348, 5.9962226 , ...,  63.32571184, 11.20955446]))

    Expand the template field (see TmpltInspiral help for how to create
    `templates` from an hdf bank file):

    >>> coincs.fieldnames
        ('event_id', 'ranking_stat', 'ifar', 'ifap', 'template_id', 'H1', 'L1')
    >>> templates.fieldnames
        ('mass1', 'mass2', 'spin1z', 'spin2z', 'template_hash', 'template_id')
    >>> coincs = coincs.expand_templates(templates)
    >>> coincs.fieldnames
        ('event_id', 'ranking_stat', 'ifar', 'ifap', 'template_id', 'H1', 'L1',
         'mass1', 'mass2', 'spin1z', 'spin2z', 'template_hash')
    """
    default_name = 'coinc_event'
    # add detectors as a persistent attribute
    __persistent_attributes__ = ['detectors'] + \
        SnglEvent.__persistent_attributes__

    @classmethod
    def default_fields(cls, ranking_stat_alias='new_snr',
            detectors=['detector1', 'detector2']):
        """
        Both the ranking stat alias and the maximum number of single-detector
        fields can be set.
        """
        fields = {
            'ifar': float,
            'ifar_exc': float,
            'ifap': float,
            'ifap_exc': float,
            (ranking_stat_alias, 'ranking_stat'): float,
            'snr': float,
            }
        sngls = {
            det: {
                'event_id': int,
                'end_time_s': int,
                'end_time_ns': int,
                }.items()
            for det in detectors
            }
        # we'll inherit SnglEvent's _static_fields
        return dict(cls._static_fields.items() + fields.items() + \
            sngls.items())

    def __new__(cls, shape, name=None, detectors=['detector1', 'detector2'],
            ranking_stat_alias='new_snr', **kwargs):
        """
        Adds nsngls and ranking_stat_alias to initialization.
        """
        field_args = {'ranking_stat_alias': ranking_stat_alias,
            'detectors': detectors}
        # add the detectors to the requested names if not already
        if 'names' in kwargs:
            names = kwargs.pop('names')
            if isinstance(names, str) or isinstance(names, unicode):
                names = [names]
            else:
                names = list(names)
            names += [det for det in detectors if det not in names]
            kwargs['names'] = names
        # Note: we need to call _LSCArrayWithDefaults directly, as using
        # super will lead to SnglEvent's __new__, which sets its own field args
        obj = _LSCArrayWithDefaults.__new__(cls, shape, name=name,
            field_args=field_args, **kwargs)
        # set the detectors attribute
        obj.addattr('detectors', tuple(sorted(detectors)))
        return obj

    @property
    def far(self):
        return 1./self.ifar

    @property
    def false_alarm_rate(self):
        return self.far

    @property
    def far_exc(self):
        return 1./self.ifar_exc

    @property
    def false_alarm_rate_exc(self):
        return self.far_exc

    @property
    def fap(self):
        return 1./self.ifap

    @property
    def false_alarm_probability(self):
        return self.fap

    @property
    def fap_exc(self):
        return 1./self.ifap_exc

    @property
    def false_alarm_probability_exc(self):
        return self.fap_exc

    @property
    def detected_in(self):
        """
        Returns the names of the detectors that contributed to each coinc
        event. This is found by returning all of the detectors for which
        ``self[detector].event_id != lscarrays.ID_NOT_SET``.
        """
        detectors = numpy.array(self.detectors)
        mask = numpy.vstack([
            self[det]['event_id'] != ID_NOT_SET \
            for det in detectors]).T
        return numpy.array([','.join(detectors[numpy.where(mask[ii,:])]) \
            for ii in range(self.size)])

    def expand_sngls(self, sngl_event_array, get_fields=None,
            sngls_map_field='event_id'):
        """
        Given an array of singles, adds sub-arrays of the single-detector
        trigger data named by ifo to self. For each detector in that is in
        `self.detectors`, the data is retrieved such that:
        ``self[detector]['event_id'] == sngl_event_array[where('detector' == detector)]['event_id']``.
        If the `sngl_event_array` only has data of a subset of the detectors
        in `self.detectors`, fields will only be added for those detectors.
        For example, if ``self.detectors == ['H1', 'L1']`` but
        `sngl_event_array` only has H1 data, then only `self['H1']` will be
        added.

        Parameters
        ----------
        sngl_event_array : (any subclass of) LSCArray
            The array of singles with additional fields to add.
        get_fields : {None | (list of) strings}
            The names of the fields to get from the sngl_event_array.
            If `None`, all fields will be retrieved.
        sngls_map_field : {'event_id' | string}
            The name of the field in the `sngl_event_array` to match.
            Default is `'event_id'`.

        Returns
        -------
        newarray : new instance of this array
            A copy of this array with the single detector fields added.

        Notes
        -----
        * The `sngl_event_array` must have a `detector` field.
        * The `sngl_event_array` must have one and only one row in it for a
          given ``detector, event_id`` pair in self. If `sngl_event_array`
          has a `detector` that is in self, but does not have an `event_id`
          for that detector, a KeyError is raised.
        * If some events in this array were not found in all of the detectors
          listed in `detectors` --- indicated by ``event_id == ID_NOT_SET``
          for a detector --- then a zeroed entry will be added to that detector
          for that row. For example, if ``detectors = ['H1', 'L1', 'V1']``,
          but row `i` has ``V1.event_id == ID_NOT_SET``, then the added fields
          for `V1` in that row will be zeroed. Likewise, the attribute
          `detected_in` will return `'H1,L1'` for that row.
        """
        #
        # Note: since this join has to do deal with possibly missing ifos, it
        # is a bit more complicated then the join that's carried out in
        # join_arrays; thus we do the join manually here, instead of calling
        # join_arrays.
        #
        # cycle over the ifos in self.sngl_ifos.ifos that are in the
        # sngl_event_array
        others_dets = numpy.unique(sngl_event_array['detector'])
        if not any(others_dets):
            raise ValueError("sngl_event_array's detector field is not " +
                "populated!")
        # if getting all fields, exclude the detector field, since that will
        # be the sub-array's name
        if get_fields is None:
            get_fields = [name for name in sngl_event_array.fieldnames \
                if name != 'detector']
        new_self = self
        for det in self.detectors:
            if det not in others_dets:
                # just skip
                continue
            if others_dets.size == 1:
                # we can just look at the whole array
                other_array = sngl_event_array
            else:
                # pull out the rows with detector == det
                other_array = sngl_event_array[numpy.where(
                    sngl_event_array['detector'] == det)]
            # figure out what rows in self have an entry for this detector
            mask = self[det]['event_id'] != ID_NOT_SET
            if mask.all():
                # every row has an entry, just map all
                map_indices = None
            else:
                map_indices = numpy.where(mask)
            # join
            new_self = new_self.join(other_array,
                '%s.event_id' % det, det,
                sngls_map_field, get_fields=get_fields,
                map_indices=map_indices)
        return new_self


# we'll vectorize TimeDelayFromEarthCenter for faster processing of end times
def _time_delay_from_center(geocent_end_time, detector, ra, dec):
    return geocent_end_time + lal.TimeDelayFromEarthCenter(detector.location,
        ra, dec, geocent_end_time)
time_delay_from_center = numpy.vectorize(_time_delay_from_center)

class SimInspiral(Waveform):
    """
    Subclasses Waveform, adding default fields for sky-location and detector-
    specific parameters for injections. Also has default field `simulation_id` 
    for indexing.
    
    Similar to `CoincEvent`, `SimInspiral` takes a `detectors` keyword at
    initialization. This adds fields named by the provided detectors for
    storing site-specific information. By default, each detector field has
    sub-fields `eff_dist`, and `sigma`. Detector end times are not stored;
    instead, they are calucated on the fly using the `end_time` method. If a
    detector name is provided to ``end_time()``, the end times at that
    detector are returned using the geocentric end times (which are stored).
    See `end_time` for details.

    The `recovered` field is used to store information about any recovered
    (coinc or single) events that the injections were associated with. By
    default, this is a sub-array with one field, `'event_id'`, which is the
    `event_id`(s) of the event(s). If given a `CoincEvent` or `SnglEvent` array
    with the events, the `recovered` field can be expanded to include all of
    recovered information; see `expand_recovered` for details. Only injections
    for which `recovered.event_id != ID_NOT_SET` are exanded; i.e., the
    recovered `event_id` field has to be set in order for an injection to be
    considered recovered. The attribute `isrecovered` checks this, returning
    a boolean array indicating whether each event was recovered or not.

    A single injection may be mapped to multiple recovered events. The maximum
    number that it may be mapped to is set by the `nrecovered` keyword at
    initialization. By default this is 1. If greater than 1, each row's
    recovered field will be a subarray with size == `nrecovered`. If an
    injection is mapped to fewer events than `nrecovered`, the extra entries
    will be left empty.

    Notes
    -----
    * When a new array is initialized with a ``recovered`` field, all of the
      injections are marked as not recovered (i.e.,
      ``arr['recovered']['isfound'] = False``).
    * By default `(Coinc|Sngl)Event` arrays do not contain information about
      the template they were found with, only the `template_id`. To get the
      parameters of the recovered waveforms, run `expand_templates` on the
      event array *first*, then run `expand_recovered` on the expanded event
      array.

    Examples
    --------
    Intialize an empty SimInspiral array:

    >>> sims = SimInspiral(10)
    >>> sorted(sims.all_fieldnames)
    ['amp_order', 'approximant', 'argument_periapsis', 'dec',
     'detector1.eff_dist', 'detector1.sigma', 'detector2.eff_dist',
     'detector2.sigma', 'distance', 'duration', 'eccentricity', 'f_max',
     'f_min', 'f_ref', 'frame_axis', 'geocent_end_time_ns',
     'geocent_end_time_s', 'inclination', 'lambda1', 'lambda2', 'mass1',
     'mass2', 'min_vol', 'modes_choice', 'phase_order', 'phi_ref',
     'polarization', 'process_id', 'quadparam1', 'quadparam2', 'ra',
     'recovered.event_id', 'sample_rate', 'segment_length', 'simulation_id',
     'spin1x', 'spin1y', 'spin1z', 'spin2x', 'spin2y', 'spin2z', 'spin_order',
     'taper', 'tidal_order', 'volume_weight']

         
    Create a SimInspiral array from an hdfinjfind file:

    .. code-block:: python

        hdfinjfind = h5py.File('BNS1INJ_coinc/H1L1-HDFINJFIND_BNS1INJ_INJ_INJ-1117400416-928800.hdf', 'r')
        hdfinj = hdfinjfind['injections']
        det1 = hdfinjfind.attrs['detector_1']
        det2 = hdfinjfind.attrs['detector_2']
        sims = SimInspiral(len(hdfinj['end_time']), detectors=[det1, det2], names=['simulation_id', 'geocent_end_time_s', 'geocent_end_time_ns', 'distance', 'mass1', 'mass2', 'ra', 'dec', 'recovered'])
        sims['simulation_id'] = numpy.arange(sims.size)
        sims['geocent_end_time_s'] = hdfinj['end_time'].value.astype(int)
        sims['geocent_end_time_ns'] = (hdfinj['end_time'].value % 1 * 1e9).astype(int)
        sims['distance'] = hdfinj['distance']
        sims['mass1'] = hdfinj['mass1']
        sims['mass2'] = hdfinj['mass2']
        sims['ra'] = hdfinj['longitude']
        sims['dec'] = hdfinj['latitude']
        sims['recovered']['event_id'][hdfinjfind['found']['injection_index']] = numpy.arange(len(hdfinjfind['found']['fap']))
        sims[det1]['eff_dist'] = hdfinj['eff_dist_%s' %(det1[0].lower())]
        sims[det2]['eff_dist'] = hdfinj['eff_dist_%s' %(det2[0].lower())]

    Add recovered information (see `CoincEvent` on how to create `coincs` from
    an hdf file of found injections, and `TmpltInspiral` on how to create
    `templates` from a bank hdf file). Note that the property `far` is copied
    from the `CoincEvent` array to the `SimInspiral` array:

    >>> coincs = coincs.expand_templates(templates)
    >>> sorted(coincs.all_fieldnames)
        ['H1.end_time_ns', 'H1.end_time_s', 'H1.event_id', 'L1.end_time_ns',
         'L1.end_time_s', 'L1.event_id', 'event_id', 'ifap', 'ifar', 'mass1',
         'mass2', 'ranking_stat', 'spin1z', 'spin2z', 'template_hash',
         'template_id']
    >>> sorted(sims.all_fieldnames)
        ['H1.eff_dist', 'H1.sigma', 'L1.eff_dist', 'L1.sigma', 'dec',
         'distance', 'geocent_end_time_ns', 'geocent_end_time_s', 'mass1',
         'mass2', 'ra', 'recovered.event_id', 'simulation_id']
    >>> 'far' in dir(coincs)
        True
    >>> 'far' in dir(sims)
        False
    >>> sims = sims.expand_recovered(coincs)
    >>> sorted(sims.all_fieldnames)
        ['H1.eff_dist', 'H1.sigma', 'L1.eff_dist', 'L1.sigma', 'dec',
         'distance', 'geocent_end_time_ns', 'geocent_end_time_s', 'mass1',
         'mass2', 'ra', 'recovered.H1.end_time_ns', 'recovered.H1.end_time_s',
         'recovered.H1.event_id', 'recovered.L1.end_time_ns',
         'recovered.L1.end_time_s', 'recovered.L1.event_id',
         'recovered.event_id', 'recovered.ifap', 'recovered.ifar',
         'recovered.mass1', 'recovered.mass2', 'recovered.ranking_stat',
         'recovered.spin1z', 'recovered.spin2z', 'recovered.template_hash',
         'recovered.template_id', 'simulation_id']
    >>> 'far' in dir(sims)
        True
    >>> sims.isrecovered, sims.recovered.far
    (array([False, False, False, ..., False, False, False], dtype=bool),
     array([ inf,  inf,  inf, ...,  inf,  inf,  inf]))
    >>> recsims = sims.get_recovered()
    >>> recsims.mchirp, recsims.recovered.mchirp, recsims.recovered.far
    (array([ 1.81795762,  1.74482645,  ...,  1.1793251 ,  1.41826222]),
     array([ 1.81705785,  1.74960184,  ...,  1.17973995,  1.41618061], dtype=float32),
     array([ 5.31908815e-05, 7.11640778e+00,  ...,  5.31908815e-05,  5.31908815e-05]))
    """
    default_name = 'sim_inspiral'
    # add detectors as a persistent attribute
    __persistent_attributes__ = ['detectors'] + \
        Waveform.__persistent_attributes__

    @classmethod
    def default_fields(cls, detectors=['detector1', 'detector2'],
            nrecovered=1):
        """
        The number of ifos stored in site params and the maximum number of
        events that can be stored in the recovered fields can be set.
        """
        if detectors is None:
            detectors = []
        fields = {
            # ids
            'simulation_id': int,
            'process_id': int,
            # location params
            'geocent_end_time_s': 'int_4s',
            'geocent_end_time_ns': 'int_4s',
            'distance': float,
            ('right_ascension', 'ra'): float,
            ('declination', 'dec'): float,
            'polarization': float,
            # recovered params
            'recovered': ({
                'event_id': int,
                }.items(), nrecovered),
            # distribution params
            'min_vol': float,
            'volume_weight': float,
            # XXX: How to store information about mass and spin
            # distributions?
            }
        # site params
        site_params = {
            det: {
                'eff_dist': float,
                'sigma': float,
                }.items()
            for det in detectors
            }
        # we'll inherit static fields from Waveform
        return dict(cls._static_fields.items() + fields.items() + \
            site_params.items())

    def __new__(cls, shape, name=None, detectors=['detector1', 'detector2'],
            nrecovered=None, **kwargs):
        """
        Adds detectors and nrecovered to initialization.
        """
        field_args = {'nrecovered': nrecovered,
            'detectors': detectors}
        # add the detectors to the requested names if not already
        if 'names' in kwargs:
            names = kwargs.pop('names')
            if isinstance(names, str) or isinstance(names, unicode):
                names = [names]
            else:
                names = list(names)
            names += [det for det in detectors if det not in names]
            kwargs['names'] = names
        obj = super(SimInspiral, cls).__new__(cls, shape, name=name,
            field_args=field_args, **kwargs)
        # set the detectors attribute
        obj.addattr('detectors', tuple(sorted(detectors)))
        return obj

    @property
    def geocent_end_time(self):
        return self.geocent_end_time_s + 1e-9*self.geocent_end_time_ns

    def end_time(self, detector=None):
        """
        Returns the end time in the given detector. If detector is None,
        returns the geocentric end time.
        """
        geocent_end_time = self.geocent_end_time
        if detector is None:
            return geocent_end_time
        else:
            detector = lalsim.DetectorPrefixToLALDetector(detector)
            return time_delay_from_center(geocent_end_time, detector, self.ra,
                self.dec)

    @property
    def isrecovered(self):
        """Returns boolean array indicating rows for which the recovered
        `event_id`s are not equal to `ID_NOT_SET`.
        """
        return self['recovered']['event_id'] != ID_NOT_SET

    def expand_recovered(self, event_array, get_fields=None,
            selfs_map_field='recovered.event_id',
            events_map_field='event_id'):
        """
        Given an array of (coinc) events, replaces the recovered field with a
        sub-array of the recovered data. This is done by getting
        all rows in the ``event_array`` such that:
        ``self[where(self.isrecovered)][selfs_map_field] == event_array[events_map_field]``.

        Parameters
        ----------
        event_array : {SimEvent | CoincEvent}
            The array of events with additional fields to add.
        get_fields : {None | (list of) strings}
            The names of the fields to get from the ``event_array``.
            If ``None``, all fields will be retrieved.
        selfs_map_field : {'recovered.event_id' | string}
            The name of the field in self to use to match to events in
            ``event_array``. Default is ``recovered.event_id``.
        events_map_field : {'event_id' | string}
            The name of the field in the ``event_array`` to match.
            Default is ``event_id``.

        Returns
        -------
        new array : new instances of this array
            A copy of this array with the ``recovered`` sub-array containing
            all of the fields specified by ``get_fields``.
        
        Notes
        -----
        * Only elements for which ``self.isrecovered == True`` will be
          expanded. All other elements will have a zeroed row in ``recovered``
          field of the output array.

        * If the injection is mapped to multiple events, all of the events will
          be expanded.

        * To get the recovered parameters, expand the event array's template
          field first, then expand the recovered here. See
          (Coinc|Sngl)Event.expand_templates for details.

        * The coincs themselves contain subarrays of the single-detector
          triggers. To get access to the single-detector recovered information,
          add the single-detector information to the coinc event array first
          (see ``CoincEvent.exapand_sngls``), then expand the coincs here.

        * If the injections are mapped to single events rather than coinc
          events, you can map to the single events by passing a SnglEvent array
          as the ``event_array``.
        """
        return self.join(event_array, selfs_map_field, 'recovered',
            other_map_field=events_map_field, get_fields=get_fields,
            map_indices=numpy.where(self.isrecovered))

    @property
    def optimal_snr(self):
        """
        Gives the maximum SNR that the injections can have.
        """
        return numpy.sqrt(numpy.array([self[det]['sigma']**2
            for det in self.detectors]).sum())/self['distance']

    @property
    def detected_in(self):
        """
        Returns the names of the detectors that the injections was detected in.
        This is done by returning all of the detectors in self's detectors for
        which ``self['recovered'][detector].event_id != lscarrays.ID_NOT_SET``.

        Note: ``expand_recovered`` must have been run first on the list of
        events.
        """
        detectors = numpy.array([det for det in self.detectors \
            if det in self['recovered'].fieldnames])
        if detectors.size == 0:
            raise ValueError("No detectors found in this array's recovered " +
                "field. Did you run expand_recovered?")
        mask = numpy.vstack([
            self['recovered'][det]['event_id'] != ID_NOT_SET \
            for det in detectors]).T
        return numpy.array([','.join(detectors[numpy.where(mask[ii,:])]) \
            for ii in range(self.size)])

    @property
    def recovered_idx(self):
        """
        Returns the indices in self that were recovered.
        """
        return numpy.where(self.isrecovered)[0]

    def get_recovered(self):
        """
        Returns the elements in self that were recovered.
        """
        return self[self.recovered_idx]

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
This modules provides definitions of, and helper functions for, FieldArrays.
FieldArrays are wrappers of numpy recarrays with additional functionality
useful for storing and retrieving data created by a search for gravitationa
waves.
"""

import os, sys, types, re, copy, numpy, inspect
from pycbc_glue.ligolw import types as ligolw_types
from pycbc import coordinates, conversions, cosmology
from pycbc.detector import Detector
from pycbc.waveform import parameters

# what functions are given to the eval in FieldArray's __getitem__:
_numpy_function_lib = {_x: _y for _x,_y in numpy.__dict__.items()
                       if isinstance(_y, (numpy.ufunc, float))}

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
VIRTUALFIELD_DTYPE = 'VIRTUAL'

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
    When FieldArrays is first loaded, the default is True.

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
    >>> from pycbc.io import FieldArrays
    >>> FieldArrays.lstring_as_obj()
        True
    >>> FieldArrays.FieldArray.from_arrays([numpy.zeros(10)], dtype=[('foo', 'lstring')])
    FieldArray([(0.0,), (0.0,), (0.0,), (0.0,), (0.0,), (0.0,), (0.0,), (0.0,),
           (0.0,), (0.0,)], 
          dtype=[('foo', 'O')])
    >>> FieldArrays.lstring_as_obj(False)
        False
    >>> FieldArrays.FieldArray.from_arrays([numpy.zeros(10)], dtype=[('foo', 'lstring')])
    FieldArray([('0.0',), ('0.0',), ('0.0',), ('0.0',), ('0.0',), ('0.0',),
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
        ilwd_as_int(_default_types_status['ilwd_as_int'])
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
    """Given a python string, gets FieldArray field names used in it. This
    differs from get_vars_from_arg in that any identifier with a '.' in it
    will be treated as one identifier. For example, if
    ``arg = '3*narf/foo.bar'``, this will return ``set(['narf', 'foo.bar'])``.
    """
    return set(_fieldparser.findall(arg))

# this parser looks for fields inside a class method function. This is done by
# looking for variables that start with self.{x} or self["{x}"]; e.g.,
# self.a.b*3 + self.c, self['a.b']*3 + self.c, self.a.b*3 + self["c"], all
# return set('a.b', 'c').
_instfieldparser = re.compile(
    r'''self(?:\.|(?:\[['"]))(?P<identifier>[\w_][.\w\d_]*)''')
def get_instance_fields_from_arg(arg):
    """Given a python string definining a method function on an instance of an
    FieldArray, returns the field names used in it. This differs from
    get_fields_from_arg in that it looks for variables that start with 'self'.
    """
    return set(_instfieldparser.findall(arg))

def get_needed_fieldnames(arr, names):
    """Given a FieldArray-like array and a list of names, determines what
    fields are needed from the array so that using the names does not result
    in an error.

    Parameters
    ----------
    arr : instance of a FieldArray or similar
        The array from which to determine what fields to get.
    names : (list of) strings
        A list of the names that are desired. The names may be either a field,
        a virtualfield, a property, a method of ``arr``, or any function of
        these. If a virtualfield/property or a method, the source code of that
        property/method will be analyzed to pull out what fields are used in
        it.

    Returns
    -------
    set
        The set of the fields needed to evaluate the names.
    """
    fieldnames = set([])
    # we'll need the class that the array is an instance of to evaluate some 
    # things
    cls = arr.__class__
    if isinstance(names, (str, unicode)):
        names = [names]
    # parse names for variables, incase some of them are functions of fields
    parsed_names = set([])
    for name in names:
        parsed_names.update(get_fields_from_arg(name))
    # only include things that are in the array's namespace
    names = list(parsed_names & (set(dir(arr)) | set(arr.fieldnames)))
    for name in names:
        if name in arr.fieldnames:
            # is a field, just add the name
            fieldnames.update([name])
        else:
            # the name is either a virtualfield, a method, or some other
            # property; we need to evaluate the source code to figure out what
            # fields we need
            try:
                # the underlying functions of properties need to be retrieved
                # using their fget attribute
                func = getattr(cls, name).fget
            except AttributeError:
                # no fget attribute, assume is an instance method
                func = getattr(arr, name)
            # evaluate the source code of the function
            try:
                sourcecode = inspect.getsource(func)
            except TypeError:
                # not a function, just pass
                continue
            # evaluate the source code for the fields
            possible_fields = get_instance_fields_from_arg(sourcecode)
            # some of the variables returned by possible fields may themselves
            # be methods/properties that depend on other fields. For instance,
            # mchirp relies on eta and mtotal, which each use mass1 and mass2;
            # we therefore need to anayze each of the possible fields
            fieldnames.update(get_needed_fieldnames(arr, possible_fields))
    return fieldnames


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
    new_arr = merge_list[0].__class__(merge_list[0].shape, dtype=new_dt)
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
        if isinstance(names, (str, unicode)):
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


#
# =============================================================================
#
#                           Base FieldArray definitions
#
# =============================================================================
#

# We'll include functions in various pycbc modules in FieldArray's function
# library. All modules used must have an __all__ list defined.
_modules_for_functionlib = [conversions, coordinates, cosmology]
_fieldarray_functionlib = {_funcname : getattr(_mod, _funcname)
                              for _mod in _modules_for_functionlib
                              for _funcname in getattr(_mod, '__all__')}

class FieldArray(numpy.recarray):
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
    To initialize an FieldArray from an already existing arrays, use the
    ``FieldArray.from_arrays`` class method. To initialize from a list of
    tuples, use ``FieldArray.from_records``. See the docstring for those methods
    for details. For more information on initalizing an empty array, see
    ``numpy.recarray`` help.

    Some additional features:

    * **Arbitrary functions**:

      You can retrive functions on fields in the same manner that you access
      individual fields. For example, if you have a FieldArray ``x`` with
      fields ``a`` and ``b``, you can access each field with
      ``x['a'], x['b']``.  You can also do ``x['a*b/(a+b)**2.']``,
      ``x[cos(a)*sin(b)]``, etc. Boolean operations are also possible, e.g.,
      ``x['(a < 3) & (b < 2)']``. Syntax for functions is python. Any numpy
      ufunc, as well as all functions listed in the functionlib attribute, may
      be used. Note that while fields may be accessed as attributes (e.g,
      field ``a`` can be accessed via ``x['a']`` or ``x.a``), functions on
      multiple fields may not (``x.a+b`` does not work, for obvious reasons).

    * **Subfields and '.' indexing**:
      Structured arrays, which are the base class for recarrays and, by
      inheritance, FieldArrays, allows for fields to themselves have fields. For
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
      FieldArray, those can be accessed in the same way as fields are. For
      example, define ``Foo`` as:

    .. code-block:: python

        class Foo(FieldArray):
            @property
            def bar(self):
                return self['a']**2.

            def narf(self, y):
                return self['a'] + y

    Then if we have an instance: ``foo = Foo(100, dtype=[('a', float)])``.
    The ``bar`` and ``narf`` attributes may be accessed via field notation:
    ``foo.bar``, ``foo['bar']``, ``foo.narf(10)`` and ``foo['narf(10)']``.

    * **Virtual fields**:
      Virtual fields are methods wrapped as properties that operate on one or
      more fields, thus returning an array of values. To outside code virtual
      fields look the same as fields, and can be called similarily. Internally,
      no additional data is stored; the operation is performed on the fly when
      the virtual field is called. Virtual fields can be added to an array
      instance with the add_virtualfields method. Alternatively, virtual fields
      can be defined by sub-classing FieldArray:

    .. code-block:: python

        class Foo(FieldArray):
            _virtualfields = ['bar']
            @property
            def bar(self):
                return self['a']**2.

    The fields property returns the names of both fields and virtual fields.

    .. note::
    
        It can happen that a field, virtual field, or function in the
        functionlib have that same name. In that case, precedence is: field,
        virtual field, function. For example, if a function called 'foo' is in
        the function library, and a virtual field is added call 'foo', then
        `a['foo']` will return the virtual field rather than the function.
        Likewise, if the array is initialized with a field called `foo`, or a
        field with that name is added, `a['foo']` will return that field
        rather than the virtual field and/or the function.

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
        array classes defined below. 

    Create an empty array with four rows and two fields called `foo` and
    `bar`, both of which are floats:

    >>> x = FieldArray(4, dtype=[('foo', float), ('bar', float)])

    Set/retrieve a fields using index or attribute syntax:

    >>> x['foo'] = [1.,2.,3.,4.]
    >>> x['bar'] = [5.,6.,7.,8.]
    >>> x
    FieldArray([(1.0, 5.0), (2.0, 6.0), (3.0, 7.0), (4.0, 8.0)], 
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

    Add a virtual field:

    >>> def c(self):
    ...     return self['a'] + self['b']
    ...
    >>> x = x.add_virtualfields('c', c)
    >>> x.fields
    ('a', 'b', 'c')
    >>> x['c']
    array([  6.,   8.,  10.,  12.])

    Create an array with subfields:

    >>> x = FieldArray(4, dtype=[('foo', [('cat', float), ('hat', int)]), ('bar', float)])
    >>> x.fieldnames
        ['foo.cat', 'foo.hat', 'bar']

    Load from a list of arrays (in this case, from an hdf5 file):

    >>> bankhdf = h5py.File('bank/H1L1-BANK2HDF-1117400416-928800.hdf')
    >>> bankhdf.keys()
        [u'mass1', u'mass2', u'spin1z', u'spin2z', u'template_hash']
    >>> templates = FieldArray.from_arrays(bankhdf.values(), names=bankhdf.keys())
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
    >>> sim_array = FieldArray.from_ligolw_table(sim_table)
    >>> sim_array.mass1
    array([ 2.27440691,  1.85058105,  1.61507106, ...,  2.0504961 ,
            2.33554196,  2.02732205], dtype=float32)
    >>> sim_array.waveform
    array([u'SpinTaylorT2', u'SpinTaylorT2', u'SpinTaylorT2', ...,
           u'SpinTaylorT2', u'SpinTaylorT2', u'SpinTaylorT2'], dtype=object)
    
    >>> sim_array = FieldArray.from_ligolw_table(sim_table, columns=['simulation_id', 'mass1', 'mass2'])
    >>> sim_array
    FieldArray([(0, 2.274406909942627, 2.6340370178222656),
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

    """
    _virtualfields = []
    _functionlib = _fieldarray_functionlib
    __persistent_attributes__ = ['name', '_virtualfields', '_functionlib']

    def __new__(cls, shape, name=None, zero=True, **kwargs):
        """Initializes a new empty array.
        """
        obj = super(FieldArray, cls).__new__(cls, shape, **kwargs).view(
            type=cls)
        obj.name = name
        obj.__persistent_attributes__ = [a
            for a in cls.__persistent_attributes__]
        obj._functionlib = {f: func for f,func in cls._functionlib.items()}
        obj._virtualfields = [f for f in cls._virtualfields]
        # zero out the array if desired
        if zero:
            default = default_empty(1, dtype=obj.dtype)
            obj[:] = default
        return obj

    def __array_finalize__(self, obj):
        """Default values are set here.

        See <https://docs.scipy.org/doc/numpy/user/basics.subclassing.html> for
        details.
        """
        if obj is None:
            return
        # copy persistent attributes
        try:
            obj.__copy_attributes__(self)
        except AttributeError:
            pass
        # numpy has some issues with dtype field names that are unicode,
        # so we'll force them to strings here
        if self.dtype.names is not None and \
                any(isinstance(name, unicode) for name in obj.dtype.names):
            self.dtype.names = map(str, self.dtype.names)

    def __copy_attributes__(self, other, default=None):
        """Copies the values of all of the attributes listed in
        `self.__persistent_attributes__` to other.
        """
        [setattr(other, attr, copy.deepcopy(getattr(self, attr, default))) \
            for attr in self.__persistent_attributes__]

    def __getattribute__(self, attr):
        """Allows fields to be accessed as attributes.
        """
        # first try to get the attribute
        try:
            return numpy.ndarray.__getattribute__(self, attr)
        except AttributeError as e:
            # might be a field, try to retrive it using getitem
            if attr in self.fields:
                return self.__getitem__(attr)
            # otherwise, unrecognized
            raise AttributeError(e)

    def __setitem__(self, item, values):
        """Wrap's recarray's setitem to allow attribute-like indexing when
        setting values.
        """
        try:
            return super(FieldArray, self).__setitem__(item, values)
        except ValueError:
            # we'll get a ValueError if a subarray is being referenced using
            # '.'; so we'll try to parse it out here
            fields = item.split('.')
            if len(fields) > 1:
                for field in fields[:-1]:
                    self = self[field]
                item = fields[-1]
            # now try again
            return super(FieldArray, self).__setitem__(item, values)

    def __getbaseitem__(self, item):
        """Gets an item assuming item is either an index or a fieldname.
        """
        # We cast to ndarray to avoid calling array_finalize, which can be slow
        out = self.view(numpy.ndarray)[item]
        # cast back to an instance of self if there are more than one field and
        # element
        if out.size > 1 and out.dtype.fields is not None:
            out = out.view(type(self))
        return out

    def __getsubitem__(self, item):
        """Gets a subfield using `field.subfield` notation.
        """
        try:
            return self.__getbaseitem__(item)
        except ValueError as err:
            subitems = item.split('.')
            if len(subitems) > 1:
                return self.__getbaseitem__(subitems[0]
                    ).__getsubitem__('.'.join(subitems[1:]))
            else:
                raise ValueError(err)

    def __getitem__(self, item):
        """Wraps recarray's  `__getitem__` so that math functions on fields and
        attributes can be retrieved. Any function in numpy's library may be
        used.
        """
        try:
            return self.__getsubitem__(item)
        except ValueError:
            #
            #   arg isn't a simple argument of row, so we'll have to eval it
            #
            # get the function library
            item_dict = dict(_numpy_function_lib.items())
            item_dict.update(self._functionlib)
            # get the set of fields & attributes we will need
            itemvars = get_vars_from_arg(item)
            # pull out any other needed attributes
            itemvars = get_fields_from_arg(item)
            d = {attr: getattr(self, attr)
                for attr in set(dir(self)).intersection(itemvars)}
            item_dict.update(d)
            # pull out the fields: note, by getting the parent fields, we
            # also get the sub fields name
            item_dict.update({fn: self.__getbaseitem__(fn) \
                for fn in set(self.fieldnames).intersection(itemvars)})
            # add any aliases
            item_dict.update({alias: item_dict[name]
                              for alias,name in self.aliases.items()
                              if name in item_dict})
            return eval(item, {"__builtins__": None}, item_dict)

    def __contains__(self, field):
        """Returns True if the given field name is in self's fields."""
        return field in self.fields

    def sort(self, axis=-1, kind='quicksort', order=None):
        """Sort an array, in-place. 

        This function extends the standard numpy record array in-place sort
        to allow the basic use of Field array virtual fields. Only a single
        field is currently supported when referencing a virtual field.

        Parameters
        ----------
        axis : int, optional
            Axis along which to sort. Default is -1, which means sort along the
            last axis.
        kind : {'quicksort', 'mergesort', 'heapsort'}, optional
            Sorting algorithm. Default is 'quicksort'.
        order : list, optional
            When `a` is an array with fields defined, this argument specifies
            which fields to compare first, second, etc.  Not all fields need be
            specified.
        """
        try:
            numpy.recarray.sort(self, axis=axis, kind=kind, order=order)
        except ValueError:
            if isinstance(order, list):
                raise ValueError("Cannot process more than one order field")
            self[:] = self[numpy.argsort(self[order])]

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
        if isinstance(names, (str, unicode)):
            names = [names]
            methods = [methods]
        for name,method in zip(names, methods):
            setattr(self, name, types.MethodType(method, self))
        
    def add_properties(self, names, methods):
        """Returns a view of self with the given methods added as properties.

        From: <http://stackoverflow.com/a/2954373/1366472>.
        """
        cls = type(self)
        cls = type(cls.__name__, (cls,), dict(cls.__dict__))
        if isinstance(names, (str, unicode)):
            names = [names]
            methods = [methods]
        for name,method in zip(names, methods):
            setattr(cls, name, property(method))
        return self.view(type=cls)

    def add_virtualfields(self, names, methods):
        """Returns a view of this array with the given methods added as virtual
        fields. Specifically, the given methods are added using add_properties
        and their names are added to the list of virtual fields. Virtual fields
        are properties that are assumed to operate on one or more of self's
        fields, thus returning an array of values.
        """
        if isinstance(names, (str, unicode)):
            names = [names]
            methods = [methods]
        out = self.add_properties(names, methods)
        if out._virtualfields is None:
            out._virtualfields = []
        out._virtualfields.extend(names)
        return out

    def add_functions(self, names, functions):
        """Adds the given functions to the function library.

        Functions are added to this instance of the array; all copies of
        and slices of this array will also have the new functions included.

        Parameters
        ----------
        names : (list of) string(s)
            Name or list of names of the functions.
        functions : (list of) function(s)
            The function(s) to call.
        """
        if isinstance(names, (str, unicode)):
            names = [names]
            functions = [functions]
        if len(functions) != len(names):
            raise ValueError("number of provided names must be same as number "
                             "of functions")
        self._functionlib.update(dict(zip(names, functions)))

    def del_functions(self, names):
        """Removes the specified function names from the function library.

        Functions are removed from this instance of the array; all copies
        and slices of this array will also have the functions removed.

        Parameters
        ----------
        names : (list of) string(s)
            Name or list of names of the functions to remove.
        """
        if isinstance(names, (str, unicode)):
            names = [names]
        for name in names:
            self._functionlib.pop(name)

    @classmethod
    def from_arrays(cls, arrays, name=None, **kwargs):
        """Creates a new instance of self from the given (list of) array(s).
        This is done by calling numpy.rec.fromarrays on the given arrays with
        the given kwargs. The type of the returned array is cast to this
        class, and the name (if provided) is set.

        Parameters
        ----------
        arrays : (list of) numpy array(s)
            A list of the arrays to create the FieldArray from.
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
        """Creates a new instance of self from the given (list of) record(s).

        A "record" is a tuple in which each element is the value of one field
        in the resulting record array. This is done by calling
        `numpy.rec.fromrecords` on the given records with the given kwargs.
        The type of the returned array is cast to this class, and the name
        (if provided) is set.

        Parameters
        ----------
        records : (list of) tuple(s)
            A list of the tuples to create the FieldArray from.
        name : {None|str}
            What the output array should be named.

        Other Parameters
        ----------------
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
    def from_kwargs(cls, **kwargs):
        """Creates a new instance of self from the given keyword arguments.
        Each argument will correspond to a field in the returned array, with
        the name of the field given by the keyword, and the value(s) whatever
        the keyword was set to. Each keyword may be set to a single value or
        a list of values. The number of values that each argument is set to
        must be the same; this will be the size of the returned array.

        Examples
        --------
        Create an array with fields 'mass1' and 'mass2':
        >>> a = FieldArray.from_kwargs(mass1=[1.1, 3.], mass2=[2., 3.])
        >>> a.fieldnames
        ('mass1', 'mass2')
        >>> a.mass1, a.mass2
        (array([ 1.1,  3. ]), array([ 2.,  3.]))

        Create an array with only a single element in it:
        >>> a = FieldArray.from_kwargs(mass1=1.1, mass2=2.)
        >>> a.mass1, a.mass2
        (array([ 1.1]), array([ 2.]))
        """
        arrays = []
        names = []
        for p,vals in kwargs.items():
            if not isinstance(vals, numpy.ndarray):
                if not isinstance(vals, list):
                    vals = [vals]
                vals = numpy.array(vals)
            arrays.append(vals)
            names.append(p)
        return cls.from_arrays(arrays, names=names)


    @classmethod
    def from_ligolw_table(cls, table, columns=None, cast_to_dtypes=None):
        """Converts the given ligolw table into an FieldArray. The `tableName`
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
        array : FieldArray
            The input table as an FieldArray.
        """
        name = table.tableName.split(':')[0]
        if columns is None:
            # get all the columns
            columns = table.validcolumns
        else:
            # note: this will raise a KeyError if one or more columns is
            # not in the table's validcolumns
            new_columns = {}
            for col in columns:
                new_columns[col] = table.validcolumns[col]
            columns = new_columns
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

    def to_array(self, fields=None, axis=0):
        """Returns an `numpy.ndarray` of self in which the fields are included
        as an extra dimension.

        Parameters
        ----------
        fields : {None, (list of) strings}
            The fields to get. All of the fields must have the same datatype.
            If None, will try to return all of the fields.
        axis : {0, int}
            Which dimension to put the fields in in the returned array. For
            example, if `self` has shape `(l,m,n)` and `k` fields, the
            returned array will have shape `(k,l,m,n)` if `axis=0`, `(l,k,m,n)`
            if `axis=1`, etc. Setting `axis=-1` will put the fields in the
            last dimension. Default is 0.
        
        Returns
        -------
        numpy.ndarray
            The desired fields as a numpy array.
        """
        if fields is None:
            fields = self.fieldnames
        if isinstance(fields, (str, unicode)):
            fields = [fields]
        return numpy.stack([self[f] for f in fields], axis=axis)

    @property
    def fieldnames(self):
        """Returns a tuple listing the field names in self. Equivalent to
        `array.dtype.names`, where `array` is self.
        """
        return self.dtype.names

    @property
    def virtualfields(self):
        """Returns a tuple listing the names of virtual fields in self.
        """
        if self._virtualfields is None:
            vfs = tuple()
        else:
            vfs = tuple(self._virtualfields)
        return vfs

    @property
    def functionlib(self):
        """Returns the library of functions that are available when calling
        items.
        """
        return self._functionlib

    @property
    def fields(self):
        """Returns a tuple listing the names of fields and virtual fields in
        self."""
        return tuple(list(self.fieldnames) + list(self.virtualfields))

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

    def add_fields(self, arrays, names=None, assubarray=False):
        """
        Adds the given arrays as new fields to self.  Returns a new instance
        with the new fields added. Note: this array does not change; the
        returned array is a new copy. 

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

    def parse_boolargs(self, args):
        """Returns an array populated by given values, with the indices of
        those values dependent on given boolen tests on self.
        
        The given `args` should be a list of tuples, with the first element the
        return value and the second argument a string that evaluates to either
        True or False for each element in self.

        Each boolean argument is evaluated on elements for which every prior
        boolean argument was False. For example, if array `foo` has a field
        `bar`, and `args = [(1, 'bar < 10'), (2, 'bar < 20'), (3, 'bar < 30')]`,
        then the returned array will have `1`s at the indices for
        which `foo.bar < 10`, `2`s where `foo.bar < 20 and not foo.bar < 10`,
        and `3`s where `foo.bar < 30 and not (foo.bar < 10 or foo.bar < 20)`.
        
        The last argument in the list may have "else", an empty string, None,
        or simply list a return value. In any of these cases, any element not
        yet populated will be assigned the last return value.
        
        Parameters
        ----------
        args : {(list of) tuples, value}
            One or more return values and boolean argument determining where
            they should go.

        Returns
        -------
        return_values : array
            An array with length equal to self, with values populated with the
            return values.
        leftover_indices : array
            An array of indices that evaluated to False for all arguments.
            These indices will not have been popluated with any value,
            defaulting to whatever numpy uses for a zero for the return
            values' dtype. If there are no leftovers, an empty array is
            returned.

        Examples
        --------
        Given the following array:

        >>> arr = FieldArray(5, dtype=[('mtotal', float)])
        >>> arr['mtotal'] = numpy.array([3., 5., 2., 1., 4.])

        Return `"TaylorF2"` for all elements with `mtotal < 4` (note that the
        elements 1 and 4 are leftover):

        >>> arr.parse_boolargs(('TaylorF2', 'mtotal<4'))
            (array(['TaylorF2', '', 'TaylorF2', 'TaylorF2', ''], 
            dtype='|S8'),
            array([1, 4]))

        Return `"TaylorF2"` for all elements with `mtotal < 4`,
        `"SEOBNR_ROM_DoubleSpin"` otherwise:

        >>> arr.parse_boolargs([('TaylorF2', 'mtotal<4'), ('SEOBNRv2_ROM_DoubleSpin', 'else')])
            (array(['TaylorF2', 'SEOBNRv2_ROM_DoubleSpin', 'TaylorF2', 'TaylorF2',
            'SEOBNRv2_ROM_DoubleSpin'], 
            dtype='|S23'),
            array([], dtype=int64))
        
        The following will also return the same:

        >>> arr.parse_boolargs([('TaylorF2', 'mtotal<4'), ('SEOBNRv2_ROM_DoubleSpin',)])
        >>> arr.parse_boolargs([('TaylorF2', 'mtotal<4'), ('SEOBNRv2_ROM_DoubleSpin', '')])
        >>> arr.parse_boolargs([('TaylorF2', 'mtotal<4'), 'SEOBNRv2_ROM_DoubleSpin'])

        Return `"TaylorF2"` for all elements with `mtotal < 3`, `"IMRPhenomD"`
        for all elements with `3 <= mtotal < 4`, `"SEOBNRv2_ROM_DoubleSpin"`
        otherwise:

        >>> arr.parse_boolargs([('TaylorF2', 'mtotal<3'), ('IMRPhenomD', 'mtotal<4'), 'SEOBNRv2_ROM_DoubleSpin'])
            (array(['IMRPhenomD', 'SEOBNRv2_ROM_DoubleSpin', 'TaylorF2', 'TaylorF2',
            'SEOBNRv2_ROM_DoubleSpin'], 
            dtype='|S23'),
            array([], dtype=int64))

        Just return `"TaylorF2"` for all elements:

        >>> arr.parse_boolargs('TaylorF2')
            (array(['TaylorF2', 'TaylorF2', 'TaylorF2', 'TaylorF2', 'TaylorF2'], 
            dtype='|S8'),
            array([], dtype=int64))

        """
        if not isinstance(args, list):
            args = [args]
        # format the arguments
        return_vals = []
        bool_args = []
        for arg in args:
            if not isinstance(arg, tuple):
                return_val = arg
                bool_arg = None
            elif len(arg) == 1:
                return_val = arg[0]
                bool_arg = None
            elif len(arg) == 2:
                return_val, bool_arg = arg
            else:
                raise ValueError("argument not formatted correctly")
            return_vals.append(return_val)
            bool_args.append(bool_arg)
        # get the output dtype
        outdtype = numpy.array(return_vals).dtype
        out = numpy.zeros(self.size, dtype=outdtype)
        mask = numpy.zeros(self.size, dtype=bool)
        leftovers = numpy.ones(self.size, dtype=bool)
        for ii,(boolarg,val) in enumerate(zip(bool_args, return_vals)):
            if boolarg is None or boolarg == '' or boolarg.lower() == 'else': 
                if ii+1 != len(bool_args):
                    raise ValueError("only the last item may not provide "
                        "any boolean arguments")
                mask = leftovers
            else:
                mask = leftovers & self[boolarg] 
            out[mask] = val
            leftovers &= ~mask
        return out, numpy.where(leftovers)[0]

    def append(self, other):
        """Appends another array to this array.

        The returned array will have all of the class methods and virutal
        fields of this array, including any that were added using `add_method`
        or `add_virtualfield`. If this array and other array have one or more
        string fields, the dtype for those fields are updated to a string
        length that can encompass the longest string in both arrays.

        .. note::
            Increasing the length of strings only works for fields, not
            sub-fields.

        Parameters
        ----------
        other : array
            The array to append values from. It must have the same fields and
            dtype as this array, modulo the length of strings. If the other
            array does not have the same dtype, a TypeError is raised.

        Returns
        -------
        array
            An array with others values appended to this array's values. The
            returned array is an instance of the same class as this array,
            including all methods and virtual fields.
        """
        try:
            return numpy.append(self, other).view(type=self.__class__)
        except TypeError:
            # see if the dtype error was due to string fields having different
            # lengths; if so, we'll make the joint field the larger of the
            # two
            str_fields = [name for name in self.fieldnames
                          if _isstring(self.dtype[name])]
            # get the larger of the two
            new_strlens = dict(
                [[name,
                  max(self.dtype[name].itemsize, other.dtype[name].itemsize)]
                 for name in str_fields]
            )
            # cast both to the new string lengths
            new_dt = []
            for dt in self.dtype.descr:
                name = dt[0]
                if name in new_strlens:
                    dt = (name, self.dtype[name].type, new_strlens[name])
                new_dt.append(dt)
            new_dt = numpy.dtype(new_dt)
            return numpy.append(
                self.astype(new_dt),
                other.astype(new_dt)
                ).view(type=self.__class__)

    @classmethod
    def parse_parameters(cls, parameters, possible_fields):
        """Parses a list of parameters to get the list of fields needed in
        order to evaluate those parameters.

        Parameters
        ----------
        parameters : (list of) string(s)
            The list of desired parameters. These can be (functions of) fields
            or virtual fields.
        possible_fields : (list of) string(s)
            The list of possible fields.

        Returns
        -------
        list :
            The list of names of the fields that are needed in order to
            evaluate the given parameters.
        """
        if isinstance(possible_fields, (str, unicode)):
            possible_fields = [possible_fields]
        possible_fields = map(str, possible_fields)
        # we'll just use float as the dtype, as we just need this for names
        arr = cls(1, dtype=zip(possible_fields,
                               len(possible_fields)*[float]))
        # try to perserve order
        return list(get_needed_fieldnames(arr, parameters))

def _isstring(dtype):
    """Given a numpy dtype, determines whether it is a string. Returns True
    if the dtype is string or unicode.
    """
    return dtype.type == numpy.unicode_ or dtype.type == numpy.string_


def aliases_from_fields(fields):
    """Given a dictionary of fields, will return a dictionary mapping the
    aliases to the names.
    """
    return dict(c for c in fields if isinstance(c, tuple))


def fields_from_names(fields, names=None):
    """Given a dictionary of fields and a list of names, will return a
    dictionary consisting of the fields specified by names. Names can be
    either the names of fields, or their aliases.
    """

    if names is None:
        return fields
    if isinstance(names, (str, unicode)):
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


#
# =============================================================================
#
#                           FieldArrays with default fields
#
# =============================================================================
#

class _FieldArrayWithDefaults(FieldArray):
    """
    Subclasses FieldArray, adding class attribute ``_staticfields``, and
    class method ``default_fields``. The ``_staticfields`` should be a
    dictionary that defines some field names and corresponding dtype. The
    ``default_fields`` method returns a dictionary of the static fields
    and any default virtualfields that were added. A field array can then
    be initialized in one of 3 ways:

     1. With just a shape. In this case, the returned array will have all
     of the default fields.

     2. With a shape and a list of names, given by the ``names`` keyword
     argument. The names may be default fields, virtual fields, a method or
     property of the class, or any python function of these things. If a
     virtual field, method, or property is in the names, the needed underlying
     fields will be included in the return array. For example, if the class
     has a virtual field called 'mchirp', which is a function of fields called
     'mass1' and 'mass2', then 'mchirp' or any function of 'mchirp' may be
     included in the list of names (e.g., names=['mchirp**(5/6)']). If so, the
     returned array will have fields 'mass1' and 'mass2' even if these were
     not specified in names, so that 'mchirp' may be used without error.
     names must be names of either default fields or virtualfields, else a
     KeyError is raised.

     3. With a shape and a dtype. Any field specified by the dtype will be
     used. The fields need not be in the list of default fields, and/or the
     dtype can be different than that specified by the default fields.

    If additional fields are desired beyond the default fields, these can
    be specified using the ``additional_fields`` keyword argument; these should
    be provided in the same way as ``dtype``; i.e, as a list of (name, dtype)
    tuples.

    This class does not define any static fields, and ``default_fields`` just
    returns an empty dictionary. This class is mostly meant to be subclassed
    by other classes, so they can add their own defaults.
    """

    _staticfields = {}
    @classmethod
    def default_fields(cls, include_virtual=True, **kwargs):
        """The default fields and their dtypes. By default, this returns
        whatever the class's ``_staticfields`` and ``_virtualfields`` is set
        to as a dictionary of fieldname, dtype (the dtype of virtualfields is
        given by VIRTUALFIELD_DTYPE). This function should be overridden by
        subclasses to add dynamic fields; i.e., fields that require some input
        parameters at initialization. Keyword arguments can be passed to this
        to set such dynamic fields.
        """
        add_fields = {}
        if include_virtual:
            add_fields.update(dict([[name, VIRTUALFIELD_DTYPE]
                for name in cls._virtualfields]))
        return dict(cls._staticfields.items() + add_fields.items())
        

    def __new__(cls, shape, name=None, additional_fields=None,
                field_kwargs=None, **kwargs):
        """The ``additional_fields`` should be specified in the same way as
        ``dtype`` is normally given to FieldArray. The ``field_kwargs`` are
        passed to the class's default_fields method as keyword arguments.
        """
        if field_kwargs is None:
            field_kwargs = {}
        if 'names' in kwargs and 'dtype' in kwargs:
            raise ValueError("Please provide names or dtype, not both")
        default_fields = cls.default_fields(include_virtual=False,
            **field_kwargs)
        if 'names' in kwargs:
            names = kwargs.pop('names')
            if isinstance(names, (str, unicode)):
                names = [names]
            # evaluate the names to figure out what base fields are needed
            # to do this, we'll create a small default instance of self (since
            # no names are specified in the following initialization, this
            # block of code is skipped)
            arr = cls(1, field_kwargs=field_kwargs)
            # try to perserve order
            sortdict = dict([[nm, ii] for ii,nm in enumerate(names)])
            names = list(get_needed_fieldnames(arr, names))
            names.sort(key=lambda x: sortdict[x] if x in sortdict
                else len(names))
            # add the fields as the dtype argument for initializing 
            kwargs['dtype'] = [(fld, default_fields[fld]) for fld in names]
        if 'dtype' not in kwargs:
            kwargs['dtype'] = default_fields.items()
        # add the additional fields
        if additional_fields is not None:
            if not isinstance(additional_fields, list):
                additional_fields = [additional_fields]
            if not isinstance(kwargs['dtype'], list):
                kwargs['dtype'] = [kwargs['dtype']]
            kwargs['dtype'] += additional_fields
        return super(_FieldArrayWithDefaults, cls).__new__(cls, shape,
            name=name, **kwargs)

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
        if isinstance(names, (str, unicode)):
            names = [names]
        default_fields = self.default_fields(include_virtual=False, **kwargs)
        # parse out any virtual fields
        arr = self.__class__(1, field_kwargs=kwargs)
        # try to perserve order
        sortdict = dict([[nm, ii] for ii,nm in enumerate(names)])
        names = list(get_needed_fieldnames(arr, names))
        names.sort(key=lambda x: sortdict[x] if x in sortdict
            else len(names))
        fields = [(name, default_fields[name]) for name in names]
        arrays = []
        names = []
        for name,dt in fields:
            arrays.append(default_empty(self.size, dtype=[(name, dt)])) 
            names.append(name)
        return self.add_fields(arrays, names)

    @classmethod
    def parse_parameters(cls, parameters, possible_fields=None):
        """Parses a list of parameters to get the list of fields needed in
        order to evaluate those parameters.

        Parameters
        ----------
        parameters : (list of) strings
            The list of desired parameters. These can be (functions of) fields
            or virtual fields.
        possible_fields : {None, dict}
            Specify the list of possible fields. Must be a dictionary given
            the names, and dtype of each possible field. If None, will use this
            class's `_staticfields`.

        Returns
        -------
        list :
            The list of names of the fields that are needed in order to
            evaluate the given parameters.
        """
        if possible_fields is not None:
            # make sure field names are strings and not unicode
            possible_fields = dict([[f, dt]
                for f,dt in possible_fields.items()])
            class ModifiedArray(cls):
                _staticfields = possible_fields
            cls = ModifiedArray
        return cls(1, names=parameters).fieldnames

#
# =============================================================================
#
#                           WaveformArray
#
# =============================================================================
#

class WaveformArray(_FieldArrayWithDefaults):
    """
    A FieldArray with some default fields and properties commonly used
    by CBC waveforms. This may be initialized in one of 3 ways:
    
    1. With just the size of the array. In this case, the returned array will
    have all of the default field names. Example:

    >>> warr = WaveformArray(10)
    >>> warr.fieldnames
        ('distance',
         'spin2x',
         'mass1',
         'mass2',
         'lambda1',
         'polarization',
         'spin2y',
         'spin2z',
         'spin1y',
         'spin1x',
         'spin1z',
         'inclination',
         'coa_phase',
         'dec',
         'tc',
         'lambda2',
         'ra')

    2. With some subset of the default field names. Example:

    >>> warr = WaveformArray(10, names=['mass1', 'mass2'])
    >>> warr.fieldnames
        ('mass1', 'mass2')

    The list of names may include virtual fields, and methods, as well as
    functions of these. If one or more virtual fields or methods are specified,
    the source code is analyzed to pull out whatever underlying fields are
    needed. Example:

    >>> warr = WaveformArray(10, names=['mchirp**(5/6)', 'chi_eff', 'cos(coa_phase)'])
    >>> warr.fieldnames
        ('spin2z', 'mass1', 'mass2', 'coa_phase', 'spin1z')

    3. By specifying a dtype. In this case, only the provided fields will
    be used, even if they are not in the defaults. Example:

    >>> warr = WaveformArray(10, dtype=[('foo', float)])
    >>> warr.fieldnames
        ('foo',)    

    Additional fields can also be specified using the additional_fields
    keyword argument. Example:

    >>> warr = WaveformArray(10, names=['mass1', 'mass2'], additional_fields=[('bar', float)])
    >>> warr.fieldnames
        ('mass1', 'mass2', 'bar')

    .. note::
        If an array is initialized with all of the default fields (case 1,
        above), then the names come from waveform.parameters; i.e., they
        are actually Parameter instances, not just strings. This means that the
        field names carry all of the metadata that a Parameter has. For
        example:

        >>> warr = WaveformArray(10)
        >>> warr.fields[0]
            'distance'
        >>> warr.fields[0].description
            'Luminosity distance to the binary (in Mpc).'
        >>> warr.fields[0].label
            '$d_L$ (Mpc)'

    """
    _staticfields = (parameters.cbc_intrinsic_params +
                     parameters.extrinsic_params).dtype_dict

    _virtualfields = [
        parameters.mchirp, parameters.eta, parameters.mtotal,
        parameters.q, parameters.primary_mass, parameters.secondary_mass,
        parameters.chi_eff,
        parameters.spin_px, parameters.spin_py, parameters.spin_pz,
        parameters.spin_sx, parameters.spin_sy, parameters.spin_sz,
        parameters.spin1_a, parameters.spin1_azimuthal, parameters.spin1_polar,
        parameters.spin2_a, parameters.spin2_azimuthal, parameters.spin2_polar]

    @property
    def primary_mass(self):
        """Returns the larger of self.mass1 and self.mass2."""
        return conversions.primary_mass(self.mass1, self.mass2)

    @property
    def secondary_mass(self):
        """Returns the smaller of self.mass1 and self.mass2."""
        return conversions.secondary_mass(self.mass1, self.mass)

    @property
    def mtotal(self):
        """Returns the total mass."""
        return conversions.mtotal_from_mass1_mass2(self.mass1, self.mass2)

    @property
    def q(self):
        """Returns the mass ratio m1/m2, where m1 >= m2."""
        return conversions.q_from_mass1_mass2(self.mass1, self.mass2)

    @property
    def eta(self):
        """Returns the symmetric mass ratio."""
        return conversions.eta_from_mass1_mass2(self.mass1, self.mass2)

    @property
    def mchirp(self):
        """Returns the chirp mass."""
        return conversions.mchirp_from_mass1_mass2(self.mass1, self.mass2)

    @property
    def chi_eff(self):
        """Returns the effective spin."""
        return conversions.chi_eff(self.mass1, self.mass2, self.spin1z,
                                   self.spin2z)

    @property
    def spin_px(self):
        """Returns the x-component of the spin of the primary mass."""
        return conversions.primary_spin(self.mass1, self.mass2, self.spin1x,
                                        self.spin2x)

    @property
    def spin_py(self):
        """Returns the y-component of the spin of the primary mass."""
        return conversions.primary_spin(self.mass1, self.mass2, self.spin1y,
                                        self.spin2y)

    @property
    def spin_pz(self):
        """Returns the z-component of the spin of the primary mass."""
        return conversions.primary_spin(self.mass1, self.mass2, self.spin1z,
                                        self.spin2z)

    @property
    def spin_sx(self):
        """Returns the x-component of the spin of the secondary mass."""
        return conversions.secondary_spin(self.mass1, self.mass2, self.spin1x,
                                        self.spin2x)

    @property
    def spin_sy(self):
        """Returns the y-component of the spin of the secondary mass."""
        return conversions.secondary_spin(self.mass1, self.mass2, self.spin1y,
                                        self.spin2y)

    @property
    def spin_sz(self):
        """Returns the z-component of the spin of the secondary mass."""
        return conversions.secondary_spin(self.mass1, self.mass2, self.spin1z,
                                        self.spin2z)

    @property
    def spin1_a(self):
        """Returns the dimensionless spin magnitude of mass 1."""
        return coordinates.cartesian_to_spherical_rho(
                                    self.spin1x, self.spin1y, self.spin1z)

    @property
    def spin1_azimuthal(self):
        """Returns the azimuthal spin angle of mass 1."""
        return coordinates.cartesian_to_spherical_azimuthal(
                                     self.spin1x, self.spin1y)

    @property
    def spin1_polar(self):
        """Returns the polar spin angle of mass 1."""
        return coordinates.cartesian_to_spherical_polar(
                                     self.spin1x, self.spin1y, self.spin1z)

    @property
    def spin2_a(self):
        """Returns the dimensionless spin magnitude of mass 2."""
        return coordinates.cartesian_to_spherical_rho(
                                    self.spin1x, self.spin1y, self.spin1z)

    @property
    def spin2_azimuthal(self):
        """Returns the azimuthal spin angle of mass 2."""
        return coordinates.cartesian_to_spherical_azimuthal(
                                     self.spin2x, self.spin2y)

    @property
    def spin2_polar(self):
        """Returns the polar spin angle of mass 2."""
        return coordinates.cartesian_to_spherical_polar(
                                     self.spin2x, self.spin2y, self.spin2z)


__all__ = ['FieldArray', 'WaveformArray']

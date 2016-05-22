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
FieldArrays are wrappers of numpy recarrays with additional functionality useful
for storing and retrieving data created by a search for gravitational waves.
"""

import os, sys, types, re, copy, numpy
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


#
# =============================================================================
#
#                           Base FieldArray definitions
#
# =============================================================================
#
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
    individual fields. For example, if you have an FieldArray ``x`` with fields
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
    FieldArray, those can be accessed in the same way as fields are. For example,
    define ``Foo`` as:

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

    Create an array with subfields:

    >>> x = FieldArray(4, dtype=[('foo', [('cat', float), ('hat', int)]), ('bar', float)])
    >>> x.all_fieldnames
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
    __persistent_attributes__ = ['name', 'id_maps']

    def __new__(cls, shape, name=None, zero=True, **kwargs):
        """Initializes a new empty array.
        """
        obj = super(FieldArray, cls).__new__(cls, shape, **kwargs).view(
            type=cls)
        obj.name = name
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


    def __getsubitem__(self, item):
        """Gets a subfield using `field.subfield` notation.
        """
        try:
            return super(FieldArray, self).__getitem__(item)
        except ValueError as err:
            subitems = item.split('.')
            if len(subitems) > 1:
                return super(FieldArray, self).__getitem__(subitems[0]
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
                super(FieldArray, self).__getitem__(fn)] \
                for fn in set(self.fieldnames).intersection(itemvars)])

            # add any aliases
            for alias, name in self.aliases.items():
                if name in item_dict:
                    item_dict[alias] = item_dict[name]

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
            A list of the tuples to create the FieldArray from.
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

__all__ = ['FieldArray']

# Copyright (C) 2006,2012  Kipp Cannon
#
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
#                                    ILWDs
#
# =============================================================================
#


"""
The ilwd:char type is used to store ID strings for objects within LIGO
Light-Weight XML files.  This module and its associated C extention module
_ilwd provide a class for memory-efficient storage of ilwd:char strings.

LIGO Light Weight XML "ilwd:char" IDs are strings of the form
"table:column:integer", for example "process:process_id:10".  Large complex
documents can have many millions of these strings, and their storage
represents a significant RAM burden.  However, while there can be millions
of ID strings in a document there might be only a small number (e.g., 10 or
fewer) unique ID prefixes in a document (the table name and column name
part).  The amount of RAM required to load a document can be significantly
reduced if the small number of unique string prefixes are stored separately
and reused.  This module provides the machinery used to do this.

The ilwdchar class in this module converts a string or unicode object
containing an ilwd:char ID into a more memory efficient representation.

Example:

>>> x = ilwdchar("process:process_id:10")
>>> print x
process:process_id:10

Like strings, the object resulting from this is immutable.  It provides two
read-only attributes, "table_name" and "column_name", that can be used to
access the table and column parts of the original ID string.  The integer
suffix can be retrieved by converting the object to an integer.

Example:

>>> x.table_name
'process'
>>> int(x)
10

The object also provides the read-only attribute "index_offset", giving the
length of the string preceding the interger suffix.

Example:

>>> x.index_offset
19

The objects support some arithmetic operations.

Example:

>>> y = x + 5
>>> str(y)
'process:process_id:15'
>>> int(y - x)
5

The objects are pickle-able.

Example:

>>> import pickle
>>> x == pickle.loads(pickle.dumps(x))
True

To simplify interaction with documents that do not contain fully-populated
columns, None is allowed as an input value and is not converted.

Example:

>>> print ilwdchar(None)
None


Implementation details
======================

Memory is reduced by storing the table_name, column_name, and index_offset
values as class attributes, so only one copy is present in memory and is
shared across all instances of the class.  This means that each unique
table_name and column_name pair requires its own class.  These classes are
created on the fly as new IDs are processed, and get added to this module's
name space.  They are all subclasses of _ilwd.ilwdchar, which implements
the low-level machinery.  After a new class is created it can be accessed
as a symbol in this module, but each of those symbols does not exist until
at least one corresponding ID string has been processed.

Example:

>>> import ilwd
>>> "foo_bar_class" in ilwd.__dict__
False
>>> x = ilwd.ilwdchar("foo:bar:0")
>>> type(x)
<class 'pycbc.ligolw.ilwd.foo_bar_class'>
>>> "foo_bar_class" in ilwd.__dict__
True
>>> print ilwd.foo_bar_class(10)
foo:bar:10

The ilwdchar class itself is never instantiated, its .__new__() method
parses the ID string parameter and creates an instance of the appropriate
subclass of _ilwd.ilwdchar, creating a new subclass before doing so if
neccessary.
"""


import copy_reg


from pycbc import version
from . import _ilwd


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % version.version
__date__ = version.date


#
# =============================================================================
#
#                                Cached Classes
#
# =============================================================================
#


#
# Function for retrieving ilwdchar subclasses.
#


def get_ilwdchar_class(tbl_name, col_name, namespace = globals()):
	"""
	Searches this module's namespace for a subclass of _ilwd.ilwdchar
	whose table_name and column_name attributes match those provided.
	If a matching subclass is found it is returned; otherwise a new
	class is defined, added to this module's namespace, and returned.

	Example:

	>>> process_id = get_ilwdchar_class("process", "process_id")
	>>> x = process_id(10)
	>>> str(type(x))
	"<class 'pycbc.ligolw.ilwd.process_process_id_class'>"
	>>> str(x)
	'process:process_id:10'

	Retrieving and storing the class provides a convenient mechanism
	for quickly constructing new ID objects.

	Example:

	>>> for i in range(10):
	...	print str(process_id(i))
	...
	process:process_id:0
	process:process_id:1
	process:process_id:2
	process:process_id:3
	process:process_id:4
	process:process_id:5
	process:process_id:6
	process:process_id:7
	process:process_id:8
	process:process_id:9
	"""
	#
	# if the class already exists, retrieve and return it
	#

	key = (str(tbl_name), str(col_name))
	cls_name = "%s_%s_class" % key
	assert cls_name != "get_ilwdchar_class"
	try:
		return namespace[cls_name]
	except KeyError:
		pass

	#
	# otherwise define a new class, and add it to the cache
	#

	class new_class(_ilwd.ilwdchar):
		__slots__ = ()
		table_name, column_name = key
		index_offset = len("%s:%s:" % key)

	new_class.__name__ = cls_name

	namespace[cls_name] = new_class

	#
	# pickle support
	#

	copy_reg.pickle(new_class, lambda x: (ilwdchar, (unicode(x),)))

	#
	# return the new class
	#

	return new_class


#
# Metaclass to redirect instantiation to the correct subclass for
# _ilwd.ilwdchar
#


class ilwdchar(object):
	"""
	Metaclass wrapper of pycbc.ligolw._ilwd.ilwdchar class.
	Instantiating this class constructs and returns an instance of a
	subclass of pycbc.ligolw._ilwd.ilwdchar.
	"""
	def __new__(cls, s):
		"""
		Convert an ilwd:char-formated string into an instance of
		the matching subclass of _ilwd.ilwdchar.  If the input is
		None then the return value is None.

		Example:

		>>> x = ilwdchar("process:process_id:10")
		>>> str(x)
		'process:process_id:10'
		>>> x.table_name
		'process'
		>>> x.column_name
		'process_id'
		>>> int(x)
		10
		>>> x.index_offset
		19
		>>> str(x)[x.index_offset:]
		'10'
		>>> print ilwdchar(None)
		None
		"""
		#
		# None is no-op
		#

		if s is None:
			return None

		#
		# try parsing the string as an ilwd:char formated string
		#

		try:
			table_name, column_name, i = s.strip().split(":")
		except (ValueError, AttributeError):
			raise ValueError("invalid ilwd:char '%s'" % repr(s))

		#
		# retrieve the matching class from the ID class cache, and
		# return an instance initialized to the desired value
		#

		return get_ilwdchar_class(table_name, column_name)(int(i))

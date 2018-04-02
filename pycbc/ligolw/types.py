# Copyright (C) 2006--2011  Kipp Cannon
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
#                                   Preamble
#
# =============================================================================
#


"""
Definitions of type strings found in LIGO Light Weight XML files.

Notes.  To guarantee that a double-precision floating-point number can be
reconstructed exactly from its representation as a decimal number, one must
use 17 decimal digits;  for single-precision, the number is 9.  Python uses
only double-precision numbers, but LIGO Light Weight XML allows for
single-precision values, so I provide distinct format specifiers for those
cases here.  In both cases, I have elected to use 1 fewer digits than are
required to uniquely reconstruct the number:  the XML written by this
library is lossy.  I made this choice to reduce the file size, for example

>>> "%.17g" % 0.1
'0.10000000000000001'

while

>>> "%.16g" % 0.1
'0.1'

In this worst case, storing full precision increases the size of the XML by
more than an order of magnitude.  If you wish to make a different choice
for your files, for example if you wish your XML files to be lossless,
simply include the lines

	pycbc.ligolw.types.FormatFunc.update({
		"real_4": u"%.9g".__mod__,
		"real_8": u"%.17g".__mod__,
		"float": u"%.9g".__mod__,
		"double": u"%.17g".__mod__,
		u"complex_8": pycbc.ligolw.types.mk_complex_format_func(u"%.9g"),
		u"complex_16": pycbc.ligolw.types.mk_complex_format_func(u"%.17g")
	})

anywhere in your code, but before you write the document to a file.

References:

	- http://docs.sun.com/source/806-3568/ncg_goldberg.html
"""


import base64


from pycbc import version
from . import ilwd


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % version.version
__date__ = version.date


#
# =============================================================================
#
#                               Type Categories
#
# =============================================================================
#


IDTypes = set([u"ilwd:char", u"ilwd:char_u"])
"""LIGO Light-Weight XML type strings for ID-like data."""

BlobTypes = set([u"blob", u"ilwd:char_u"])
"""LIGO Light-Weight XML type strings for binary blob-like data."""

StringTypes = set([u"char_s", u"char_v", u"lstring", u"string", u"ilwd:char"])
"""LIGO Light-Weight XML type strings for string-like data."""

IntTypes = set([u"int_2s", u"int_2u", u"int_4s", u"int_4u", u"int_8s", u"int_8u", u"int"])
"""LIGO Light-Weight XML type strings for integer-like data."""

FloatTypes = set([u"real_4", u"real_8", u"float", u"double"])
"""LIGO Light-Weight XML type strings for floating-point-like data."""

ComplexTypes = set([u"complex_8", u"complex_16"])
"""LIGO Light-Weight XML type strings for complex-like data."""

NumericTypes = IntTypes | FloatTypes | ComplexTypes
"""LIGO Light-Weight XML type strings for number-like data."""

TimeTypes = set([u"GPS", u"Unix", u"ISO-8601"])
"""LIGO Light-Weight XML type strings for time-like data."""

Types = BlobTypes | StringTypes | NumericTypes | TimeTypes
"""All valid LIGO Light-Weight XML type strings."""


#
# =============================================================================
#
#                         Output Format Look-up Table
#
# =============================================================================
#


def string_format_func(s):
	"""
	Function used internally to format string data for output to XML.
	Escapes back-slashes and quotes, and wraps the resulting string in
	quotes.
	"""
	return u"\"%s\"" % unicode(s).replace(u"\\", u"\\\\").replace(u"\"", u"\\\"")


def blob_format_func(b):
	"""
	Function used internally to format binary data.  Base64-encodes the
	data and wraps the resulting string in quotes.
	"""
	return u"\"%s\"" % base64.standard_b64encode(b)


def mk_complex_format_func(fmt):
	"""
	Function used internally to generate functions to format complex
	valued data.
	"""
	fmt = fmt + u"+i" + fmt
	def complex_format_func(z):
		return fmt % (z.real, z.imag)
	return complex_format_func


FormatFunc = {
	u"char_s": string_format_func,
	u"char_v": string_format_func,
	u"ilwd:char": u"\"%s\"".__mod__,
	u"ilwd:char_u": blob_format_func,
	u"blob": blob_format_func,
	u"lstring": string_format_func,
	u"string": string_format_func,
	u"int_2s": u"%d".__mod__,
	u"int_2u": u"%u".__mod__,
	u"int_4s": u"%d".__mod__,
	u"int_4u": u"%u".__mod__,
	u"int_8s": u"%d".__mod__,
	u"int_8u": u"%u".__mod__,
	u"int": u"%d".__mod__,
	u"real_4": u"%.8g".__mod__,
	u"real_8": u"%.16g".__mod__,
	u"float": u"%.8g".__mod__,
	u"double": u"%.16g".__mod__,
	u"complex_8": mk_complex_format_func(u"%.8g"),
	u"complex_16": mk_complex_format_func(u"%.16g")
}
"""
Look-up table mapping LIGO Light-Weight XML data type strings to functions
for formating Python data for output.  This table is used universally by
glue.ligolw XML writing codes.
"""


#
# =============================================================================
#
#                  Conversion To And From Native Python Types
#
# =============================================================================
#


ToPyType = {
	u"char_s": unicode,
	u"char_v": unicode,
	u"ilwd:char": ilwd.ilwdchar,
	u"ilwd:char_u": lambda s: buffer(base64.b64decode(s)),
	u"blob": lambda s: buffer(base64.b64decode(s)),
	u"lstring": unicode,
	u"string": unicode,
	u"int_2s": int,
	u"int_2u": int,
	u"int_4s": int,
	u"int_4u": int,
	u"int_8s": int,
	u"int_8u": int,
	u"int": int,
	u"real_4": float,
	u"real_8": float,
	u"float": float,
	u"double": float,
	u"complex_8": lambda s: complex(*map(float, s.split(u"+i"))),
	u"complex_16": lambda s: complex(*map(float, s.split(u"+i")))
}
"""
Look-up table mapping LIGO Light-Weight XML data type strings to functions
for parsing Python data from input.  This table is used universally by
glue.ligolw XML parsing codes.
"""


FromPyType = {
	ilwd.ilwdchar: u"ilwd:char",
	buffer: u"blob",
	str: u"lstring",
	unicode: u"lstring",
	bool: u"int_4s",
	int: u"int_8s",
	long: u"int_8s",
	float: u"real_8",
	complex: u"complex_16"
}
"""
Look-up table used to guess LIGO Light-Weight XML data type strings from
Python types.  This table is used when auto-generating XML from Python
objects.
"""


#
# =============================================================================
#
#                  Conversion To and From Native Numpy Types
#
# =============================================================================
#


ToNumPyType = {
	u"int_2s": "int16",
	u"int_2u": "uint16",
	u"int_4s": "int32",
	u"int_4u": "uint32",
	u"int_8s": "int64",
	u"int_8u": "uint64",
	u"int": "int32",
	u"real_4": "float32",
	u"real_8": "float64",
	u"float": "float32",
	u"double": "float64",
	u"complex_8": "complex64",
	u"complex_16": "complex128"
}
"""
Look-up table mapping LIGO Light-Weight XML data type strings to numpy
array type strings.  Used by glue.ligolw array reading codes.
"""


FromNumPyType = {
	"int16": u"int_2s",
	"uint16": u"int_2u",
	"int32": u"int_4s",
	"uint32": u"int_4u",
	"int64": u"int_8s",
	"uint64": u"int_8u",
	"float32": u"real_4",
	"float64": u"real_8",
	"complex64": u"complex_8",
	"complex128": u"complex_16"
}
"""
Look-up table mapping numpy array type strings to LIGO Light-Weight XML
data type strings.  Uesd by glue.ligolw array writing codes.
"""


#
# =============================================================================
#
#                 Conversion To and From Native Database Types
#
# =============================================================================
#


#
# SQL does not support complex numbers.  Documents containing
# complex-valued table columns cannot be stored in SQL databases at this
# time.
#


ToMySQLType = {
	u"char_s": "CHAR(20)",
	u"char_v": "VARCHAR(64)",
	u"ilwd:char": "VARCHAR(64)",
	u"ilwd:char_u": "BLOB",
	u"blob": "BLOB",
	u"lstring": "VARCHAR(255)",
	u"string": "VARCHAR(255)",
	u"int_2s": "SMALLINT",
	u"int_2u": "SMALLINT",
	u"int_4s": "INTEGER",
	u"int_4u": "INTEGER",
	u"int_8s": "BIGINT",
	u"int_8u": "BIGINT",
	u"int": "INTEGER",
	u"real_4": "FLOAT",
	u"real_8": "DOUBLE",
	u"float": "FLOAT",
	u"double": "DOUBLE"
}
"""
Look-up table mapping LIGO Light-Weight XML data type strings to MySQL
column types.  Used by XML --> MySQL conversion codes.
"""


ToSQLiteType = {
	u"char_s": "TEXT",
	u"char_v": "TEXT",
	u"ilwd:char": "TEXT",
	u"ilwd:char_u": "BLOB",
	u"blob": "BLOB",
	u"lstring": "TEXT",
	u"string": "TEXT",
	u"int_2s": "INTEGER",
	u"int_2u": "INTEGER",
	u"int_4s": "INTEGER",
	u"int_4u": "INTEGER",
	u"int_8s": "INTEGER",
	u"int_8u": "INTEGER",
	u"int": "INTEGER",
	u"real_4": "REAL",
	u"real_8": "REAL",
	u"float": "REAL",
	u"double": "REAL"
}
"""
Look-up table mapping LIGO Light-Weight XML data type strings to SQLite
column types.  Used by XML --> SQLite conversion codes.
"""


FromSQLiteType = {
	"BLOB": u"blob",
	"TEXT": u"lstring",
	"STRING": u"lstring",
	"INTEGER": u"int_4s",
	"REAL": u"real_8"
}
"""
Look-up table used to guess LIGO Light-Weight XML data type strings from
SQLite column types.  Used when auto-generating XML from the contents of an
SQLite database.
"""

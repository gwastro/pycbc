# Copyright (C) 2006--2015  Kipp Cannon
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
While the ligolw module provides classes and parser support for reading and
writing LIGO Light Weight XML documents, this module supplements that code
with classes and parsers that add intelligence to the in-RAM document
representation.

In particular, the document tree associated with a Table element is
enhanced.  During parsing, the Stream element in this module converts the
character data contained within it into a list of objects.  The list
contains one object for each row of the table, and the objects' attributes
are the names of the table's columns.  When the document is written out
again, the Stream element serializes the row objects back into character
data.

The Table element exports a list-like interface to the rows.  The Column
elements also provide list-like access to the values in the corresponding
columns of the table.
"""


import copy
import itertools
import re
import sys
import warnings
from xml.sax.saxutils import escape as xmlescape
from xml.sax.xmlreader import AttributesImpl


from pycbc import version
from . import ligolw
from . import tokenizer
from . import types as ligolwtypes


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % version.version
__date__ = version.date


#
# =============================================================================
#
#                           Column Name Manipulation
#
# =============================================================================
#


# Regular expression to extract the significant part of a column name
# according to the LIGO LW naming conventions.

# FIXME: the pattern should be
#
# r"(?:\A[a-z0-9_]+:|\A)(?P<FullName>(?:[a-z0-9_]+:|\A)(?P<Name>[a-z0-9_]+))\Z"
#
# but people are putting upper case letters in names!!!!!  Someone is going
# to get the beats.  There is a reason for requiring names to be all lower
# case:  SQL table and column names are case insensitive, therefore (i)
# when converting a document to SQL the columns "Rho" and "rho" would
# become indistinguishable and so it would be impossible to convert a
# document with columns having names like this into an SQL database;  and
# (ii) even if that degeneracy is not encountered the case cannot be
# preserved and so when converting back to XML the correct capitalization
# is lost.  Requiring names to be all lower-case creates the same
# degeneracies in XML representations that exist in SQL representations
# ensuring compatibility and defines the correct case to restore the names
# to when converting to XML.  Other rules can be imagined that would work
# as well, this is the one that got chosen.


ColumnPattern = re.compile(r"(?:\A\w+:|\A)(?P<FullName>(?:(?P<Table>\w+):|\A)(?P<Name>\w+))\Z")


def StripColumnName(name):
	"""
	Return the significant portion of a column name according to LIGO
	LW naming conventions.

	Example:

	>>> StripColumnName("process_params_group:process_params:program")
	'program'
	>>> StripColumnName("process_params:program")
	'program'
	>>> StripColumnName("program")
	'program'
	"""
	try:
		return ColumnPattern.search(name).group("Name")
	except AttributeError:
		return name


def CompareColumnNames(name1, name2):
	"""
	Convenience function to compare two column names according to LIGO
	LW naming conventions:  StripColumnName() is applied to both names
	and the results compared.

	Example:

	>>> # compare as equal
	>>> CompareColumnNames("process_params:program", "process:program")
	0
	>>> # not equal
	>>> CompareColumnNames("program", "start_time")
	-1

	Note that "process_params:program", "process:program" compare as
	equal because both columns are named "program" although they are
	from different tables.
	"""
	return cmp(StripColumnName(name1), StripColumnName(name2))


def getColumnsByName(elem, name):
	"""
	Return a list of Column elements named name under elem.  The name
	comparison is done with CompareColumnNames().
	"""
	name = StripColumnName(name)
	return elem.getElements(lambda e: (e.tagName == ligolw.Column.tagName) and (e.Name == name))


#
# =============================================================================
#
#                           Table Name Manipulation
#
# =============================================================================
#


# Regular expression used to extract the signifcant portion of a table or
# stream name, according to LIGO LW naming conventions.


TablePattern = re.compile(r"(?:\A[a-z0-9_]+:|\A)(?P<Name>[a-z0-9_]+):table\Z")


def StripTableName(name):
	"""
	Return the significant portion of a table name according to LIGO LW
	naming conventions.

	Example:

	>>> StripTableName("sngl_burst_group:sngl_burst:table")
	'sngl_burst'
	>>> StripTableName("sngl_burst:table")
	'sngl_burst'
	>>> StripTableName("sngl_burst")
	'sngl_burst'
	"""
	if name.lower() != name:
		warnings.warn("table name \"%s\" is not lower case" % name)
	try:
		return TablePattern.search(name).group("Name")
	except AttributeError:
		return name


def CompareTableNames(name1, name2):
	"""
	Convenience function to compare two table names according to LIGO
	LW naming conventions.  StripTableName() is applied to both names
	and the results compared.

	Example:

	>>> # compare as equal
	>>> CompareTableNames("sngl_inspiral:table", "sngl_inspiral")
	0
	>>> # not equal
	>>> CompareTableNames("sngl_burst_group:sngl_burst:table", "sngl_inspiral")
	-1
	"""
	return cmp(StripTableName(name1), StripTableName(name2))


def getTablesByName(elem, name):
	"""
	Return a list of Table elements named name under elem.  The name
	comparison is done using CompareTableNames().
	"""
	name = StripTableName(name)
	return elem.getElements(lambda e: (e.tagName == ligolw.Table.tagName) and (e.Name == name))


#
# =============================================================================
#
#                                  Utilities
#
# =============================================================================
#


def new_from_template(template):
	"""
	Deprecated legacy wrapper of .copy() method of Table instances.
	"""
	import warnings
	warnings.warn("pycbc.ligolw.table.new_from_template() is deprecated.  Use .copy() method of Table instances instead.", DeprecationWarning)
	return template.copy()


def get_table(xmldoc, name):
	"""
	Scan xmldoc for a Table element named name.  The comparison is done
	using CompareTableNames().  Raises ValueError if not exactly 1 such
	table is found.

	NOTE:  if a Table sub-class has its .tableName attribute set, then
	its .get_table() class method can be used instead.  This is true
	for all Table classes in the pycbc.ligolw.lsctables module, and it
	is recommended to always use the .get_table() class method of those
	classes to retrieve those standard tables instead of calling this
	function and passing the .tableName attribute.  The example below
	shows both techniques.

	Example:

	>>> import ligolw
	>>> import lsctables
	>>> xmldoc = ligolw.Document()
	>>> xmldoc.appendChild(ligolw.LIGO_LW()).appendChild(lsctables.New(lsctables.SnglInspiralTable))
	[]
	>>> # find table with this function
	>>> sngl_inspiral_table = get_table(xmldoc, lsctables.SnglInspiralTable.tableName)
	>>> # find table with .get_table() class method (preferred)
	>>> sngl_inspiral_table = lsctables.SnglInspiralTable.get_table(xmldoc)

	See also the .get_table() class method of the Table class.
	"""
	tables = getTablesByName(xmldoc, name)
	if len(tables) != 1:
		raise ValueError("document must contain exactly one %s table" % StripTableName(name))
	return tables[0]


def reassign_ids(elem):
	"""
	Recurses over all Table elements below elem whose next_id
	attributes are not None, and uses the .get_next_id() method of each
	of those Tables to generate and assign new IDs to their rows.  The
	modifications are recorded, and finally all ID attributes in all
	rows of all tables are updated to fix cross references to the
	modified IDs.

	This function is used by ligolw_add to assign new IDs to rows when
	merging documents in order to make sure there are no ID collisions.
	Using this function in this way requires the .get_next_id() methods
	of all Table elements to yield unused IDs, otherwise collisions
	will result anyway.  See the .sync_next_id() method of the Table
	class for a way to initialize the .next_id attributes so that
	collisions will not occur.

	Example:

	>>> import ligolw
	>>> import lsctables
	>>> xmldoc = ligolw.Document()
	>>> xmldoc.appendChild(ligolw.LIGO_LW()).appendChild(lsctables.New(lsctables.SnglInspiralTable))
	[]
	>>> reassign_ids(xmldoc)
	"""
	mapping = {}
	for tbl in elem.getElementsByTagName(ligolw.Table.tagName):
		if tbl.next_id is not None:
			tbl.updateKeyMapping(mapping)
	for tbl in elem.getElementsByTagName(ligolw.Table.tagName):
		tbl.applyKeyMapping(mapping)


def reset_next_ids(classes):
	"""
	For each class in the list, if the .next_id attribute is not None
	(meaning the table has an ID generator associated with it), set
	.next_id to 0.  This has the effect of reseting the ID generators,
	and is useful in applications that process multiple documents and
	add new rows to tables in those documents.  Calling this function
	between documents prevents new row IDs from growing continuously
	from document to document.  There is no need to do this, it's
	purpose is merely aesthetic, but it can be confusing to open a
	document and find process ID 300 in the process table and wonder
	what happened to the other 299 processes.

	Example:

	>>> import lsctables
	>>> reset_next_ids(lsctables.TableByName.values())
	"""
	for cls in classes:
		if cls.next_id is not None:
			cls.set_next_id(type(cls.next_id)(0))


#
# =============================================================================
#
#                                Column Element
#
# =============================================================================
#


class Column(ligolw.Column):
	"""
	High-level column element that provides list-like access to the
	values in a column.

	Example:

	>>> from xml.sax.xmlreader import AttributesImpl
	>>> import sys
	>>> tbl = Table(AttributesImpl({u"Name": u"test"}))
	>>> col = tbl.appendChild(Column(AttributesImpl({u"Name": u"test:snr", u"Type": u"real_8"})))
	>>> tbl.appendChild(TableStream(AttributesImpl({u"Name": u"test"})))	# doctest: +ELLIPSIS
	<pycbc.ligolw.table.TableStream object at ...>
	>>> tbl._update_column_info()
	>>> col.Name
	u'snr'
	>>> col.Type
	u'real_8'
	>>> # append 3 rows (with nothing in them)
	>>> tbl.append(tbl.RowType())
	>>> tbl.append(tbl.RowType())
	>>> tbl.append(tbl.RowType())
	>>> # assign values to the rows, in order, in this column
	>>> col[:] = [8.0, 10.0, 12.0]
	>>> col[:]
	[8.0, 10.0, 12.0]
	>>> col.asarray()
	array([  8.,  10.,  12.])
	>>> tbl.write(sys.stdout)	# doctest: +NORMALIZE_WHITESPACE
	<Table Name="test">
		<Column Type="real_8" Name="test:snr"/>
		<Stream Name="test">
			8,
			10,
			12
		</Stream>
	</Table>
	>>> col.index(10)
	1
	>>> col[0] = 9.
	>>> col[1] = 9.
	>>> col[2] = 9.
	>>> tbl.write(sys.stdout)		# doctest: +NORMALIZE_WHITESPACE
	<Table Name="test">
		<Column Type="real_8" Name="test:snr"/>
		<Stream Name="test">
			9,
			9,
			9
		</Stream>
	</Table>
	>>> col.count(9)
	3

	NOTE:  the .Name attribute returns the stripped "Name" attribute of
	the element, e.g. as would be obtained with StripColumnName(), but
	when assigning to the .Name attribute the value provided is stored
	without modification, i.e. there is no attempt to reattach the
	table's name to the string.  The calling code is responsible for
	doing the correct manipulations.  Therefore, the assignment
	operation below

	>>> col.Name, col.getAttribute("Name")
	(u'snr', u'test:snr')
	>>> col.Name = col.Name
	>>> col.Name, col.getAttribute("Name")
	(u'snr', u'snr')

	does not preserve the value of the "Name" attribute (though it does
	preserve the stripped form reported by the .Name property).  This
	asymmetry is necessary because the correct table name string to
	reattach to the attribute's value cannot always be known, e.g., if
	the Column object is not part of an XML tree and does not have a
	parent node.
	"""
	Name = ligolw.attributeproxy(u"Name", dec = StripColumnName)

	def __len__(self):
		"""
		The number of values in this column.
		"""
		return len(self.parentNode)

	def __getitem__(self, i):
		"""
		Retrieve the value in this column in row i.
		"""
		if isinstance(i, slice):
			return map(lambda r: getattr(r, self.Name), self.parentNode[i])
		else:
			return getattr(self.parentNode[i], self.Name)

	def __setitem__(self, i, value):
		"""
		Set the value in this column in row i.  i may be a slice.

		NOTE:  Unlike normal Python lists, the length of the Column
		cannot be changed as it is tied to the number of rows in
		the Table.  Therefore, if i is a slice, value should be an
		iterable with exactly the correct number of items.  No
		check is performed to ensure that this is true:  if value
		contains too many items the extras will be ignored, and if
		value contains too few items only as many rows will be
		updated as there are items.
		"""
		if isinstance(i, slice):
			for r, val in itertools.izip(self.parentNode[i], value):
				setattr(r, self.Name, val)
		else:
			setattr(self.parentNode[i], self.Name, value)

	def __delitem__(self, *args):
		raise NotImplementedError

	def __iter__(self):
		"""
		Return an iterator object for iterating over values in this
		column.
		"""
		for row in self.parentNode:
			yield getattr(row, self.Name)

	def count(self, value):
		"""
		Return the number of rows with this column equal to value.
		"""
		return sum(getattr(row, self.Name) == value for row in self.parentNode)

	def index(self, value):
		"""
		Return the smallest index of the row(s) with this column
		equal to value.
		"""
		for i in xrange(len(self.parentNode)):
			if getattr(self.parentNode[i], self.Name) == value:
				return i
		raise ValueError(value)

	def __contains__(self, value):
		"""
		Returns True or False if there is or is not, respectively,
		a row containing val in this column.
		"""
		for i in xrange(len(self.parentNode)):
			if getattr(self.parentNode[i], self.Name) == value:
				return True
		return False

	def asarray(self):
		"""
		Construct a numpy array from this column.  Note that this
		creates a copy of the data, so modifications made to the
		array will *not* be recorded in the original document.
		"""
		# most codes don't use this feature, this is the only place
		# numpy is used here, and importing numpy can be
		# time-consuming, so we derfer the import until needed.
		import numpy
		try:
			dtype = ligolwtypes.ToNumPyType[self.Type]
		except KeyError as e:
			raise TypeError("cannot determine numpy dtype for Column '%s': %s" % (self.getAttribute("Name"), e))
		return numpy.fromiter(self, dtype = dtype)


#
# =============================================================================
#
#                                Stream Element
#
# =============================================================================
#


#
# A subclass of tokenizer.RowBuilder that interns strings.
#


class InterningRowBuilder(tokenizer.RowBuilder):
	"""
	This subclass of the tokenizer.RowBuilder class respects the
	"interning" hints provided by table definitions, and attempts to
	replace the values of row attributes associated with interned
	columns with references to shared instances of those values.  This
	results in a reduction in memory use which is small for most
	documents, but can be subtantial when dealing with tables
	containing large volumes of repeated information.

	Example:

	>>> class Row(object):
	...	pass
	...
	>>> # 3rd arg is optional list of attributes to intern
	>>> rows = InterningRowBuilder(Row, ["name", "age"], ("name",))
	>>> l = list(rows.append(["Dick", 20., "Jane", 75., "Dick", 22.]))
	>>> l[0].name
	'Dick'
	>>> l[2].name
	'Dick'
	>>> l[2].name is l[0].name
	True

	Note that Python naturally interns short strings, so this example
	would return True regardless;  it is intended only to demonstrate
	the use of the class.

	The values are stored in a dictionary that is shared between all
	instances of this class, and which survives forever.  Nothing is
	ever naturally "uninterned", so the string dictionary grows without
	bound as more documents are processed.  This can be a problem in
	some use cases, and the work-around is to run

	>>> InterningRowBuilder.strings.clear()

	to reset the dictionary at appropriate points in the application.
	Typically this would be done immediately after each document is
	loaded.
	"""
	strings = {}
	def append(self, tokens):
		interns = self.interns
		setdefault = self.strings.setdefault
		for row in super(InterningRowBuilder, self).append(tokens):
			for col in interns:
				val = getattr(row, col)
				setattr(row, col, setdefault(val, val))
			yield row


#
# Select the RowBuilder class to use when parsing tables.
#


RowBuilder = tokenizer.RowBuilder


#
# Stream class
#


class TableStream(ligolw.Stream):
	"""
	High-level Stream element for use inside Tables.  This element
	knows how to parse the delimited character stream into row objects
	that it appends into the list-like parent element, and knows how to
	turn the parent's rows back into a character stream.
	"""
	def __init__(self, *args):
		super(TableStream, self).__init__(*args)
		self._tokenizer = tokenizer.Tokenizer(self.Delimiter)
		self._rowbuilder = None

	def config(self, parentNode):
		# some initialization that requires access to the
		# parentNode, and so cannot be done inside the __init__()
		# function.
		loadcolumns = set(parentNode.columnnames)
		if parentNode.loadcolumns is not None:
			# FIXME:  convert loadcolumns attributes to sets to
			# avoid the conversion.
			loadcolumns &= set(parentNode.loadcolumns)
		self._tokenizer.set_types([(pytype if colname in loadcolumns else None) for pytype, colname in zip(parentNode.columnpytypes, parentNode.columnnames)])
		columnnames = [name for name in parentNode.columnnames if name in loadcolumns]
		# FIXME:  convert interncolumns attributes to sets to
		# simplify computing the intersection
		interncolumns = [name for name in (parentNode.interncolumns or set()) if name in columnnames]
		self._rowbuilder = RowBuilder(parentNode.RowType, columnnames, interncolumns)
		return self

	def appendData(self, content):
		# tokenize buffer, pack into row objects, and append to
		# table
		appendfunc = self.parentNode.append
		for row in self._rowbuilder.append(self._tokenizer.append(content)):
			appendfunc(row)

	def unlink(self):
		"""
		Break internal references within the document tree rooted
		on this element to promote garbage collection.
		"""
		self._tokenizer = None
		self._rowbuilder = None
		super(TableStream, self).unlink()

	def endElement(self):
		# stream tokenizer uses delimiter to identify end of each
		# token, so add a final delimiter to induce the last token
		# to get parsed but only if there's something other than
		# whitespace left in the tokenizer's buffer.
		if not self._tokenizer.data.isspace():
			self.appendData(self.Delimiter)
		# call parent's _end_of_rows() hook.
		self.parentNode._end_of_rows()

	def write(self, fileobj = sys.stdout, indent = u""):
		# retrieve the .write() method of the file object to avoid
		# doing the attribute lookup in loops
		w = fileobj.write
		# loop over parent's rows.  This is complicated because we
		# need to not put a delimiter at the end of the last row
		# unless it ends with a null token
		w(self.start_tag(indent))
		rowdumper = tokenizer.RowDumper(self.parentNode.columnnames, [ligolwtypes.FormatFunc[coltype] for coltype in self.parentNode.columntypes], self.Delimiter)
		rowdumper.dump(self.parentNode)
		try:
			line = rowdumper.next()
		except StopIteration:
			# table is empty
			pass
		else:
			# write first row
			newline = u"\n" + indent + ligolw.Indent
			w(newline)
			# the xmlescape() call replaces things like "<"
			# with "&lt;" so that the string will not confuse
			# an XML parser when the file is read.  turning
			# "&lt;" back into "<" during file reading is
			# handled by the XML parser, so there is no code
			# in Glue related to that.
			w(xmlescape(line))
			# now add delimiter and write the remaining rows
			newline = rowdumper.delimiter + newline
			for line in rowdumper:
				w(newline)
				w(xmlescape(line))
			if rowdumper.tokens and rowdumper.tokens[-1] == u"":
				# the last token of the last row was null:
				# add a final delimiter to indicate that a
				# token is present
				w(rowdumper.delimiter)
		w(u"\n" + self.end_tag(indent) + u"\n")


#
# =============================================================================
#
#                                Table Element
#
# =============================================================================
#


class TableRow(object):
	"""
	Helpful parent class for row objects.  Also used as the default row
	class by Table instances.
	"""
	pass


class Table(ligolw.Table, list):
	"""
	High-level Table element that knows about its columns and rows.
	"""
	validcolumns = None
	loadcolumns = None
	interncolumns = None
	constraints = None
	how_to_index = None
	RowType = TableRow
	next_id = None

	def __init__(self, *args):
		"""
		Initialize
		"""
		super(Table, self).__init__(*args)
		self.columnnames = []
		self.columntypes = []
		self.columnpytypes = []

	Name = ligolw.attributeproxy(u"Name", enc = (lambda name: u"%s:table" % name), dec = StripTableName)


	#
	# Table retrieval
	#


	@classmethod
	def get_table(cls, xmldoc):
		"""
		Equivalent to the module-level function get_table(), but
		uses the .tableName attribute of this class to provide the
		name of the table to search for.  The Table parent class
		does not provide a .tableName attribute, but sub-classes,
		especially those in lsctables.py, do provide a value for
		that attribute, and in those cases this class method
		provides a cleaner way to retrieve them.

		Example:

		>>> import ligolw
		>>> import lsctables
		>>> xmldoc = ligolw.Document()
		>>> xmldoc.appendChild(ligolw.LIGO_LW()).appendChild(lsctables.New(lsctables.SnglInspiralTable))
		[]
		>>> sngl_inspiral_table = lsctables.SnglInspiralTable.get_table(xmldoc)
		"""
		return get_table(xmldoc, cls.tableName)

	def copy(self):
		"""
		Construct and return a new Table document subtree whose
		structure is the same as this table, that is it has the
		same columns etc..  The rows are not copied.  Note that a
		fair amount of metadata is shared between the original and
		new tables.  In particular, a copy of the Table object
		itself is created (but with no rows), and copies of the
		child nodes are created.  All other object references are
		shared between the two instances, such as the RowType
		attribute on the Table object.
		"""
		new = copy.copy(self)
		new.childNodes = map(copy.copy, self.childNodes)
		for child in new.childNodes:
			child.parentNode = new
		del new[:]
		new._end_of_columns()
		new._end_of_rows()
		return new


	@classmethod
	def CheckElement(cls, elem):
		"""
		Return True if element is a Table element whose Name
		attribute matches the .tableName attribute of this class
		according to CompareTableNames();  return False otherwise.
		See also .CheckProperties().
		"""
		return cls.CheckProperties(elem.tagName, elem.attributes)


	@classmethod
	def CheckProperties(cls, tagname, attrs):
		"""
		Return True if tagname and attrs are the XML tag name and
		element attributes, respectively, of a Table element whose
		Name attribute matches the .tableName attribute of this
		class according to CompareTableNames();  return False
		otherwise.  The Table parent class does not provide a
		.tableName attribute, but sub-classes, especially those in
		lsctables.py, do provide a value for that attribute.  See
		also .CheckElement()

		Example:

		>>> import lsctables
		>>> lsctables.ProcessTable.CheckProperties(u"Table", {u"Name": u"process:table"})
		True
		"""
		return tagname == cls.tagName and not CompareTableNames(attrs[u"Name"], cls.tableName)


	#
	# Column access
	#

	def getColumnByName(self, name):
		"""
		Retrieve and return the Column child element named name.
		The comparison is done using CompareColumnNames().  Raises
		KeyError if this table has no column by that name.

		Example:

		>>> import lsctables
		>>> tbl = lsctables.New(lsctables.SnglInspiralTable)
		>>> col = tbl.getColumnByName("mass1")
		"""
		try:
			col, = getColumnsByName(self, name)
		except ValueError:
			# did not find exactly 1 matching child
			raise KeyError(name)
		return col


	def appendColumn(self, name):
		"""
		Append a Column element named "name" to the table.  Returns
		the new child.  Raises ValueError if the table already has
		a column by that name, and KeyError if the validcolumns
		attribute of this table does not contain an entry for a
		column by that name.

		Note that the name string is assumed to be "pre-stripped",
		that is it is the significant portion of the elements Name
		attribute.  The Column element's Name attribute will be
		constructed by pre-pending the stripped Table element's
		name and a colon.

		Example:

		>>> import lsctables
		>>> process_table = lsctables.New(lsctables.ProcessTable, [])
		>>> col = process_table.appendColumn("program")
		>>> col.getAttribute("Name")
		'process:program'
		>>> col.Name
		'program'
		"""
		try:
			self.getColumnByName(name)
			# if we get here the table already has that column
			raise ValueError("duplicate Column '%s'" % name)
		except KeyError:
			pass
		column = Column(AttributesImpl({u"Name": "%s:%s" % (StripTableName(self.tableName), name), u"Type": self.validcolumns[name]}))
		streams = self.getElementsByTagName(ligolw.Stream.tagName)
		if streams:
			self.insertBefore(column, streams[0])
		else:
			self.appendChild(column)
		return column


	#
	# Row access
	#

	def appendRow(self, *args, **kwargs):
		"""
		Create and append a new row to this table, then return it

		All positional and keyword arguments are passed to the RowType
		constructor for this table.
		"""
		row = self.RowType(*args, **kwargs)
		self.append(row)
		return row


	#
	# Element methods
	#

	def _update_column_info(self):
		"""
		Used for validation during parsing, and additional
		book-keeping.  For internal use only.
		"""
		del self.columnnames[:]
		del self.columntypes[:]
		del self.columnpytypes[:]
		for child in self.getElementsByTagName(ligolw.Column.tagName):
			if self.validcolumns is not None:
				try:
					if self.validcolumns[child.Name] != child.Type:
						raise ligolw.ElementError("invalid type '%s' for Column '%s' in Table '%s', expected type '%s'" % (child.Type, child.getAttribute("Name"), self.getAttribute("Name"), self.validcolumns[child.Name]))
				except KeyError:
					raise ligolw.ElementError("invalid Column '%s' for Table '%s'" % (child.getAttribute("Name"), self.getAttribute("Name")))
			if child.Name in self.columnnames:
				raise ligolw.ElementError("duplicate Column '%s' in Table '%s'" % (child.getAttribute("Name"), self.getAttribute("Name")))
			self.columnnames.append(child.Name)
			self.columntypes.append(child.Type)
			try:
				self.columnpytypes.append(ligolwtypes.ToPyType[child.Type])
			except KeyError:
				raise ligolw.ElementError("unrecognized Type '%s' for Column '%s' in Table '%s'" % (child.Type, child.getAttribute("Name"), self.getAttribute("Name")))

	def _verifyChildren(self, i):
		"""
		Used for validation during parsing, and additional
		book-keeping.  For internal use only.
		"""
		super(Table, self)._verifyChildren(i)
		child = self.childNodes[i]
		if child.tagName == ligolw.Column.tagName:
			self._update_column_info()
		elif child.tagName == ligolw.Stream.tagName:
			# require agreement of non-stripped strings
			if child.getAttribute("Name") != self.getAttribute("Name"):
				raise ligolw.ElementError("Stream name '%s' does not match Table name '%s'" % (child.getAttribute("Name"), self.getAttribute("Name")))

	def _end_of_columns(self):
		"""
		Called during parsing to indicate that the last Column
		child element has been added.  Subclasses can override this
		to perform any special action that should occur following
		the addition of the last Column element.
		"""
		pass

	def _end_of_rows(self):
		"""
		Called during parsing to indicate that the last row has
		been added.  Subclasses can override this to perform any
		special action that should occur following the addition of
		the last row.
		"""
		pass

	def removeChild(self, child):
		"""
		Remove a child from this element.  The child element is
		returned, and it's parentNode element is reset.
		"""
		super(Table, self).removeChild(child)
		if child.tagName == ligolw.Column.tagName:
			self._update_column_info()
		return child

	def unlink(self):
		"""
		Break internal references within the document tree rooted
		on this element to promote garbage collection.
		"""
		super(Table, self).unlink()
		del self[:]

	def endElement(self):
		# Table elements are allowed to contain 0 Stream children,
		# but _end_of_columns() and _end_of_rows() hooks must be
		# called regardless, so we do that here if needed.
		if self.childNodes[-1].tagName != ligolw.Stream.tagName:
			self._end_of_columns()
			self._end_of_rows()

	#
	# Row ID manipulation
	#

	@classmethod
	def get_next_id(cls):
		"""
		Returns the current value of the next_id class attribute,
		and increments the next_id class attribute by 1.  Raises
		ValueError if the table does not have an ID generator
		associated with it.
		"""
		# = None if no ID generator
		id = cls.next_id
		cls.next_id += 1
		return id

	@classmethod
	def set_next_id(cls, id):
		"""
		Sets the value of the next_id class attribute.  This is a
		convenience function to help prevent accidentally assigning
		a value to an instance attribute instead of the class
		attribute.
		"""
		cls.next_id = id

	def sync_next_id(self):
		"""
		Determines the highest-numbered ID in this table, and sets
		the table's .next_id attribute to the next highest ID in
		sequence.  If the .next_id attribute is already set to a
		value greater than the highest value found, then it is left
		unmodified.  The return value is the ID identified by this
		method.  If the table's .next_id attribute is None, then
		this function is a no-op.

		Note that tables of the same name typically share a common
		.next_id attribute (it is a class attribute, not an
		attribute of each instance) so that IDs can be generated
		that are unique across all tables in the document.  Running
		sync_next_id() on all the tables in a document that are of
		the same type will have the effect of setting the ID to the
		next ID higher than any ID in any of those tables.

		Example:

		>>> import lsctables
		>>> tbl = lsctables.New(lsctables.ProcessTable)
		>>> print tbl.sync_next_id()
		process:process_id:0
		"""
		if self.next_id is not None:
			if len(self):
				n = max(self.getColumnByName(self.next_id.column_name)) + 1
			else:
				n = type(self.next_id)(0)
			if n > self.next_id:
				self.set_next_id(n)
		return self.next_id

	def updateKeyMapping(self, mapping):
		"""
		Used as the first half of the row key reassignment
		algorithm.  Accepts a dictionary mapping old key --> new
		key.  Iterates over the rows in this table, using the
		table's next_id attribute to assign a new ID to each row,
		recording the changes in the mapping.  Returns the mapping.
		Raises ValueError if the table's next_id attribute is None.
		"""
		if self.next_id is None:
			raise ValueError(self)
		try:
			column = self.getColumnByName(self.next_id.column_name)
		except KeyError:
			# table is missing its ID column, this is a no-op
			return mapping
		for i, old in enumerate(column):
			if old is None:
				raise ValueError("null row ID encountered in Table '%s', row %d" % (self.getAttribute("Name"), i))
			if old in mapping:
				column[i] = mapping[old]
			else:
				column[i] = mapping[old] = self.get_next_id()
		return mapping

	def applyKeyMapping(self, mapping):
		"""
		Used as the second half of the key reassignment algorithm.
		Loops over each row in the table, replacing references to
		old row keys with the new values from the mapping.
		"""
		for coltype, colname in zip(self.columntypes, self.columnnames):
			if coltype in ligolwtypes.IDTypes and (self.next_id is None or colname != self.next_id.column_name):
				column = self.getColumnByName(colname)
				for i, old in enumerate(column):
					try:
						column[i] = mapping[old]
					except KeyError:
						pass


#
# =============================================================================
#
#                               Content Handler
#
# =============================================================================
#


#
# Override portions of a ligolw.LIGOLWContentHandler class
#


def use_in(ContentHandler):
	"""
	Modify ContentHandler, a sub-class of
	pycbc.ligolw.LIGOLWContentHandler, to cause it to use the Table,
	Column, and Stream classes defined in this module when parsing XML
	documents.

	Example:

	>>> from pycbc.ligolw import ligolw
	>>> class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
	...	pass
	...
	>>> use_in(LIGOLWContentHandler)
	<class 'pycbc.ligolw.table.LIGOLWContentHandler'>
	"""
	def startColumn(self, parent, attrs):
		return Column(attrs)

	def startStream(self, parent, attrs, __orig_startStream = ContentHandler.startStream):
		if parent.tagName == ligolw.Table.tagName:
			parent._end_of_columns()
			return TableStream(attrs).config(parent)
		return __orig_startStream(self, parent, attrs)

	def startTable(self, parent, attrs):
		return Table(attrs)

	ContentHandler.startColumn = startColumn
	ContentHandler.startStream = startStream
	ContentHandler.startTable = startTable

	return ContentHandler

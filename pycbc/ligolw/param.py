# Copyright (C) 2006--2009,2012--2015  Kipp Cannon
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
High-level support for Param elements.
"""


import pickle
import re
import sys
from xml.sax.saxutils import escape as xmlescape
from xml.sax.xmlreader import AttributesImpl as Attributes
try:
	import yaml
except ImportError:
	# yaml serialization is optional
	pass


from pycbc import version
from . import ligolw
from . import types as ligolwtypes


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % version.version
__date__ = version.date


#
# =============================================================================
#
#                           Param Name Manipulation
#
# =============================================================================
#


#
# Regular expression used to extract the signifcant portion of a param
# name, according to LIGO LW naming conventions.
#


ParamPattern = re.compile(r"(?P<Name>[a-z0-9_:]+):param\Z")


def StripParamName(name):
	"""
	Return the significant portion of a param name according to LIGO LW
	naming conventions.
	"""
	try:
		return ParamPattern.search(name).group("Name")
	except AttributeError:
		return name


def CompareParamNames(name1, name2):
	"""
	Convenience function to compare two param names according to LIGO
	LW naming conventions.
	"""
	return cmp(StripParamName(name1), StripParamName(name2))


def getParamsByName(elem, name):
	"""
	Return a list of params with name name under elem.
	"""
	name = StripParamName(name)
	return elem.getElements(lambda e: (e.tagName == ligolw.Param.tagName) and (e.Name == name))


#
# =============================================================================
#
#                                  Utilities
#
# =============================================================================
#


def new_param(name, type, value, start = None, scale = None, unit = None, dataunit = None, comment = None):
	"""
	Construct a LIGO Light Weight XML Param document subtree.  FIXME:
	document keyword arguments.
	"""
	elem = Param()
	elem.Name = name
	elem.Type = type
	elem.pcdata = value
	# FIXME:  I have no idea how most of the attributes should be
	# encoded, I don't even know what they're supposed to be.
	if dataunit is not None:
		elem.DataUnit = dataunit
	if scale is not None:
		elem.Scale = scale
	if start is not None:
		elem.Start = start
	if unit is not None:
		elem.Unit = unit
	if comment is not None:
		elem.appendChild(ligolw.Comment())
		elem.childNodes[-1].pcdata = comment
	return elem


def get_param(xmldoc, name):
	"""
	Scan xmldoc for a param named name.  Raises ValueError if not
	exactly 1 such param is found.
	"""
	params = getParamsByName(xmldoc, name)
	if len(params) != 1:
		raise ValueError("document must contain exactly one %s param" % StripParamName(name))
	return params[0]


def from_pyvalue(name, value, **kwargs):
	"""
	Convenience wrapper for new_param() that constructs a Param element
	from an instance of a Python builtin type.  See new_param() for a
	description of the valid keyword arguments.
	"""
	return new_param(name, ligolwtypes.FromPyType[type(value)], value, **kwargs)


def get_pyvalue(xml, name):
	"""
	Convenience wrapper for get_param() that recovers an instance of a
	Python builtin type from a Param element.
	"""
	# Note:  the Param is automatically parsed into the correct Python
	# type, so this function is mostly a no-op.
	return get_param(xml, name).pcdata


#
# =============================================================================
#
#                        (De-)Serialization via Params
#
# =============================================================================
#


def pickle_to_param(obj, name):
	"""
	Return the top-level element of a document sub-tree containing the
	pickled serialization of a Python object.
	"""
	return from_pyvalue(u"pickle:%s" % name, unicode(pickle.dumps(obj)))


def pickle_from_param(elem, name):
	"""
	Retrieve a pickled Python object from the document tree rooted at
	elem.
	"""
	return pickle.loads(str(get_pyvalue(elem, u"pickle:%s" % name)))


def yaml_to_param(obj, name):
	"""
	Return the top-level element of a document sub-tree containing the
	YAML serialization of a Python object.
	"""
	return from_pyvalue(u"yaml:%s" % name, unicode(yaml.dump(obj)))


def yaml_from_param(elem, name):
	"""
	Retrieve a YAMLed Python object from the document tree rooted at
	elem.
	"""
	return yaml.load(get_pyvalue(elem, u"yaml:%s" % name))


try:
	yaml
except NameError:
	# yaml module not loaded, disable (de-)serializers
	def yaml_to_param(obj, name):
		"""
		Not available.  Install yaml to use.
		"""
		raise NameError("yaml not installed")
	def yaml_from_param(elem, name):
		"""
		Not available.  Install yaml to use.
		"""
		raise NameError("yaml not installed")


#
# =============================================================================
#
#                               Element Classes
#
# =============================================================================
#


#
# FIXME: params of type string should be quoted in order to correctly
# delimit their extent.  If that were done, then the pcdata in a Param
# element could be parsed using the Stream tokenizer (i.e., as though it
# were a single-token stream), which would guarantee that Stream data and
# Param data is parsed using the exact same rules.  Unfortunately, common
# practice is to not quote Param string values, so we parse things
# differently here.  In particular, we strip whitespace from the start and
# stop of all Param pcdata.  If this causes your string Param values to be
# corrupted (because you need leading and trailing white space preserved),
# then you need to make everyone switch to quoting their string Param
# values, and once that is done then this code will be changed.  Perhaps a
# warning should be emitted for non-quoted strings to encourage a
# transition?
#


class Param(ligolw.Param):
	"""
	High-level Param element.  The value is stored in the pcdata
	attribute as the native Python type rather than as a string.
	"""
	def __init__(self, *args):
		"""
		Initialize a new Param element.
		"""
		super(Param, self).__init__(*args)
		self.pytype = ligolwtypes.ToPyType[self.Type]

	def endElement(self):
		if self.pcdata is not None:
			# convert pcdata from string to native Python type
			self.pcdata = self.pytype(self.pcdata.strip())

	def write(self, fileobj = sys.stdout, indent = u""):
		fileobj.write(self.start_tag(indent))
		for c in self.childNodes:
			if c.tagName not in self.validchildren:
				raise ligolw.ElementError("invalid child %s for %s" % (c.tagName, self.tagName))
			c.write(fileobj, indent + ligolw.Indent)
		if self.pcdata is not None:
			# we have to strip quote characters from string
			# formats (see comment above)
			fileobj.write(xmlescape(ligolwtypes.FormatFunc[self.Type](self.pcdata).strip(u"\"")))
		fileobj.write(self.end_tag(u"") + u"\n")

	Name = ligolw.attributeproxy(u"Name", enc = (lambda name: u"%s:param" % name), dec = StripParamName)
	Scale = ligolw.attributeproxy(u"Scale", enc = ligolwtypes.FormatFunc[u"real_8"], dec = ligolwtypes.ToPyType[u"real_8"])
	Type = ligolw.attributeproxy(u"Type", default = u"lstring")


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
	pycbc.ligolw.LIGOLWContentHandler, to cause it to use the Param
	class defined in this module when parsing XML documents.

	Example:

	>>> from pycbc.ligolw import ligolw
	>>> def MyContentHandler(ligolw.LIGOLWContentHandler):
	...	pass
	...
	>>> use_in(MyContentHandler)
	"""
	def startParam(self, parent, attrs):
		return Param(attrs)

	ContentHandler.startParam = startParam

	return ContentHandler

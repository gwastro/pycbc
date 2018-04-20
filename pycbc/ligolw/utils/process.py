# Copyright (C) 2006--2013,2015  Kipp Cannon
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


"""
A collection of utilities to assist applications in manipulating the
process and process_params tables in LIGO Light-Weight XML documents.
"""


import os
import socket
import StringIO
import time


from pycbc import version
from .. import ligolw
from .. import lsctables
from .. import types as ligolwtypes


try:
	from lal import UTCToGPS as _UTCToGPS
except ImportError:
	# lal is optional
	# FIXME:  make it not optional
	from glue import gpstime
	_UTCToGPS = lambda utc: int(gpstime.GpsSecondsFromPyUTC(time.mktime(utc)))


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>, Larne Pekowsky <lppekows@physics.syr.edu>"
__version__ = "git id %s" % version.version
__date__ = version.date


#
# =============================================================================
#
#                               Process Metadata
#
# =============================================================================
#


def get_username():
	"""
	Try to retrieve the username from a variety of sources.  First the
	environment variable LOGNAME is tried, if that is not set the
	environment variable USERNAME is tried, if that is not set the
	password database is consulted (only on Unix systems, if the import
	of the pwd module succeeds), finally if that fails KeyError is
	raised.
	"""
	try:
		return os.environ["LOGNAME"]
	except KeyError:
		pass
	try:
		return os.environ["USERNAME"]
	except KeyError:
		pass
	try:
		import pwd
		return pwd.getpwuid(os.getuid())[0]
	except (ImportError, KeyError):
		raise KeyError


def append_process(xmldoc, program = None, version = None, cvs_repository = None, cvs_entry_time = None, comment = None, is_online = False, jobid = 0, domain = None, ifos = None):
	"""
	Add an entry to the process table in xmldoc.  program, version,
	cvs_repository, comment, and domain should all be strings or
	unicodes.  cvs_entry_time should be a string or unicode in the
	format "YYYY/MM/DD HH:MM:SS".  is_online should be a boolean, jobid
	an integer.  ifos should be an iterable (set, tuple, etc.) of
	instrument names.

	See also register_to_xmldoc().
	"""
	try:
		proctable = lsctables.ProcessTable.get_table(xmldoc)
	except ValueError:
		proctable = lsctables.New(lsctables.ProcessTable)
		xmldoc.childNodes[0].appendChild(proctable)

	proctable.sync_next_id()

	process = proctable.RowType()
	process.program = program
	process.version = version
	process.cvs_repository = cvs_repository
	# FIXME:  remove the "" case when the git versioning business is
	# sorted out
	if cvs_entry_time is not None and cvs_entry_time != "":
		try:
			# try the git_version format first
			process.cvs_entry_time = _UTCToGPS(time.strptime(cvs_entry_time, "%Y-%m-%d %H:%M:%S +0000"))
		except ValueError:
			# fall back to the old cvs format
			process.cvs_entry_time = _UTCToGPS(time.strptime(cvs_entry_time, "%Y/%m/%d %H:%M:%S"))
	else:
		process.cvs_entry_time = None
	process.comment = comment
	process.is_online = int(is_online)
	process.node = socket.gethostname()
	try:
		process.username = get_username()
	except KeyError:
		process.username = None
	process.unix_procid = os.getpid()
	process.start_time = _UTCToGPS(time.gmtime())
	process.end_time = None
	process.jobid = jobid
	process.domain = domain
	process.instruments = ifos
	process.process_id = proctable.get_next_id()
	proctable.append(process)
	return process


def set_process_end_time(process):
	"""
	Set the end time in a row in a process table to the current time.
	"""
	process.end_time = _UTCToGPS(time.gmtime())
	return process


def append_process_params(xmldoc, process, params):
	"""
	xmldoc is an XML document tree, process is the row in the process
	table for which these are the parameters, and params is a list of
	(name, type, value) tuples one for each parameter.

	See also process_params_from_dict(), register_to_xmldoc().
	"""
	try:
		paramtable = lsctables.ProcessParamsTable.get_table(xmldoc)
	except ValueError:
		paramtable = lsctables.New(lsctables.ProcessParamsTable)
		xmldoc.childNodes[0].appendChild(paramtable)

	for name, typ, value in params:
		row = paramtable.RowType()
		row.program = process.program
		row.process_id = process.process_id
		row.param = unicode(name)
		if typ is not None:
			row.type = unicode(typ)
			if row.type not in ligolwtypes.Types:
				raise ValueError("invalid type '%s' for parameter '%s'" % (row.type, row.param))
		else:
			row.type = None
		if value is not None:
			row.value = unicode(value)
		else:
			row.value = None
		paramtable.append(row)
	return process


def get_process_params(xmldoc, program, param, require_unique_program = True):
	"""
	Return a list of the values stored in the process_params table for
	params named param for the program(s) named program.  The values
	are returned as Python native types, not as the strings appearing
	in the XML document.  If require_unique_program is True (default),
	then the document must contain exactly one program with the
	requested name, otherwise ValueError is raised.  If
	require_unique_program is not True, then there must be at least one
	program with the requested name otherwise ValueError is raised.
	"""
	process_ids = lsctables.ProcessTable.get_table(xmldoc).get_ids_by_program(program)
	if len(process_ids) < 1:
		raise ValueError("process table must contain at least one program named '%s'" % program)
	elif require_unique_program and len(process_ids) != 1:
		raise ValueError("process table must contain exactly one program named '%s'" % program)
	return [row.pyvalue for row in lsctables.ProcessParamsTable.get_table(xmldoc) if (row.process_id in process_ids) and (row.param == param)]


def doc_includes_process(xmldoc, program):
	"""
	Return True if the process table in xmldoc includes entries for a
	program named program.
	"""
	return program in lsctables.ProcessTable.get_table(xmldoc).getColumnByName(u"program")


def process_params_from_dict(paramdict):
	"""
	Generator function yields (name, type, value) tuples constructed
	from a dictionary of name/value pairs.  The tuples are suitable for
	input to append_process_params().  This is intended as a
	convenience for converting command-line options into process_params
	rows.  The name values in the output have "--" prepended to them
	and all "_" characters replaced with "-".  The type strings are
	guessed from the Python types of the values.  If a value is a
	Python list (or instance of a subclass thereof), then one tuple is
	produced for each of the items in the list.

	Example:

	>>> list(process_params_from_dict({"verbose": True, "window": 4.0, "include": ["/tmp", "/var/tmp"]}))
	[(u'--window', u'real_8', 4.0), (u'--verbose', None, None), (u'--include', u'lstring', '/tmp'), (u'--include', u'lstring', '/var/tmp')]
	"""
	for name, values in paramdict.items():
		# change the name back to the form it had on the command line
		name = u"--%s" % name.replace("_", "-")

		if values is True or values is False:
			yield (name, None, None)
		elif values is not None:
			if not isinstance(values, list):
				values = [values]
			for value in values:
				yield (name, ligolwtypes.FromPyType[type(value)], value)


def register_to_xmldoc(xmldoc, program, paramdict, **kwargs):
	"""
	Register the current process and params to an XML document.
	program is the name of the program.  paramdict is a dictionary of
	name/value pairs that will be used to populate the process_params
	table;  see process_params_from_dict() for information on how these
	name/value pairs are interpreted.  Any additional keyword arguments
	are passed to append_process().  Returns the new row from the
	process table.
	"""
	process = append_process(xmldoc, program = program, **kwargs)
	append_process_params(xmldoc, process, process_params_from_dict(paramdict))
	return process


# The tables in the segment database declare most fields "NOT NULL", so provide stub values
def register_to_ldbd(client, program, paramdict, version = u'0', cvs_repository = u'-', cvs_entry_time = 0, comment = u'-', is_online = False, jobid = 0, domain = None, ifos = u'-'):
	"""
	Register the current process and params to a database via a
	LDBDClient.  The program and paramdict arguments and any additional
	keyword arguments are the same as those for register_to_xmldoc().
	Returns the new row from the process table.
	"""
	xmldoc = ligolw.Document()
	xmldoc.appendChild(ligolw.LIGO_LW())
	process = register_to_xmldoc(xmldoc, program, paramdict, version = version, cvs_repository = cvs_repository, cvs_entry_time = cvs_entry_time, comment = comment, is_online = is_online, jobid = jobid, domain = domain, ifos = ifos)

	fake_file = StringIO.StringIO()
	xmldoc.write(fake_file)
	client.insert(fake_file.getvalue())

	return process

# Copyright (C) 2007-2015  Kipp Cannon
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
This module provides an implementation of the Table element that uses a
database engine for storage.  On top of that it then re-implements a number
of the tables from the lsctables module to provide versions of their
methods that work against the SQL database.
"""


import itertools
import operator
import os
import re
import shutil
import signal
import sys
import tempfile
import threading
from xml.sax.xmlreader import AttributesImpl
import warnings


from pycbc import version
from glue import offsetvector
from glue import segments
from . import ilwd
from . import ligolw
from . import table
from . import lsctables
from . import types as ligolwtypes


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % version.version
__date__ = version.date


#
# =============================================================================
#
#                                  Connection
#
# =============================================================================
#


def connection_db_type(connection):
	"""
	A totally broken attempt to determine what type of database a
	connection object is attached to.  Don't use this.

	The input is a DB API 2.0 compliant connection object, the return
	value is one of the strings "sqlite3" or "mysql".  Raises TypeError
	when the database type cannot be determined.
	"""
	if "sqlite" in repr(connection):
		return "sqlite"
	if "mysql" in repr(connection):
		return "mysql"
	raise TypeError(connection)


#
# Module-level variable used to hold references to
# tempfile.NamedTemporaryFiles objects to prevent them from being deleted
# while in use.  NOT MEANT FOR USE BY CODE OUTSIDE OF THIS MODULE!
#


temporary_files = {}
temporary_files_lock = threading.Lock()


#
# Module-level variable to hold the signal handlers that have been
# overridden as part of the clean-up-scratch-files-on-signal feature.  NOT
# MEANT FOR USE BY CODE OUTSIDE OF THIS MODULE!
#


origactions = {}


def install_signal_trap(signums = (signal.SIGTERM, signal.SIGTSTP), retval = 1):
	"""
	Installs a signal handler to erase temporary scratch files when a
	signal is received.  This can be used to help ensure scratch files
	are erased when jobs are evicted by Condor.  signums is a squence
	of the signals to trap, the default value is a list of the signals
	used by Condor to kill and/or evict jobs.

	The logic is as follows.  If the current signal handler is
	signal.SIG_IGN, i.e. the signal is being ignored, then the signal
	handler is not modified since the reception of that signal would
	not normally cause a scratch file to be leaked.  Otherwise a signal
	handler is installed that erases the scratch files.  If the
	original signal handler was a Python callable, then after the
	scratch files are erased the original signal handler will be
	invoked.  If program control returns from that handler, i.e.  that
	handler does not cause the interpreter to exit, then sys.exit() is
	invoked and retval is returned to the shell as the exit code.

	Note:  by invoking sys.exit(), the signal handler causes the Python
	interpreter to do a normal shutdown.  That means it invokes
	atexit() handlers, and does other garbage collection tasks that it
	normally would not do when killed by a signal.

	Note:  this function will not replace a signal handler more than
	once, that is if it has already been used to set a handler
	on a signal then it will be a no-op when called again for that
	signal until uninstall_signal_trap() is used to remove the handler
	from that signal.

	Note:  this function is called by get_connection_filename()
	whenever it creates a scratch file.
	"""
	# NOTE:  this must be called with the temporary_files_lock held.
	# ignore signums we've already replaced
	signums = set(signums) - set(origactions)

	def temporary_file_cleanup_on_signal(signum, frame):
		with temporary_files_lock:
			temporary_files.clear()
		if callable(origactions[signum]):
			# original action is callable, chain to it
			return origactions[signum](signum, frame)
		# original action was not callable or the callable
		# returned.  invoke sys.exit() with retval as exit code
		sys.exit(retval)

	for signum in signums:
		origactions[signum] = signal.getsignal(signum)
		if origactions[signum] != signal.SIG_IGN:
			# signal is not being ignored, so install our
			# handler
			signal.signal(signum, temporary_file_cleanup_on_signal)


def uninstall_signal_trap(signums = None):
	"""
	Undo the effects of install_signal_trap().  Restores the original
	signal handlers.  If signums is a sequence of signal numbers only
	the signal handlers for those signals will be restored (KeyError
	will be raised if one of them is not one that install_signal_trap()
	installed a handler for, in which case some undefined number of
	handlers will have been restored).  If signums is None (the
	default) then all signals that have been modified by previous calls
	to install_signal_trap() are restored.

	Note:  this function is called by put_connection_filename() and
	discard_connection_filename() whenever they remove a scratch file
	and there are then no more scrach files in use.
	"""
	# NOTE:  this must be called with the temporary_files_lock held.
	if signums is None:
		signums = origactions.keys()
	for signum in signums:
		signal.signal(signum, origactions.pop(signum))


#
# Functions to work with database files in scratch space
#


def get_connection_filename(filename, tmp_path = None, replace_file = False, verbose = False):
	"""
	Utility code for moving database files to a (presumably local)
	working location for improved performance and reduced fileserver
	load.
	"""
	def mktmp(path, suffix = ".sqlite", verbose = False):
		with temporary_files_lock:
			# make sure the clean-up signal traps are installed
			install_signal_trap()
			# create the remporary file and replace it's
			# unlink() function
			temporary_file = tempfile.NamedTemporaryFile(suffix = suffix, dir = path if path != "_CONDOR_SCRATCH_DIR" else os.getenv("_CONDOR_SCRATCH_DIR"))
			def new_unlink(self, orig_unlink = temporary_file.unlink):
				# also remove a -journal partner, ignore all errors
				try:
					orig_unlink("%s-journal" % self)
				except:
					pass
				orig_unlink(self)
			temporary_file.unlink = new_unlink
			filename = temporary_file.name
			# hang onto reference to prevent its removal
			temporary_files[filename] = temporary_file
		if verbose:
			print >>sys.stderr, "using '%s' as workspace" % filename
		# mkstemp() ignores umask, creates all files accessible
		# only by owner;  we should respect umask.  note that
		# os.umask() sets it, too, so we have to set it back after
		# we know what it is
		umsk = os.umask(0777)
		os.umask(umsk)
		os.chmod(filename, 0666 & ~umsk)
		return filename

	def truncate(filename, verbose = False):
		if verbose:
			print >>sys.stderr, "'%s' exists, truncating ..." % filename,
		try:
			fd = os.open(filename, os.O_WRONLY | os.O_TRUNC)
		except Exception as e:
			if verbose:
				print >>sys.stderr, "cannot truncate '%s': %s" % (filename, str(e))
			return
		os.close(fd)
		if verbose:
			print >>sys.stderr, "done."

	def cpy(srcname, dstname, verbose = False):
		if verbose:
			print >>sys.stderr, "copying '%s' to '%s' ..." % (srcname, dstname),
		shutil.copy2(srcname, dstname)
		if verbose:
			print >>sys.stderr, "done."
		try:
			# try to preserve permission bits.  according to
			# the documentation, copy() and copy2() are
			# supposed preserve them but don't.  maybe they
			# don't preserve them if the destination file
			# already exists?
			shutil.copystat(srcname, dstname)
		except Exception as e:
			if verbose:
				print >>sys.stderr, "warning: ignoring failure to copy permission bits from '%s' to '%s': %s" % (filename, target, str(e))

	database_exists = os.access(filename, os.F_OK)

	if tmp_path is not None:
		# for suffix, can't use splitext() because it only keeps
		# the last bit, e.g. won't give ".xml.gz" but just ".gz"
		target = mktmp(tmp_path, suffix = ".".join(os.path.split(filename)[-1].split(".")[1:]), verbose = verbose)
		if database_exists:
			if replace_file:
				# truncate database so that if this job
				# fails the user won't think the database
				# file is valid
				truncate(filename, verbose = verbose)
			else:
				# need to copy existing database to work
				# space for modifications
				i = 1
				while True:
					try:
						cpy(filename, target, verbose = verbose)
					except IOError as e:
						import errno
						import time
						if e.errno not in (errno.EPERM, errno.ENOSPC):
							# anything other
							# than out-of-space
							# is a real error
							raise
						if i < 5:
							if verbose:
								print >>sys.stderr, "warning: attempt %d: %s, sleeping and trying again ..." % (i, errno.errorcode[e.errno])
							time.sleep(10)
							i += 1
							continue
						if verbose:
							print >>sys.stderr, "warning: attempt %d: %s: working with original file '%s'" % (i, errno.errorcode[e.errno], filename)
						with temporary_files_lock:
							del temporary_files[target]
						target = filename
					break
	else:
		with temporary_files_lock:
			if filename in temporary_files:
				raise ValueError("file '%s' appears to be in use already as a temporary database file and is to be deleted" % filename)
		target = filename
		if database_exists and replace_file:
			truncate(target, verbose = verbose)

	del mktmp
	del truncate
	del cpy

	return target


def set_temp_store_directory(connection, temp_store_directory, verbose = False):
	"""
	Sets the temp_store_directory parameter in sqlite.
	"""
	if temp_store_directory == "_CONDOR_SCRATCH_DIR":
		temp_store_directory = os.getenv("_CONDOR_SCRATCH_DIR")
	if verbose:
		print >>sys.stderr, "setting the temp_store_directory to %s ..." % temp_store_directory,
	cursor = connection.cursor()
	cursor.execute("PRAGMA temp_store_directory = '%s'" % temp_store_directory)
	cursor.close()
	if verbose:
		print >>sys.stderr, "done"


def put_connection_filename(filename, working_filename, verbose = False):
	"""
	This function reverses the effect of a previous call to
	get_connection_filename(), restoring the working copy to its
	original location if the two are different.  This function should
	always be called after calling get_connection_filename() when the
	file is no longer in use.

	During the move operation, this function traps the signals used by
	Condor to evict jobs.  This reduces the risk of corrupting a
	document by the job terminating part-way through the restoration of
	the file to its original location.  When the move operation is
	concluded, the original signal handlers are restored and if any
	signals were trapped they are resent to the current process in
	order.  Typically this will result in the signal handlers installed
	by the install_signal_trap() function being invoked, meaning any
	other scratch files that might be in use get deleted and the
	current process is terminated.
	"""
	if working_filename != filename:
		# initialize SIGTERM and SIGTSTP trap
		deferred_signals = []
		def newsigterm(signum, frame):
			deferred_signals.append(signum)
		oldhandlers = {}
		for sig in (signal.SIGTERM, signal.SIGTSTP):
			oldhandlers[sig] = signal.getsignal(sig)
			signal.signal(sig, newsigterm)

		# replace document
		if verbose:
			print >>sys.stderr, "moving '%s' to '%s' ..." % (working_filename, filename),
		shutil.move(working_filename, filename)
		if verbose:
			print >>sys.stderr, "done."

		# remove reference to tempfile.TemporaryFile object.
		# because we've just deleted the file above, this would
		# produce an annoying but harmless message about an ignored
		# OSError, so we create a dummy file for the TemporaryFile
		# to delete.  ignore any errors that occur when trying to
		# make the dummy file.  FIXME: this is stupid, find a
		# better way to shut TemporaryFile up
		try:
			open(working_filename, "w").close()
		except:
			pass
		with temporary_files_lock:
			del temporary_files[working_filename]

		# restore original handlers, and send ourselves any trapped signals
		# in order
		for sig, oldhandler in oldhandlers.iteritems():
			signal.signal(sig, oldhandler)
		while deferred_signals:
			os.kill(os.getpid(), deferred_signals.pop(0))

		# if there are no more temporary files in place, remove the
		# temporary-file signal traps
		with temporary_files_lock:
			if not temporary_files:
				uninstall_signal_trap()


def discard_connection_filename(filename, working_filename, verbose = False):
	"""
	Like put_connection_filename(), but the working copy is simply
	deleted instead of being copied back to its original location.
	This is a useful performance boost if it is known that no
	modifications were made to the file, for example if queries were
	performed but no updates.

	Note that the file is not deleted if the working copy and original
	file are the same, so it is always safe to call this function after
	a call to get_connection_filename() even if a separate working copy
	is not created.
	"""
	if working_filename == filename:
		return
	with temporary_files_lock:
		if verbose:
			print >>sys.stderr, "removing '%s' ..." % working_filename,
		# remove reference to tempfile.TemporaryFile object
		del temporary_files[working_filename]
		if verbose:
			print >>sys.stderr, "done."
		# if there are no more temporary files in place, remove the
		# temporary-file signal traps
		if not temporary_files:
			uninstall_signal_trap()


#
# =============================================================================
#
#                                  ID Mapping
#
# =============================================================================
#


def idmap_create(connection):
	"""
	Create the _idmap_ table.  This table has columns "old" and "new"
	containing text strings mapping old IDs to new IDs.  The old column
	is a primary key (is indexed and must contain unique entries).  The
	table is created as a temporary table, so it will be automatically
	dropped when the database connection is closed.

	This function is for internal use, it forms part of the code used
	to re-map row IDs when merging multiple documents.
	"""
	connection.cursor().execute("CREATE TEMPORARY TABLE _idmap_ (old TEXT PRIMARY KEY NOT NULL, new TEXT NOT NULL)")


def idmap_reset(connection):
	"""
	Erase the contents of the _idmap_ table, but leave the table in
	place.

	This function is for internal use, it forms part of the code used
	to re-map row IDs when merging multiple documents.
	"""
	connection.cursor().execute("DELETE FROM _idmap_")


def idmap_sync(connection):
	"""
	Iterate over the tables in the database, ensure that there exists a
	custom DBTable class for each, and synchronize that table's ID
	generator to the ID values in the database.
	"""
	xmldoc = get_xml(connection)
	for tbl in xmldoc.getElementsByTagName(DBTable.tagName):
		tbl.sync_next_id()
	xmldoc.unlink()


def idmap_get_new(connection, old, tbl):
	"""
	From the old ID string, obtain a replacement ID string by either
	grabbing it from the _idmap_ table if one has already been assigned
	to the old ID, or by using the current value of the Table
	instance's next_id class attribute.  In the latter case, the new ID
	is recorded in the _idmap_ table, and the class attribute
	incremented by 1.

	This function is for internal use, it forms part of the code used
	to re-map row IDs when merging multiple documents.
	"""
	cursor = connection.cursor()
	cursor.execute("SELECT new FROM _idmap_ WHERE old == ?", (old,))
	new = cursor.fetchone()
	if new is not None:
		# a new ID has already been created for this old ID
		return ilwd.ilwdchar(new[0])
	# this ID was not found in _idmap_ table, assign a new ID and
	# record it
	new = tbl.get_next_id()
	cursor.execute("INSERT INTO _idmap_ VALUES (?, ?)", (old, new))
	return new


def idmap_get_max_id(connection, id_class):
	"""
	Given an ilwd:char ID class, return the highest ID from the table
	for whose IDs that is the class.

	Example:

	>>> event_id = ilwd.ilwdchar("sngl_burst:event_id:0")
	>>> print event_id
	sngl_inspiral:event_id:0
	>>> max_id = get_max_id(connection, type(event_id))
	>>> print max_id
	sngl_inspiral:event_id:1054
	"""
	cursor = connection.cursor()
	cursor.execute("SELECT MAX(CAST(SUBSTR(%s, %d, 10) AS INTEGER)) FROM %s" % (id_class.column_name, id_class.index_offset + 1, id_class.table_name))
	maxid = cursor.fetchone()[0]
	cursor.close()
	if maxid is None:
		return None
	return id_class(maxid)


#
# =============================================================================
#
#                             Database Information
#
# =============================================================================
#


#
# SQL parsing
#


_sql_create_table_pattern = re.compile(r"CREATE\s+TABLE\s+(?P<name>\w+)\s*\((?P<coldefs>.*)\)", re.IGNORECASE)
_sql_coldef_pattern = re.compile(r"\s*(?P<name>\w+)\s+(?P<type>\w+)[^,]*")


#
# Database info extraction utils
#


def get_table_names(connection):
	"""
	Return a list of the table names in the database.
	"""
	cursor = connection.cursor()
	cursor.execute("SELECT name FROM sqlite_master WHERE type == 'table'")
	return [name for (name,) in cursor]


def get_column_info(connection, table_name):
	"""
	Return an in order list of (name, type) tuples describing the
	columns in the given table.
	"""
	cursor = connection.cursor()
	cursor.execute("SELECT sql FROM sqlite_master WHERE type == 'table' AND name == ?", (table_name,))
	statement, = cursor.fetchone()
	coldefs = re.match(_sql_create_table_pattern, statement).groupdict()["coldefs"]
	return [(coldef.groupdict()["name"], coldef.groupdict()["type"]) for coldef in re.finditer(_sql_coldef_pattern, coldefs) if coldef.groupdict()["name"].upper() not in ("PRIMARY", "UNIQUE", "CHECK")]


def get_xml(connection, table_names = None):
	"""
	Construct an XML document tree wrapping around the contents of the
	database.  On success the return value is a ligolw.LIGO_LW element
	containing the tables as children.  Arguments are a connection to
	to a database, and an optional list of table names to dump.  If
	table_names is not provided the set is obtained from get_table_names()
	"""
	ligo_lw = ligolw.LIGO_LW()

	if table_names is None:
		table_names = get_table_names(connection)

	for table_name in table_names:
		# build the table document tree.  copied from
		# lsctables.New()
		try:
			cls = TableByName[table_name]
		except KeyError:
			cls = DBTable
		table_elem = cls(AttributesImpl({u"Name": u"%s:table" % table_name}), connection = connection)
		for column_name, column_type in get_column_info(connection, table_elem.Name):
			if table_elem.validcolumns is not None:
				# use the pre-defined column type
				column_type = table_elem.validcolumns[column_name]
			else:
				# guess the column type
				column_type = ligolwtypes.FromSQLiteType[column_type]
			table_elem.appendChild(table.Column(AttributesImpl({u"Name": u"%s:%s" % (table_name, column_name), u"Type": column_type})))
		table_elem._end_of_columns()
		table_elem.appendChild(table.TableStream(AttributesImpl({u"Name": u"%s:table" % table_name, u"Delimiter": table.TableStream.Delimiter.default, u"Type": table.TableStream.Type.default})))
		ligo_lw.appendChild(table_elem)
	return ligo_lw


#
# =============================================================================
#
#                            DBTable Element Class
#
# =============================================================================
#


class DBTable(table.Table):
	"""
	A special version of the Table class using an SQL database for
	storage.  Many of the features of the Table class are not available
	here, but instead the user can use SQL to query the table's
	contents.

	The constraints attribute can be set to a text string that will be
	added to the table's CREATE statement where constraints go, for
	example you might wish to set this to "PRIMARY KEY (event_id)" for
	a table with an event_id column.

	Note:  because the table is stored in an SQL database, the use of
	this class imposes the restriction that table names be unique
	within a document.

	Also note that at the present time there is really only proper
	support for the pre-defined tables in the lsctables module.  It is
	possible to load unrecognized tables into a database from LIGO
	Light Weight XML files, but without developer intervention there is
	no way to indicate the constraints that should be imposed on the
	columns, for example which columns should be used as primary keys
	and so on.  This can result in poor query performance.  It is also
	possible to extract a database' contents to a LIGO Light Weight XML
	file even when the database contains unrecognized tables, but
	without developer intervention the column types will be guessed
	using a generic mapping of SQL types to LIGO Light Weight types.

	Each instance of this class must be connected to a database.  The
	(Python DBAPI 2.0 compatible) connection object is passed to the
	class via the connection parameter at instance creation time.

	Example:

	>>> import sqlite3
	>>> connection = sqlite3.connection()
	>>> tbl = dbtables.DBTable(AttributesImpl({u"Name": u"process:table"}), connection = connection)

	A custom content handler must be created in order to pass the
	connection keyword argument to the DBTable class when instances are
	created, since the default content handler does not do this.  See
	the use_in() function defined in this module for information on how
	to create such a content handler

	If a custom pycbc.ligolw.Table subclass is defined in
	pycbc.ligolw.lsctables whose name matches the name of the DBTable
	being constructed, the lsctables class is added to the list of
	parent classes.  This allows the lsctables class' methods to be
	used with the DBTable instances but not all of the methods will
	necessarily work with the database-backed version of the class.
	Your mileage may vary.
	"""
	def __new__(cls, *args, **kwargs):
		# does this class already have table-specific metadata?
		if not hasattr(cls, "tableName"):
			# no, try to retrieve it from lsctables
			attrs, = args
			name = table.StripTableName(attrs[u"Name"])
			if name in lsctables.TableByName:
				# found metadata in lsctables, construct
				# custom subclass.  the class from
				# lsctables is added as a parent class to
				# allow methods from that class to be used
				# with this class, however there is no
				# guarantee that all parent class methods
				# will be appropriate for use with the
				# DB-backend object.
				lsccls = lsctables.TableByName[name]
				class CustomDBTable(cls, lsccls):
					tableName = lsccls.tableName
					validcolumns = lsccls.validcolumns
					loadcolumns = lsccls.loadcolumns
					constraints = lsccls.constraints
					next_id = lsccls.next_id
					RowType = lsccls.RowType
					how_to_index = lsccls.how_to_index

				# save for re-use (required for ID
				# remapping across multiple documents in
				# ligolw_sqlite)
				TableByName[name] = CustomDBTable

				# replace input argument with new class
				cls = CustomDBTable
		return table.Table.__new__(cls, *args)

	def __init__(self, *args, **kwargs):
		# chain to parent class
		table.Table.__init__(self, *args)

		# retrieve connection object from kwargs
		self.connection = kwargs.pop("connection")

		# pre-allocate a cursor for internal queries
		self.cursor = self.connection.cursor()

	def copy(self, *args, **kwargs):
		"""
		This method is not implemented.  See
		pycbc.ligolw.table.Table for more information.
		"""
		raise NotImplemented

	def _end_of_columns(self):
		table.Table._end_of_columns(self)
		# dbcolumnnames and types have the "not loaded" columns
		# removed
		if self.loadcolumns is not None:
			self.dbcolumnnames = [name for name in self.columnnames if name in self.loadcolumns]
			self.dbcolumntypes = [name for i, name in enumerate(self.columntypes) if self.columnnames[i] in self.loadcolumns]
		else:
			self.dbcolumnnames = self.columnnames
			self.dbcolumntypes = self.columntypes

		# create the table
		ToSQLType = {
			"sqlite": ligolwtypes.ToSQLiteType,
			"mysql": ligolwtypes.ToMySQLType
		}[connection_db_type(self.connection)]
		try:
			statement = "CREATE TABLE IF NOT EXISTS " + self.Name + " (" + ", ".join(map(lambda n, t: "%s %s" % (n, ToSQLType[t]), self.dbcolumnnames, self.dbcolumntypes))
		except KeyError as e:
			raise ValueError("column type '%s' not supported" % str(e))
		if self.constraints is not None:
			statement += ", " + self.constraints
		statement += ")"
		self.cursor.execute(statement)

		# record the highest internal row ID
		self.last_maxrowid = self.maxrowid() or 0

		# construct the SQL to be used to insert new rows
		params = {
			"sqlite": ",".join("?" * len(self.dbcolumnnames)),
			"mysql": ",".join(["%s"] * len(self.dbcolumnnames))
		}[connection_db_type(self.connection)]
		self.append_statement = "INSERT INTO %s (%s) VALUES (%s)" % (self.Name, ",".join(self.dbcolumnnames), params)
		self.append_attrgetter = operator.attrgetter(*self.dbcolumnnames)

	def _end_of_rows(self):
		# FIXME:  is this needed?
		table.Table._end_of_rows(self)
		self.connection.commit()

	def sync_next_id(self):
		if self.next_id is not None:
			max_id = idmap_get_max_id(self.connection, type(self.next_id))
			if max_id is None:
				self.set_next_id(type(self.next_id)(0))
			else:
				self.set_next_id(max_id + 1)
		return self.next_id

	def maxrowid(self):
		self.cursor.execute("SELECT MAX(ROWID) FROM %s" % self.Name)
		return self.cursor.fetchone()[0]

	def __len__(self):
		self.cursor.execute("SELECT COUNT(*) FROM %s" % self.Name)
		return self.cursor.fetchone()[0]

	def __iter__(self):
		cursor = self.connection.cursor()
		cursor.execute("SELECT * FROM %s" % self.Name)
		for values in cursor:
			yield self.row_from_cols(values)

	# FIXME:  is adding this a good idea?
	#def __delslice__(self, i, j):
	#	# sqlite numbers rows starting from 1:  [0:10] becomes
	#	# "rowid between 1 and 10" which means 1 <= rowid <= 10,
	#	# which is the intended range
	#	self.cursor.execute("DELETE FROM %s WHERE ROWID BETWEEN %d AND %d" % (self.Name, i + 1, j))

	def _append(self, row):
		"""
		Standard .append() method.  This method is for intended for
		internal use only.
		"""
		self.cursor.execute(self.append_statement, self.append_attrgetter(row))

	def _remapping_append(self, row):
		"""
		Replacement for the standard .append() method.  This
		version performs on the fly row ID reassignment, and so
		also performs the function of the updateKeyMapping()
		method.  SQLite does not permit the PRIMARY KEY of a row to
		be modified, so it needs to be done prior to insertion.
		This method is intended for internal use only.
		"""
		if self.next_id is not None:
			# assign (and record) a new ID before inserting the
			# row to avoid collisions with existing rows
			setattr(row, self.next_id.column_name, idmap_get_new(self.connection, getattr(row, self.next_id.column_name), self))
		self._append(row)

	append = _append

	def row_from_cols(self, values):
		"""
		Given an iterable of values in the order of columns in the
		database, construct and return a row object.  This is a
		convenience function for turning the results of database
		queries into Python objects.
		"""
		row = self.RowType()
		for c, t, v in zip(self.dbcolumnnames, self.dbcolumntypes, values):
			if t in ligolwtypes.IDTypes:
				v = ilwd.ilwdchar(v)
			setattr(row, c, v)
		return row
	# backwards compatibility
	_row_from_cols = row_from_cols

	def unlink(self):
		table.Table.unlink(self)
		self.connection = None
		self.cursor = None

	def applyKeyMapping(self):
		"""
		Used as the second half of the key reassignment algorithm.
		Loops over each row in the table, replacing references to
		old row keys with the new values from the _idmap_ table.
		"""
		assignments = ", ".join("%s = (SELECT new FROM _idmap_ WHERE old == %s)" % (colname, colname) for coltype, colname in zip(self.dbcolumntypes, self.dbcolumnnames) if coltype in ligolwtypes.IDTypes and (self.next_id is None or colname != self.next_id.column_name))
		if assignments:
			# SQLite documentation says ROWID is monotonically
			# increasing starting at 1 for the first row unless
			# it ever wraps around, then it is randomly
			# assigned.  ROWID is a 64 bit integer, so the only
			# way it will wrap is if somebody sets it to a very
			# high number manually.  This library does not do
			# that, so I don't bother checking.
			self.cursor.execute("UPDATE %s SET %s WHERE ROWID > %d" % (self.Name, assignments, self.last_maxrowid))
			self.last_maxrowid = self.maxrowid() or 0


#
# =============================================================================
#
#                                  LSC Tables
#
# =============================================================================
#


class ProcessParamsTable(DBTable):
	tableName = lsctables.ProcessParamsTable.tableName
	validcolumns = lsctables.ProcessParamsTable.validcolumns
	constraints = lsctables.ProcessParamsTable.constraints
	next_id = lsctables.ProcessParamsTable.next_id
	RowType = lsctables.ProcessParamsTable.RowType
	how_to_index = lsctables.ProcessParamsTable.how_to_index

	def append(self, row):
		if row.type is not None and row.type not in ligolwtypes.Types:
			raise ligolw.ElementError("unrecognized type '%s'" % row.type)
		DBTable.append(self, row)


class TimeSlideTable(DBTable):
	tableName = lsctables.TimeSlideTable.tableName
	validcolumns = lsctables.TimeSlideTable.validcolumns
	constraints = lsctables.TimeSlideTable.constraints
	next_id = lsctables.TimeSlideTable.next_id
	RowType = lsctables.TimeSlideTable.RowType
	how_to_index = lsctables.TimeSlideTable.how_to_index

	def as_dict(self):
		"""
		Return a ditionary mapping time slide IDs to offset
		dictionaries.
		"""
		return dict((ilwd.ilwdchar(id), offsetvector.offsetvector((instrument, offset) for id, instrument, offset in values)) for id, values in itertools.groupby(self.cursor.execute("SELECT time_slide_id, instrument, offset FROM time_slide ORDER BY time_slide_id"), lambda (id, instrument, offset): id))

	def get_time_slide_id(self, offsetdict, create_new = None, superset_ok = False, nonunique_ok = False):
		"""
		Return the time_slide_id corresponding to the offset vector
		described by offsetdict, a dictionary of instrument/offset
		pairs.

		If the optional create_new argument is None (the default),
		then the table must contain a matching offset vector.  The
		return value is the ID of that vector.  If the table does
		not contain a matching offset vector then KeyError is
		raised.

		If the optional create_new argument is set to a Process
		object (or any other object with a process_id attribute),
		then if the table does not contain a matching offset vector
		a new one will be added to the table and marked as having
		been created by the given process.  The return value is the
		ID of the (possibly newly created) matching offset vector.

		If the optional superset_ok argument is False (the default)
		then an offset vector in the table is considered to "match"
		the requested offset vector only if they contain the exact
		same set of instruments.  If the superset_ok argument is
		True, then an offset vector in the table is considered to
		match the requested offset vector as long as it provides
		the same offsets for the same instruments as the requested
		vector, even if it provides offsets for other instruments
		as well.

		More than one offset vector in the table might match the
		requested vector.  If the optional nonunique_ok argument is
		False (the default), then KeyError will be raised if more
		than one offset vector in the table is found to match the
		requested vector.  If the optional nonunique_ok is True
		then the return value is the ID of one of the matching
		offset vectors selected at random.
		"""
		# look for matching offset vectors
		if superset_ok:
			ids = [id for id, slide in self.as_dict().items() if offsetdict == dict((instrument, offset) for instrument, offset in slide.items() if instrument in offsetdict)]
		else:
			ids = [id for id, slide in self.as_dict().items() if offsetdict == slide]
		if len(ids) > 1:
			# found more than one
			if nonunique_ok:
				# and that's OK
				return ids[0]
			# and that's not OK
			raise KeyError(offsetdict)
		if len(ids) == 1:
			# found one
			return ids[0]
		# offset vector not found in table
		if create_new is None:
			# and that's not OK
			raise KeyError(offsetdict)
		# that's OK, create new vector
		id = self.get_next_id()
		for instrument, offset in offsetdict.items():
			row = self.RowType()
			row.process_id = create_new.process_id
			row.time_slide_id = id
			row.instrument = instrument
			row.offset = offset
			self.append(row)

		# return new ID
		return id


#
# =============================================================================
#
#                                Table Metadata
#
# =============================================================================
#


def build_indexes(connection, verbose = False):
	"""
	Using the how_to_index annotations in the table class definitions,
	construct a set of indexes for the database at the given
	connection.
	"""
	cursor = connection.cursor()
	for table_name in get_table_names(connection):
		# FIXME:  figure out how to do this extensibly
		if table_name in TableByName:
			how_to_index = TableByName[table_name].how_to_index
		elif table_name in lsctables.TableByName:
			how_to_index = lsctables.TableByName[table_name].how_to_index
		else:
			continue
		if how_to_index is not None:
			if verbose:
				print >>sys.stderr, "indexing %s table ..." % table_name
			for index_name, cols in how_to_index.iteritems():
				cursor.execute("CREATE INDEX IF NOT EXISTS %s ON %s (%s)" % (index_name, table_name, ",".join(cols)))
	connection.commit()


#
# =============================================================================
#
#                                Table Metadata
#
# =============================================================================
#


#
# Table name ---> table type mapping.
#


TableByName = {
	table.StripTableName(ProcessParamsTable.tableName): ProcessParamsTable,
	table.StripTableName(TimeSlideTable.tableName): TimeSlideTable
}


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
	pycbc.ligolw.LIGOLWContentHandler, to cause it to use the DBTable
	class defined in this module when parsing XML documents.  Instances
	of the class must provide a connection attribute.  When a document
	is parsed, the value of this attribute will be passed to the
	DBTable class' .__init__() method as each table object is created,
	and thus sets the database connection for all table objects in the
	document.

	Example:

	>>> import sqlite3
	>>> from pycbc.ligolw import ligolw
	>>> class MyContentHandler(ligolw.LIGOLWContentHandler):
	...	def __init__(self, *args):
	...		super(MyContentHandler, self).__init__(*args)
	...		self.connection = sqlite3.connection()
	...
	>>> use_in(MyContentHandler)

	Multiple database files can be in use at once by creating a content
	handler class for each one.
	"""
	ContentHandler = lsctables.use_in(ContentHandler)

	def startTable(self, parent, attrs):
		name = table.StripTableName(attrs[u"Name"])
		if name in TableByName:
			return TableByName[name](attrs, connection = self.connection)
		return DBTable(attrs, connection = self.connection)

	ContentHandler.startTable = startTable

	return ContentHandler

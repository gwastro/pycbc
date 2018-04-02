/*
 * Copyright (C) 2007-2009,2011  Kipp C. Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


/*
 * ============================================================================
 *
 *                         tokenizer.RowDumper Class
 *
 * ============================================================================
 */


#include <Python.h>
#include <structmember.h>
#include <stdlib.h>
#include <tokenizer.h>

/* Gain access to 64-bit addressing where possible
 * http://www.python.org/dev/peps/pep-0353/#conversion-guidelines */
#if PY_VERSION_HEX < 0x02050000 && !defined(PY_SSIZE_T_MIN)
typedef int Py_ssize_t;
#define PY_SSIZE_T_MAX INT_MAX
#define PY_SSIZE_T_MIN INT_MIN
#endif

/*
 * ============================================================================
 *
 *                              Row Dumper Type
 *
 * ============================================================================
 */


/*
 * Structure
 */


typedef struct {
	PyObject_HEAD
	/* delimiter character to be used in row construction */
	PyObject *delimiter;
	/* tuple of attribute names as Python strings */
	PyObject *attributes;
	/* tuple of row element format functions */
	PyObject *formats;
	/* the source of row objects to be turned to unicode strings */
	PyObject *iter;
	/* number of rows converted so far */
	Py_ssize_t rows_converted;
	/* tuple of unicode tokens from most recently converted row */
	PyObject *tokens;
} ligolw_RowDumper;


/*
 * __del__() method
 */


static void __del__(PyObject *self)
{
	ligolw_RowDumper *rowdumper = (ligolw_RowDumper *) self;

	Py_XDECREF(rowdumper->delimiter);
	Py_XDECREF(rowdumper->attributes);
	Py_XDECREF(rowdumper->formats);
	Py_XDECREF(rowdumper->iter);
	Py_XDECREF(rowdumper->tokens);

	self->ob_type->tp_free(self);
}


/*
 * __init__() method
 */


static int __init__(PyObject *self, PyObject *args, PyObject *kwds)
{
	ligolw_RowDumper *rowdumper = (ligolw_RowDumper *) self;
	Py_UNICODE default_delimiter = ',';

	rowdumper->delimiter = NULL;
	if(!PyArg_ParseTuple(args, "OO|U", &rowdumper->attributes, &rowdumper->formats, &rowdumper->delimiter))
		return -1;

	if(rowdumper->delimiter)
		Py_INCREF(rowdumper->delimiter);
	else
		rowdumper->delimiter = PyUnicode_FromUnicode(&default_delimiter, 1);
	rowdumper->attributes = llwtokenizer_build_attributes(rowdumper->attributes);
	rowdumper->formats = llwtokenizer_build_formats(rowdumper->formats);
	if(!rowdumper->delimiter || !rowdumper->attributes || !rowdumper->formats)
		/* memory clean-up happens in __del__() */
		return -1;

	if(PyTuple_GET_SIZE(rowdumper->attributes) != PyTuple_GET_SIZE(rowdumper->formats)) {
		/* memory clean-up happens in __del__() */
		PyErr_SetString(PyExc_ValueError, "len(attributes) != len(formats)");
		return -1;
	}

	rowdumper->rows_converted = 0;
	rowdumper->iter = Py_None;
	Py_INCREF(rowdumper->iter);
	rowdumper->tokens = Py_None;
	Py_INCREF(rowdumper->tokens);

	return 0;
}


/*
 * dump() method
 */


static PyObject *dump(PyObject *self, PyObject *iterable)
{
	ligolw_RowDumper *rowdumper = (ligolw_RowDumper *) self;
	PyObject *iter = PyObject_GetIter(iterable);

	if(!iter)
		return NULL;

	Py_DECREF(rowdumper->iter);
	rowdumper->iter = iter;

	Py_INCREF(self);
	return self;
}


/*
 * __iter__() method
 */


static PyObject *__iter__(PyObject *self)
{
	Py_INCREF(self);
	return self;
}


/*
 * next() method
 */


static PyObject *next(PyObject *self)
{
	ligolw_RowDumper *rowdumper = (ligolw_RowDumper *) self;
	const Py_ssize_t n = PyTuple_GET_SIZE(rowdumper->attributes);
	PyObject *row;
	PyObject *result;
	Py_ssize_t i;

	if(!PyIter_Check(rowdumper->iter)) {
		PyErr_SetObject(PyExc_TypeError, rowdumper->iter);
		return NULL;
	}
	row = PyIter_Next(rowdumper->iter);
	if(!row) {
		if(!PyErr_Occurred()) {
			Py_DECREF(rowdumper->iter);
			rowdumper->iter = Py_None;
			Py_INCREF(rowdumper->iter);
			PyErr_SetNone(PyExc_StopIteration);
		}
		return NULL;
	}

	Py_DECREF(rowdumper->tokens);
	rowdumper->tokens = PyTuple_New(n);
	if(!rowdumper->tokens) {
		rowdumper->tokens = Py_None;
		Py_INCREF(rowdumper->tokens);
		Py_DECREF(row);
		return NULL;
	}

	for(i = 0; i < n; i++) {
		PyObject *val = PyObject_GetAttr(row, PyTuple_GET_ITEM(rowdumper->attributes, i));
		PyObject *token;

		if(!val) {
			Py_DECREF(rowdumper->tokens);
			rowdumper->tokens = Py_None;
			Py_INCREF(rowdumper->tokens);
			Py_DECREF(row);
			return NULL;
		}

		if(val == Py_None)
			token = PyUnicode_FromUnicode(NULL, 0); /* u"" */
		else
			token = PyObject_CallFunctionObjArgs(PyTuple_GET_ITEM(rowdumper->formats, i), val, NULL);
		Py_DECREF(val);

		if(!token) {
			Py_DECREF(rowdumper->tokens);
			rowdumper->tokens = Py_None;
			Py_INCREF(rowdumper->tokens);
			Py_DECREF(row);
			return NULL;
		}

		PyTuple_SET_ITEM(rowdumper->tokens, i, token);
	}

	Py_DECREF(row);

	result = PyUnicode_Join(rowdumper->delimiter, rowdumper->tokens);

	rowdumper->rows_converted += result != NULL;

	return result;
}


/*
 * Type information
 */


static struct PyMemberDef members[] = {
	{"delimiter", T_OBJECT, offsetof(ligolw_RowDumper, delimiter), READONLY, "The delimiter as a unicode string."},
	{"attributes", T_OBJECT, offsetof(ligolw_RowDumper, attributes), READONLY, "In-order tuple of attribute names as strings."},
	{"formats", T_OBJECT, offsetof(ligolw_RowDumper, formats), READONLY, "In-order tuple of row element format functions."},
	{"iter", T_OBJECT, offsetof(ligolw_RowDumper, iter), 0, "The iterator being used to provide rows for conversion."},
	{"rows_converted", T_LONG, offsetof(ligolw_RowDumper, rows_converted), 0, "Count of rows converted."},
	{"tokens", T_OBJECT, offsetof(ligolw_RowDumper, tokens), READONLY, "In-order tuple of unicode tokens from most recently converted row."},
	{NULL,}
};


static struct PyMethodDef methods[] = {
	{"dump", dump, METH_O, "Set the Python iterable from which row objects will be retrieved for dumping."},
	{NULL,}
};


PyTypeObject ligolw_RowDumper_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(ligolw_RowDumper),
	.tp_dealloc = __del__,
	.tp_doc =
"An iterator for converting row objects into string tokens.\n" \
"\n" \
"Example:\n" \
"\n" \
">>> class Row(object):\n" \
"...     pass\n" \
"... \n" \
">>> rows = [Row(), Row(), Row()]\n" \
">>> rows[0].snr = 10.1\n" \
">>> rows[1].snr = 15.2\n" \
">>> rows[2].snr = 20.3\n" \
">>> rows[0].status = \"bad\"\n" \
">>> rows[1].status = \"bad\"\n" \
">>> rows[2].status = \"good\"\n" \
">>> rowdumper = RowDumper((\"snr\", \"status\"), (\"%.16g\".__mod__, \"\\\"%s\\\"\".__mod__))\n" \
">>> for line in rowdumper.dump(rows):\n" \
"...     print line\n" \
"... \n" \
"10.1,\"bad\"\n" \
"15.2,\"bad\"\n" \
"20.3,\"good\"\n" \
"\n" \
"An instance of RowDumper is initialized with two arguments and an optional\n" \
"third argument.  The first argument is a sequence of attribute names.  The\n" \
"second argument is a sequence of Python format strings.  The third, optional,\n" \
"argument is the unicode string to use as the delimiter between tokens (the\n" \
"default is u\",\").  The row dumper is started by calling the .dump() method\n" \
"which takes a Python iterable as its single argument.  After the .dump()\n" \
"method has been called, when a RowDumper instance is iterated over it\n" \
"retrieves objects, one-by-one, from the iterable passed to the .dump() method\n" \
"and yields a sequence of unicode strings containing the delimited string\n" \
"representations of the values of the attributes of those objects.  The\n" \
"attribute values are printed in the order specified when the RowDumper was\n" \
"created, and using the formats specified.  An attribute whose value is None\n" \
"is printed as an empty string regardless of the requested format.",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
	.tp_init = __init__,
	.tp_iter = __iter__,
	.tp_iternext = next,
	.tp_members = members,
	.tp_methods = methods,
	.tp_name = MODULE_NAME ".RowDumper",
	.tp_new = PyType_GenericNew,
};

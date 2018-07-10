/*
 * Copyright (C) 2007-2009  Kipp C. Cannon
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
 *                         tokenizer.RowBuilder Class
 *
 * ============================================================================
 */


#include <Python.h>
#include <structmember.h>
#include <stdlib.h>
#include <tokenizer.h>


/*
 * ============================================================================
 *
 *                              Row Builder Type
 *
 * ============================================================================
 */


/*
 * Structure
 */


typedef struct {
	PyObject_HEAD
	/* class to be instantiated for new rows */
	PyTypeObject *rowtype;
	/* tuple of attribute names */
	PyObject *attributes;
	/* tuple of internable attributes */
	PyObject *interns;
	/* current row */
	PyObject *row;
	/* current attribute index */
	int i;
	/* the iterable passed to append() */
	PyObject *iter;
} ligolw_RowBuilder;


/*
 * append() method
 */


static PyObject *append(PyObject *self, PyObject *iter)
{
	ligolw_RowBuilder *rowbuilder = (ligolw_RowBuilder *) self;

	Py_XDECREF(rowbuilder->iter);
	rowbuilder->iter = PyObject_GetIter(iter);
	if(!rowbuilder->iter)
		return NULL;

	Py_INCREF(self);
	return self;
}


/*
 * __del__() method
 */


static void __del__(PyObject *self)
{
	ligolw_RowBuilder *rowbuilder = (ligolw_RowBuilder *) self;

	Py_XDECREF(rowbuilder->rowtype);
	Py_XDECREF(rowbuilder->attributes);
	Py_XDECREF(rowbuilder->interns);
	Py_XDECREF(rowbuilder->row);
	Py_XDECREF(rowbuilder->iter);

	self->ob_type->tp_free(self);
}


/*
 * __init__() method
 */


static int __init__(PyObject *self, PyObject *args, PyObject *kwds)
{
	ligolw_RowBuilder *rowbuilder = (ligolw_RowBuilder *) self;

	rowbuilder->interns = NULL;
	if(!PyArg_ParseTuple(args, "OO|O", &rowbuilder->rowtype, &rowbuilder->attributes, &rowbuilder->interns))
		return -1;

	Py_INCREF(rowbuilder->rowtype);

	rowbuilder->attributes = llwtokenizer_build_attributes(rowbuilder->attributes);
	if(rowbuilder->interns)
		rowbuilder->interns = PySequence_Tuple(rowbuilder->interns);
	else
		rowbuilder->interns = PyTuple_New(0);

	if(!rowbuilder->attributes || !rowbuilder->interns)
		return -1;

	rowbuilder->row = Py_None;
	Py_INCREF(rowbuilder->row);
	rowbuilder->i = 0;
	rowbuilder->iter = NULL;

	return 0;
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
	ligolw_RowBuilder *rowbuilder = (ligolw_RowBuilder *) self;
	PyObject *item;

	if(!rowbuilder->iter) {
		PyErr_SetNone(PyExc_StopIteration);
		return NULL;
	}

	while((item = PyIter_Next(rowbuilder->iter))) {
		int result;
		if(rowbuilder->row == Py_None) {
			rowbuilder->row = PyType_GenericNew(rowbuilder->rowtype, NULL, NULL);
			if(!rowbuilder->row) {
				rowbuilder->row = Py_None;
				return NULL;
			}
			Py_DECREF(Py_None);
		}
		result = PyObject_SetAttr(rowbuilder->row, PyTuple_GET_ITEM(rowbuilder->attributes, rowbuilder->i), item);
		Py_DECREF(item);
		if(result < 0)
			return NULL;
		if(++rowbuilder->i >= PyTuple_GET_SIZE(rowbuilder->attributes)) {
			PyObject *row = rowbuilder->row;
			rowbuilder->row = Py_None;
			Py_INCREF(rowbuilder->row);
			rowbuilder->i = 0;
			return row;
		}
	}

	if(!PyErr_Occurred()) {
		PyErr_SetNone(PyExc_StopIteration);
		Py_DECREF(rowbuilder->iter);
		rowbuilder->iter = NULL;
	}

	return NULL;
}


/*
 * Type information
 */


static struct PyMemberDef members[] = {
	{"rowtype", T_OBJECT, offsetof(ligolw_RowBuilder, rowtype), 0, "row class"},
	{"attributes", T_OBJECT, offsetof(ligolw_RowBuilder, attributes), READONLY, "in-order tuple of attribute names"},
	{"interns", T_OBJECT, offsetof(ligolw_RowBuilder, interns), 0, "names of attributes suitable for interning"},
	{"row", T_OBJECT, offsetof(ligolw_RowBuilder, row), 0, "current row object"},
	{"i", T_INT, offsetof(ligolw_RowBuilder, i), 0, "current attribute index"},
	{NULL,}
};


static struct PyMethodDef methods[] = {
	{"append", append, METH_O,
"Append a sequence of tokens to the row builder, returning an iterator for\n"\
"generating a sequence of new row instances.  The tokens argument should be\n"\
"an iterable, producing a sequence of token objects.  If fewer tokens are\n"\
"yielded from the iterable than are required to construct a complete row,\n"\
"then the row is stored in its partially-populated state and its\n"\
"construction will continue upon the next invocation.  Note that it is\n"\
"possible that a call to this method will yield no new rows at all.\n"\
"\n"\
"Example:\n"\
"\n"\
">>> for row in rows.append([10, 6.8, 15, 29.1]):\n"\
"...     print row.snr\n"\
"..."
	},
	{NULL,}
};


PyTypeObject ligolw_RowBuilder_Type = {
	PyObject_HEAD_INIT(NULL)
	.tp_basicsize = sizeof(ligolw_RowBuilder),
	.tp_dealloc = __del__,
	.tp_doc =
"This class provides the logic required to transform a sequence of of tokens\n"\
"parsed out of the delimited text of a Stream element into a sequence of row\n"\
"objects for insertion into a Table element.  An instance of this class is\n"\
"initialized with a Python class to be instantiated to form row objects,\n"\
"and an iterable providing the names of the row class' attributes to which\n"\
"tokens will be assigned in order.\n"\
"\n"\
"Example:\n"\
"\n"\
">>> import tokenizer\n"\
">>> class Row(object):\n"\
"...     pass\n"\
"...\n"\
">>> t = tokenizer.Tokenizer(u\",\")\n"\
">>> t.set_types([int, float])\n"\
">>> rows = tokenizer.RowBuilder(Row, [\"time\", \"snr\"])\n"\
">>> l = list(rows.append(t.append(u\"10,6.8,15,29.1,\")))\n"\
">>> l[0].snr\n"\
"6.8\n"\
">>> l[1].time\n"\
"15\n"\
"\n"\
"Hint:  If you wish to try to save memory by \"interning\" the values in\n"\
"certain columns of a table, try sub-classing this and replacing the append\n"\
"method with your own.\n"\
"\n"\
"Example:\n"\
"\n"\
">>> strings = {}\n"\
">>> OldRowBuilder = RowBuilder\n"\
">>> class MyRowBuilder(RowBuilder):\n"\
"...     def append(self, tokens):\n"\
"...             for row in OldRowBuilder.append(self, tokens):\n"\
"...                     if hasattr(row, \"channel\"):\n"\
"...                             row.channel = strings.setdefault(row.channel, row.channel)\n"\
"...                     yield row\n"\
"...\n"\
">>> RowBuilder = MyRowBuilder\n"\
"\n"\
"This will significantly slow down table parsing, but for now this approach\n"\
"of allowing individual applications to override row construction on an\n"\
"as-desired basis seems to be the best way to implement the feature without\n"\
"adding a great deal of complexity.  Note that when initialized the\n"\
"RowBuilder class is passed the interns argument which is an iterable of\n"\
"attribute names that should be considered for interning.  These names come\n"\
"from hints stored in the Table class definitions, and will have been\n"\
"filtered so that only names corresponding to columns actually in the table\n"\
"will be listed.",
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
	.tp_init = __init__,
	.tp_iter = __iter__,
	.tp_iternext = next,
	.tp_members = members,
	.tp_methods = methods,
	.tp_name = MODULE_NAME ".RowBuilder",
	.tp_new = PyType_GenericNew,
};

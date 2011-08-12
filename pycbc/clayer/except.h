// Copyright (C) 2011 Karsten Wiesner
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
// Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


//
// =============================================================================
//
//                                   Preamble
//
// =============================================================================
//
// headerfile to extend C by exceptions

#include <setjmp.h>

extern jmp_buf pycbc_clayer_exception_buffer;
extern int pycbc_clayer_exception_status;
extern char pycbc_clayer_exception_message[2048];

#define try if ((pycbc_clayer_exception_status = \
setjmp(pycbc_clayer_exception_buffer)) == 0)

#define catch(val) else if (pycbc_clayer_exception_status == val)

#define throw(val, msg) \
{snprintf( pycbc_clayer_exception_message, \
 sizeof(pycbc_clayer_exception_message)/\
 sizeof(*pycbc_clayer_exception_message), msg); \
 longjmp(pycbc_clayer_exception_buffer,val);}

#define finally else

// Exception codes (according to SWIG_exception())
#define IOError             1
#define RuntimeError        2
#define IndexError          3
#define TypeError           4
#define DivisionByZero      5
#define OverflowError       6
#define SyntaxError         7
#define ValueError          8
#define SystemError         9
#define AttributeError     10
#define MemoryError        11
#define NullReferenceError 12

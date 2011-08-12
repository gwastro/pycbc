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
// swigfile to extend C by exceptions

%include "exception.i"

%exception {
    try {
        $action
    } catch(IOError) {
        SWIG_exception(SWIG_IOError, pycbc_clayer_exception_message);
    } catch(RuntimeError) {
        SWIG_exception(SWIG_RuntimeError, pycbc_clayer_exception_message);
    } catch(IndexError) {
        SWIG_exception(SWIG_IndexError, pycbc_clayer_exception_message);
    } catch(TypeError) {
        SWIG_exception(SWIG_TypeError, pycbc_clayer_exception_message);
    } catch(DivisionByZero) {
        SWIG_exception(SWIG_DivisionByZero, pycbc_clayer_exception_message);
    } catch(OverflowError) {
        SWIG_exception(SWIG_OverflowError, pycbc_clayer_exception_message);
    } catch(SyntaxError) {
        SWIG_exception(SWIG_SyntaxError, pycbc_clayer_exception_message);
    } catch(ValueError) {
        SWIG_exception(SWIG_ValueError, pycbc_clayer_exception_message);
    } catch(SystemError) {
        SWIG_exception(SWIG_SystemError, pycbc_clayer_exception_message);
    } catch(AttributeError) {
        SWIG_exception(SWIG_AttributeError, pycbc_clayer_exception_message);
    } catch(MemoryError) {
        SWIG_exception(SWIG_MemoryError, pycbc_clayer_exception_message);
    } catch(NullReferenceError) {
        SWIG_exception(SWIG_NullReferenceError, pycbc_clayer_exception_message);
    } finally {
        SWIG_exception(SWIG_UnknownError, pycbc_clayer_exception_message);
    }
}

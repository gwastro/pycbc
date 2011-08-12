%include "exception.i"

%exception {
    try {
        $action
    } catch(IOError) {
        SWIG_exception(SWIG_IOError, "IO Error");
    } catch(RuntimeError) {
        SWIG_exception(SWIG_RuntimeError, "Runtime Error");
    } catch(IndexError) {
        SWIG_exception(SWIG_IndexError, "Index Error");
    } catch(TypeError) {
        SWIG_exception(SWIG_TypeError, "Type Error");
    } catch(DivisionByZero) {
        SWIG_exception(SWIG_DivisionByZero, "Division By Zero");
    } catch(OverflowError) {
        SWIG_exception(SWIG_OverflowError, "Overflow Error");
    } catch(SyntaxError) {
        SWIG_exception(SWIG_SyntaxError, "Syntax Error");
    } catch(ValueError) {
        SWIG_exception(SWIG_ValueError, "Value Error");
    } catch(SystemError) {
        SWIG_exception(SWIG_SystemError, "System Error");
    } catch(AttributeError) {
        SWIG_exception(SWIG_AttributeError, "Attribute Error");
    } catch(MemoryError) {
        SWIG_exception(SWIG_MemoryError, "Memory Error");
    } catch(NullReferenceError) {
        SWIG_exception(SWIG_NullReferenceError, "Null Reference Error");
    } finally {
        SWIG_exception(SWIG_UnknownError, "Unknown error");
    }
}

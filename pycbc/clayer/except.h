

#ifndef PYCBC_EXCEPT_H
#define PYCBC_EXCEPT_H

/* error codes */
#define PYCBC_NO_ERROR            0        /* no error */
#define PYCBC_ATTRIBUTE_ERROR     1        /* invalid attribute accessed */
#define PYCBC_EOF_ERROR           2        /* end of file on io */
#define PYCBC_IO_ERROR            3        /* io error */
#define PYCBC_INDEX_ERROR         4        /* subscript out of range */
#define PYCBC_TYPE_ERROR          5        /* invalid type */
#define PYCBC_VALUE_ERROR         6        /* right type, wrong value */
#define PYCBC_MEMORY_ERROR        7        /* out of memory */
#define PYCBC_NAME_ERROR          8        /* name not found */
#define PYCBC_OVERFLOW_ERROR      9        /* artithmetic operation too large */
#define PYCBC_ZERO_DIVISION_ERROR 10       /* divide by zero */
#define PYCBC_RUNTIME_ERROR       11       /* everything else */
#define PYCBC_STOP_ITERATION       12      

#define EXCEPTION_MESSAGE_SIZE 255

void pycbc_throw_exception(int num,...);
void pycbc_throw_exception_bare(int num, char* msg);
void pycbc_clear_exception(void);
int  pycbc_check_exception(void);
char* pycbc_get_error_message(void);

#endif

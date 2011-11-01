#include <stdio.h>
#include <stdarg.h>
#include <pycbc/clayer/except.h>

static char error_message[EXCEPTION_MESSAGE_SIZE+1];
static int error_status = 0;

void pycbc_throw_exception_bare(int num, char* msg)
{
  strncpy(error_message,msg,EXCEPTION_MESSAGE_SIZE);
  error_status = num;
}

void pycbc_throw_exception(int error_num, ...)
{
	va_list argp;
	va_start(argp, error_num);
	char* fmt=va_arg(argp,char*);
	
    // If the first argument is a c-string format with printf
	if (fmt != (char*)0)
	  vsprintf(error_message, fmt, argp);
	  
	va_end(argp); 
  error_status = error_num;
}

void pycbc_clear_exception()
{
  error_status = PYCBC_NO_ERROR;
}

int pycbc_check_exception() 
{
  return error_status;
}

char* pycbc_get_error_message()
{
  return error_message;
}

#include <Python.h>
#include <pycbc/clayer/log.h>
#include <stdio.h>
#include <stdarg.h> 

char logger_name[LOGGER_NAME_SIZE+1];

void pycbc_set_logger(char* name){
    strncpy(logger_name,name,LOGGER_NAME_SIZE);
}

void pycbc_clear_logger(){
    logger_name[0]='\0';
}

void pycbc_log(int level, char* fmt,...){
    char log_message[LOG_MESSAGE_SIZE+1];

    //format the message string
    va_list argp;
	va_start(argp, fmt);
	vsprintf(log_message, fmt, argp);  
	va_end(argp); 

    //import the logger
    PyObject* module = PyImport_ImportModule("logging");
     
    if (logger_name){
        //Use the given logger
        PyObject_CallMethod(module,"basicConfig","");
        PyObject* logger = PyObject_CallMethod(module,"getLogger","s",logger_name);   
        PyObject_CallMethod(logger,"log","is",level,log_message);
        Py_DECREF(logger);
    }
    else{
        //Use the root logger by default
        PyObject_CallMethod(module,"log","is",level,log_message);
    }


    //Deallocate the python objects 
    Py_DECREF(module);

}

#ifndef PYCBC_LOGGING_H
#define PYCBC_LOGGING_H

#define NOTSET    0
#define DEBUG    10
#define INFO     20
#define WARNING  30
#define ERROR    40
#define CRITICAL 50

#define LOG_MESSAGE_SIZE 255
#define LOGGER_NAME_SIZE 255

/* error codes */

void pycbcset_logger(char* name);
void pycbc_log(int level,char* fmt, ...);
void pycbc_clear_logger(void);

#endif

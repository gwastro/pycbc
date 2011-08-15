/*#######################################
#
# Various utilities for OpenCL programming
#  o Logging and log levels
#  o Error checking and exit options
#  o Platform, context and device initialisation
#
#           -- Gergely Debreczeni
#
#########################################*/

#ifndef _GCLUTILS_H
#define _GCLUTILS_H

#include <stdio.h>


// All OpenCL headers
#if defined (__APPLE__) || defined(MACOSX)
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif 


#include <string.h>

// Definition of global variables

#define gclMaxNumberOfErrorCodes  255
#define gclShortErrorMessageSize  100
#define gclMaxLogMessageSize      1000
#define gclMaxPlatforms           100
#define gclMaxPrograms            100

//
#ifndef GCL_DEVICE_TYPE
#define GCL_DEVICE_TYPE  CL_DEVICE_TYPE_GPU
#endif

// Variables for platform, context and device

cl_int                 gcl_err             = CL_SUCCESS;
cl_uint                gcl_numPlatforms;
cl_platform_id         gcl_platform        = NULL;
cl_platform_id*        gcl_platforms       = new cl_platform_id[gclMaxPlatforms];
size_t                 gcl_devicenumber;
cl_device_id*          gcl_devices;
cl_device_id           gcl_available_device= NULL;
cl_context_properties  gcl_cps[3]          = {CL_CONTEXT_PLATFORM, (cl_context_properties) gcl_platform, 0};
cl_context_properties* gcl_cprops;
cl_context             gcl_context;
cl_program             gcl_programs[gclMaxPrograms]; 

typedef enum gclContextDevices
{
    gclAllGPUs,
    gclFirstAvailable
} gclContextDeviceMode;

// Variables for logging and error messages

enum gclLogLevels  {gclALL = 6, gclDEBUG = 5, gclINFO = 4, gclWARN =3, gclERR = 2, gclCRIT = 1, gclNONE = 0 };
int  gclLogLevel = gclINFO;

char * gclErrorCodes[gclMaxNumberOfErrorCodes];
char * gcllogmessage;
char * gclLogLevelNames[7];

void gclInitErrorMessages() {
    for (int i = 0; i < gclMaxNumberOfErrorCodes; i++) {
        gclErrorCodes[i] = (char *) malloc(gclShortErrorMessageSize);
    }
    strcpy(gclErrorCodes[0],"CL_SUCCESS");
    strcpy(gclErrorCodes[1],"CL_DEVICE_NOT_FOUND");
    strcpy(gclErrorCodes[2],"CL_DEVICE_NOT_AVAILABLE");
    strcpy(gclErrorCodes[3],"CL_COMPILER_NOT_AVAILABLE");
    strcpy(gclErrorCodes[4],"CL_MEM_OBJECT_ALLOCATION_FAILURE");
    strcpy(gclErrorCodes[5],"CL_OUT_OF_RESOURCES");
    strcpy(gclErrorCodes[6],"CL_OUT_OF_HOST_MEMORY");
    strcpy(gclErrorCodes[7],"CL_PROFILING_INFO_NOT_AVAILABLE");
    strcpy(gclErrorCodes[8],"CL_MEM_COPY_OVERLAP");
    strcpy(gclErrorCodes[9],"CL_IMAGE_FORMAT_MISMATCH");
    strcpy(gclErrorCodes[10],"CL_IMAGE_FORMAT_NOT_SUPPORTED");
    strcpy(gclErrorCodes[11],"CL_BUILD_PROGRAM_FAILURE");
    strcpy(gclErrorCodes[12],"CL_MAP_FAILURE");
    strcpy(gclErrorCodes[13],"CL_MISALIGNED_SUB_BUFFER_OFFSET");
    strcpy(gclErrorCodes[14],"CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST");
    strcpy(gclErrorCodes[30],"CL_INVALID_VALUE");
    strcpy(gclErrorCodes[31],"CL_INVALID_DEVICE_TYPE");
    strcpy(gclErrorCodes[32],"CL_INVALID_PLATFORM");
    strcpy(gclErrorCodes[33],"CL_INVALID_DEVICE");
    strcpy(gclErrorCodes[34],"CL_INVALID_CONTEXT");
    strcpy(gclErrorCodes[35],"CL_INVALID_QUEUE_PROPERTIES");
    strcpy(gclErrorCodes[36],"CL_INVALID_COMMAND_QUEUE");
    strcpy(gclErrorCodes[37],"CL_INVALID_HOST_PTR");
    strcpy(gclErrorCodes[38],"CL_INVALID_MEM_OBJECT");
    strcpy(gclErrorCodes[39],"CL_INVALID_IMAGE_FORMAT_DESCRIPTOR");
    strcpy(gclErrorCodes[40],"CL_INVALID_IMAGE_SIZE");
    strcpy(gclErrorCodes[41],"CL_INVALID_SAMPLER");
    strcpy(gclErrorCodes[42],"CL_INVALID_BINARY");
    strcpy(gclErrorCodes[43],"CL_INVALID_BUILD_OPTIONS");
    strcpy(gclErrorCodes[44],"CL_INVALID_PROGRAM");
    strcpy(gclErrorCodes[45],"CL_INVALID_PROGRAM_EXECUTABLE");
    strcpy(gclErrorCodes[46],"CL_INVALID_KERNEL_NAME");
    strcpy(gclErrorCodes[47],"CL_INVALID_KERNEL_DEFINITION");
    strcpy(gclErrorCodes[48],"CL_INVALID_KERNEL");
    strcpy(gclErrorCodes[49],"CL_INVALID_ARG_INDEX");
    strcpy(gclErrorCodes[50],"CL_INVALID_ARG_VALUE");
    strcpy(gclErrorCodes[51],"CL_INVALID_ARG_SIZE");
    strcpy(gclErrorCodes[52],"CL_INVALID_KERNEL_ARGS");
    strcpy(gclErrorCodes[53],"CL_INVALID_WORK_DIMENSION");
    strcpy(gclErrorCodes[54],"CL_INVALID_WORK_GROUP_SIZE");
    strcpy(gclErrorCodes[55],"CL_INVALID_WORK_ITEM_SIZE");
    strcpy(gclErrorCodes[56],"CL_INVALID_GLOBAL_OFFSET");
    strcpy(gclErrorCodes[57],"CL_INVALID_EVENT_WAIT_LIST");
    strcpy(gclErrorCodes[58],"CL_INVALID_EVENT");
    strcpy(gclErrorCodes[59],"CL_INVALID_OPERATION");
    strcpy(gclErrorCodes[60],"CL_INVALID_GL_OBJECT");
    strcpy(gclErrorCodes[61],"CL_INVALID_BUFFER_SIZE");
    strcpy(gclErrorCodes[62],"CL_INVALID_MIP_LEVEL");
    strcpy(gclErrorCodes[63],"CL_INVALID_GLOBAL_WORK_SIZE");
}

void gclInitLogMessages() {
    gcllogmessage = (char *) malloc(gclMaxLogMessageSize);
    for (int i = 0 ; i< 7; i ++ ) {
        gclLogLevelNames[i] = (char *) malloc(10);
    }
    strcpy(gclLogLevelNames[0],"NONE");
    strcpy(gclLogLevelNames[1],"CRIT");
    strcpy(gclLogLevelNames[2],"ERROR");
    strcpy(gclLogLevelNames[3],"WARN");
    strcpy(gclLogLevelNames[4],"INFO");
    strcpy(gclLogLevelNames[5],"DEBUG");
    strcpy(gclLogLevelNames[6],"ALL");
}

void gclLog(int debuglevel, char * message) {
    if (debuglevel <= gclLogLevel) fprintf(stderr,"%s: %s\n",gclLogLevelNames[debuglevel], message);
}

void gclCheckError(int err, char * desc, bool exitOnError = false) {
    int i;

    sprintf(gcllogmessage,"Checking error code: %s (%d)", desc, err);
    gclLog(gclDEBUG,gcllogmessage);

    for (i = 0; i < 13; i++) {
        if ((abs(err) == i ) && ( i != CL_SUCCESS)) {
            fprintf(stderr,"ERROR: %s %s (%d) \n",desc, gclErrorCodes[i], err);
        }
    }
    for (i = 30; i < 64; i++) {
        if ((abs(err) == i ) && ( i != CL_SUCCESS)) {
            fprintf(stderr,"ERROR: %s %s (%d) \n",desc, gclErrorCodes[i], err);
        }
    }

    if ( exitOnError && err != CL_SUCCESS ) {
        sprintf(gcllogmessage,"--- Exiting.\n");
        gclLog(gclERR,gcllogmessage);
        exit (-1);
    }

}


cl_int gclInitGPU(gclContextDeviceMode mode = gclAllGPUs) {

    sprintf(gcllogmessage,"Entering to gpuInit.");
    gclLog(gclDEBUG,gcllogmessage);

    // Start of platform layer

    gcl_err = clGetPlatformIDs(0, NULL, &gcl_numPlatforms);
    gclCheckError(gcl_err, "clGetPlatformIDs(gcl_numPlatforms)");

    if (0 < gcl_numPlatforms)
    {
        gcl_err = clGetPlatformIDs(gcl_numPlatforms, gcl_platforms, NULL);
        gclCheckError(gcl_err, "clGetPlatformIDs(platforms)");
        gcl_platform = gcl_platforms[0];
    } else {return -1;}

    //Start of context layer and device assignment

    gcl_cps[0] = CL_CONTEXT_PLATFORM;
    gcl_cps[1] = (cl_context_properties) gcl_platform;
    gcl_cps[2] = 0;
    gcl_cprops = (NULL == gcl_platform) ? NULL : gcl_cps;

    if (mode == gclAllGPUs)
    {
        gcl_context = clCreateContextFromType(gcl_cprops, GCL_DEVICE_TYPE, NULL, NULL, &gcl_err);
        gclCheckError(gcl_err, "clCreateContextFromType()");
        
        clGetContextInfo(gcl_context, CL_CONTEXT_DEVICES, 0, NULL, &gcl_devicenumber);
        gcl_devices = (cl_device_id*)malloc(gcl_devicenumber);
        clGetContextInfo(gcl_context, CL_CONTEXT_DEVICES, gcl_devicenumber, gcl_devices, NULL);
        gclCheckError( gcl_devicenumber != 0 ? CL_SUCCESS : -1, "devicenumber <= 0");
    } 
    else //if (mode == gclFirstAvailable)
    {
        cl_uint numDevices;
        gcl_err = clGetDeviceIDs(gcl_platform, CL_DEVICE_TYPE_GPU, 0, NULL, &numDevices);
        gclCheckError(gcl_err, "cgGetDeviceIDs()");
        
        gcl_devices = (cl_device_id*)malloc(sizeof(cl_device_id) * numDevices);
        gcl_err = clGetDeviceIDs(gcl_platform, CL_DEVICE_TYPE_GPU, numDevices, gcl_devices, NULL);
        char name[200];
        for (unsigned int i = 0; i < numDevices; ++i)
        {
            cl_bool available;
            gcl_err = clGetDeviceInfo(gcl_devices[i], CL_DEVICE_AVAILABLE, sizeof(cl_bool), &available, NULL);
            gclCheckError(gcl_err, "clGetDeviceInfo()");
            
            if(available)
            {
                gcl_available_device = gcl_devices[i];
                gcl_err = clGetDeviceInfo(gcl_devices[i], CL_DEVICE_NAME, sizeof(name), name, NULL);
                sprintf(gcllogmessage,"Found available device: %s", name);
                gclLog(gclINFO, gcllogmessage);
                break;
            } 
            else
            {
                gcl_err = clGetDeviceInfo(gcl_devices[i], CL_DEVICE_NAME, sizeof(name), name, NULL);
                if(gcl_err == CL_SUCCESS)
                    sprintf(gcllogmessage, "Device %s not available for compute.", name);
                else
                    sprintf(gcllogmessage, "Device #%d not available for compute.", i);
                gclLog(gclINFO, gcllogmessage);
            }
        }
        if (gcl_available_device == NULL)
        {
            gclLog(gclERR, "None of the devices were available for compute. The program will abort...");
            exit(-1);
        }
        gcl_context = clCreateContext(gcl_cprops, 1, &gcl_available_device, NULL, NULL, &gcl_err);
        gclCheckError(gcl_err, "clCreateContext()", true);
    }
    

    // Exiting
    sprintf(gcllogmessage,"Leaving from gpuInit.");
    gclLog(gclDEBUG,gcllogmessage);

    return gcl_err;
}


cl_int gclCheckBuild(cl_int err, uint prognum) {

    sprintf(gcllogmessage,"Entering to gclCheckBuild");
    gclLog(gclDEBUG,gcllogmessage);
    gclCheckError(err, "clBuildProgram()", false);
    if (err != CL_BUILD_SUCCESS) {
        /*          size_t        deviceNum;
        cl_device_id  *cdDevices;
        gcl_err = clGetContextInfo(gcl_context, CL_CONTEXT_DEVICES, 0, NULL, &deviceNum);
        gclCheckError(gcl_err, "clGetContextInfo3");

        cdDevices = (cl_device_id *) malloc(deviceNum);
        gcl_err  = clGetContextInfo(gcl_context, CL_CONTEXT_DEVICES, deviceNum, cdDevices, NULL);
        gclCheckError(gcl_err, "clGetContextInfo4");
*/        
        size_t LogSize;

        gcl_err = clGetProgramBuildInfo(gcl_programs[prognum], gcl_devices[0], CL_PROGRAM_BUILD_LOG , 0 , NULL, &LogSize);
        gclCheckError(err, "clGetProgramBuildInfo1");

        char *  LogMessage =  (char *) malloc(LogSize);

        gcl_err = clGetProgramBuildInfo(gcl_programs[prognum], gcl_devices[0], CL_PROGRAM_BUILD_LOG , LogSize, LogMessage, NULL);
        gclCheckError(err, "clGetProgramBuildInfo2()");

        sprintf(gcllogmessage,"The build log message is: \n");
        gclLog(gclERR,gcllogmessage);
        gclLog(gclERR, LogMessage);
        return -1;
    }
    sprintf(gcllogmessage,"Leaving from gclCheckBuild");
    gclLog(gclDEBUG,gcllogmessage);
    return err;
}

#endif


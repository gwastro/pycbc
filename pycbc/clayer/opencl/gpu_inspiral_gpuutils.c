#include <stdio.h>
#include <string.h>
#include <CL/opencl.h>
#include "pycbcopencl_types.h"

cl_int gpuinsp_InitGPU(cl_context_t* c, unsigned device_id)
{
  int err;
  cl_uint 		 numPlatforms;
  cl_platform_id         Platform             = 0;
  cl_platform_id*	 Platforms            = new cl_platform_id[10];
  cl_context_properties  cps[3]               = {CL_CONTEXT_PLATFORM, (cl_context_properties) Platform, 0};
  cl_context_properties* cprops;
  cl_uint                numDevices;
  cl_device_id*          Devices;
  cl_device_id           AvailableDevice      = NULL;
  char                   DeviceName[200];

//Getting the platforms
  err = clGetPlatformIDs(0, NULL, &numPlatforms);
  printf("-1Err: %d\n",err);
  printf("BlaBla\n");
  if (0 < numPlatforms)
    {
        err = clGetPlatformIDs(numPlatforms, Platforms, NULL);
        printf("0Err: %d\n",err);
        c->platform = Platforms[0];
    } else {return -1;}

//Getting the context layer
  cps[0] = CL_CONTEXT_PLATFORM;
  cps[1] = (cl_context_properties) c->platform;
  cps[2] = 0;
  cprops = (NULL == c->platform) ? NULL : cps;

  err = clGetDeviceIDs(c->platform, CL_DEVICE_TYPE_GPU, 0, NULL, &numDevices);
  printf("1Err: %d\n",err);
  Devices = (cl_device_id*) malloc(sizeof(cl_device_id) * numDevices);
  err = clGetDeviceIDs(c->platform, CL_DEVICE_TYPE_GPU, numDevices, Devices, NULL);
  printf("2Err: %d\n",err);

//Getting the first available device
  for (unsigned int i = 0; i < numDevices; ++i)
      {
        cl_bool available;
        err = clGetDeviceInfo(Devices[i], CL_DEVICE_AVAILABLE, sizeof(cl_bool), &available, NULL);
        printf("3Err: %d\n",err);
        if (available)
            {
                AvailableDevice = Devices[i];
                err = clGetDeviceInfo(Devices[i], CL_DEVICE_NAME, sizeof(DeviceName), DeviceName, NULL);
                c->device = AvailableDevice;
                break;
            }
            else
            {
                err = clGetDeviceInfo(Devices[i], CL_DEVICE_NAME, sizeof(DeviceName), DeviceName, NULL);
                if (err == CL_SUCCESS)
                    printf("Device %s available for compute.", DeviceName);
                else
                    printf("Device #%d not available for compute.", i);
            }
        if (AvailableDevice == NULL)
        {
          return(-1);
        }
      }
  c->context = (cl_context) clCreateContext(cprops, 1, &AvailableDevice, NULL, NULL, &err);

//Creating the kernel and command queues

  c->kernel_queue = clCreateCommandQueue(c->context, c->device, 0, &err);
  printf("4Err: %d\n",err);
  c->io_queue     = clCreateCommandQueue(c->context, c->device, 0, &err);
  printf("5Err: %d\n",err);

  return(err);
}

cl_int gpuinsp_InitGPU() (cl_context_t* c, unsigned device_id)
{
  int err;
  cl_uint 		 numPlatforms;
  cl_platform_id*	 Platforms            = new cl_platform_id[10];
  cl_context_properties  cps[3]               = {CL_CONTEXT_PLATFORM, (cl_context_properties) gcl_platform, 0};
  cl_context_properties* cprops;
  cl_uint                numDevices;
  cl_device_id*          Devices;
  char                   DeviceName[200];
  cl_device_id           AvailableDevice      = NULL;

//Getting the platforms
  err = clGetPlatformIDs(0, NULL, numPlatforms);
  if (0 < numPlatforms)
    {
        err = clGetPlatformIDs(numPlatforms, Platforms, NULL);
        c->platform = Platforms[0];
    } else {return -1;}

//Getting the context layer
  cps[0] = CL_CONTEXT_PLATFORM;
  cps[1] = (cl_context_properties) c->platform;
  cps[2] = 0;
  cprops = (NULL == c->platform) ? NULL : cps;

  err = clGetDeviceIDs(c->platform, CL_DEVICE_TYPE_GPU, 0, NULL, &numDevices);
  Devices = (cl_device_id*) malloc(sizeof(cl_device_id) * numDevices);
  err = clGetDeviceIDs(gcl_platform, CL_DEVICE_TYPE_GPU, numDevices, Devices, NULL);
  for (unsigned int i = 0; i < numDevices; ++i)
      {
        cl_bool available;
        err = clGetDeviceInfo(Devices[i], CL_DEVICE_AVAILABLE, sizeof(cl_bool), &available, NULL);
        if (available)
            {
                AvailableDevice = Devices[i];
                err = clGetDeviceInfo(Devices[i], CL_DEVICE_NAME, sizeof(DeviceName), DeviceName, NULL);
                break;
            }
            else
            {
                gcl_err = clGetDeviceInfo(Devices[i], CL_DEVICE_NAME, sizeof(DeviceName), DeviceName, NULL);
                if (gcl_err == CL_SUCCESS)
                    printf("Device %s not available for compute.", DeviceName);
                else
                    printf("Device #%d not available for compute.", DeviceName);
            }
        }
        if (AvailableDevice == NULL)
        {
          return(-1);
        }
        c->device = AvailableDevice;
      }

    c->context = clCreateContext(cprops, 1, &AvailableDevice, NULL, NULL, &err);

}

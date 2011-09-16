#define complexMul(a,b) ((float2)(mad(-(a).y, (b).y, (a).x * (b).x), mad((a).y, (b).x, (a).x * (b).y)))
#define GROUPSIZE 256

__kernel void gpuSnrProduct(unsigned long int buffer_length,
                            __global float* stilde_real, __global float* stilde_imag,
                            __global float* htilde_real, __global float* htilde_imag,
                            __global float* snr)
{
    unsigned int gid   = get_global_id(0);
    unsigned int gsize = get_global_size(0);

    #pragma unroll 2

    for(unsigned int i = gid; i < buffer_length; i+=2*gsize)
    {
     float tr = htilde_real[i];
     float ti = htilde_imag[i];
     float sr = stilde_real[i];
     float si = stilde_imag[i];
     snr[i] = sqrt(pow(sr*tr+si*ti,2)+pow(si*tr-sr*ti,2));
    }
}

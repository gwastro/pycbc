#define complexMul(a,b) ((float2)(mad(-(a).y, (b).y, (a).x * (b).x), mad((a).y, (b).x, (a).x * (b).y)))
#define GROUPSIZE 256

__kernel void gpuSnrProduct(unsigned int buffer_length, unsigned int templ_offset,
                            __global float* templ_real, __global float* templ_imag,
                            __global float* noise_real, __global float* noise_imag,
                            __global float* psd, __global float* out_real,
                            __global float* out_imag)
{
    unsigned int gid = get_global_id(0);
    unsigned int gsize = get_global_size(0);
    float2 cTempl  = (float2)(templ_real[templ_offset + gid], templ_imag[templ_offset + gid]);
    cTempl.y = -cTempl.y;
    float psd_val  = psd[gid];
    float psd_val_sq = pown(psd_val, 2);

    #pragma unroll 2
    for(unsigned int i = gid; i < buffer_length; i+=2*gsize)
    {
        float2 cNoise  = (float2)(noise_real[i], noise_imag[i]);
        float2 cMul = complexMul(cTempl, cNoise);
        out_real[i] = (cMul.x * psd_val) / psd_val_sq;
        out_imag[i] = (cMul.y * psd_val) / psd_val_sq;
        out_real[i+gsize] = 0.0;
        out_imag[i+gsize] = 0.0;
    }
}


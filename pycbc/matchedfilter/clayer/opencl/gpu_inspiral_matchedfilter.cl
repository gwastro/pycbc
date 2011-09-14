#define complexMul(a,b) ((float2)(mad(-(a).y, (b).y, (a).x * (b).x), mad((a).y, (b).x, (a).x * (b).y)))
#define GROUPSIZE 256

#include "../include/SI.h"

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

__kernel void gpuSnrNormalize(unsigned int buffer_length, unsigned int filter_offset, 
                              __global float* snr_real, __global float* snr_imag, 
                              __constant float* sqrt_of_variances, 
                              unsigned int templates_are_normalized,
                              unsigned int template_index, 
                              float normalizing_factor, __global float* filtered)
{
    unsigned int gid = get_global_id(0);
    unsigned int gsize = get_global_size(0);
    float sqrt_of_variance;
    if(!templates_are_normalized)
        sqrt_of_variance = sqrt_of_variances[template_index];
    else
        sqrt_of_variance = 1.0;

    #pragma unroll 2
    for(unsigned int i = gid; i < buffer_length; i += 2*gsize)
    {
     	float2 snr = (float2)(snr_real[i], snr_imag[i]);
        filtered[filter_offset + i] = sqrt(pown(normalizing_factor*snr.x, 2) + pown(normalizing_factor*snr.y, 2)) / sqrt_of_variance;
    }
}


__kernel void gpuTemplateProduct(unsigned int buffer_length, __global float* templ_real, __global float* templ_imag,
                     __global float* psd, __global float* out_real)
{
    unsigned int gid = get_global_id(0);
    unsigned int gsize = get_global_size(0);
    float psd_val = psd[gid];

    #pragma unroll 2
    for(unsigned int i = gid; i < buffer_length; i+=2*gsize)
    {
     	float2 templ_val  = (float2)(templ_real[i], templ_imag[i]);
        out_real[i] = (templ_val.x*templ_val.x + templ_val.y*templ_val.y) / psd_val;
    }
}

__kernel void gpuTemplateProductChi2() 
{
}

__kernel void gpuInitBufferReal(__global float* buffer, unsigned int buffer_length, unsigned int offset)
{
    unsigned int gid = get_global_id(0);
    unsigned int gsize = get_global_size(0);

    #pragma unroll 2
    for(unsigned int i = gid + offset; i < buffer_length; i+=gsize)
    {
     	buffer[i] = 0.0;
    }
}

__kernel void gpuInitPeakNumBuffers(__global int* num_of_peaks, __global int* num_of_local_maxes)
{
    unsigned int gid = get_global_id(0);
    num_of_peaks[gid] = 0;
    num_of_local_maxes[gid] = 0;
}


__kernel void gpuSumReduceVariance(__global float* in_real, unsigned int vector_length, 
                                   unsigned int segment_length, unsigned int segment_valid_length, 
                                   unsigned int groups_per_slice, __global float* partResult_real)
{
    __local float shared_real[GROUPSIZE];

    unsigned int tid = get_local_id(0);
    unsigned int gid = get_global_id(0);
    unsigned int slice_id = (2*gid) / segment_valid_length;
    unsigned int group_id = get_group_id(0);
    unsigned int i = slice_id*segment_length + (group_id % groups_per_slice)*GROUPSIZE*2 + tid;

    shared_real[tid] = in_real[i];
    if (i + GROUPSIZE < vector_length)
    {
        shared_real[tid] += in_real[i+GROUPSIZE];
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    // do reduction in shared mem
    #if (GROUPSIZE >= 512)
        if (tid < 256)
        {
            shared_real[tid] += shared_real[tid + 256];
        }
	barrier(CLK_LOCAL_MEM_FENCE);
    #endif
    #if (GROUPSIZE >= 256)
        if (tid < 128)
        {
            shared_real[tid] += shared_real[tid + 128];
        }
	barrier(CLK_LOCAL_MEM_FENCE);
    #endif
    #if (GROUPSIZE >= 128)
        if (tid <  64)
        {
            shared_real[tid] += shared_real[tid +  64];
        }
	barrier(CLK_LOCAL_MEM_FENCE);
    #endif
    #if (GROUPSIZE >=  64)
        if (tid < 32)
        {
            shared_real[tid] += shared_real[tid + 32];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    #endif
    #if (GROUPSIZE >=  32)
        if(tid < 16)
        {
            shared_real[tid] += shared_real[tid + 16];
        }
    #endif
    #if (GROUPSIZE >=  16)
        if(tid < 8)
        {
            shared_real[tid] += shared_real[tid +  8];
        }
    #endif
    #if (GROUPSIZE >=   8)
        if(tid < 4)
        {
            shared_real[tid] += shared_real[tid +  4];
	}
    #endif
    #if (GROUPSIZE >=   4)
        if(tid < 2)
        {
            shared_real[tid] += shared_real[tid +  2];
        }
    #endif
    #if (GROUPSIZE >=   2)
        if(tid < 1)
        {
            shared_real[tid] += shared_real[tid +  1];
        }
    #endif
    // write result for this block to global mem
    if (tid == 0)
    {
        partResult_real[group_id] = shared_real[0];
    }
}


/*
 * Used for the final step of the reduction phase.
 * The result will be written to the lower bound of the `result_real'
 * and `result_imag' arrays.
*/
__kernel void gpuSumReduceVarianceFinal(__global float* in_real, unsigned int slice_valid_length, float normalizing_factor,
                             __global float* result_real, __local float* shared_real)
{
    unsigned int tid = get_local_id(0);
    unsigned int localSize = get_local_size(0);
    unsigned int group_id = get_group_id(0);
    unsigned int i = group_id*slice_valid_length + tid;

    shared_real[tid] = in_real[i];
    shared_real[tid] += in_real[i+localSize];
    barrier(CLK_LOCAL_MEM_FENCE);

    if(localSize >= 512 && tid < 256)
    {
        shared_real[tid] += shared_real[tid + 256];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    if(localSize >= 256 && tid < 128)
    {
        shared_real[tid] += shared_real[tid + 128];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    if(localSize >= 128 && tid < 64)
    {
        shared_real[tid] += shared_real[tid + 64];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
   if(localSize >= 64 && tid < 32)
    {
        shared_real[tid] += shared_real[tid + 32];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    if(localSize >= 32 && tid < 16)
    {
       	shared_real[tid] += shared_real[tid + 16];
    }
    if(localSize >= 16 && tid < 8)
    {
       	shared_real[tid] += shared_real[tid + 8];
    }
    if(localSize >= 8 && tid < 4)
    {
        shared_real[tid] += shared_real[tid + 4];
    }
    if(localSize >= 4 && tid < 2)
    {
     	shared_real[tid] += shared_real[tid + 2];
    }
    if(localSize >= 2 && tid < 1)
    {
     	shared_real[tid] += shared_real[tid + 1];
    }

    // write back final result to the first element of the result array
    if (tid == 0)
    {
        result_real[group_id] = sqrt(normalizing_factor*shared_real[0]);
    }
}

__kernel void gpuSumReduceChi2() 
{
}

__kernel void gpuSumReduceChi2Final() 
{
}


// This kernel is broken for the moment

__kernel void gpuPeakSearch(unsigned int buffer_length, __constant float4* filtered, 
                            float threshold, volatile __global int* num_of_peaks, 
                            __global int* peak_starts, volatile __global int* num_of_local_maxes, 
                           __global int* max_indexes, __global float* max_values)
{
}
/*

    unsigned int gid = get_global_id(0);
    unsigned int gsize = get_global_size(0);
    unsigned int tid = get_local_id(0);

    #pragma unroll 2
    for(unsigned int i = gid; i < buffer_length / 4; i += gsize)
    {
        __local float2 f_values[GROUPSIZE];
        int peak_id;
        int max_id;
        float4 tval = filtered[i];
        float wg_prec_val;

        #pragma unroll 
        for(unsigned int j = 0; j < 3; j++)
        {
         if((tval[j] <= threshold) && (tval[j+1] > threshold))
     	    {
               	peak_id = atomic_inc(num_of_peaks);
                peak_starts[peak_id] = 4*i+j+1;
        
            }
        }
   }
}
        #pragma unroll
       	for(unsigned int j = 1; j <= 2; j++)
        {
            if(tval[j] > threshold && tval[j] > tval[j-1] && tval[j] > tval[j+1])
            {
                max_id = atomic_inc(num_of_local_maxes);
               	max_values[max_id] = tval[j];
                max_indexes[max_id] = 4*i+j;
            }
        }

        f_values[tid][0] = tval[0];
        f_values[tid][1] = tval[3];

        if(i > 0 && tid == 0)
            wg_prec_val = filtered[i-1][3];

        barrier(CLK_LOCAL_MEM_FENCE);

        if( tval[0] > threshold && ( tid > 0 && f_values[tid-1][1] < threshold || i == 0 || tid == 0 && wg_prec_val <= threshold ))
        {
            peak_id = atomic_inc(num_of_peaks);
            peak_starts[peak_id] = 4*i;
        }

        if(tval[0] > threshold && tval[1] < tval[0] && ((tid == 0 && wg_prec_val < tval[0]) || (tid > 0 && f_values[tid-1][1] < tval[0])))
        {
            max_id = atomic_inc(num_of_local_maxes);
            max_vales[max_id] = tval[0];
            max_indexes[max_id] = 4*i;
        }
        if(tval[3] > threshold && tval[2] < tval[3] && (tid == GROUPSIZE - 1 || f_values[tid+1][0] < tval[3]))
        {
            max_id = atomic_inc(num_of_local_maxes);
     	    max_vales[max_id] = tval[3];
            max_indexes[max_id] = 4*i+3;
        }
    }
}

*/
__kernel void gpuGenerateTemplates(int n, int batch_number, float flow, __global float * masses,
                                   int templ_length, __global float * templ_real,__global float * templ_imag)

/*
   This routine calculates one frequency bin of one template.
   When called global_size must be equal signal_batch_size*templ_length/2,
   in order to work properly.

   n            - number of templates to generate
   batch_number - serial number of batch
   flow         - lower cutoff frequency
   masses       - array of mass parameters
   templ_length - length of the template
   templ_real   - real part of the template buffer
   templ_imag   - imag part of the template buffer

*/

{
 int gid  = get_global_id(0);

 int templ_num = (int)(gid*2/templ_length);
 int frek_num  = (gid-templ_num*templ_length)/2;

 float m1    = masses[(batch_number*n+templ_num)*2];
 float m2    = masses[(batch_number*n+templ_num)*2+1];

 float Mtot  = m1 + m2;
 float mu    = m1*m2/Mtot;
 float eta   = mu/Mtot;

 float SI_c3 = SI_c*SI_c*SI_c;

 float df    = 1./templ_length;
 float fisco = SI_c3/(pow(6,1.5)*PI*SI_G*Mtot);

 float c1    = (15293365.0/508032.0+27145./504.*eta+3085/72*eta*eta);
 float c2    = 16*PI;
 float c3    = (3715./756.+55./9.*eta);

 float phase = 0;
float dist  = 1;

 float f     = frek_num * df;

 if (( f >= flow ) && ( f <= fisco )) {

 float amp   = pow(5*PI/24,0.5)*(SI_G*MSun/SI_c/SI_c/SI_pc)*
              pow(PI*SI_G*MSun/SI_c3,-1./6.)*pow(mu/MSun,0.5)*
              pow(Mtot/MSun,0.5)*pow(f,-7./6.);
 float v     = pow((SI_G*Mtot/SI_c3*PI*f),1./3.);
 float psi   = phase + 3/(128*eta)*(pow(v,-5.)+c3*pow(v,-3.)+c2*pow(v,-2)+c1*pow(v,-1));

 float realhk = cos(psi)*amp*dist;
 float imaghk = sin(psi)*amp*dist;

 templ_real[templ_num * templ_length + frek_num] = realhk;
 templ_imag[templ_num * templ_length + frek_num] = imaghk;

 templ_real[(templ_num + 1) * templ_length - frek_num] = realhk;
 templ_imag[(templ_num + 1) * templ_length - frek_num] = imaghk;

 }
}

__kernel void gpuHammingWeight(int N, int Nidx, __global float * source, __global float * dest)
/*
   This kernel weights the data with the hamming window values.
   Other, better window functions could be choosen, as well, however
   for performance tests this does the job.

   N      = length of the segment
   Nidx   = number of segment
   source = source data
   dest   = weighted data
*/

{
 int   gid  = get_global_id(0);
 int   indx = Nidx * N + gid;
 float  hamm = 0.54 - 0.46*cos(2*PI*gid/(N-1));
 dest[gid] = source[indx]*hamm;
}

__kernel void gpuAvrgPSD(__global float * source_real,
                         __global float * source_imag, __global float * target_psd)

/*
 This is the kernel averages the partial PSDs
*/
{
 int gid  = get_global_id(0);
 target_psd[gid] = pow(source_real[gid],2) + pow(source_imag[gid],2);
}



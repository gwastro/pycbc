#include <stdio.h>
#include <datavectorcpu.h>
#include <except.h> /* For PyCBC exceptions */
#include <string.h> /* For memcpy */
#include "straindatacpu.h"
#include "straindatacpu_private.h"

/*
  Temporary functions to copy the contents of a subvector of one datavector to
  another.  To be eventually replaced by slice arguments to datavectors.
 */

int copy_subvector_complex_double_cpu(complex_vector_double_cpu_t *src,
				      complex_vector_double_cpu_t *dst,
				      unsigned long offset,
				      unsigned long length)
{
  if (dst->meta_data.vector_length < length){
    pycbc_throw_exception(PYCBC_INDEX_ERROR,"Destination datavector too small to hold source subvector\n");
    return PYCBC_INDEX_ERROR;
  }
  if (src->meta_data.vector_length < offset + length){
    pycbc_throw_exception(PYCBC_INDEX_ERROR,"Source subvector ends after source vector\n");
    return PYCBC_INDEX_ERROR;
  }

  memcpy((void *) dst->data, (void *) &(src->data[offset]), length*sizeof(*(src->data)));
  return PYCBC_NO_ERROR;
}

int copy_subvector_complex_single_cpu(complex_vector_single_cpu_t *src,
				      complex_vector_single_cpu_t *dst,
				      unsigned long offset,
				      unsigned long length)
{
  if (dst->meta_data.vector_length < length){
    pycbc_throw_exception(PYCBC_INDEX_ERROR,"Destination datavector too small to hold source subvector\n");
    return PYCBC_INDEX_ERROR;
  }
  if (src->meta_data.vector_length < offset + length){
    pycbc_throw_exception(PYCBC_INDEX_ERROR,"Source subvector ends after source vector\n");
    return PYCBC_INDEX_ERROR;
  }

  memcpy((void *) dst->data, (void *) &(src->data[offset]), length*sizeof(*(src->data)));
  return PYCBC_NO_ERROR;
}

int copy_subvector_real_double_cpu(real_vector_double_cpu_t *src,
				   real_vector_double_cpu_t *dst,
				   unsigned long offset,
				   unsigned long length)
{
  if (dst->meta_data.vector_length < length){
    pycbc_throw_exception(PYCBC_INDEX_ERROR,"Destination datavector too small to hold source subvector\n");
    return PYCBC_INDEX_ERROR;
  }
  if (src->meta_data.vector_length < offset + length){
    pycbc_throw_exception(PYCBC_INDEX_ERROR,"Source subvector ends after source vector\n");
    return PYCBC_INDEX_ERROR;
  }

  memcpy((void *) dst->data, (void *) &(src->data[offset]), length*sizeof(*(src->data)));
  return PYCBC_NO_ERROR;
}

int copy_subvector_real_single_cpu(real_vector_single_cpu_t *src,
				   real_vector_single_cpu_t *dst,
				   unsigned long offset,
				   unsigned long length)
{
  if (dst->meta_data.vector_length < length){
    pycbc_throw_exception(PYCBC_INDEX_ERROR,"Destination datavector too small to hold source subvector\n");
    return PYCBC_INDEX_ERROR;
  }
  if (src->meta_data.vector_length < offset + length){
    pycbc_throw_exception(PYCBC_INDEX_ERROR,"Source subvector ends after source vector\n");
    return PYCBC_INDEX_ERROR;
  }

  memcpy((void *) dst->data, (void *) &(src->data[offset]), length*sizeof(*(src->data)));
  return PYCBC_NO_ERROR;
}

//
// Prototype implementation of the ffts for segmenting straindata
//

void* fftw_generate_plan(unsigned long length, real_vector_single_cpu_t* in_tmp,
                    complex_vector_single_cpu_t* out_tmp, char* sign, char* style)
{
    void* plan;

    // for testing 3 * 4 byte buffer as "plan" prototype object
    plan = calloc( 3 , 4 );

    //printf("fftw_generate_plan: length= %ld, in_tmp = %p, out_tmp = %p, sign= %s, style= %s ==> plan= %p\n",
    //       length, in_tmp, out_tmp, sign, style, plan);

    return plan;

}

int fftw_transform_segments(void* plan, real_vector_single_cpu_t* in_buf,
                            unsigned long input_buf_offset,
                            complex_vector_single_cpu_t* out_buf)
{

    //printf("fftw_transform_segments: plan= %p, in_buf + offset = %p, out_buf = %p\n",
    //       plan, in_buf->data + input_buf_offset, out_buf->data);

    // free(plan); // just a test makes no sense for frequently calls
    // The plan lives in python as:
    // <Swig Object of type 'void *' at 0x1004ebfc0>
    // pointers a correctly wrapped by swig in this  way
    //
    // To ensure that the plan is constructed and first of all destructed
    // correctly by the fftw destructor
    // we should define the plan as a regular swigged object (like datavector)
    // with a new_ and a delete_ function in C
    // so the lifetime is determined by the owner object
    // (FftSegmentsImplementation) in python

    return 0;
}

int frame_cpp_read_frames(real_vector_double_cpu_t* out_buf, char* channel_name,
                          unsigned long gps_start_time, unsigned long gps_end_time,
                          char* cache_url)
{

    //printf("frame_cpp_read_frames: out_buf: %p, channel_name= %s, gps_start_time = %ld, gps_end_time = %ld, cache_url= %s, vector_length: %d\n",
    //       out_buf, channel_name, gps_start_time, gps_end_time, cache_url, out_buf->meta_data.vector_length);

    return 0;
}



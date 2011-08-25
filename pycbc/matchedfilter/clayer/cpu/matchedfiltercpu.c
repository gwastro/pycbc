#include <stdio.h>
#include "../../../datavector/clayer/cpu/datavectorcpu_types.h"
#include "matchedfiltercpu_prototypes.h"
#include "matchedfiltercpu_types.h"

void gen_snr_cpu(cpu_context_t* context,
		 real_vector_single_t* snr,
		 complex_vector_single_t* stilde, 
		 complex_vector_single_t* htilde,
		 complex_vector_single_t* q,
		 complex_vector_single_t* qtilde,
		 /*complex_fft_plan_t* plan,*/
		 double f_min,
		 double sigma_sq)
{    

    static unsigned cnt=0;
    int j;
    const int N = 2 * (stilde->meta_data.vector_length - 1);
    double norm = 4.0 / ((double) N * (double) N * stilde->meta_data.delta_x); 
    norm *= norm;
    
    printf("%d: called gen_snr_cpu with context: %p, snr: %p s: %p h: %p q: %p qtilde %p\n", 
           cnt++, context, snr, stilde, htilde, q, qtilde);
    
    /* perform the correlation */
    correlate_complex_freq_vectors( qtilde, stilde, htilde, f_min );
    /* complex_vec_fft( q, qtilde, plan->theplan ); */
  
    /* normalize the snr */
    for ( j = 0; j < q->meta_data.vector_length; ++j )
        snr->data[j] = (norm/sigma_sq) * 
        (__real__ q->data[j] * __real__ q->data[j] + 
         __imag__ q->data[j] * __imag__ q->data[j]);

    snr->meta_data.delta_x = q->meta_data.delta_x;
  
    return;
}


/* freq domain correlation of two vectors (suitably whitened) */
void correlate_complex_freq_vectors( complex_vector_single_t* out,
				     complex_vector_single_t* x, 
				     complex_vector_single_t* y, 
				     double f_min)
{
    int k;
    const double df = x->meta_data.delta_x;
    const int bin_min = f_min / df > 1 ? f_min / df : 1;
    const int len = x->meta_data.vector_length;

    for ( k = bin_min; k < len-1; k++ )
    {
        __real__ out->data[k] = __real__ x->data[k] * __real__ y->data[k] - 
                                __imag__ x->data[k] * 
                                ( 0.0 - __imag__ y->data[k] );
        __imag__ out->data[k] = __real__ x->data[k] * 
                                (0.0 - __imag__ y->data[k]) + 
                                __imag__ x->data[k] * __real__ y->data[k];
    }

    out->meta_data.delta_x = df;

    return;
}

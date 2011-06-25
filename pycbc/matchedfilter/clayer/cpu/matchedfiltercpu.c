#include <stdio.h>
#include "../../../datavector/clayer/cpu/datavectorcpu_types.h"

int gen_snr_cpu(complex_vector_single_t* stilde, 
                complex_vector_single_t* htilde,
                real_vector_single_t* snr)
{
    
    printf("in gen_snr_cpu. complex_vector_single_t* stilde= %p\n", stilde);
    printf("in gen_snr_cpu. complex_vector_single_t* htilde= %p\n", htilde);
    printf("in gen_snr_cpu. real_vector_single_t* snr= %p\n", snr);
        
    printf("stilde->data = %p \n", stilde->data);
    
    return 0;  // return 0 or error code

}

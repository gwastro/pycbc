#include <stdio.h>
#include "../../../datavector/clayer/cpu/datavectorcpu_types.h"
//#include "matchedfiltercpu_types.h"
//#include "matchedfiltercpu_prototypes.h"


int gen_snr_cpu(real_vector_single_t* stilde, 
                real_vector_single_t* htilde,
                real_vector_single_t* snr)
{
    
    printf("in gen_snr_cpu. real_vector_single_t* stilde= %p\n", stilde);
    printf("in gen_snr_cpu. real_vector_single_t* htilde= %p\n", htilde);
    printf("in gen_snr_cpu. real_vector_single_t* snr= %p\n", snr);
    
    printf("stilde->meta_data.t_start = %ld \n", stilde->meta_data.t_start);
    
    printf("stilde->data = %p \n", stilde->data);
    
    return 1;  // return true or error code

}

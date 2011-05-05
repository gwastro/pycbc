#include <stdio.h>
#include "../../../datavector/clayer/cpu/datavectorcpu_types.h"
//#include "matchedfiltercpu_types.h"
//#include "matchedfiltercpu_prototypes.h"


real_vector_t* gen_snr_cpu(real_vector_t* stilde, real_vector_t* htilde)
{
    
    printf("in gen_snr_cpu. real_vector_t* stilde= %p\n", stilde);
    printf("in gen_snr_cpu. real_vector_t* htilde= %p\n", htilde);
    
    printf("stilde->t_start = %d \n", stilde->t_start);
    
    printf("stilde->data = %p \n", stilde->data);
    
    return stilde;

}

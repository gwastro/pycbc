#include <stdio.h>
#include "../../../datavector/clayer/cpu/datavectorcpu_types.h"
#include "../../../clayer/cpu/pycbccpu_types.h"


void gen_snr_cpu(cpu_context_t* context,
                real_vector_single_t* snr,
                complex_vector_single_t* stilde, 
                complex_vector_single_t* htilde)
{    
    static unsigned cnt=0;
    
    printf("%d: called gen_snr_cpu with context: %p, snr: %p s: %p h: %p\n", 
           cnt++, context, snr, stilde, htilde);
    
    
    return;

}

#include <stdio.h>
#include <memory.h>
#include <datavectorcpu.h>
#include <math.h>
#include <complex.h>
//#include <lal/LALConstants.h>
#define LAL_PI        3.1415926535897932384626433832795029  /**< pi */
#define LAL_GAMMA     0.5772156649015328606065120900824024  /**< gamma */
#define LAL_MTSUN_SI  4.92549095e-6   /**< Geometrized solar mass, s */

/*
 * internal structure and prototype
 */

typedef struct {
    double phN;
    double ph2;
    double ph3;
    double ph4;
    double ph5;
    double ph5l;
    double ph6;
    double ph7;
    double ph8;
} pycbc_templatebank_pnconstants;

static void sp_phase_engine(
    complex_vector_single_cpu_t* exp_psi,
    const double M,
    const pycbc_templatebank_pnconstants pn_consts,
    double f_min,
    double f_max,
    real_vector_single_cpu_t* minus_one_by_three
    );

/* 
 * functions
 */


void new_kfac_vec(
    real_vector_single_cpu_t* vec,
    const double kfac
    )
{
    int i;

    for (i = 0; i < vec->meta_data.vector_length; i++)
    {
        vec->data[i] = pow(i * vec->meta_data.delta_x, kfac);
    }

    return;
}


void precondition_factor(
    real_vector_single_cpu_t* vec
    )
{
    new_kfac_vec(vec, -7./6.);
    return;
}


void compute_template_phasing(
    complex_vector_single_cpu_t* exp_psi,
    double M,
    double eta,
    int order,
    double f_min,
    double f_max,
    real_vector_single_cpu_t* minus_one_by_three
    )
{
    pycbc_templatebank_pnconstants pn_consts;
    memset(&pn_consts, 0, sizeof(pycbc_templatebank_pnconstants));

    switch (order)
    {
        case 8:
            /* pseudo 4.0 pN term */
            pn_consts.ph8 = 3923.0;
        case 7:
            /* 3.5 pN term */
            pn_consts.ph7 =
                LAL_PI*(77096675.0/254016.0 + 378515.0*eta/1512.0 
                - 74045.0*eta*eta/756.0);
        case 6:
            /* 3.0 pN term */
            /* According to 
               Blanchet et. al., Phys. Rev. D., 65, 061501
               the sign of the 47324.0/63.0 term should be +
               This also puts this code in agreement with LAL
            */
            pn_consts.ph6 = (11583231236531.0/4694215680.0 
                - 640.0*LAL_PI*LAL_PI/3.0 - 6848*LAL_GAMMA/21.0)
                + eta*(-15335597827.0/3048192.0 
                    + 2255.0*LAL_PI*LAL_PI/12.0 + 47324.0/63.0 - 7948.0/9.0)
                + eta*eta*76055.0/1728.0
                - eta*eta*eta*127825.0/1296.0;
        case 5:
            /* contribution to 2.5 pN term which is multiplied by (1+3ln(v/v0)) */
            pn_consts.ph5 = LAL_PI*(38645.0/756.0 - 65.0*eta/9.0);
            /* contribution to 2.5 pN term which is multiplied by ln(v/v0)      */
            pn_consts.ph5l = 0;
        case 4:
            /* 2.0 pN term */
            pn_consts.ph4 = 15293365.0/508032.0 
                + eta*(27145.0/504.0 + eta*3085.0/72.0);
        case 3:
            /* 1.5 pN term */
            pn_consts.ph3 = -4.0*(4.0*LAL_PI);
        case 2:
            /* 1.0 pN term */
            pn_consts.ph2 = 3715.0/756.0 + eta*55.0/9.0;
        case 1:
        case 0:
            /* 0.0 pN term */
            pn_consts.phN = 3.0/(eta*128.0);
            break;
        default:
            abort();
    }

    sp_phase_engine(exp_psi, M, pn_consts, f_min, f_max, minus_one_by_three);
    return;
}

static void sp_phase_engine(
    complex_vector_single_cpu_t* exp_psi,
    const double M,
    const pycbc_templatebank_pnconstants pn_consts,
    double f_min,
    double f_max,
    real_vector_single_cpu_t* minus_one_by_three
    )
{
    int k, k_min, k_max;
    double x, x1, psi, psi0, psi1, psi2;

    /* length of frequency vector */
    const int N = exp_psi->meta_data.vector_length;

    /* frequency interval */
    const double df = exp_psi->meta_data.delta_x;

    /* initial velocity */
    const double v0 = cbrt( LAL_PI * M * LAL_MTSUN_SI * f_min );

    /* post-Newtonian constants */
    const double c0  = pn_consts.phN;
    const double c10 = pn_consts.ph2;
    const double c15 = pn_consts.ph3;
    const double c20 = pn_consts.ph4;
    const double c25_a = pn_consts.ph5;
    const double c25_b = pn_consts.ph5l;
    const double c30 = pn_consts.ph6;
    const double c35 = pn_consts.ph7;
    const double c40 = pn_consts.ph8;


    /* chebychev coefficents for expansion of sin and cos */
    const double s2 = -0.16605;
    const double s4 =  0.00761;
    const double c2 = -0.49670;
    const double c4 =  0.03705;

    memset( exp_psi->data, 0, exp_psi->meta_data.vector_length * sizeof(complex float) );

    /* compute cutoff indices */
    k_min = f_min / df > 1 ? f_min / df : 1;
    k_max = f_max / df < N ? f_max / df : N;

    x1 = 1. / cbrt(LAL_PI * M * LAL_MTSUN_SI * df);
    x = x1 * minus_one_by_three->data[k_min];
    psi = c0 * x * ( c20 + x * ( c15 + x * (c10 + x * x ) ) );
    psi += c0 * (c25_a*(1.0 + 3.0*log( 1.0/(x*v0) )) + c25_b*log( 1.0/(x*v0)) );
    if ( c30 )
        psi += c0 * (c30 - 6848.0*log(4.0/x)/21.0) * (1.0/x);
    psi += c0 * c35 * (1.0/(x*x));
    psi += c0 * c40 * -1.0 * log(x) * (1.0/(x*x*x));
    psi0 = -2 * LAL_PI * ( floor ( 0.5 * psi / LAL_PI ) );

    for ( k = k_min; k < k_max ; ++k )
    {
        x = x1 * minus_one_by_three->data[k];
        psi = c0 * x * ( c20 + x * ( c15 + x * (c10 + x * x ) ) );
        psi += c0 * (c25_a*(1.0 + 3.0*log( 1.0/(x*v0) )) + c25_b*log( 1.0/(x*v0)) );
        if ( c30 )
            psi += c0 * (c30 - 6848.0*log(4.0/x)/21.0) * (1.0/x);
        psi += c0 * c35 * (1.0/(x*x));
        psi += c0 * c40 * -1.0 * log(x) * (1.0/(x*x*x));
        psi1 = psi + psi0;

        /* range reduction of psi1 */
        while ( psi1 < -LAL_PI )
        {
            psi1 += 2 * LAL_PI;
            psi0 += 2 * LAL_PI;
        }
        while ( psi1 > LAL_PI )
        {
            psi1 -= 2 * LAL_PI;
            psi0 -= 2 * LAL_PI;
        }

        /* compute approximate sine and cosine */
        /* FIXME change the g99 __real__ and __imag__ syntax to something more portable */
        if ( psi1 < -LAL_PI/2 )
        {
            psi1 = -LAL_PI - psi1;
            psi2 = psi1 * psi1;
            /* XXX note minus sign on imag part is due to LAL sign convention vs PN literature sign convention */
            __imag__ exp_psi->data[k] = - psi1 * ( 1 + psi2 * ( s2 + psi2 * s4 ) );
            __real__ exp_psi->data[k] = -1 - psi2 * ( c2 + psi2 * c4 );
        }
        else if ( psi1 > LAL_PI/2 )
        {
            psi1 = LAL_PI - psi1;
            psi2 = psi1 * psi1;
            /* XXX note minus sign on imag part is due to LAL sign convention vs PN literature sign convention */
            __real__ exp_psi->data[k] = - psi1 * ( 1 + psi2 * ( s2 + psi2 * s4 ) );
            __imag__ exp_psi->data[k] = -1 - psi2 * ( c2 + psi2 * c4 );
        }
        else
        {
            psi2 = psi1 * psi1;
            /* XXX note minus sign on imag part is due to LAL sign convention vs PN literature sign convention */
            __real__ exp_psi->data[k] = - psi1 * ( 1 + psi2 * ( s2 + psi2 * s4 ) );
            __imag__ exp_psi->data[k] = 1 + psi2 * ( c2 + psi2 * c4 );
        }
    }

    return;
}

TITLE Poisson Shot Noise (Correct Model for High Variance Conductance)



NEURON {

    THREADSAFE

    POINT_PROCESS ShotNoiseProcess

    RANGE g, E, g0, std, tau, rate, q

    NONSPECIFIC_CURRENT i

    

    BBCOREPOINTER rng

    

    : Misc

    RANGE synapseID

}



UNITS {

    (nA) = (nanoamp)

    (mV) = (millivolt)

    (uS) = (microsiemens)

}



PARAMETER {

    E     = 0      (mV)

    g0    = 0.0121 (uS)

    std   = 0.0030 (uS)

    tau   = 2.728  (ms)



    : Misc

    synapseID = 0

}



VERBATIM

#include <stdlib.h>

#include <stdio.h>

#include <math.h>

#include <string.h>



#ifndef NRN_VERSION_GTEQ_8_2_0

#include "nrnran123.h"



#ifndef CORENEURON_BUILD

extern int ifarg(int iarg);

extern void* vector_arg(int iarg);

extern double* vector_vec(void* vv);

extern int vector_capacity(void* vv);

#endif



extern double nrn_random_pick(void* r);

extern void* nrn_random_arg(int argpos);

#define RANDCAST

#else

#define RANDCAST (Rand*)

#endif

ENDVERBATIM



STATE {

    g (uS)        : MUST be a state variable if it has g'

}



ASSIGNED {

    v       (mV)

    i       (nA)

    rate    (1/ms)

    q       (uS)

    usingR123

    rng

}



INITIAL {

    g = 0



    if (g0 > 0 && tau > 0) {

        : q = 2 * variance / mean

        q = 2.0 * (std*std) / g0



        if (q > 0) {

            : rate = mean / (q * tau)

            rate = g0 / (q * tau)

        } else {

            rate = 0

        }

    } else {

        q = 0

        rate = 0

    }



    VERBATIM

    if( usingR123 ) {

        nrnran123_setseq((nrnran123_State*)_p_rng, 0, 0);

    }

    ENDVERBATIM



    if (rate > 0) {

        LOCAL delay

        delay = invl(rate)

        : Force positive and ensure minimum delay for CoreNEURON compatibility

        delay = sqrt(delay * delay)

        if (delay < 0.1) { delay = 0.1 }

        net_send(delay, 1)

    }

}



BREAKPOINT {

    SOLVE state METHOD cnexp

    i = g * (v - E)

}



DERIVATIVE state {

    g' = -g / tau

}



NET_RECEIVE (w) {

    LOCAL delay

    g = g + q

    delay = invl(rate)

    : Force positive and ensure minimum delay for CoreNEURON compatibility

    delay = sqrt(delay * delay)

    if (delay < 0.1) { delay = 0.1 }

    net_send(delay, 1)

}





FUNCTION invl(mean_rate (1/ms)) (ms) {

    VERBATIM

    double r = 0.5;

    double result = 1.0;

    double mrate = _lmean_rate;

    

    // Get random number if RNG pointer is available

    if (_p_rng) {

        #ifdef CORENEURON_BUILD

        // In CoreNEURON, RNG is always R123 type

        r = nrnran123_dblpick((nrnran123_State*)_p_rng);

        #else

        // In NEURON, check which RNG type

        if (usingR123) {

            r = nrnran123_dblpick((nrnran123_State*)_p_rng);

        } else {

            r = nrn_random_pick(RANDCAST _p_rng);

        }

        #endif

    }

    // else keep default r = 0.5

    

    // Safety clamp r to (0, 1)

    if (r != r) r = 0.5; // NaN check

    if (r <= 0.0) r = 1e-9;

    if (r >= 1.0) r = 0.999;

    

    // Calculate interval

    if (mrate <= 1e-9) {

        result = 100.0;

    } else {

        result = -log(r) / mrate;

    }

    

    // Force positive

    if (result < 0.0) result = -result;

    

    // Minimum bound

    if (result < 0.001) result = 0.001;

    

    // Final NaN check

    if (result != result) result = 1.0;

    

    _linvl = result;

    ENDVERBATIM

}









PROCEDURE setRNG() {

VERBATIM

    #ifndef CORENEURON_BUILD

    // For compatibility, allow for either MCellRan4 or Random123

    usingR123 = 0;

    

    if( ifarg(1) && hoc_is_double_arg(1) ) {

        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);

        uint32_t a2 = 0;

        uint32_t a3 = 0;



        if (*pv) {

            nrnran123_deletestream(*pv);

            *pv = (nrnran123_State*)0;

        }

        if (ifarg(2)) {

            a2 = (uint32_t)*getarg(2);

        }

        if (ifarg(3)) {

            a3 = (uint32_t)*getarg(3);

        }

        *pv = nrnran123_newstream3((uint32_t)*getarg(1), a2, a3);

        usingR123 = 1;

    } else if( ifarg(1) ) {  

        // Legacy Random object

        void** pv = (void**)(&_p_rng);

        *pv = nrn_random_arg(1);

    } else {  

        void** pv = (void**)(&_p_rng);

        *pv = (void*)0;

    }

    #endif

ENDVERBATIM

}



PROCEDURE clearRNG() {

VERBATIM

    #ifndef CORENEURON_BUILD

    if (usingR123) {

        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);

        if (*pv) {

            nrnran123_deletestream(*pv);

            *pv = (nrnran123_State*)0;

        }

    } else {

        void** pv = (void**)(&_p_rng);

        if (*pv) {

            *pv = (void*)0;

        }

    }

    #endif

ENDVERBATIM

}



FUNCTION urand() {

    VERBATIM

    double value = 0.5; // Default safe value

    

    if (usingR123 && _p_rng) {

         value = nrnran123_dblpick((nrnran123_State*)_p_rng);

    } else if (_p_rng) {

         #ifndef CORENEURON_BUILD

         value = nrn_random_pick(RANDCAST _p_rng);

         #endif

    }

    // else keep default 0.5

    

    // Check for NaN or infinity

    if (value != value) { value = 0.5; } // NaN check

    if (value < 0.0 || value > 1.0e308) { value = 0.5; } // Catch negatives and inf

    

    // Strict clamp to (0, 1) open interval

    if (value <= 0.0) { value = 1e-9; }

    if (value >= 1.0) { value = 0.999999; }

    

    _lurand = value;

    ENDVERBATIM

}



VERBATIM

static void bbcore_write(double* x, int* d, int* x_offset, int* d_offset, _threadargsproto_) {

  if (d) {

    uint32_t* di = ((uint32_t*)d) + *d_offset;

    

    // Get address of rng pointer

    nrnran123_State** pv = (nrnran123_State**)(&_p_rng);

    

    // Check if pointer is valid before accessing

    // Note: The reference implementation assumes it is valid if we are calling this?

    // Actually ProbAMPANMDA_EMS does not check for null *pv, implying it assumes initialized.

    // But we should be safe.

    if (*pv) {

        nrnran123_getids3(*pv, di, di+1, di+2);

        char which;

        nrnran123_getseq(*pv, di+3, &which);

        di[4] = (int)which;

    } else {

        di[0]=0; di[1]=0; di[2]=0; di[4]=0;

    }

  }

  // Reserve 5 ints

  *d_offset += 5;

}



static void bbcore_read(double* x, int* d, int* x_offset, int* d_offset, _threadargsproto_) {

  uint32_t* di = ((uint32_t*)d) + *d_offset;

  

  // If seeds are valid

  if (di[0] != 0 || di[1] != 0 || di[2] != 0) {

      nrnran123_State** pv = (nrnran123_State**)(&_p_rng);

      

      #if !NRNBBCORE

      if(*pv) { nrnran123_deletestream(*pv); }

      #endif

      

      *pv = nrnran123_newstream3(di[0], di[1], di[2]);

      char which = (char)di[4];

      nrnran123_setseq(*pv, di[3], which);

      // Note: DO NOT set usingR123 here - it's set by setRNG in NEURON

  }



  *d_offset += 5;

}

ENDVERBATIM
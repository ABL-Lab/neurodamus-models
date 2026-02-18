: from https://senselab.med.yale.edu/ModelDB/ShowModel.cshtml?model=168148&file=/stadler2014_layerV/kBK.mod
TITLE large-conductance calcium-activated potassium channel (BK)
	:Mechanism according to Gong et al 2001 and Womack&Khodakakhah 2002,
	:adapted for Layer V cells on the basis of Benhassine&Berger 2005.
	:NB: concentrations in mM

	:Modified by Dhuruva Priyan (gmdhuruva@gmail.com)
	:Added Q10 correction

	:FITTED PARAMETERS:
    :=================================================
    : caPh    = 1.7090e-05 mM  (half-max Ca for P0)
    : caPk    = 0.351      (Hill coeff for Ca-P0)
    : caPmax  = 0.916      (max P0)
    : caPmin  = 0.500      (min P0)
    : caVhh   = 3.0647e-03 mM  (half-max Ca for Vh)
    : caVhk   = -1.098      (Hill coeff for Ca-Vh)
    : caVhmax = 150.00 mV    (max Vh)
    : caVhmin = -38.84 mV    (min Vh)
    : k       = 28.90 mV      (voltage sensitivity)
	
NEURON {
	SUFFIX kBK
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE gpeak, gkact, caPh, caPk, caPmax, caPmin
    RANGE p, pinf, tau
	RANGE caVhh, CaVhk, caVhmax, caVhmin, k, tau
        GLOBAL pinfmin : cutoff - if pinf < pinfmin, set to 0.; by default cutoff not used (pinfmin==0)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) 	= (millimolar)
}



PARAMETER {
		:maximum conductance (Benhassine 05)
	gpeak   = 268e-4	(mho/cm2) <0, 1e9>
	
	                                    : Calcium dependence of opening probability (Gong 2001)
	caPh    = 1.7090e-05     (mM)             : conc. with half maximum open probaility
	caPk    = 0.351                      : Steepness of calcium dependence curve
	caPmax  = 0.916                       : max and
	caPmin  = 0.500                        : min open probability
		
	                                    : Calcium dependence of Vh shift (Womack 2002)
	caVhh   = 3.0647e-03    (mM)              : Conc. for half of the Vh shift
	caVhk   = -1.098                 : Steepness of the Vh-calcium dependence curve
	caVhmax = 150.00 (mV)               : max and
	caVhmin = -38.84 (mV)               : min Vh
	
	                                    : Voltage dependence of open probability (Gong 2001)
	                                    : must not be zero
	k       = 28.90	(mV)
	
	                                    : Timeconstant of channel kinetics
	                                    : no data for a description of a calcium&voltage dependence
	                                    : some points (room temp) in Behassine 05 & Womack 02
	tau     = 1 (ms) <1e-12, 1e9>
	scale   = 1                    : scaling to incorporate higher ca conc near ca channels
        
    pinfmin = 0.0                       : cutoff for pinf - less than that set pinf to 0.0
	tau_act = 1 (ms)
    tau_deact = 3 (ms)
	
} 	


ASSIGNED {
	v 		(mV)
	ek		(mV)
	ik 		(mA/cm2)
    	cai  		(mM)
	caiScaled	(mM)
	pinf		(1)
	qt (1)
}


STATE {
        p
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gpeak*p* (v - ek)
}

DERIVATIVE states {     
        rate(v, cai)
        :p' =  (pinf - p) / (tau)

		if (p < pinf) {
			p' = (pinf - p) / tau_act  : activation
		} else {
			p' = (pinf - p) / tau_deact  : deactivation
		}
}

INITIAL {     
        rate(v, cai)
        p = pinf
}

PROCEDURE rate(v(mV), ca(mM))  {
		LOCAL q10
		q10 = 3.0
		qt = q10^((celsius - 22)/10)

        caiScaled = ca*scale
        pinf = P0ca(caiScaled) / ( 1 + exp( (Vhca(caiScaled)-v)/k ) )
        if(pinf < pinfmin) { pinf = 0.0 }
}

FUNCTION P0ca(ca(mM)) (1) {
		
	if (ca < 1E-18) { 		:check for division by zero		
	P0ca = caPmin
	} else {
	P0ca = caPmin + ( (caPmax - caPmin) / ( 1 + (caPh/ca)^caPk ))
	}
}

FUNCTION Vhca(ca(mM)) (mV) {
		
	if (ca < 1E-18) {		:check for division by zero
	Vhca = caVhmax
	} else {
	Vhca = caVhmin + ( (caVhmax - caVhmin ) / ( 1 + ((caVhh/ca)^caVhk)) )
	}
}	


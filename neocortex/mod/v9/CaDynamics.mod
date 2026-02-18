TITLE decay of internal calcium concentration
:
: Internal calcium concentration calculated from calcium currents
: and buffered by endogenous buffer and extrusion mechanism.
:
: Uses differential equations from Helmchen 1996
:dCa/dt = (dCa_T delta_t - (gamma*(dCa - Ca_rest)))/kb
: or dCa/dt = (dCa_T delta_t)/kb - (dCa - Ca_rest)/taur
: with  taur = kb/gamma
:
: to add exogenous buffer kb = 1+kendo+kexo
: for OGB-1 kexo = concOGB1/kd = 200uM/0.2uM => kb=1020
: for OGB-6 kexo = concOGB6/kd = 200uM/3uM   => kb=80
: for Fluo-4 kexo = concFluo4/kd = 300uM/0.335uM => kb=895

:
: mod file was modified from original version (Destexhe 92)
: use diam/4 instead of depth to calculate [Ca]
: Units checked using "modlunit" -> factor 10000 needed in ca entry
:
: Written by B Kampa May 2006


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cad
	USEION ca READ ica, cai WRITE cai
	GLOBAL cainf
	RANGE cai_OGB1_dye, cai_OGB6_dye, cai_Fluo4_dye
    RANGE tau, tau_OGB1_dye, tau_OGB6_dye, tau_Fluo4_dye
}

UNITS {
	(molar) = (1/liter)			: moles do not appear in units
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)	= (ms mM)
	FARADAY = (faraday) (coulomb)
}


PARAMETER {
	diam		(um)
    
    gamma_0 = 0.24 : Intrinsic Extrusion Rate | Obtained from Cornelisse et al., 2007
    
    kE = 62 : Endogenous buffer capacity of dendrite | Obtained from Cornelisse et al., 2007

	kexo_OGB1 = 1000 : Exogenous buffer capacity for OGB-1 | OGB-1 kexo = concOGB1/kd = 200uM/0.2uM | Methods section of Kampa and Stuart 2006
	kexo_OGB6 = 66.67 : Exogenous buffer capacity for OGB-6 | OGB-6 kexo = concOGB6/kd = 200uM/3uM | Methods section of Kampa and Stuart 2006
	kexo_Fluo4 = 895 : Exogenous buffer capacity for Fluo-4 | Fluo-4 kexo = concFluo4/kd = 300uM/0.345uM | https://www.thermofisher.com/order/catalog/product/F14200

	cainf = 6.5e-5 (mM) : baseline calcium | CA1: (Sabatini 2002) 65 nM
	
	cai		(mM)
	cai_OGB1_dye	(mM)
	cai_OGB6_dye	(mM)
	cai_Fluo4_dye	(mM)

	tau		(ms) : time constant for calcium decay
	tau_OGB1_dye	(ms) : time constant for OGB-1 dye decay
	tau_OGB6_dye	(ms) : time constant for OGB-6 dye decay
	tau_Fluo4_dye	(ms) : time constant for Fluo-4 dye decay
}	

STATE {
	ca		(mM) <1e-5>
	ca_OGB1_dye	(mM) <1e-5>
	ca_OGB6_dye	(mM) <1e-5>
	ca_Fluo4_dye	(mM) <1e-5>
}

INITIAL {
	ca = cainf
	cai = ca
	cai_OGB1_dye = ca_OGB1_dye
	cai_OGB6_dye = ca_OGB6_dye
	cai_Fluo4_dye = ca_Fluo4_dye
}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
}

BREAKPOINT {
	SOLVE state METHOD euler
}

DERIVATIVE state {
	LOCAL SVR, gamma

	SVR = 4/diam
	drive_channel =  - (10000) * ica * SVR / (2 * FARADAY)
    
	if (drive_channel <= 0.) { drive_channel = 0. }	: cannot pump inward
	
    gamma = gamma_0*SVR
    
    tau = (1+kE)/gamma
    ca' = (drive_channel/(1+kE)) + ((cainf-ca)/tau)
	cai = ca
    
    tau_OGB1_dye = (1+kE+kexo_OGB1)/gamma 
	ca_OGB1_dye' = (drive_channel/(1+kE+kexo_OGB1)) + ((cainf-ca_OGB1_dye)/tau_OGB1_dye)
	cai_OGB1_dye = ca_OGB1_dye

	tau_OGB6_dye = (1+kE+kexo_OGB6)/gamma
	ca_OGB6_dye' = (drive_channel/(1+kE+kexo_OGB6)) + ((cainf-ca_OGB6_dye)/tau_OGB6_dye)
	cai_OGB6_dye = ca_OGB6_dye

	tau_Fluo4_dye = (1+kE+kexo_Fluo4)/gamma
	ca_Fluo4_dye' = (drive_channel/(1+kE+kexo_Fluo4)) + ((cainf-ca_Fluo4_dye)/tau_Fluo4_dye)
	cai_Fluo4_dye = ca_Fluo4_dye
}
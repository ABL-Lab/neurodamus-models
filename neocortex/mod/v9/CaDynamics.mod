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
:
: 2026-07-23: dye-kinetics states (OGB-1, OGB-6, Fluo-4, Fluo-5F) removed from
: the hot path -- they were never reported in production runs and added 5 extra
: STATE ODEs per segment across every apical/basal compartment of pyramidal
: cells, causing ~7x runtime regression vs CaDynamics_DC0 (1 state) at
: identical dt/tstop. See git history / DEES_cell_packages for the full
: multi-dye version if imaging-fit work needs it again:
:
: to add exogenous buffer kb = 1+kendo+kexo
: for OGB-1 kexo = concOGB1/kd = 200uM/0.2uM => kb=1020
: for OGB-6 kexo = concOGB6/kd = 200uM/3uM   => kb=80
: for Fluo-4 kexo = concFluo4/kd = 300uM/0.335uM => kb=895
: Fluo-5F (Nevian 2007 / Antic 2009): 200 uM dye, Kd=2300 nM,
:   kon=250 mM-1ms-1, koff=0.575 ms-1, explicit 2-state kon/koff kinetics


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cad
	USEION ca READ ica, cai WRITE cai
	GLOBAL cainf
	RANGE tau
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

	cainf = 6.5e-5 (mM) : baseline calcium | CA1: (Sabatini 2002) 65 nM

	cai		(mM)

	tau		(ms) : time constant for calcium decay
}

STATE {
	ca		(mM) <1e-5>
}

INITIAL {
	ca = cainf
	cai = ca
}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
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
}
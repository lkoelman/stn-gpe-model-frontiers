TITLE  Ion channels for STN cells

:
: Na+, K, Ca_T, Leakage and Ca diffusion
: 

COMMENT
Model described in Rubin, Terman (2002) & (2004)
Parameters are for the 2004 model, 2002 values are tagged with 'RT2002' where different.

Next to the HH currents 
  - I_Na - fast sodium current
  - I_K - delayed rectifier K+ current
  - I_l - leak current
, the GP neuron contains four additional currents:
  - I_T - T-type Ca2+ currents, low threshold
  - I_Ca - Calcium current
  - I_AHP - Ca2+ activated K+ channels (repolarization component of plateu potentials)
  
ENDCOMMENT


NEURON {
	SUFFIX stnRT
	USEION ca READ eca, cai WRITE ica, cai
	USEION k READ ek WRITE ik
	USEION na READ ena WRITE ina
	
	:############ Nonspecific (no ion) ############
	NONSPECIFIC_CURRENT ibias
	NONSPECIFIC_CURRENT il
	RANGE il, ibias, bias_amp
	RANGE gl, el

	:############ Na ion ############
	:Na fast
	RANGE ina, ena0
	RANGE gnabar
	RANGE m_inf, h_inf, tau_h, phi_h

	:############ K ion ############
	RANGE ik, ek0
	:K delayed rectifier
	RANGE ikDR
	RANGE gkdrbar
	RANGE n_inf, tau_n, phi_n
	:K AHP (Ca-dependent)
	RANGE ikAHP
	RANGE gkcabar, k1

	:############ Ca ion ############
	RANGE ica, eca0
	:Regular Ca current
	RANGE icaS
	RANGE gcabar
	RANGE s_inf
	:T-type ca current
	RANGE icaT
	RANGE gcatbar
	RANGE a_inf, r_inf, phi_r, b_inf, tau_r
	: (export T-current parameters for modification)
	RANGE theta_a, sigma_a, theta_b, sigma_b
	RANGE theta_r, sigma_r, phi_r, tau_r0, tau_r1
	: Ca ion dynamics
	RANGE kca, epsilon : cavol, caarea
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S)  = (siemens)
	(nS) = (nanosiemens)
	(molar) = (1/liter)
	(mM)	= (millimolar)
	FARADAY = (faraday) (coulomb)  :units are really coulombs/mole
	PI	= (pi) (1)
}

PARAMETER {
	R = 8.31441 (Gas constant)
	T 			(Absolute temp)
	celsius		(degC)

::: Na ion :::
:Na fast
	ena0 = 55 (mV) : vna in .ode
	gnabar = 37.5 (S/cm2) : gna in .ode
	: steady state values for rate constants
	theta_m = -30 : thetam in .ode (INVERTED)
	sigma_m = 15  : sm in .ode
	theta_h = -39 : thetah in .ode
	sigma_h = 3.1 : sh in .ode
	: time constants for rate constants with first-order kinetics
	tau_h0 = 1 (ms) : tauh0 in .ode
	tau_h1 = 500 (ms) : tauh1 in .ode
	theta_tauh = 57 :thh in .ode
	sigma_tauh = 3 : sigmah in .ode
	phi_h = 0.75 : phi in .ode
	
::: K ion :::
:K delayed rectifier
	ek0 = -80 (mV) : vkg in .ode
	gkdrbar = 45 (S/cm2) : gkg in .ode
	: steady state values for rate constants
	theta_n = -32 : thetan in .ode
	sigma_n = -8 : sn in .ode
	: time constants for rate constants with first-order kinetics
	tau_n0 = 1 (ms) : taun0 in .ode
	tau_n1 = 100 (ms) : taun1 in .ode
	theta_taun = 80 : thn in .ode
	sigma_taun = 26 : sigman in .ode
	phi_n = 0.75 : phi in .ode
:K AHP (Ca-dependent)
	gkcabar = 9 (S/cm2) : gahp in .ode
	k1 = 15.0 : k1 in .ode

::: Ca ion :::
:Regular Ca current
	eca0 = 140 (mV) : vca in .ode
	gcabar = 0.5 (S/cm2) : gca in .ode
	: steady state values for rate constants
	theta_s = -39 : thetas in .ode (INVERTED)
	sigma_s = 8   : ss in .ode
:T-type Ca current
	gcatbar = 0.5 (S/cm2) : gt in .ode
	: steady state values for rate constants
	theta_a = -63 : thetat in .ode
	sigma_a = -7.8 : kt in .ode
	theta_b = 0.25 : rth in .ode (0.4 in RT2002)
	sigma_b = -0.07 : rsig in .ode (-0.1 in RT2002)
	theta_r = -67 : thetar in .ode (-67 in RT2002)
	sigma_r = -2 : kr in .ode (INVERTED) (-2.0 in RT2002)
	: time constants for r
	tau_r0 = 7.1 : taur0 in .ode (40 in RT2002)
	tau_r1 = 17.5 : taur1 in .ode
	theta_taur = -68 : thr in .ode
	sigma_taur = 2.2 : sigmar in .ode
	phi_r = 0.5 : phir in .ode (0.2 in RT2002)
: Ca ion dynamics
	epsilon = 5e-05 (1/ms) : eps in .ode (3.75e-5 in RT2002)
	kca = 22.5 : kca in .ode
	: area and volume are for cylindrical section with diam=60 and L=60 micron
	: caarea = 11309.73: pi*diam*L, provide in (micron^2)
	: cavol = 1.6964600329e-10 (L) : pi*R^2*L*conversion factor (micron^2) to (L)

::: Nonspecific (no ion) :::
: Leak current
	gl = 2.25 (S/cm2) : gl in .ode
	el = -60 (mV) : vl in .ode
	bias_amp = 25 (mA/cm2) : inp in .ode (0 in RT2002)
}

: Variables that are
: 	- declared and assigned by mechanism
:		- (i.e. variables assigned that are not parameters or state variables)
: 	- declared by NEURON and assigned OR used by mechanism (e.g. v, ena, ek, eca)
ASSIGNED {
	v	(mV)
	
: Nonspecific (no ion)/leak
	il	(mA/cm2)
	ibias	(mA/cm2)

::: Na ion :::
	ina	(mA/cm2)
	ena     (mV)   := 60
:Na fast
	m_inf
	h_inf
	tau_m	(ms)
	tau_h	(ms)
	
::: K ion :::
	ik	(mA/cm2) : sum of K currents
	ek      (mV) := -90
:K delayed rectifier
	ikDR	(mA/cm2)  
	n_inf
	tau_n	(ms)
:K AHP (Ca-dependent)
	ikAHP	(mA/cm2)

::: Ca ion :::
	ica		(mA/cm2) : sum of Ca currents
	eca 	(mV)
:Regular Ca current
	icaS 	(mA/cm2)
	s_inf
:T-type ca current
	icaT 	(mA/cm2)
	a_inf
	b_inf
	r_inf
	tau_r
}

STATE {
	n h r : rate constants governed by first order kinetics

	: including ion concentrations in STATE allows us to assign independent initial values
	cai (mM)
}


BREAKPOINT {
	: compute rate constants governed by first-order kinetics
	SOLVE states METHOD cnexp

	: compute steady-state rate constants dependent on other rate constants
	set_ss_rates_rdep(r)

	: compute individual currents
	il = gl * (v - el)
	ibias = -bias_amp
	ina = gnabar * m_inf*m_inf*m_inf*h * (v - ena0)
	ikDR = gkdrbar * n*n*n*n * (v - ek0)
	ikAHP = gkcabar * (v - ek0) * (cai/(cai+k1)) : Ca-dependent K-current
	icaS = gcabar * s_inf*s_inf * (v - eca0)
	icaT = gcatbar * a_inf*a_inf*a_inf*b_inf*b_inf*(v - eca0) : T-type Ca current
	
	: Combine ionic currents
	ica = icaT + icaS
	ik = ikDR + ikAHP
	
}

DERIVATIVE states {
	: compute voltage-dependent steady-state rates and time constants
	set_timeconsts(v)
	set_ss_rates(v)

	: compute rate constants governed by first-order kinetics
	n' = phi_n*(n_inf - n)/tau_n
	h' = phi_h*(h_inf - h)/tau_h
	r' = phi_r*(r_inf - r)/tau_r

	: update intracellular Ca2+ ion concentration (all Ca currents contribute)
	cai' = epsilon*(-icaS - icaT - kca*cai)
	: cai' = epsilon*(-ica*caarea*1e-11/(2*FARADAY*cavol) - kca*cai)
	: cai' = epsilon*(-ica - kca*cai)
    
    : Units for the ica equation are mol/(L*sec)
    : (mM/msec) = (Ica mA/cm2)*(area um2)*(1e-8 cm2/um2)*(1e-3 A/mA)*(1/(2*F) mol/C)*(1e-3 sec/msec)*(1e3 mMol/mol)(1/volume 1/L)
    : - you multiply current with area because currents are actually current densities
    : - Faraday constant has units [C/mol]

	
}

:UNITSOFF

INITIAL {
	: set v-dependent time constants
	set_timeconsts(v)
	: set v-dependent rate constants
	set_ss_rates(v)
	n = n_inf
	h = h_inf
	r = r_inf
	set_ss_rates_rdep(r)

}

FUNCTION sigmshift(v, tau_x0, tau_x1, mid, slope) {
	: shifted sigmoidal function with minimum tau_x0 and max. tau_x0+tau_x1
	sigmshift = tau_x0 + tau_x1/(1.0+exp(-(v-mid)/slope))
}

FUNCTION sigm(v, mid, slope) {
	: sigmoid (logistic function)
	sigm = 1.0/(1.0+exp(-(v-mid)/slope))
}

FUNCTION sigmmir(v, mid, slope) {
	: mirrored sigmoid (logistic function)
	sigmmir = 1.0/(1.0+exp((v-mid)/slope))
}

PROCEDURE set_timeconsts(v(mV)) {
	: set voltage-dependent time constants for rates with first-order kinetics
	tau_n = sigmshift(-v, tau_n0, tau_n1, theta_taun, sigma_taun) : taun in .ode
	tau_h = sigmshift(-v, tau_h0, tau_h1, theta_tauh, sigma_tauh) : tauh in .ode
	tau_r = sigmshift(-v, tau_r0, tau_r1, theta_taur, sigma_taur) : taur in .ode
}

PROCEDURE set_ss_rates(v(mV)) { 
	: set voltage-dependent steady-state rates for rates with first-order kinetics
	n_inf = sigmmir(v, theta_n, sigma_n) : ninf in .ode
	h_inf = sigmmir(v, theta_h, sigma_h) : hinf in .ode
	r_inf = sigm(v, theta_r, sigma_r) : rinf in .ode

	: set voltage-dependent steady-state rates for instantaneous rates
	m_inf = sigm(v, theta_m, sigma_m) : minf in .ode
	a_inf = sigmmir(v, theta_a, sigma_a) : tinf in .ode
	s_inf = sigm(v, theta_s, sigma_s) : sinf in .ode
}

PROCEDURE set_ss_rates_rdep(r) {
	: set steady state rates that are dependent on other rates/gating variables
	: set r-dependent rate constant for T-type Ca current
	b_inf = sigmmir(r, theta_b, sigma_b) - sigmmir(0, theta_b, sigma_b) : rnew in .ode
}

:UNITSON
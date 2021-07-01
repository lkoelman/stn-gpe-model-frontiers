TITLE transient sodium channel

COMMENT
Sodium channel for transient sodium current including voltage-dependent
slow inactivation and recovery from inactivation (Kuo and Bean 1994) 

Modified by Lucas Koelman, UCD Neuromuscular Systems group (November 2016)
	- Technical improvements:
		- Q10 correction of maximum conductance
		- use tables for rates
	- Mechanism:
		- see parameters in analyze_stn_channels.py to implement a specific channel model

Modified by Akemann and Knoepfel, J.Neurosci. 26 (2006) 4602
	- Q10 correction of all rate constants
	- changed kinetic parameters to model the transient Na current with
	  slow inactivation + recovery instead of resurgent current
		a) epsilon = 1e-12 1/ms (from epsilon = 1.75 1/ms in Narsg)
			- prevents transition to 'open channel block' state 'OB'
		b) Oon = 2.3 1/ms (from Oon = 0.75 1/ms in Narsg)
			- rapid inactivation of transient Na current
		c) gbar = 0.008 mho/cm2 (from 0.015 mho/cm2)

Khaliq & Raman (2003), J.Neurosci. 23 (2003) 4899
	- implementation of model of resurgent Na current and publication to ModelDB
	- updated kinetic parameters

Raman & Bean (2001), Biophys.J. 80 (2001) 729 
	- addition of blocking particle to model resurgent Na current

Kuo & Bean (1994)
	- initial formulation of Markov-type gating model for sodium current
	  with slow inactivation and recovery/deinactivation
	  	- see fig. 6

ENDCOMMENT

NEURON {
  SUFFIX NaTB : NOTE: change Natrans->Na and gbar->gna for interoperability with default Na mechanism
  USEION na READ ena WRITE ina
  RANGE g, gbar, ina
  GLOBAL activate_Q10_rates, activate_Q10_gbar

  : channel state info for recording
  : O is open state
  : B is blocked state
  : I6 is open+inactivated state
  : Ctot = C1+C2+C3+C4+C5 = sum of closed states
  : Itot = I1+I2+I3+I4+I5 = sum of closed+inactivated states
  RANGE O, I6, B, Ctot, Itot
}

UNITS { 
	(mV) = (millivolt)
	(S) = (siemens)
}

CONSTANT {
	ratesQ10 = 3
	gmaxQ10 = 1.980105147e+00
}

PARAMETER {
	activate_Q10_rates = 1
	activate_Q10_gbar = 0

	gbar = 0.014 (S/cm2) 
	celsius (degC)			

	: On and Off rates for inactivating particle
	Con = 0.004			(1/ms)		: closed -> inactivated transitions
	Coff = 9.0			(1/ms)		: inactivated -> closed transitions
	Oon = 0.8			(1/ms)		: open -> Ineg transition
	Ooff = 0.004		(1/ms)		: Ineg -> open transition

	: Transition rates between weakly/tightly bound states
	alpha = 140			(1/ms)		: activation
	beta = 10			(1/ms)		: deactivation
	gamma = 150			(1/ms)		: opening
	delta = 40			(1/ms)		: closing, greater than BEAN/KUO = 0.2
	epsilon = 1e-12		(1/ms)		: open -> Iplus for tau = 0.3 ms at +30 with x5
	zeta = 0.03			(1/ms)		: Iplus -> open for tau = 25 ms at -30 with x6
	: NOTE: by setting epsilon so low, transition to the 'B' state (resurgent current)
	:       is disabled in this mechanism (set this greater to implement resurgent current).
	:       This makes the state diagram

	: Slopes for voltage dependence of transition rates
	x1 = 27				(mV)								: Vdep of activation (alpha)
	x2 = -27			(mV)								: Vdep of deactivation (beta)
	: BELOW no V dependece in Taddese & Bean 2002: set so exp(x)~=1.0
	x3 = 1e12			(mV)								: Vdep of opening (gamma)
	x4 = 1e12			(mV)								: Vdep of closing (delta)
	x5 = 1e12			(mV)								: Vdep into Ipos (epsilon)
	x6 = 1e12			(mV)								: Vdep out of Ipos (zeta)
}

ASSIGNED {
	alfac   				: microscopic reversibility factors
	btfac

	: rates
	f01  		(/ms)
	f02  		(/ms)
	f03 		(/ms)
	f04			(/ms)
	f0O 		(/ms)
	fip 		(/ms)
	f11 		(/ms)
	f12 		(/ms)
	f13 		(/ms)
	f14 		(/ms)
	f1n 		(/ms)
	fi1 		(/ms)
	fi2 		(/ms)
	fi3 		(/ms)
	fi4 		(/ms)
	fi5 		(/ms)
	fin 		(/ms)

	b01 		(/ms)
	b02 		(/ms)
	b03 		(/ms)
	b04		(/ms)
	b0O 		(/ms)
	bip 		(/ms)
	b11  		(/ms)
	b12 		(/ms)
	b13 		(/ms)
	b14 		(/ms)
	b1n 		(/ms)
	bi1 		(/ms)
	bi2 		(/ms)
	bi3 		(/ms)
	bi4 		(/ms)
	bi5 		(/ms)
	bin 		(/ms)
	
	v					(mV)
 	ena					(mV)
	ina 					(milliamp/cm2)
	g					(S/cm2)
	qt
	gqt

	: channel states for recording
	Ctot : sum of closed state fractions
	Itot : sum of inactivated state fractions
}

STATE {
	C1 FROM 0 TO 1
	C2 FROM 0 TO 1
	C3 FROM 0 TO 1
	C4 FROM 0 TO 1
	C5 FROM 0 TO 1
	I1 FROM 0 TO 1
	I2 FROM 0 TO 1
	I3 FROM 0 TO 1
	I4 FROM 0 TO 1
	I5 FROM 0 TO 1
	O FROM 0 TO 1
	B FROM 0 TO 1
	I6 FROM 0 TO 1
}

BREAKPOINT {
	SOLVE activation METHOD sparse
 	g = gbar * gqt * O
 	ina = g * (v - ena)

 	: channel state info for recording
 	Ctot = C1 + C2 + C3 + C4 + C5 : total closed
 	Itot = I1 + I2 + I3 + I4 + I5 : total closed and inactivated
}

INITIAL {
	LOCAL Con1,Coff1,Con5,Coff5
	if (activate_Q10_rates) {
		qt = ratesQ10^((celsius-22 (degC))/10 (degC))
	} else {
		qt = 1
	}
	if (activate_Q10_gbar) {
		gqt = gmaxQ10^((celsius-23 (degC))/10 (degC))
	} else {
		gqt = 1
	}

	Con1 = 0.004 : inactivation rate at tightest level (1)
	Coff1 = 9.0 : recovery rate at tightest level (1)
	Con5 = 0.05 : inactivation rate at weakest bound level (5)
	Coff5 = 0.00025  : recovery rate at weakest bound level (5)
	alfac = (Con5/Con1)^(1.0/4.0) : alfac = 1.88 at room temperature
	btfac = 1.0/((Coff1/Coff5)^(1.0/4.0)) : btfac = 1./13.8 at room temperature
	
	rates(v)

	SOLVE seqinitial
}

KINETIC activation
{
	rates(v)
	~ C1 <-> C2					(f01,b01) : first is A->B, second is A<-B
	~ C2 <-> C3					(f02,b02)
	~ C3 <-> C4					(f03,b03)
	~ C4 <-> C5					(f04,b04)
	~ C5 <-> O					(f0O,b0O)
	~ O <-> B					(fip,bip)
	~ O <-> I6					(fin,bin)
	~ I1 <-> I2					(f11,b11)
	~ I2 <-> I3					(f12,b12)
	~ I3 <-> I4					(f13,b13)
	~ I4 <-> I5					(f14,b14)
	~ I5 <-> I6					(f1n,b1n)
	~ C1 <-> I1					(fi1,bi1)
	~ C2 <-> I2					(fi2,bi2)
	~ C3 <-> I3					(fi3,bi3)
 	~ C4 <-> I4					(fi4,bi4)
 	~ C5 <-> I5					(fi5,bi5)

CONSERVE C1 + C2 + C3 + C4 + C5 + O + B + I1 + I2 + I3 + I4 + I5 + I6 = 1
}

LINEAR seqinitial { : sets initial equilibrium
 ~          I1*bi1 + C2*b01 - C1*(    fi1+f01) = 0 : equilibrium for C1
 ~ C1*f01 + I2*bi2 + C3*b02 - C2*(b01+fi2+f02) = 0 : equilibrium for C2
 ~ C2*f02 + I3*bi3 + C4*b03 - C3*(b02+fi3+f03) = 0 : equilibrium for C3
 ~ C3*f03 + I4*bi4 + C5*b04 - C4*(b03+fi4+f04) = 0 : equilibrium for C4
 ~ C4*f04 + I5*bi5 + O*b0O - C5*(b04+fi5+f0O) = 0  : equilibrium for C5
 ~ C5*f0O + B*bip + I6*bin - O*(b0O+fip+fin) = 0   : equilibrium for O
 ~ O*fip + B*bip = 0 :                               equilibrium for B

 ~          C1*fi1 + I2*b11 - I1*(    bi1+f11) = 0 : equilibrium for I1
 ~ I1*f11 + C2*fi2 + I3*b12 - I2*(b11+bi2+f12) = 0 : equilibrium for I2
 ~ I2*f12 + C3*fi3 + I4*bi3 - I3*(b12+bi3+f13) = 0 : equilibrium for I3
 ~ I3*f13 + C4*fi4 + I5*b14 - I4*(b13+bi4+f14) = 0 : equilibrium for I4
 ~ I4*f14 + C5*fi5 + I6*b1n - I5*(b14+bi5+f1n) = 0 : equilibrium for I5
 
 ~ C1 + C2 + C3 + C4 + C5 + O + B + I1 + I2 + I3 + I4 + I5 + I6 = 1
}

PROCEDURE rates(v(mV) )
{
 TABLE f01,f02,f03,f04,f0O,fip,f11,f12,f13,f14,f1n,fi1,fi2,fi3,fi4,fi5,fin,b01,b02,b03,b04,b0O,bip,b11,b12,b13,b14,b1n,bi1,bi2,bi3,bi4,bi5,bin DEPEND celsius FROM -100 TO 100 WITH 400
 
 : Top row, right direction: C(i) -> C(i+1) (weaker binding)
 f01 = 4 * alpha * exp(v/x1) * qt
 f02 = 3 * alpha * exp(v/x1) * qt
 f03 = 2 * alpha * exp(v/x1) * qt
 f04 = 1 * alpha * exp(v/x1) * qt
 f0O = gamma * exp(v/x3) * qt
 fip = epsilon * exp(v/x5) * qt : Open -> Blocked

 : Top row, left direction: C(i) -> C(i-1) (tighter binding)
 b01 = 1 * beta * exp(v/x2) * qt
 b02 = 2 * beta * exp(v/x2) * qt
 b03 = 3 * beta * exp(v/x2) * qt
 b04 = 4 * beta * exp(v/x2) * qt
 b0O = delta * exp(v/x4) * qt
 bip = zeta * exp(v/x6) * qt

 : Bottom row, right direction: I(i) -> I(i+1) (weaker binding)
 f11 = 4 * alpha * alfac * exp(v/x1) * qt
 f12 = 3 * alpha * alfac * exp(v/x1) * qt
 f13 = 2 * alpha * alfac * exp(v/x1) * qt
 f14 = 1 * alpha * alfac * exp(v/x1) * qt
 f1n = gamma * exp(v/x3) * qt

 : Bottom row, left direction: I(i) -> I(i-1) (tighter binding)
 b11 = 1 * beta * btfac * exp(v/x2) * qt
 b12 = 2 * beta * btfac * exp(v/x2) * qt
 b13 = 3 * beta * btfac * exp(v/x2) * qt
 b14 = 4 * beta * btfac * exp(v/x2) * qt
 b1n = delta * exp(v/x4) * qt

 : Vertical down (inactivation)
 fi1 = Con * qt
 fi2 = Con * alfac * qt
 fi3 = Con * alfac^2 * qt
 fi4 = Con * alfac^3 * qt
 fi5 = Con * alfac^4 * qt
 fin = Oon * qt

 : Vertical up (deinactivation)
 bi1 = Coff * qt
 bi2 = Coff * btfac * qt
 bi3 = Coff * btfac^2 * qt
 bi4 = Coff * btfac^3 * qt
 bi5 = Coff * btfac^4 * qt
 bin = Ooff * qt
}
TITLE potassium delayed rectifier membrane channels for STh

COMMENT 

 Delayed rectifier potassium from pyramidal, Traub 1991. He
 based them on Sah (1988) data, which were at 22-24degC, but he scaled
 them so that they "were fast enough"?  

 How the q10 works: There is a q10 for the rates (alpha and beta's)
 called Q10 and a Q10 for the maximum conductance called gmaxQ10.  The
 q10s should have been measured at specific temperatures temp1 and
 temp2 (that are 10degC apart). Ideally, as Q10 is temperature
 dependant, we should know these two temperatures.  We are going to
 follow the more formal Arrhenius derived Q10 approach.  The
 temperature at which this channel's kinetics were recorded is tempb
 (base temperature).  What we then need to calculate is the desired
 rate scale for now working at temperature celsius (rate_k).  This is
 given by the empirical Arrhenius equation, using the Q10. 
ENDCOMMENT

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX KDR
	USEION k READ ki,ek WRITE ik
	RANGE gk
	GLOBAL rest,activate_Q10,Q10,gmaxQ10,rate_k,gmax_k,temp1,temp2,tempb
}

PARAMETER {
        v (mV)
	dt (ms)
	gk    = 3.842910637e-03 (mho/cm2)
	rest  = -60.0 (mV) : for conversion from Traub
	ek
	ki
	celsius
	
	activate_Q10 = 1
	Q10 = 1.200000603e+00
	gmaxQ10 = 1.200000603e+00
	temp1 = 19.0 (degC)
	temp2 = 29.0 (degC)
	tempb = 23.0 (degC)
}

STATE {
        n 
}

ASSIGNED { 
	ik (mA/cm2)
	alphan (/ms)
	betan (/ms)
	rate_k
	gmax_k
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik   = (gk*gmax_k)*n*(v-ek)
}

UNITSOFF

INITIAL {
	LOCAL ktemp,ktempb,ktemp1,ktemp2
	if (activate_Q10>0) {
	  ktemp  = celsius+273.0
	  ktempb = tempb+273.0
	  ktemp1 = temp1+273.0
	  ktemp2 = temp2+273.0
	  rate_k = exp( log(Q10)*((1/ktempb)-(1/ktemp))/((1/ktemp1)-(1/ktemp2)) )
	  gmax_k = exp( log(gmaxQ10)*((1/ktempb)-(1/ktemp))/((1/ktemp1)-(1/ktemp2)) )
	}else{
	  : Note, its not 1.0, as we have rescaled the kinetics
          :  (reverting the scaleing Traub did), the original is
          :  acheived using this rate
	  rate_k = 1.60
	  gmax_k = 1.60
	}
        settables(v)
	n = alphan/(alphan+betan)
}

DERIVATIVE states {
	settables(v)      :Computes state variables at the current v and dt.
	n' = alphan * (1-n) - betan * n
}

PROCEDURE settables(v) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
                          :Voltage shift (for temp effects) of 0.60650122.
        LOCAL vadj
        TABLE alphan, betan DEPEND rest,celsius FROM -100 TO 100 WITH 400
	vadj  = v - rest + 0.60650122

	        :"n" potassium activation system
	alphan = rate_k * 0.01 * vtrap((35.1-vadj),5.0)
        betan = rate_k * 0.156 * exp((20.0-vadj)/40.0)
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON


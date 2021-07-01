TITLE potassium Kv3.1 membrane channel for STh

COMMENT
 Potassium Kv3.1 membrane channel for STh.  Based on derived kinetics
 from Wigmore & Lacey 2000.  Their primary experiments were performed
 at 32degC

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
	SUFFIX Kv31
	USEION k READ ki,ek WRITE ik
	RANGE gk
	GLOBAL activate_Q10,Q10,gmaxQ10,rate_k,gmax_k,temp1,temp2,tempb
}

PARAMETER {
        v (mV)
	dt (ms)
	gk   = 0.015 (mho/cm2)
	ek
	ki
	celsius

	activate_Q10 = 1
	Q10 = 1.700025939e+00
	gmaxQ10 = 1.700025939e+00
	temp1 = 20.0 (degC)
	temp2 = 30.0 (degC)
	tempb = 32.0 (degC)
}

STATE {
        p 
}

ASSIGNED { 
	ik (mA/cm2)
	pinf
	ptau (ms)
	rate_k
	gmax_k
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik   = (gk*gmax_k)*p*(v-ek)
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
	  rate_k = 1.0
	  gmax_k = 1.0
	}
        settables(v)
	p = pinf
}

DERIVATIVE states {  
        settables(v)
	p' = (pinf-p)/ptau
}

PROCEDURE settables(v) {:Computes rate and other constants at current v.
                        :Call once from HOC to initialize inf at resting v.
			:Voltage shift (for temp effects) of -5.08
	TABLE pinf, ptau DEPEND celsius FROM -100 TO 100 WITH 400

	pinf = 1.0/(1.0+exp((v + -0.083699749)/ -9.0))
	ptau = ((7.3/(exp((v + 32.9163003)/-14.0)+exp((v + 2.91630025)/16.0)))+1) / rate_k

}

UNITSON


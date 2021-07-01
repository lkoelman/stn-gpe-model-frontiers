TITLE sodium membrane channels for STh

COMMENT
 Sodium from pyramidal, Traub 91.  He based them on Sah (1988) data,
 which were at 22degC, but he scaled them so that they "were fast
 enough"?  The Q10 measured from Sah (at least for the peak
 conductance) was 1.5

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
	SUFFIX Na
	USEION na READ nai,ena WRITE ina
	RANGE gna, m, h
	GLOBAL rest,activate_Q10,Q10,gmaxQ10,rate_k,gmax_k,temp1,temp2,tempb
}

PARAMETER {
        v (mV)
	dt (ms)
	gna   = 1.0e-7 (mho/cm2)
	rest  = -60.0 (mV) : for conversion from Traub WAS -60
	ena
	nai
	celsius
	
	activate_Q10 = 1
	Q10 = 1.980105147e+00
	gmaxQ10 = 1.980105147e+00
	temp1 = 19.0 (degC)
	temp2 = 29.0 (degC)
	tempb = 23.0 (degC)
}

STATE {
        m h  
}

ASSIGNED { 
        ina (mA/cm2)
	alpham (/ms)
	betam (/ms)
	alphah (/ms)
	betah (/ms)
	rate_k
	gmax_k
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina  = (gna*gmax_k)*m*m*h*(v-ena)
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
        m = alpham/(alpham+betam)
        h = alphah/(alphah+betah)

}

DERIVATIVE states {
	settables(v)      :Computes state variables at the current v and dt.
	m' = alpham * (1-m) - betam * m
	h' = alphah * (1-h) - betah * h
}

PROCEDURE settables(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL vadj
        TABLE alpham, betam, alphah, betah DEPEND rest,celsius FROM -100 TO 100 WITH 400
	vadj  = v - rest

		:"m" sodium activation system
	alpham = rate_k * 0.2 * vtrap((13.1-vadj),4.0)
        betam =  rate_k * 0.175 * vtrap((vadj-40.1),1.0)

                :"h" sodium inactivation system
        alphah = rate_k * 0.08 * exp((17.0-vadj)/18.0)
        betah = rate_k * 2.5 / (exp((40.0-vadj)/5.0) + 1)
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON


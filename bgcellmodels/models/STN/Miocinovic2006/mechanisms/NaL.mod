TITLE sodium leak

COMMENT 

 Sodium leak.  This TTX sensitive sodium current is active
 between spikes and has some voltage dependence (although no
 inactivation that I know of) (Do & Bean 2003).  I based the
 voltage dependence of this channel on the Do & Bean data, but found a
 pure leak had the same effect, so I stuck to the pure leak..

 How the q10 works: There is a q10 for the rates (alpha and beta's)
 called Q10 and a Q10 for the maximum conductance called gmaxQ10.
 Here, we only use gmaxQ10.  The q10s should have been measured at
 specific temperatures temp1 and temp2 (that are 10degC
 apart). Ideally, as Q10 is temperature dependant, we should know
 these two temperatures.  We are going to follow the more formal
 Arrhenius derived Q10 approach.  The temperature at which this
 channel's kinetics were recorded is tempb (base temperature).  What
 we then need to calculate is the desired rate scale for now working
 at temperature celsius (rate_k).  This is given by the empirical
 Arrhenius equation, using the Q10. 
ENDCOMMENT

UNITS {
      (mv) = (millivolt)
      (mA) = (milliamp)
}

NEURON {
       SUFFIX NaL
       USEION na READ ena,nai WRITE ina
       RANGE gna,inaL
       GLOBAL activate_Q10,gmaxQ10,gmax_k,temp1,temp2,tempb
}

PARAMETER {
        v (mV)
	gna = 0.81e-5 (mho/cm2)
	inaL = 0.0 (mA/cm2)
	ena
	nai
	celsius

	activate_Q10 = 1
	gmaxQ10 = 1.5
	temp1 = 25.0 (degC)
	temp2 = 35.0 (degC)
	tempb = 23.0 (degC)	
}

ASSIGNED { 
        ina (mA/cm2)
	gmax_k
}

BREAKPOINT {
	   ina	= gna*gmax_k*(v-ena)
	   inaL = ina
}
UNITSOFF

INITIAL {
	LOCAL ktemp,ktempb,ktemp1,ktemp2
	if (activate_Q10>0) {
	  ktemp  = celsius+273.0
	  ktempb = tempb+273.0
	  ktemp1 = temp1+273.0
	  ktemp2 = temp2+273.0
	  gmax_k = exp( log(gmaxQ10)*((1/ktempb)-(1/ktemp))/((1/ktemp1)-(1/ktemp2)) )
	}else{
	  gmax_k = 1.0
	}
}	
UNITSON

TITLE Basic Potassium current
 
COMMENT
 from "Gamma Oscillation by Synaptic Inhibition in a Hippocampal Interneuronal Network Model" (Wang and Buzsaki 1996)
 Used in Role of a Striatal Slowly Inactivating Potassion Current in Short-term Facilitation of Corticostriatal Inputs" A computer Simulation Study" (Mahon et al. 2000)
Implemented by Kevin M. Biddell kevin.biddell@gmail.com
7/12/06

NOTE: 1S=1mho Neuron wants the units in mhos not millisiemens, please note the conversion!

Phi =5 and no q10 or temp adjustment according to Bruno Delord 11/13/06

ENDCOMMENT
 
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}
 
NEURON {
    SUFFIX Km
    USEION k WRITE ik
    RANGE gkmbar, gkm
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
    
    ek  = -90   (mV)
    gkmbar  = 0.006 (mho/cm2) : 6mS
    phi = 5 < 0, 1e9 > : according to Bruno delord 11/13/06
    Van = -27 :NOT the original value from wang and Buzsaki
    Kan = 1
    Vbn = -37 :NOT the original value from wang and Buzsaki
    Kbn = 80
}
 
STATE {
    n
}
 
ASSIGNED {
    v       (mV)
    ik      (mA/cm2)
    celsius (degC)
    ninf
    taun
    gkm
}
 
BREAKPOINT {
    SOLVE states METHOD cnexp
    gkm = gkmbar*n^4
    ik = gkm*(v - ek)
  
}
 
UNITSOFF
 
INITIAL {
    rates(v)
    n = ninf
}

DERIVATIVE states {
    rates(v)
    : n' = phi*(an*(1-n)-Bn*n)
    n' = (ninf - n) / taun

}
 
PROCEDURE rates(v) {
    LOCAL an, Bn
    TABLE ninf, taun FROM -100 TO 100 WITH 400

    : an = (-0.01*(v-Van)/Kan / (exp(-0.1*(v-Van)/Kan)-1))
    an = 0.1/Kan * vtrap(-0.1*(v-Van), Kan)
    Bn = 0.125*exp(-(v-Vbn)/Kbn)
    
    ninf = an / (an+Bn)
    taun = 1 / ((an+Bn) * phi)  
          
}

: Traps for 0 in denominator of rate eqns.
FUNCTION vtrap(x,y) {  
    if (fabs(x/y) < 1e-6) {
            vtrap = y*(1 - x/y/2)
    }else{
            vtrap = x/(exp(x/y) - 1)
    }
}
 
UNITSON


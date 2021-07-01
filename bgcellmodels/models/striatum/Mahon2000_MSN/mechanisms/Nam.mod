TITLE Basic sodium current
 
COMMENT
 from "Gamma Oscillation by Synaptic Inhibition in a Hippocampal Interneuronal Network Model" (Wang and Buzsaki 1996)
 Used in Role of a Striatal Slowly Inactivating Potassion Current in Short-term Facilitation of Corticostriatal Inputs" A computer Simulation Study" (Mahon et al. 2000)
Implemented by Kevin M. Biddell kevin.biddell@gmail.com
7/11/06

NOTE: 1S=1mho Neuron wants the units in mhos not millisiemens, please note the conversion!

Phi =5 and no q10 or temp adjustment according to Bruno Delord 11/13/06

ENDCOMMENT
 
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt) 
}
 
NEURON {
    SUFFIX Nam
    USEION na WRITE ina
    RANGE gnabar, gna
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
    
    ena = 55    (mV)
    gnabar  = 0.035 (mho/cm2) : 35mS
    phi = 5 < 0, 1e9 > : from delord 11/13/06 
    Vam = -28 :NOT the original value from wang and Buzsaki
    Kam = 1
    Vbm = -53 :NOT the original value from wang and Buzsaki
    Kbm = 18
    Vah = -51 :NOT the original value from wang and Buzsaki
    Kah = 20
    Vbh = -21 :NOT the original value from wang and Buzsaki
    Kbh = 1
           
}
 
STATE {
    : m is always minf
    h
}
 
ASSIGNED {
    v       (mV)
    ina     (mA/cm2)
    celsius (degC)
    minf
    hinf
    tauh
    gna
}
 
BREAKPOINT {
    SOLVE states METHOD cnexp
    gna = gnabar * minf^3 * h
    ina = gna * (v - ena)
  
}
 
UNITSOFF
 
INITIAL {
    rates(v)
    : m = minf
    h = hinf
}

DERIVATIVE states {
    rates(v)
    : h' = phi * (ah*(1-h) - Bh*h)
    h' = (hinf - h) / tauh
}
 
PROCEDURE rates(v) {
    LOCAL am, Bm, ah, Bh
    TABLE minf, hinf, tauh FROM -100 TO 100 WITH 400
    : DEPEND Vam, Kam, Vbm, Kbm, Vah, Kah, Vbh, Kbh
        
    : am = (-0.1*(v-Vam)/Kam / (exp(-0.1 * (v-Vam)/Kam) - 1))
    am = 1/Kam * vtrap(-0.1*(v-Vam), Kam)
    Bm = 4 * exp(-(v-Vbm)/Kbm)
    ah = 0.07 * exp(-(v-Vah)/Kah)
    Bh =  1 / (1 + exp(-0.1*(v-Vbh)/Kbh))
    
    minf = am / (am+Bm)
    hinf = ah / (ah+Bh)
    : tauh = 1 / (ah+Bh)
    tauh = 1 / ((ah+Bh) * phi) : phi moved here from state equation
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


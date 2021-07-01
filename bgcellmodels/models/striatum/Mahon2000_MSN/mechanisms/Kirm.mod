TITLE Inward rectifying Potassium current
 
COMMENT
  Used in Role of a Striatal Slowly Inactivating Potassion Current in Short-term 
  Facilitation of Corticostriatal Inputs" A computer Simulation Study" (Mahon et al. 2000)
  Implemented by Kevin M. Biddell kevin.biddell@gmail.com
  7/12/06
note: activation is instantaneous but we will assume for implementation purposes that the mtau is small enough to simulate this behavior on our time scale
NOTE: 1S=1mho Neuron wants the units in mhos not millisiemens, please note the conversion!
ENDCOMMENT
 
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}
 
NEURON {
    SUFFIX Kirm
    USEION k WRITE ik
    RANGE gkirmbar, gkirm, minf
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
    ek  = -90 (mV)
    gkirmbar= 0.00015 (mho/cm2) :0.15
    tau = 0.01   (ms)
    Vsm = -100
    ksm = -10   
}
 
STATE {
    m
}
 
ASSIGNED {
    v       (mV)
    ik      (mA/cm2)
    celsius (degC)
    minf
    mtau
    gkirm
}
 
BREAKPOINT {
    SOLVE states METHOD cnexp
    gkirm = gkirmbar*m
    ik = gkirm*(v - ek)
}
 
UNITSOFF
 
INITIAL {
    rates(v)
    m = minf
    mtau = tau  
}

DERIVATIVE states {  :Computes states variable m at the current v and dt.
    rates(v)       
    m' = ( minf - m ) / mtau
}
 
PROCEDURE rates(v) {  : Computes rate and other constants at current v.
                      : Call once from HOC to initialize inf at resting v.
                      : The activation is instantaneous therefore not temperature dependent
    TABLE minf FROM -100 TO 100 WITH 400
    minf = 1/(1+exp(-(v-Vsm)/ksm))   
}
 
UNITSON


TITLE A fast A-type Potassium current
 
COMMENT
  Used in Role of a Striatal Slowly Inactivating Potassion Current in Short-term 
  Facilitation of Corticostriatal Inputs" A computer Simulation Study" (Mahon et al. 2000)
  Implemented by Kevin M. Biddell kevin.biddell@gmail.com
  7/13/06
NOTE: 1S=1mho Neuron wants the units in mhos not millisiemens, please note the conversion!
ENDCOMMENT
 
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}
 
NEURON {
    SUFFIX KAfm
    USEION k WRITE ik
    RANGE gkafmbar, gkafm
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
    ek  = -73   (mV)
    gkafmbar= 0.00009 (mho/cm2) :0.09 mS
    Etemp   = 22 :delord correspondence 11/15/06
    Vsm = -33.1
    ksm = 7.5
    tom = 1
    Vsh = -70.4
    ksh = -7.6
    toh = 25         
}
 
STATE {
    m h
}
 
ASSIGNED {
    ik      (mA/cm2)
    v       (mV)
    celsius (degC)
    minf
    hinf
    mtau
    htau
    gkafm
    
}
 
BREAKPOINT {
    SOLVE states METHOD cnexp
    gkafm = gkafmbar*m*h
    ik = gkafm*(v - ek)
  
}
 
UNITSOFF
 
INITIAL {
    rates(v)
    m= minf
    h= hinf
}

DERIVATIVE states {  :Computes states variable m at the current v and dt.
    rates(v)      
    m' = ( minf - m ) / mtau
    h' = (hinf - h ) / htau
}
 
PROCEDURE rates(v) {  :Computes rate and other constants at current v. Call once from HOC to initialize inf at resting v.
    LOCAL  q10,tadj
    TABLE minf, hinf, mtau, htau DEPEND celsius FROM -100 TO 100 WITH 400
    
    q10 = 2.5
    tadj=q10^((celsius-Etemp)/10)
    
    minf=1/(1+exp(-(v-Vsm)/ksm))
    hinf=1/(1+exp(-(v-Vsh)/ksh))
    mtau=tom/tadj
    htau=toh/tadj         
}
 
UNITSON


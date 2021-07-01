TITLE A persistant Sodium current
 
COMMENT
  Used in Role of a Striatal Slowly Inactivating Potassion Current in Short-term 
  Facilitation of Corticostriatal Inputs" A computer Simulation Study" (Mahon et al. 2000)
  Implemented by Kevin M. Biddell kevin.biddell@gmail.com
  7/17/06
NOTE: 1S=1mho Neuron wants the units in mhos not millisiemens, please note the conversion!
ENDCOMMENT
 
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}
 
NEURON {
    SUFFIX NaPm
    USEION na WRITE ina
    RANGE gnapmbar, gnapm, minf, mtau
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
    
    ena = 45    (mV)
    gnapmbar= 0.00002 (mho/cm2) :0.02
    Etemp   = 22    :delord correspondence 11/15/06
    Vsm = -47.8
    ksm = 3.1
    tom = 1
    
             
}
 
STATE {
        m
}
 
ASSIGNED {
    v       (mV)
    ina     (mA/cm2)
    celsius (degC)
    minf
    mtau
    gnapm
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        gnapm = gnapmbar*m
        ina = gnapm*(v - ena)
  
}
 
UNITSOFF
 
INITIAL {
    rates(v)
    m = minf
}

DERIVATIVE states {  :Computes states variable m at the current v and dt.
    rates(v)      
    m' = ( minf - m ) / mtau
}
 
PROCEDURE rates(v) {  :Computes rate and other constants at current v. Call once from HOC to initialize inf at resting v.
    LOCAL  q10, tadj
    TABLE minf, mtau DEPEND celsius FROM -100 TO 100 WITH 400
    
    q10 = 2.5
    tadj = q10^((celsius-Etemp)/10)
    
    minf = 1/(1+exp(-(v-Vsm)/ksm))
    mtau = tom/tadj
      
}
 
UNITSON


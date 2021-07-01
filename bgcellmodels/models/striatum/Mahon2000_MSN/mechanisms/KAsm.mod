TITLE A Slow A-type Potassium current
 
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
    SUFFIX KAsm
    USEION k WRITE ik
    RANGE gkasmbar, gkasm, minf, hinf, mtau, htau
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
    
    ek  = -85   (mV)
    gkasmbar= 0.00032 (mho/cm2) :0.32mS
    Etemp   = 22  : Delord correspondence 11/15/06 and fitting

    Vsm = -25.6
    ksm = 13.3
    Vsh = -78.8
    ksh = -10.4
    tom = 131.4
    Vtm = -37.4
    ktm = 27.3
    Vth = -38.2
    kth = 28
    hint = 0.46     : Delord correspondence 11/15/06 = 0.2968

}
 
STATE {
    m h
}
 
ASSIGNED {
    v       (mV)
    ik      (mA/cm2)
    celsius (degC)
    minf
    hinf
    mtau
    htau
    gkasm
}
 
BREAKPOINT {
    SOLVE states METHOD cnexp
    gkasm = gkasmbar*m*h
    ik = gkasm*(v - ek)
  
}
 
UNITSOFF
 
INITIAL {
    rates(v)
    m= minf
    h= hint
}

DERIVATIVE states {  :Computes states variable m at the current v and dt.
        rates(v)      
       
    m' = ( minf - m ) / mtau
    h' = (hinf - h ) / htau
}
 
PROCEDURE rates(v) {  :Computes rate and other constants at current v. Call once from HOC to initialize inf at resting v.
    LOCAL  q10, tadj
    TABLE minf, hinf, mtau, htau DEPEND celsius FROM -100 TO 100 WITH 400
    
    q10 = 2.5
    tadj = q10^((celsius-Etemp)/10)
    
    minf = 1/(1+exp(-(v-Vsm)/ksm))
    hinf = 1/(1+exp(-(v-Vsh)/ksh))
    mtau = tom/(exp(-(v-Vtm)/ktm)+exp((v-Vtm)/ktm))/tadj
    htau = (1790+2930*exp(-((v-Vth)/kth)^2)*((v-Vth)/kth))/tadj
          
}
 
UNITSON


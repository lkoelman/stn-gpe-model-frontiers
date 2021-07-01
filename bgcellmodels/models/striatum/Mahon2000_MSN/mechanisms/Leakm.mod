TITLE Leakm current
:NOTE 1S=1mho Neuron wants the units in mhos not millisiemens, please note the conversion! 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Leakm
        NONSPECIFIC_CURRENT il
        RANGE  gl, el
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
       
        gl = .000075 (mho/cm2) :0.075 mS
        el = -75 (mV) : -75 mV in Mahon et al, -90 mV in Corbit et al.
}
  
ASSIGNED {
	 v (mV)
        il (mA/cm2)
}
 
BREAKPOINT {
:        SOLVE states
        il = gl*(v - el)
}
 

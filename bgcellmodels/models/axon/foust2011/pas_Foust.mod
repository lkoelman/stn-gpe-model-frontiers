TITLE Passive Leak Current for Cortical Neuron Axon Compartments

COMMENT
  
  Model Reference: 
  
  Foust, A.J., Yu, Y., Popovic, M., Zecevic, D. and McCormick, D.A., 
  2011. "Somatic membrane potential and Kv1 channels control spike 
  repolarization in cortical axon collaterals and presynaptic boutons." 
  Journal of Neuroscience, 31(43), pp.15490-15498.
  
  Implemented by John Fleming - john.fleming@ucdconnect.ie - 06/12/18
  
  Edits: 
  
ENDCOMMENT


UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
 (S) = (siemens)
}

NEURON {
	SUFFIX pas_Foust
	NONSPECIFIC_CURRENT i			: Declare ASSIGNED variables as RANGE variables so that they can be accessed outside of mod file
	RANGE i, g, e				: leak current, specific conductance and equilibrium potential
}

PARAMETER {
	g = 3.33333e-5 (S/cm2)
	e = -70 (mV)
}

ASSIGNED {
	v (mV)	
	i (mA/cm2)
}

BREAKPOINT {
	i = g*(v - e)
}
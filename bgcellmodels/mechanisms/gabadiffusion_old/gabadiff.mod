TITLE Extracellular GABA diffusion

COMMENT
TODO : check units of differential equation (changed current from mA
       which is mC/s to mM/ms which is concentration/time i.e. 
       quantity / (time * volume))

USAGE
-----

In synapse mechanism add the following:

  NEURON { USEION gaba WRITE igaba ... }
  ASSIGNED { igaba (mM/ms) ... }
  BREAKPOINT { igaba = value ... }

ENDCOMMENT


NEURON {
  SUFFIX gabadiff

  : gaba is an 'ion species' declared here and will therefore have an
  : automatically generated mechanism 'gaba_ion' with associated global
  : variables 'gabao0_gaba_ion' and 'gabai0_gaba_ion'
  USEION gaba READ gabao, igaba WRITE gabao VALENCE 0

  RANGE Dgaba, Dleak, Vmax, Km
  RANGE dshell, gabac
}

UNITS {
  (molar) = (1/liter)
  (mM)    = (millimolar)
  (um)    = (micron)
  (mA)    = (milliamp)
  FARADAY = (faraday)  (10000 coulomb)
  PI      = (pi)       (1)
}

PARAMETER {

  : article eq. 1 : diffusion constant D = 8e-6 cm^2/s
  Dgaba = 8e-1  (um^2/ms)
  
  : article p. 9515 : Km = 4 uM estimated from studies (ref. 21)
  Km = 0.004    (mM)
  
  : article p. 9515 : Vmax = 0.1 M/s = 0.1 mM/ms
  Vmax = 0.1    (mM/ms)
  
  : article p. 9516 : a leak in each compartment of 1e-8 cm^2/s
  : leak is optional, can set to zero
  Dleak = 1e-3  (um^2/ms)
  
  : Size of shell that defines volume of diffusion
  dshell = 1.0  (um)
  vshell        (um)
}

CONSTANT {
  facD = 1e5 : 1 cm2/s = 1e5 um2/ms
}

: variables defined elsewhere that are read or assigned
ASSIGNED {
  igaba     (mA/cm2) : igaba is normally mA/cm2 (defined as ion) but interpret as mM/ms
  gabao     (mM)
  release
  uptake
  leakage
}


STATE {
  gabac      (mM)
}

BREAKPOINT {
  SOLVE state METHOD sparse
}

INITIAL {
   gabac = gabao
   vshell = dshell * dshell * dshell
}

KINETIC state {
  : the COMPARTMENT <volume> multiplies the dSTATE/dt left hand side 
  : of the equivalent differential equations, converting it to an extensive 
  : quantity and making it consistent with flux terms in units of absolute 
  : quantity per time. Because of the COMPARTMENT statement, the left hand side 
  : of the differential equation is not d[ion]/dt but d(total ion in the shell)/dt.
  : This is consistent with the right hand side of the equation, which is in 
  : units of mass per time.
  COMPARTMENT dshell*dshell*dshell {gabac}
  : in LONGITUDINAL_DIFFUSION <flux>, flux is the product of the diffusion
  : constant and the area of the cross section between adjacent compartments. 
  : Units of the flux must be (um^4/ms),
  LONGITUDINAL_DIFFUSION Dgaba*dshell*dshell {gabac} : GABA diffusion in extracellular space

  : GABA release from synapses
  release = vshell * igaba                          : um3 * mM / ms

  : GABA uptake (Destexhe & Sejnowsi 1995, eq. 1)
  uptake = (-Vmax * gabac * vshell / (gabac + Km))  : um3 * mM / ms

  : GABA leakage
  leakage = (-Dleak * dshell * gabac)               : um3 * mM / ms

  : GABA total flux
  ~ gabac << (release + uptake + leakage)

  : GABA buffering (not present in article)
  : ~ gaba + buffer <-> gababuffer  (k1buf*dsqvol, k2buf*dsqvol)
  
  gabao = gabac
}
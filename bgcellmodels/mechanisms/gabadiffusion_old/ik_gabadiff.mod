TITLE K-current linked via G-protein to GABA-B receptors

COMMENT

Kinetic model of GABA-B receptors by Destexhe & Sejnowski (1995),
with GABA concentration read from external mechanism.

Extracellular GABA was defined as ion species 'gabao' in NEURON,
see mechanism gababuff.mod.

For details on GABA-B signaling cascade and its parameters,
see mechanism 'gagab3.mod' obtainable at http://cns.iaf.cnrs-gif.fr/alain_demos.html
under package "Kinetic Models of Synaptic Transmission (zip format)"

 TODO : gpcr coupling from Destexhe data
  - gating function
  - add coupling to GPCR (Destexhe: bound receptor fraction and G-protein
    concentration) - or do it in gababuff mechanism
ENDCOMMENT


NEURON {
  POINT_PROCESS ik_gabadiff

  : gaba is an 'ion species' declared here and will therefore have an
  : automatically generated mechanism 'gaba_ion' with associated global
  : variables 'gabao0_gaba_ion' and 'gabai0_gaba_ion'
  USEION gaba READ gabao

  RANGE R, G, g, gmax
  NONSPECIFIC_CURRENT i
  RANGE K1, K2, K3, K4, KD, Erev
}
UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (umho) = (micromho)
  (mM) = (milli/liter)
}

PARAMETER {
:
: From Kfit with long pulse (5ms 0.5mM)
:
  K1  = 0.52  (/ms mM)  : forward binding rate to receptor
  K2  = 0.0013 (/ms)    : backward (unbinding) rate of receptor
  K3  = 0.098 (/ms)     : rate of G-protein production
  K4  = 0.033 (/ms)     : rate of G-protein decay
  KD  = 100             : dissociation constant of K+ channel
  n   = 4               : nb of binding sites of G-protein on K+
  Erev  = -95 (mV)      : reversal potential (E_K)
  gmax = 0.001 (umho)   : maximum conductance
}


ASSIGNED {
  gabao   (mM)    : extracellular gaba concentration (replaces C)
  v       (mV)    : postsynaptic voltage
  i       (nA)    : current = g*(v - Erev)
  g       (umho)  : conductance
  Gn              : exponentiated G
}


STATE {
  R       : fraction of activated receptor
  G       : fraction of activated G-protein
}


INITIAL {
  R = 0
  G = 0
}

BREAKPOINT {
  SOLVE bindkin METHOD euler
  Gn = G^n
  g = gmax * Gn / (Gn+KD)
  i = g*(v - Erev)
}


DERIVATIVE bindkin {

  R' = K1 * gabao * (1-R) - K2 * R : = K1 * gabao - R * (K1 * gabao + K2)
  G' = K3 * R - K4 * G

}
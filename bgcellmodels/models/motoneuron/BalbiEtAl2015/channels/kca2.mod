
 COMMENT
 
 kca2.mod
 
 Calcium-dependent potassium channel
 Based on
 Pennefather (1990) -- sympathetic ganglion cells taken from
 Reuveni et al (1993) -- neocortical cells
 
 Author: Zach Mainen, Salk Institute, 1995, zach@salk.edu
 modified Jan,2000 by RKP;modified July 2005 to include a contribution from calcium
 flowing through L channels
 	
 ENDCOMMENT

 NEURON {
 	SUFFIX kca2
 	USEION k READ ek WRITE ik
 	USEION ca READ ica WRITE cai
    USEION caL READ icaL VALENCE 2
    RANGE n, g,ik,cai,ica,icaL,depth1,taur1,depth2,taur2
 	GLOBAL Ra, Rb, caix
 }
 
 UNITS {
 	(mA)    = (milliamp)
 	(mV)    = (millivolt)
 	(S)     = (siemens)
 	(um)    = (micron)
 	(molar) = (1/liter)			: moles do not appear in units
 	(mM)	= (millimolar)
 	(msM)	= (ms mM)
 	FARADAY = (faraday) (coulomb)
 } 
 
 PARAMETER {
 	g       = 0.03  (S/cm2)	
 	v 		        (mV)
 	cai  		    (mM)
 	caix    = 2	
    cainf   = 0.0001
 	depth1	= .1	(um)		: depth of shell
 	taur1	= 20	(ms)		: rate of calcium removal
 	depth2	= 10	(um)		: depth of shell
 	taur2	= 200	(ms)		: rate of calcium removal
								
    Ra   = 0.1		: max act rate  
 	Rb   = 0.1		: max deact rate 
 
 	celsius		(degC)
 } 
 
 
 ASSIGNED {
 	ik 	            (mA/cm2)
	ica             (mA/cm2)
 	icaL            (mA/cm2)
 	ek	            (mV)
 	ninf
 	ntau 		    (ms)	
    drive_channel1	(mM/ms)
    drive_channel2	(mM/ms)
 }
  
 
 STATE { 
 	n 
 	ca (mM)
	caL (mM)
}
 
 INITIAL { 
	ca  = cainf
 	caL = 0
    cai = cainf
    rates(cai)
 	n   = ninf
 }
 
 BREAKPOINT {
         SOLVE states METHOD cnexp
 		 ik = g * n * (v - ek)
 } 
 

DERIVATIVE states {  
 	drive_channel1 = - (10000) * ica/ (2 * FARADAY * depth1)
 	if (drive_channel1 <= 0.) { drive_channel1 = 0. }	: cannot pump inward
 	ca' = drive_channel1 + (cainf-ca)/taur1
 	drive_channel2 =  - (10000) * icaL/ (2 * FARADAY * depth2)
 		: conversione della corrente del calcio in flusso di ione (mM/ms)
 		: depth2 = spessore della regione intracitoplasmatica sub-membrana
 		: nella quale si disperde la concentrazione di ione
 	if (drive_channel2 <= 0.) { drive_channel2 = 0. }	: cannot pump inward
 	caL' = drive_channel2 + (cainf-caL)/taur2
 	cai = ca + caL

    rates(cai)    
    n' = (ninf-n)/ntau
}
PROCEDURE rates(cai(mM)) {  LOCAL a,b
					UNITSOFF
         a = Ra * (1e3*(cai - cainf))^caix		: rate constant depends on cai in uM
         b = Rb
         ntau = 1/(a+b)
         ninf = a*ntau
					UNITSON
 }

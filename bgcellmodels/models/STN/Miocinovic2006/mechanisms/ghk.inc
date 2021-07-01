

FUNCTION ghkg(v(mV), ci(mM), co(mM), z) (mV) {
        LOCAL nu,f,enu,fnu
              : Here we calculate an effective drive from the GHK equation
              : define
              :    f   = 10^3 RT/(zF)
              :    nu  = v/f  
              :        = z v10^-3 F / (RT) 
              : note the 10e-3 converts [mV] to [V]
              :    nu  = z V F / (RT)
              :
              :    enu = exp(nu)
              :        = exp(z V F / (RT))
              :
              :    fnu = nu/(enu-1) 
              :        = (z V F / (RT)) / (exp(z V F / (RT))-1)
              :        = (z V F / (RT))   (exp(-zV F / (RT))/(1-exp(-zV F / (RT))))
              :
              : now the effective drive is calculated as
              :
              :   ghkg = -f (1 - (ci/co)  enu) fnu
              :        = -10^3 RT/(zF)  (1 - (ci/co) exp(z V F / (RT))) *
              :         (z V F / (RT)) (exp(-zV F / (RT))/(1-exp(-zV F / (RT))))
              :        = -10^3 V (1/co) (co - ci exp(z V F / (RT))) (exp(-zV F / (RT))/(1-exp(-zV F / (RT))))
              :        = 10^3 V/co (ci - co exp(-zV F / (RT)))/(1-exp(-zV F / (RT)))
	      :
              : [note, the 10^3 converts back to mV]
              : and you can see this is the ghk equation if the relationship
              : between conductance and permeability is
              :
              :      g = rho z^2 co F^2/RT
              :
              : Then g*ghkg reduces to the GHK current equation
              :
        f   = (1.0e3/z)*R*(celsius+273.15)/FARADAY
        nu  = v/f
	enu = exp(nu)
	if (fabs(nu) < 1e-4) {
                fnu = 1 - nu/2
        }else{
                fnu = nu/(enu - 1) 
        }
        ghkg= -f*(1 - (ci/co)*enu)*fnu
}
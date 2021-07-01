"""
Hybrid of the Gillies & Willshaw (2006) STN neuron model with Na
channel mechanism from Raman & Bean (2001) accomodating resurgent
component of Na current.

@author Lucas Koelman
@date   14-11-2016
@note   must be run from script directory or .hoc files not found

"""

import numpy as np
import matplotlib.pyplot as plt

import neuron
from neuron import h
nrn = neuron
hoc = h

def sigm(x, theta, sigma):
    """ Sigmoid/logistic/smoothy heavyside function
    @param x    x-value or x-axis vector
    @param theta    midpoint of sigmoid
    @param sigma    slope of sigmoid
    """
    return 1./(1.+np.exp(-(x-theta)/sigma))

def sigmshift(x, y0, y1, theta, sigma):
    """ Vertically shifted sigmoid/logistic function
    @param x    x-value or x-axis vector
    @param y0   y for x << theta
    @param y1   y for x >> theta
    @param theta    midpoint of sigmoid
    @param sigma    slope of sigmoid
    """
    return y0 + y1/(1.0+np.exp(-(x-theta)/sigma))

def bellshift(x, y0, y1, theta_a, theta_b, sigma_a, sigma_b):
    """" Vertically shifted bell curve
    @param y0       y outside of bell region
    @param y1       y inside bell region
    @param theta_a  midpoint of first sigmoid
    @param sigma_a  slope of first sigmoid
    @param theta_b  midpoint of second sigmoid
    @param sigma_b  slope of second sigmoid (opposite sign of sigma_a)
    """
    return y0 + y1/(np.exp(-(x-theta_a)/sigma_a) + np.exp(-(x-theta_b)/sigma_b))

def rates_to_steadystates(v, alphafun, betafun):
    """ Steady state activation and time constant from rate functions """
    alpha = alphafun(v)
    beta = betafun(v)
    tau_inf = 1.0 / (alpha + beta)
    act_inf = alpha / (alpha + beta)
    return act_inf, tau_inf

def vtrap(x,y):
    """ Traps for 0 in denominator of rate eqns. """
    res = x/(np.exp(x/y) - 1)
    divzero = np.abs(x/y) < 1e-6
    res[divzero] = y*(1 - x[divzero]/y/2)
    return res

def plot_Na_vars_Gillies(celsius=25):
    """ Plot sodium channels state variables in Gillies & Willshaw model """
    celsius = 25
    v = np.arange(-100, 100, 0.1)

    ### Gillies & Willshaw ###
    # v-dependence of rate coefficients
    tempb = 23.0
    rest = -60.
    Q10 = 1.980105147
    rate_k = Q10**((celsius-tempb)/10.)
    vadj = v - rest
    # Na activation
    alpham = rate_k * 0.2 * vtrap((13.1-vadj),4.0)
    betam =  rate_k * 0.175 * vtrap((vadj-40.1),1.0)
    # Na inactivation
    alphah = rate_k * 0.08 * np.exp((17.0-vadj)/18.0)
    betah = rate_k * 2.5 / (np.exp((40.0-vadj)/5.0) + 1)

    # plot Gillies & Willshaw values
    plt.figure()
    plt.suptitle('Gillies & Willshaw transient Na rate coefficients')

    plt.subplot(2,1,1)
    plt.plot(v, alpham, label='alpha_m')
    plt.plot(v, betam, label='beta_m')
    plt.xlabel('V (mV)')
    plt.ylabel('(de)activation')
    plt.legend()

    plt.subplot(2,1,2)
    plt.plot(v, alphah, label='alpha_h')
    plt.plot(v, betah, label='beta_h')
    plt.xlabel('V (mV)')
    plt.ylabel('de(inactivation)')
    plt.legend()

def plot_Na_vars_MSM(voltage, celsius=22, version='KuoBean1994_transient'):
    """ Plot parameters of the Markov state model for the Na channel that was used
        and extended in several publication

    NOTE: parameter names in code are based on the diagram in Raman & Bean (2001).

    Different models of the Na channel that employ the same Markov
    State Model with different parameters are:

    ############################################################################
    # Kuo & Bean (1994) #

    - specify 'KuoBean1994_transient'

    - Hypocampal CA1 neurons
        - where recovery from inactivation is not accompanied by ionic current

    - formulates new Markov model that captures interaction between activation
      and inactivation variables

    - Model captures the newly discovered dynamics of development of and 
      _recovery_ from inactivation at hyperpolarized levels

    ############################################################################
    # Taddese & Bean (2002) #

    'TaddeseBean2002_transient_persistent'

    - Tuberomammillary neurons (a form of central neurons)

    - They show that the same V-dependent TTX-sensitive Na channels that generate
      the transient current also generate a persistent current that is sufficient
      to drive spontaneous firing

    - they adapt the Kuo & Bean (1994) Markov model parameters so that the channel
      also produces a persistent sodium current

    ############################################################################
    # Raman & bean (2001) #

    - specify 'RamanBean2001_transient_resurgent'

    - Purkinje neurons
        - where recovery fron inactivation is accompanied by sizeable ionic current

    - they add a new inactivation mechanism that accounts for this ionic current
      ('resurgent current') by adding extra state to Markov model (and changing parameters)
        - a different type of NaV channel that open opens transiently during recovery from inactivation

    - their model is essentially a model of the Na channel with two distinct inactivation
      mechanisms (one additional shaker/ball-and-chain blocking particle that allow for
      current to flow during unblocking)

    - figure 8 A-B show that this channel model accounts for both the transient and 
      resurgent current

    ############################################################################
    ### Do & Bean (2003) ###

    - provides NO MODEL, only experiments

    - see summary in PDF bookmarks

    - both resurgent + persistent component of Na current important for pacemaking
        - persistent drives spontaneous firing
        - resurgent promotes rapid firing by flowing immediately after spikes
          (i.e. a 're-surge' in Na current after the AP depolarizing surge)

    - all three components of Na current are strongly regulated by process of slow inactivation
        - 'slow inactivaton' = upon maintainted or repeated depolarization, a fraction of 
           Na channels enter a strong inactivation state from which recovery is slower than 
           for normal fast inactivation
        - performs role of spike-rate adaptation: promotes a constant frequency

    ############################################################################
    ### Khaliq, Raman (2003) ###

    - same Markov model as Raman & Bean (2001) except change to factor 'a' and 'b'

    - eliminated resurgent current by making mice that lacked gene expressing Nav1.6 channels

    - show that this slows high frequency firing

    - speculate that passing through blocked state facilitates fast recovery from inactivation
      (as opposed to regular inactivated state)
    """

    # v-dependence of rate coefficients
    q10 = 3. # Q10 for rates
    qt = q10**((celsius-22.)/10.) # Room temperature defined as 22 degrees Celsius

    if 'KuoBean1994_transient':
        Con = 0.004
        Coff = 4.5
        Oon = 4.0
        Ooff = 0.008

    elif 'TaddeseBean2002_transient_persistent':
        Con = 0.004 #       (/ms)                   : closed -> inactivated transitions
        Coff = 9.0  #       (/ms)                   : inactivated -> closed transitions
        Oon = 0.8  #   (/ms)                   : open -> Ineg transition
        Ooff = 0.004 #      (/ms)                   : Ineg -> open transition

    elif ('RamanBean2001_transient_resurgent' 
          or 'KhaliqRaman2003_transient_resurgent'):
        Con = 0.005 #       (/ms)                   : closed -> inactivated transitions
        Coff = 0.5  #       (/ms)                   : inactivated -> closed transitions
        Oon = 0.75  #   (/ms)                   : open -> Ineg transition
        Ooff = 0.005 #      (/ms)                   : Ineg -> open transition
    
    elif 'KhaliqRaman2003_transient_only':
        Con = 0.005
        Coff = 0.5
        Oon = 2.3 # NOTE: increased from 0.75
        Ooff = 0.005

    Con_func = lambda v: Con * qt
    Coff_func = lambda v: Coff * qt
    Oon_func = lambda v: Oon * qt
    Ooff_func = lambda v: Ooff * qt
    
    if 'KuoBean1994_transient':
        # Values at room temperature (22 degrees Celsius)
        alpha = 20.
        x1 = 40.
        alpha_func = lambda v: alpha * np.exp(v/x1) * qt
        beta = 3.
        x2 = -40.
        beta_func = lambda v: beta * np.exp(v/x2) * qt
        gamma = 50.
        x3 = 100.
        gamma_func = lambda v: gamma * np.exp(v/x3) * qt
        delta = 0.2
        x4 = -25.
        delta_func = lambda v: delta * np.exp(v/x4) * qt

    elif 'TaddeseBean2002_transient_persistent':
        alpha = 140. #       (/ms)                   : activation
        beta = 10.    #       (/ms)                  : deactivation
        x1 = 27. #           (mV)                    : Vdep of activation (alpha)
        x2 = -27. #          (mV)                    : Vdep of deactivation (beta)
        alpha_func = lambda v: alpha * np.exp(v/x1) * qt
        beta_func = lambda v: beta * np.exp(v/x2) * qt
        gamma = 150. #       (/ms)                   : opening
        delta = 40.  #       (/ms)                   : closing, greater than BEAN/KUO = 0.2
        gamma_func = lambda v: gamma * qt # no voltage dependence
        delta_func = lambda v: delta * qt # no voltage demendence

    elif ('RamanBean2001_transient_resurgent' 
          or 'KhaliqRaman2003_transient_resurgent'
          or 'KhaliqRaman2003_transient_only'):
        alpha = 150. #       (/ms)                   : activation
        x1 = 20. #           (mV)                    : Vdep of activation (alpha)
        alpha_func = lambda v: alpha * np.exp(v/x1) * qt
        beta = 3.    #       (/ms)                   : deactivation
        x2 = -20. #          (mV)                    : Vdep of deactivation (beta)
        beta_func = lambda v: beta * np.exp(v/x2) * qt
        gamma = 150. #       (/ms)                   : opening
        x3 = 1e12 # NOTE: this make exp(V) = 1 (no voltage dependence)
        gamma_func = lambda v: gamma * np.exp(v/x3) * qt
        delta = 40.  #       (/ms)                   : closing, greater than BEAN/KUO = 0.2
        x4 = -1e12 # NOTE: this make exp(V) = 1 (no voltage dependece)
        delta_func = lambda v: delta * np.exp(v/x4) * qt

    


    # Parameters for extra 'Blocked' state for resurgent component
    if 'KuoBean1994' or 'TaddeseBean2002_transient_persistent':
        # RESURGENT CURRENT NOT MODELED
        epsilon = 1e-12
        zeta = 0.03
        epsilon_func = lambda v: epsilon * qt # no voltage dependence
        zeta_func = lambda v: zeta * qt
    
    elif ('RamanBean2001_transient_resurgent' 
          or 'KhaliqRaman2003_transient_resurgent'):
        epsilon = 1.75 #    (/ms)                   : open -> blocked for tau = 0.3 ms at +30 with x5
        x5 = 1e12 #         (mV)                     : Vdep into Ipos (epsilon)
        epsilon_func = lambda v: epsilon * np.exp(v/x5) * qt
        zeta = 0.03 #       (/ms)                   : blocked -> open for tau = 25 ms at -30 with x6
        x6 = -25.#          (mV)                    : Vdep out of Ipos (zeta)
        zeta_func = lambda v: zeta * np.exp(v/x6) * qt

    elif 'KhaliqRaman2003_transient_only':
        epsilon = 1e12 # NOTE: changed to eliminate entry into the blocked state
        x5 = 1e12
        epsilon_func = lambda v: epsilon * np.exp(v/x5) * qt
        zeta = 0.03
        x6 = -25.
        zeta_func = lambda v: zeta * np.exp(v/x6) * qt
   
        
    
    # SET FINAL TRANSITION RATES
    v = voltage
    a_factor = ((Coff_func(v)/Con_func(v))/(Ooff_func(v)/Oon_func(v)))**(1./8.) # In Kuo & Bean
    if 'KuoBean1994_transient':
        alfac = a_factor
        btfac = 1./a_factor
    elif 'RamanBean2001_transient_resurgent':
        alfac = a_factor
        btfac = 1./a_factor
    elif 'TaddeseBean2002_transient_persistent':
        Con1 = 0.004 * qt # inactivation rate at tightest level (1)
        Coff1 = 9.0 * qt # recovery rate at tightest level (1)
        Con5 = 0.05 * qt # inactivation rate at weakest bound level (5)
        Coff5 = 0.00025 * qt # recovery rate at weakest bound level (5)
        alfac = (Con5/Con1)**(1./4.)
        btfac = 1/((Coff1/Coff5)**(1./4.))
        # alfac = 1.88 # at room temperature
        # btfac = 1./13.8 # at room temperature
    elif 'KhaliqRaman2003_transient_resurgent' or 'KhaliqRaman2003_transient_only':
        # change a and b factor with respect to Raman & Bean (2001)
        alfac = (Oon_func(v)/Con_func(v))**(1./4.) # factor for alpha ('a' in Raman & bean paper)
        btfac = (Ooff_func(v)/Coff_func(v))**(1./4.) # factor for beta ('b' in Raman & Bean paper)

    # Create transition matrix and indices
    nstates = 13
    rates = np.zeros((nstates, nstates))
    O = 0
    C1 = 1
    C2 = 2
    C3 = 3
    C4 = 4
    C5 = 5
    I1 = 6
    I2 = 7
    I3 = 8
    I4 = 9
    I5 = 10
    I6 = 11
    B = 12

    # Closed -> Closed (i+1) (weaker binding)
    rates[C1,C2] = 4 * alpha_func(v)
    rates[C2,C3] = 3 * alpha_func(v)
    rates[C3,C4] = 2 * alpha_func(v)
    rates[C4,C5] = 1 * alpha_func(v)

    # Closed -> Closed (i-1) (tighter binding)
    rates[C2,C1] = 1 * beta_func(v)
    rates[C3,C2] = 2 * beta_func(v)
    rates[C4,C3] = 3 * beta_func(v)
    rates[C5,C4] = 4 * beta_func(v)

    # Inactivated -> Inactivated (i+1) (weaker binding)
    rates[I1,I2] = 4 * alfac * alpha_func(v)
    rates[I2,I3] = 3 * alfac * alpha_func(v)
    rates[I3,I4] = 2 * alfac * alpha_func(v)
    rates[I4,I5] = 1 * alfac * alpha_func(v)
    rates[I5,I6] = gamma_func(v)

    # Inactivated -> Inactivated (i-1) (tighter binding)
    rates[I6,I5] = delta_func(v)
    rates[I5,I4] = 4 * btfac * beta_func(v)
    rates[I4,I3] = 3 * btfac * beta_func(v)
    rates[I3,I2] = 2 * btfac * beta_func(v)
    rates[I2,I1] = 1 * btfac * beta_func(v)

    # Open <-> Closed
    rates[C5,O] = gamma_func(v)
    rates[O,C5] = delta_func(v)

    # Open <-> Blocked
    rates[O,B] = epsilon_func(v)
    rates[B,O] = zeta_func(v)

    # Open <-> Inactivated
    rates[O,I6] = Oon_func(v)
    rates[I6,O] = Ooff_func(v)

    # Closed -> Inactivated states
    rates[C1,I1] = Con_func(v)
    rates[C2,I2] = alfac * Con_func(v)
    rates[C3,I3] = alfac**2 * Con_func(v)
    rates[C4,I4] = alfac**3 * Con_func(v)
    rates[C5,I5] = alfac**4 * Con_func(v)

    # Inactivated -> Closed states
    rates[I1,C1] = Coff_func(v)
    rates[I2,C2] = btfac * Coff_func(v)
    rates[I3,C3] = btfac**2 * Coff_func(v)
    rates[I4,C4] = btfac**3 * Coff_func(v)
    rates[I5,C5] = btfac**4 * Coff_func(v)


    # Calculate equilibrium fractions
    sinks = [C1, C2, C3, C4, C5, I1, I2, I3, I4, I5, I6, B] # order of columns and rows (all states minus one)
    states = sinks + [O] # all states
    for source in states:
        # Store the state's total out flux under its own index
        rates[source, source] = -sum([rates[source, sink] for sink in states])

    # Row i has entries X{i,j} representing the fluxes into sink i from source j
    sinks_dCdt = [[rates[source, sink] for source in states] for sink in sinks] # one row for all states minus one (sama as differential equations)
    rhs_sinks_dCdt = [0.0 for sink in sinks]

    mass_conservation = [[1.0 for source in states]]
    rhs_mass_conservation = [1.0]

    flux_mat = np.array(sinks_dCdt + mass_conservation)
    rhs_vec = np.array(rhs_sinks_dCdt + rhs_mass_conservation)

    # Calculate quilibrium values for state variables
    eq_fracs = np.linalg.solve(flux_mat, rhs_vec)

    # Report results
    closed_states = [C1, C2, C3, C4, C5]
    inact_states = [I1, I2, I3, I4, I5, I6]
    C_id = [states.index(Ci) for Ci in closed_states]
    I_id = [states.index(Ii) for Ii in inact_states]
    Ctot = sum([eq_fracs[j] for j in C_id])
    Itot = sum([eq_fracs[j] for j in I_id])
    Otot = eq_fracs[states.index(O)]

    print("Open fraction is %f" % Otot)
    print("Closed fraction is %f" % Ctot)
    print("Inactivated fraction is %f" % Itot)


    # Plot voltage dependence of rates
    # plt.figure()
    # plt.suptitle('Raman & Bean Na channel rate coefficients')

    # plt.subplot(3,1,1)
    # plt.plot(v, f0O, label='gamma')
    # plt.plot(v, b0O, label='delta')
    # plt.xlabel('V (mV)')
    # plt.ylabel('closing')
    # plt.legend()

    # plt.subplot(3,1,2)
    # plt.plot(v, fin_arr, label='Oon')
    # plt.plot(v, bin_arr, label='Ooff')
    # plt.xlabel('V (mV)')
    # plt.ylabel('inactivation')
    # plt.legend()

    # plt.subplot(3,1,3)
    # plt.plot(v, fip, label='epsilon')
    # plt.plot(v, bip, label='zeta')
    # plt.xlabel('V (mV)')
    # plt.ylabel('inactivation')
    # plt.legend()

    # plt.show(block=False)


if __name__ == '__main__':
    for v in [-80., -60., -40., -20., -10., 0., 10., 20.]:
        print "\nVoltage: {} mV".format(v)
        plot_Na_vars_MSM(v, version='RamanBean2001_transient_resurgent')
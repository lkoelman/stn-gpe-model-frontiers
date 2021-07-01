import numpy as np
from scipy.integrate import odeint
from pylab import *

#################################################
# PARAMETERS
#################################################

# Physical constants
Z = 2 # valence of Ca ion
R = 8.31441
FARADAY = 96485.33289

# Reversal potentials (fixed)
ena = 60
ek = -90
el  = -60

# Leakage current
gl  = 0.35

# Ca dynamics
kca   = 2.
alphaca = 0.000005182136

############ Fast Na channel ############
gnabar   = 49
# gating variable m (eq. dm/dt has params m_inf and tau_m of sigmoidal form)
theta_m = -40
k_m = -8.
tau_m0 = 0.2
tau_m1 = 3.
tht_m = -53.
sig_m = -0.7
# gating variable h (param tau_h in dh/dt has equation of bell-shaped form)
theta_h = -45.5
k_h = 6.4
tau_h0 = 0.
tau_h1 = 24.5
sig_h1 = -15.
sig_h2 = 16.
tht_h1 = -50.
tht_h2 = -50.

############ Delayed rectifier K ############
gkdrbar  = 57
# gating variable m (dm/dt has params n_inf and tau_n of bell-shaped form)
theta_n = -41.
k_n = -14.
tau_n0 = 0.
tau_n1 = 11.
tht_n1 = -40.
tht_n2 = -40.
sig_n1 = -40.
sig_n2 = 50.

############ T-type Ca current ############
gcatbar   = 5
# gating variable p (dp/dt has params q_inf and tau_p of bell-shaped form)
theta_p = -56.
k_p = -6.7
tau_p0 = 5.
tau_p1 = 0.33
tht_p1 = -27.
tht_p2 = -102.
sig_p1 = -10.
sig_p2 = 15.
# gating variable q (dq/dt has params q_inf and tau_q of bell-shaped form)
theta_q = -85.
k_q = 5.8
tau_q0 = 0.
tau_q1 = 400.
tht_q1 = -50.
tht_q2 = -50.
sig_q1 = -15.
sig_q2 = 16.

############ L-type Ca current ############
gcalbar   = 15
# gating variable c (dc/dt has params c_inf and tau_c of bell-shaped form)
theta_c = -30.6
k_c = -5.
tau_c0 = 45.
tau_c1 = 10.
tht_c1 = -27.
tht_c2 = -50.
sig_c1 = -20.
sig_c2 = 15.
# gating variable d1 (d1/dt has params d1_inf and tau_d1 of bell-shaped form)
theta_d1 = -60.
k_d1 = 7.5
tau_d10 = 400.
tau_d11 = 500.
tht_d11 = -40.
tht_d12 = -20.
sig_d11 = -15.
sig_d12 = 20.
# gating variable d2 (d2/dt has params d2_inf and tau_d2)
#   - d2_inf is of same form as others but [Ca]_i dependent instead of voltage-dependent
#   - tau_d2 is fixed (not voltage dependent)
theta_d2 = 0.1e-3
k_d2 = 0.02e-3
tau_d2 = 130.

############ A-type K current ############
gkabar  = 5
# ating variable a (da/dt has params a_inf and tau_a of sigmoidal form)  
theta_a = -45.
k_a = -14.7  
k_b = 7.5 
tau_a0 = 1.
tau_a1 = 1.
tht_a = -40.
sig_a = -0.5
# gating variable b (db/dt has params b_inf and tau_b of bell-shaped form)
theta_b = -90.
tau_b0 = 0.
tau_b1 = 200.
tht_b1 = -60.
tht_b2 = -40.
sig_b1 = -30.
sig_b2 = 10.

############ AHP current (Ca dependent K current) ############
gkcabar = 1
# gating variable 4 (dr/dt has fixed tau_r and [Ca]_i-dependent r_inf)
theta_r = 0.17e-3
k_r = -0.08e-3
tau_r = 2.
power_r = 2.

#################################################
# RATE FUNCTIONS
#################################################

def xinf(v, theta_x, k_x):
    """ compute voltage-dependent steady state value """
    return 1.0/(1.0+np.exp((v-theta_x)/k_x))

def taux(v, tau_x0, tau_x1, theta_taux, sigma_taux):
    """ compute voltage-dependent time constant """
    return tau_x0 + tau_x1/(1.0+np.exp(-(v-theta_taux)/sigma_taux))

def taux_bi(v, tau_x0, tau_x1, tht_x1, tht_x2, sig_x1, sig_x2):
    """ compute voltage-dependent time constant """
    return tau_x0 + tau_x1/(np.exp(-(v-tht_x1)/sig_x1) + np.exp(-(v-tht_x2)/sig_x2))

#################################################
# SIMULATION
#################################################

# simulation parameters
cm = 1.0
dt = 0.025
dur = 1000
t = np.arange(0, dur, dt)

# stimulation
bias_amp = 4.0 # 0.0 = spontaneous firing
stim = np.zeros_like(t) + bias_amp

# initial state (initial conditions)
celsius = 30
cao = 2
v0 = -58
cai0 = 5e-6

# Initial condictions: steady state values at V_init
m0 = xinf(v0, theta_m, k_m)
h0 = xinf(v0, theta_h, k_h)
n0 = xinf(v0, theta_n, k_n)
a0 = xinf(v0, theta_a, k_a)
b0 = xinf(v0, theta_b, k_b)
c0 = xinf(v0, theta_c, k_c)
d10 = xinf(v0, theta_d1, k_d1)
d20 = xinf(cai0, theta_d2, k_d2)
p0 = xinf(v0, theta_p, k_p)
q0 = xinf(v0, theta_q, k_q)
r0 = xinf(cai0, theta_r, k_r)

def initialconds_kang():
    """ Initial conditions used in Kang paper """
    cai0 = 0.00001505
    m0 = 0.094693
    h0 = 0.876822
    n0 = 0.228174
    a0 = 0.291417
    b0 = 0.013946
    c0 = 0.004102
    d10 = 0.435728
    d20 = 1.
    p0 = 0.423694
    q0 = 0.009521
    r0 = 0.000001
    globals().update(locals()) # hack: update global variables with local variables

# initialconds_kang()
state0 = [v0, cai0, m0, h0, n0, a0, b0, c0, d10, d20, p0, q0, r0]

def stn_neuron(state, t):
    """ Compute derivatives at local point <state> in phase space.
        NOTE: params not necessary because captured from module scope
    """
    # state/phase variables
    v = state[0]
    cai = state[1]
    m = state[2]
    h = state[3]
    n = state[4]
    a = state[5]
    b = state[6]
    c = state[7]
    d1 = state[8]
    d2 = state[9]
    p = state[10]
    q = state[11]
    r = state[12]
    
    ### compute currents ###
    # Na currents
    ina = gnabar * m*m*m*h * (v - ena)
    # K currents
    ikDR = gkdrbar * n*n*n*n * (v - ek)
    ikA = gkabar * a*a*b * (v - ek)
    ikCA = gkcabar * (v - ek)*np.power(r, power_r)
    ik = ikDR + ikA + ikCA
    # Ca currents
    T = celsius + 273.15
    eca = -1e3*(R*T)/(2*FARADAY)*np.log(cai/cao)
    icaT = gcatbar * p*p*q * (v - eca)
    icaL = gcalbar * c*c*d1*d2 * (v - eca)
    ica = icaT+icaL
    # Nonspecific currents
    ilk = gl * (v - el)

    ### update membrane voltage (derivative) ###
    ibias = 3 # stim[idx]
    dv = -ilk - ina - ik - ica + ibias # RHS divided by C_m = 1.0

    ### compute voltage-dependent rate parameters ###
    # element-wise operations reduce overhead/function-calls
    # h_inf, m_inf, n_inf, p_inf, q_inf, c_inf, d1_inf, a_inf, b_inf, d2_inf, r_inf = xinf(
    #     np.array([v, v, v, v, v, v, v, v, v, cai, cai]),
    #     np.array([theta_h, theta_m, theta_n, theta_p, theta_q, theta_c, theta_d1, theta_a, theta_b, theta_d2, theta_r]),
    #     np.array([k_h, k_m, k_n, k_p, k_q, k_c, k_d1, k_a, k_b, k_d2, k_r]))
    # tau_m, tau_a = taux(
    #     np.array([v, v]), np.array([tau_m0, tau_a0]), np.array([tau_m1, tau_a1]), 
    #     np.array([tht_m, tht_a]), np.array([sig_m, sig_a]))
    # tau_h, tau_n, tau_p, tau_q, tau_c, tau_d1, tau_b = taux_bi(
    #     np.array([v, v, v, v, v, v, v]),
    #     np.array([tau_h0, tau_n0, tau_p0, tau_q0, tau_c0, tau_d10, tau_b0]),
    #     np.array([tau_h1, tau_n1, tau_p1, tau_q1, tau_c1, tau_d11, tau_b1]),
    #     np.array([tht_h1, tht_n1, tht_p1, tht_q1, tht_c1, tht_d11, tht_b1]),
    #     np.array([tht_h2, tht_n2, tht_p2, tht_q2, tht_c2, tht_d12, tht_b2]),
    #     np.array([sig_h1, sig_n1, sig_p1, sig_q1, sig_c1, sig_d11, sig_b1]),
    #     np.array([sig_h2, sig_n2, sig_p2, sig_q2, sig_c2, sig_d12, sig_b2]))
    m_inf = xinf(v, theta_m, k_m)
    h_inf = xinf(v, theta_h, k_h)
    n_inf = xinf(v, theta_n, k_n)
    a_inf = xinf(v, theta_a, k_a)
    b_inf = xinf(v, theta_b, k_b)
    c_inf = xinf(v, theta_c, k_c)
    d1_inf = xinf(v, theta_d1, k_d1)
    d2_inf = xinf(cai, theta_d2, k_d2) # ca-dependent
    p_inf = xinf(v, theta_p, k_p)
    q_inf = xinf(v, theta_q, k_q)
    r_inf = xinf(cai, theta_r, k_r) # ca-dependent

    tau_m = taux(v, tau_m0, tau_m1, tht_m, sig_m)
    tau_h = taux_bi(v, tau_h0, tau_h1, tht_h1, tht_h2, sig_h1, sig_h2)
    tau_n = taux_bi(v, tau_n0, tau_n1, tht_n1, tht_n2, sig_n1, sig_n2)
    tau_a = taux(v, tau_a0, tau_a1, tht_a, sig_a)
    tau_b = taux_bi(v, tau_b0, tau_b1, tht_b1, tht_b2, sig_b1, sig_b2)
    tau_c = taux_bi(v, tau_c0, tau_c1, tht_c1, tht_c2, sig_c1, sig_c2)
    tau_p = taux_bi(v, tau_p0, tau_p1, tht_p1, tht_p2, sig_p1, sig_p2)
    tau_q = taux_bi(v, tau_q0, tau_q1, tht_q1, tht_q2, sig_q1, sig_q2)
    tau_d1 = taux_bi(v, tau_d10, tau_d11, tht_d11, tht_d12, sig_d11, sig_d12)

    ### update gating variables (derivatives) ###
    dh = (h_inf - h)/tau_h
    dm = (m_inf - m)/tau_m
    dn = (n_inf - n)/tau_n
    dp = (p_inf - p)/tau_p
    dq = (q_inf - q)/tau_q
    da = (a_inf - a)/tau_a
    db = (b_inf - b)/tau_b
    dc = (c_inf - c)/tau_c
    dd1 = (d1_inf - d1)/tau_d1
    dd2 = (d2_inf - d2)/tau_d2
    dr = (r_inf - r)/tau_r
    dcai = -alphaca*ica - kca*cai
    return [dv, dcai, dm, dh, dn, da, db, dc, dd1, dd2, dp, dq, dr]

def runsim():
    """ integrate differential equations """
    state = odeint(stn_neuron, state0, t)
    return state

#################################################
# ANALYSIS
#################################################

def zero_crossings(y_axis, window = 11):
    """ Find zero crossings """
    return np.where(np.diff(np.sign(y_axis)))[0]

def plotsimple(state):
    # Generate all signals
    v = state[:,0]
    cai = state[:,1]

    figure()
    interval = (0, 1000)

    subplot(3,1,1)
    # title('No applied current (spontaneous firing)')
    plot(t, v, label="$V_m$")
    xlim(interval[0], interval[1])
    ylim(-100, 100)
    ylabel("$V_m$ (mV)")
    xlabel("time (ms)")
    
    subplot(3,1,2)
    # title('Na current')
    plot(t, stim)
    xlim(interval[0], interval[1])
    ylabel("$I_{app}$ (mA)")
    xlabel("time (ms)")

    subplot(3,1,3)
    # title('K current (Delayed Rectifier)')
    plot(t, cai)
    xlim(interval[0], interval[1])
    ylim(0, 1.1*cai.max())
    ylabel("$[Ca]_i$")
    xlabel("time (ms)")
    
    show()

def test_gating(state):
    """ Plot voltage, currents, and gating variables """
    # Get signals
    v = state[:,0]
    cai = state[:,1]
    m = state[:,2]
    hg = state[:,3]
    n = state[:,4]
    a = state[:,5]
    b = state[:,6]
    c = state[:,7]
    d1 = state[:,8]
    d2 = state[:,9]
    p = state[:,10]
    q = state[:,11]
    r = state[:,12]
    ina = gnabar * m*m*m*hg * (v - ena)
    ikDR = gkdrbar * n*n*n*n * (v - ek)
    ikA = gkabar * a*a*b * (v - ek)
    ikCA = gkcabar * (v - ek)*np.power(r, power_r)
    T = celsius + 273.15
    eca = -1e3*(R*T)/(2*FARADAY)*np.log(cai/cao)
    icaT = gcatbar * p*p*q * (v - eca)
    icaL = gcalbar * c*c*d1*d2 * (v - eca)
    ilk = gl * (v - el)

    # Plot using pylab
    interval = (0, 1000) # Plot interval
    figure()
    subplot(6,2,1)
    plot(t, v, label="$V_m$")
    xlim(interval[0], interval[1])
    ylim(-100, 100)
    ylabel("$V_m$ (mV)")
    
    subplot(6,2,3)
    plot(t, ina)
    xlim(interval[0], interval[1])
    # ylim(-1000, 1000)
    ylabel("$I_{Na}$ (mA)")

    subplot(6,2,5)
    plot(t, ikDR)
    xlim(interval[0], interval[1])
    # ylim(-500, 1500)
    ylabel("$I_{kDR}$ (mA)")

    subplot(6,2,7)
    plot(t, ikA, label='ikA')
    plot(t, ikCA, label='ikCA')
    xlim(interval[0], interval[1])
    # ylim(-0.5, 1)
    legend()
    ylabel("$I_{kA}$,$I_{kCA}$ (mA)")

    subplot(6,2,9)
    plot(t, icaL)
    xlim(interval[0], interval[1])
    # ylim(-20, 10)
    ylabel("$I_{CaL}$ (mA)")

    subplot(6,2,11)
    plot(t, icaT)
    xlim(interval[0], interval[1])
    # ylim(-12, 1)
    ylabel("$I_{CaT}$ (mA)")

    subplot(6,2,2)
    plot(t, m, label='m')
    plot(t, hg, label='h')
    xlim(interval[0], interval[1])
    ylim(-0.1, 1.1)
    legend()
    ylabel("m,h")

    subplot(6,2,4)
    plot(t, n, label='n')
    xlim(interval[0], interval[1])
    ylim(-0.1, 1.1)
    ylabel("n")

    subplot(6,2,6)
    plot(t, a, label='a')
    plot(t, b, label='b')
    xlim(interval[0], interval[1])
    ylim(-0.1, 1.1)
    legend()
    ylabel("a,b")

    subplot(6,2,8)
    plot(t, p, label='p')
    plot(t, q, label='q')
    xlim(interval[0], interval[1])
    ylim(-0.1, 1.1)
    legend()
    ylabel("p,q")

    subplot(6,2,10)
    plot(t, c, label='c')
    plot(t, d1, label='d1')
    plot(t, d2, label='d2')
    xlim(interval[0], interval[1])
    ylim(-0.1, 1.1)
    legend()
    ylabel("c,d1,d2")

    subplot(6,2,12)
    caiv = cai
    caim = caiv.max()
    plot(t, caiv, label='$Ca_i$')
    xlim(interval[0], interval[1])
    ylabel("$[Ca]_i$")

    show()

state = odeint(stn_neuron, state0, t)
plotsimple(state)
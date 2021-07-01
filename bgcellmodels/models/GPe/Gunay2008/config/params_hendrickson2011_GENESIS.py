"""
Parameters from https://senselab.med.yale.edu/ModelDB/showmodel.cshtml?model=127728

Parameters were copied from GENESIS script and put in dictionary format using regex replace: 
    
    r"float\s(\w+)\s+\=\s+(-?[\d.]+)" -> r"\1 : \2"


GENESIS uses SI units
  - Cm : F/m^2
  - Rm : Ohm*m^2
  - Ra : Ohm*m
  - gbar : S/m^2
  - E : V

NEURON uses units:
  - Cm : uF/cm^2  == 1e-6/1e-4 * F/m^2 == 1e-2 * F/m^2    => x 1e2
  - Rm : see gbar                                         => x 1e4
  - Ra : Ohm*cm   == 1e-2 * Ohm*m                         => x 1e2
  - gbar : S/cm^2 == 1/1e-4 * S/m^2 == 1e4 * S/m^2        => x 1e-4
  - E : mV == V*1e-3                                      => x 1e3

"""

## parameters in script /articleCode/commonGPFull/GP1_defaults.g
script_gp_defaults = {
    'RA' : 1.74,     # uniform
    'CM' : 0.024,    # all unmyelinated regions
    'CM_my' : 0.00024,   # myelinated axon segments.
    'RM_sd' : 1.47,  # soma
    'RM_ax' : 1.47,   # unmyelinated axonal regions
    'RM_my' : 10,    # myelinated axon segments.
    'ELEAK_sd' : -0.060,    # soma & dend
    'ELEAK_ax' : -0.060,    # axon
    'EREST_ACT' : -0.060,
}

## parameters in script /articleCode/commonGPFull/simdefaults.g
script_sim_defaults = {
    #Sodium channel kinetics & voltage dependence
    'Vhalfm_NaF' : -0.0324,
    'Km_NaF' : 0.005,
    'taummax_NaF' : 0.000028,
    'taummin_NaF' : 0.000028,

    'V0h_NaF' : -0.048,
    'Kh_NaF' : -0.0028,
    'tauhV0_NaF' : -0.043,
    'tauhmax_NaF' : 0.004,
    'tauhmin_NaF' : 0.00025,   # 0.0002

    'V0s_NaF' : -0.040,
    'Ks_NaF' : -0.0054,
    'mins_NaF' : 0.15,
    'Ktaus1_NaF' : 0.0183,
    'Ktaus2_NaF' : 0.010,
    'tausmax_NaF' : 1,
    'tausmin_NaF' : 0.01,

    'Vhalfm_NaP' : -0.050,
    'V0h_NaP' : -0.057,
    'Kh_NaP' : -0.004,
    'hmin_NaP' : 0.154,
    'V0s_NaP' : -0.01,
    'Abeta_NaP' : 6.94,
    'Bbeta_NaP' : 0.447,

    #Kv2 properties
    'npower_Kv2' : 4,
    'Vhalfn_Kv2' : -0.018,
    'Kn_Kv2' : 0.0091,
    'taunmax_Kv2' : 0.03,
    'taunmin_Kv2' : 0.0001,
    'hmin_Kv2' : 0.2,

    #Kv3 properties
    'npower_Kv3' : 4,
    'Vhalfn_Kv3' : -0.013,    # Actual Vhalf
    'Kn_Kv3' : 0.0078,    # Yields K = 6 mV with Xpower = 4
    'hmin_Kv3' : 0.6,

    #Kv4 properties
    'V0n_Kv4' : -0.049,    # Yields Vhalf = -27 mV when Xpower = 4
    'Kn_Kv4' : 0.0125,    # Yields K = 9.6 mV when Xpower = 4
    'Ktaun1_Kv4' : 0.029,
    'Ktaun2_Kv4' : 0.029,

    'V0h_Kv4' : -0.083,    # changed from -0.072 02/17/2005 to match 
                                        # Tkatch et al
    'Kh_Kv4' : -0.01, # changed from -0.0047 02/17/2005 to match 
                                        # Tkatch et al
    'Ktauh1_Kv4' : 0.010,
    'Ktauh2_Kv4' : 0.010,

    #KCNQ properties
    'Vhalfn_KCNQ' : -0.0285,
    'Kn_KCNQ' : 0.0195,    # Yields K = 15 mV for 1st order Boltzmann
                                    #  when Xpower = 4.

    #SK channel properties
    'EC50_SK' : 0.00035, # SI unit = mM; default = 350 nM.
    'dQ10_SK' : 2,

    #CaHVA properties
    'npower_CaHVA' : 1,
    'Vhalfn_CaHVA' : -0.02,
    'Kn_CaHVA' : 0.007, 

    #Voltage-gated ion channel reversal potentials
    'ENa' : 0.050,
    'ECa' : 0.130,
    'EK' : -0.090,
    'Eh' : -0.03,

    #Calcium concentration parameters
    'B_Ca_GP_conc' : 4.0/3.0*5.2e-12,#3.6e-7 #changed on 10/15/2009 to be consistent with GPcomps.g
    'shell_thick' : 20e-9,  #  meters 
    'tau_CaClearance' : 0.001,   #  time constant for Ca2+ clearance (sec)
}

## parameters in script /articleCode/commonGPFull/actpars.g
script_actpars = {
    'dendNaF' : 40,
    #Voltage-gated ion channel densities
    'G_NaF_soma' : 2500,
    'G_NaP_soma' : 1,
    'G_Kv2_soma' : 320,
    'G_Kv3_soma' : 640,
    'G_Kv4f_soma' : 160,
    'G_Kv4s_soma' : 160*1.5,
    'G_KCNQ_soma' : 0.4,
    'G_SK_soma' : 50,
    'G_Ca_HVA_soma' : 2, 
    'G_h_HCN_soma' : 0.2,
    'G_h_HCN2_soma' : 0.2*2.5,

    'G_NaF_axon' : 5000,
    'G_NaP_axon' : 40,
    'G_Kv2_axon' : 640,
    'G_Kv3_axon' : 1280,  
    'G_Kv4f_axon' : 1600,
    'G_Kv4s_axon' : 1600*1.5  ,
    'G_KCNQ_axon' : 0.4,
    
    'G_NaF_dend' : 40,
    'G_NaP_dend' : 1,
    'G_Kv2_dend' : 64,
    'G_Kv3_dend' : 128,   
    'G_Kv4f_dend' : 160,
    'G_Kv4s_dend' : 160*1.5,
    'G_KCNQ_dend' : 0.4,
    'G_SK_dend' : 4,
    'G_Ca_HVA_dend' : 0.15,  
    'G_h_HCN_dend' : 0.2,
    'G_h_HCN2_dend' : 0.2*2.5,
}

# Calcium buffering parameters #################################################
from math import pi as PI

B_Ca_GP_conc = 4.0/3.0*5.2e-12  # 3.6e-7 #changed on 10/15/2009 to be consistent with GPcomps.g
shell_thick = 20e-9             #  meters 
tau_CaClearance = 0.001         #  time constant for Ca2+ clearance (sec)

### Soma sections
dia = 1.                        # (meters)
rad = dia/2.
rad_core = rad - shell_thick    # shell_thick is in (m) so rad too
surf = 4*PI*rad*rad
vol = 4.0/3.0*PI*rad*rad*rad
core_vol = 4.0/3.0*PI*rad_core*rad_core*rad_core
shell_vol = vol - core_vol      # (m^3)

soma_shell_vol = shell_vol
soma_B = B_Ca_GP_conc / shell_vol

### Axon hillock sections
L = 1.
dia = 1.
rad = dia / 2.
surf = 2*PI*rad*L
vol = PI*rad*rad*L
if dia > (shell_thick*2):
    rad_core = rad - shell_thick
    core_vol = PI*rad_core*rad_core*L
    shell_vol = vol - core_vol
else:
    shell_vol = vol

axon_shell_vol = shell_vol
axon_B = B_Ca_GP_conc / shell_vol

### Dendritic sections
L = 1.
dia = 1.
rad = dia / 2.
surf = 2*PI*rad*L
vol = PI*rad*rad*L
if dia > (shell_thick*2):
    rad_core = rad - shell_thick
    core_vol = PI*rad_core*rad_core*L
    shell_vol = vol - core_vol
else:
    shell_vol = vol

dend_shell_vol = shell_vol
dend_B = B_Ca_GP_conc / shell_vol

ca_buffering_params = {
    'B_Ca_GP_conc': B_Ca_GP_conc,
    'shell_thick': shell_thick,
    'tau_CaClearance': tau_CaClearance,
    'soma_shell_vol': soma_shell_vol,
    'axon_shell_vol': axon_shell_vol,
    'dend_shell_vol': dend_shell_vol,
}

################################################################################

all_params = {}
for param_dict in script_sim_defaults, script_gp_defaults, script_actpars, ca_buffering_params:
    all_params.update(param_dict)


def params_to_json(json_filename):
    """
    Write parameters dict to JSON file.
    """
    import json
    with open(json_filename, 'w') as outfile:
        json.dump(all_params, outfile, indent=4, sort_keys=True)
    print("Wrote parameters to json file {}".format(json_filename))
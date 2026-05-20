
# -------------------------------------------------------------------------
# Functions - from MuSICA code
#

def q_air_eair(e_air, pressure):
    """
       Computes air water vapour mixing ratio (kg kg-1) from
       air water vapour pressure (Pa) and air pressure (Pa).

       INPUT
       - air water vapour pressure (Pa)
       - air pressure (Pa)

       RETURNED VALUE
       - air water vapour mixing ratio (kg kg-1)
    """
    # molar weight of dry air (kg mol-1)
    mol_weight_air = 28.966e-03
    # molar weight of water (kg mol-1)
    mol_weight_h2o = 18.016e-03

    c0 = mol_weight_h2o/mol_weight_air
    c1 = 0.
    q_air = c0*e_air/(pressure-c1*e_air)

    return q_air


def e_air_sat(tk_air):
    """
       Computes saturated air water vapour pressure (Pa) from air temperature (K).
       This empirical formula is accurate at 0.1% between -25 degC and +35 degC
       The main advantage is that it it is invertible.

       INPUT:
       - air temperature (K)

       RETURNED VALUE:
       - saturated air water vapour pressure (Pa)

       REFERENCES
       Empirical formula derived for the AWIPS 4.2 software developped
         as part of the COMET program (Cooperative Program for
         Operational Meteorology, Education and Training).
       Website: http://meted.ucar.edu/awips/validate/
    """
    eta  = 6.113e+2
    alfa = 4.137e+1
    beta = 1.420e+0
    psi  = 3.212e-2
    mu   = 3.186e-4

    tc_air = tk_air - 273.15
    e_sat = eta + (alfa + (beta + (psi + mu*tc_air)*tc_air)*tc_air)*tc_air

    return e_sat


def eair_q_air(q_air, pressure):
    """
       Computes air water vapour pressure (Pa) from
       air water vapour mixing ratio (kg kg-1) and air pressure (Pa).

       INPUT
       - air water vapour mixing ratio (kg kg-1)
       - air pressure (Pa)

       RETURNED VALUE
       - air water vapour pressure (Pa)
    """
    # molar weight of dry air (kg mol-1)
    mol_weight_air = 28.966e-03
    # molar weight of water (kg mol-1)
    mol_weight_h2o = 18.016e-03

    c0 = mol_weight_air/mol_weight_h2o
    q_air_mol = c0*q_air
    e_air = c0 * q_air * pressure

    return e_air

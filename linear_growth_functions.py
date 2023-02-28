# linear_growth_funcitons.py
# Callable functions to calculate the linear growth rate based on expressions
#   in various papers
# All use SI units

import scipy.constants as const
import numpy as np

def H76(k_perp, k_z, L_n, L_T, T_e, nu_e):
    """
    Hudson and Kelley, 1976

    k_perp [m^-1] - perpendicular wave number
    k_z [m^-1] - parallel wave number
    L_n [m] - density gradient scale length
    L_T [m] - temperature gradient scale length
    T_e [K] - electron temperature
    nu_e [s^-1] - total electron collision frequency (nu_ei + nu_en)
    """

    # Magnetic field (B*c) in units of eV/C*km
    B = 3.2e-5  # T

    omega_n = k_perp*(const.k*T_e)/(const.e*B)/L_n
    omega_T = k_perp*(const.k*T_e)/(const.e*B)/L_T

    # Defined in Hudson & Kelley, 1976
    a_e = np.sqrt(2*const.k*T_e/const.m_e)  # Thermal speed
    lambda_e = a_e/nu_e

    # U = omega_l/k_z; U = current driven electron drift??
    omega_l = 0.

    # Equation 6c
    gamma = 1.03/(k_z*lambda_e)*omega_n/(k_z*a_e)*(omega_l - 3/2*omega_T)

    # print('gamma = ', gamma)
    # print('tau = ', 1./gamma)

    return gamma





def K04(ky, T, L_n, L_T):
    """
    Keskinen et al., 2004
    
    ky [m^-1] - perpendicular wave number
    T [K] - electron temperature
    L_n [m] - density gradiet scale length
    L_T [m] - temperature gradient scale length
    """

    B = 3.2e-5      # T

    gamma = ky*const.k*T/(B*const.e)/np.sqrt(L_n*L_T)

    # print('gamma = ', gamma)
    # print('tau = ', 1./gamma)

    return gamma




if __name__=='__main__':

    k_perp = 1/800.  # m^-1
    k_z = 3.21e-6   # m^-1
    L_n = 1000.   # m
    L_T = 1000.   # m
    T_e = 0.1*const.value('electron volt-kelvin relationship')   # K
    nu_e = 8.2  # s^-1

    print(H76(k_perp, k_z, L_n, L_T, T_e, nu_e))
    print(K04(k_perp, T_e, L_n, L_T))

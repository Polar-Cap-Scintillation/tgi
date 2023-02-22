# hudson_linear_growth.py
# Check linear growth rate calculated from the equations adn parameters given
#   in Hudson & Kelley, 1976

import scipy.constants as const
import numpy as np

def hudson_exact():
    # Table 1 parameters
    # # lambda_perp = 40 m
    # k_z = 3.21e-3   # km^-1
    # L_n = 150   # km
    # L_T = 580   # km
    # omega_n = 3.5e-3    # s^-1
    # omega_T = 9.0e-4    # s^-1
    # T_e = 0.1   # eV
    # nu_e = 8.2  # s^-1
    # # gamma = 1.2e-4

    # lambda_perp = 800 m
    k_z = 3.21e-3   # km^-1
    L_n = 150   # km
    L_T = 580   # km
    omega_n = 1.75e-4    # s^-1
    omega_T = 4.5e-5    # s^-1
    T_e = 0.1   # eV
    nu_e = 8.2  # s^-1
    # gamma = 2.4e-7


    # Electron mass in eV*s^2/km^2
    m = const.value('electron mass energy equivalent in MeV')*1e6/const.c**2*(1e3)**2

    # Defined in Hudson & Kelley, 1976
    a_e = np.sqrt(2*T_e/m)  # Thermal speed
    lambda_e = a_e/nu_e

    # U = omega_l/k_z; U = current driven electron drift??
    omega_l = 0.

    # omega_l = gamma/(1.03/(k_z*lambda_e)*omega_n/(k_z*a_e))+3/2*omega_T
    # print(omega_l/k_z)

    # Equation 6c
    gamma = 1.03/(k_z*lambda_e)*omega_n/(k_z*a_e)*(omega_l - 3/2*omega_T)

    print('gamma = ', gamma)
    print('tau = ', 1./gamma)

def hudson_formulation():

    k_perp = 1/0.8  # km^-1
    k_z = 3.21e-3   # km^-1
    L_n = 1.   # km
    L_T = 1.   # km
    T_e = 0.1   # eV
    nu_e = 8.2  # s^-1

    # Magnetic field (B*c) in units of eV/C*km
    B = 3.2e-5*const.value('joule-electron volt relationship')*const.c*1e3
    c = const.c/1e3

    omega_n = k_perp*(c*T_e)/(const.e*B)/L_n
    omega_T = k_perp*(c*T_e)/(const.e*B)/L_T


    # Electron mass in eV*s^2/km^2
    m = const.value('electron mass energy equivalent in MeV')*1e6/const.c**2*(1e3)**2

    # Defined in Hudson & Kelley, 1976
    a_e = np.sqrt(2*T_e/m)  # Thermal speed
    lambda_e = a_e/nu_e

    # U = omega_l/k_z; U = current driven electron drift??
    omega_l = 0.

    # Equation 6c
    gamma = 1.03/(k_z*lambda_e)*omega_n/(k_z*a_e)*(omega_l - 3/2*omega_T)

    print('gamma = ', gamma)
    print('tau = ', 1./gamma)

def keskinen_formulation():

    # Units don't work???
    gamma = c*ky*T/(B*n*e)*np.sqrt(dndx*dTdx)

if __name__=='__main__':
    # hudson_exact()
    # hudson_formulation()

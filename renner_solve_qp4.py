# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 12:18:07 2015

@author: Fort
"""

import numpy as np
from cmath import sqrt as csqrt, phase, pi
from scipy import special
from scipy.optimize import newton_krylov

def QPphase(per, ri, del_qp, PH_qp1):
    #per = period in seconds
    #ri = well radius in meters
    #del_qp = P/Q in MPa(m^3/s)
    #PH_qp1 = phase lag in fraction of cycle
    omega = 2*pi/per          # radial freq 
    func = lambda D1: float(PH_qp1- np.angle(csqrt(omega*1j/D1)* \
        special.kv(1, np.sqrt(omega*1j/D1)*ri)/\
        special.kv(0, np.sqrt(omega*1j/D1)*ri))/2/pi) #phase function for optirmizer

    D_init = 1e-6 #Diffusivity initial guess
    D_sol = newton_krylov(func, D_init, method='lgmres', verbose=1)
    #Solve for diffusivity
    eta = np.sqrt(omega*1j/D_sol)
    K_0_ri = special.kv(0, eta*ri)
    K_1_ri = special.kv(1, eta*ri)
    T = 1/del_qp/(2*pi*ri/10000*abs(eta*K_1_ri/K_0_ri))/1e6
    #del_qp = P/Q MPa(m^3/s)
    # amplitude from q 1*pi*T*ri/10000*
    print(per, " ", D_sol, " ",  T, " ", T/D_sol )
    return func(0.1)
       

    #QPphase(per, ri, del_qp, PH_qp1)
if __name__ == "__main__":
    QPphase(7200, 0.0486, 9.52, .037)
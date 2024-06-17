# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 14:19:15 2024

@author: nurekeye
"""
import numpy as np
import matplotlib.pyplot as plt
from sympy import *
init_printing(use_latex=True)

x,r0,rxn,D,sigma = symbols('x r0 rxn D sigma')
rxn=5
D=5.6*10**-4
sigma=4.64
x=6
Omega0 = 1 - (rxn * erfc((r0 - rxn) / (4*D*x)) / r0)

Omega = integrate (Omega0 * exp(-0.5*r0**2/sigma**2) * 4*pi*r0**2 / sqrt(8*pi**3*sigma**6), (r0, rxn, oo) )

Omega0 = 1 - ( rxn * special.erfc( (r0 - rxn) / np.sqrt(4*D*(i-x0)) ) / r0)
# Omega0p = amp_diff * Omega0 * np.exp(- grid_r0**2 / (2*sigma**2)) * 4 * np.pi * grid_r0**2 / np.sqrt(8 * np.pi**3 * sigma**6)
# Omega.append( np.trapz(y = list(Omega0p), x=list(grid_r0)) )


x0=-100
y0=0
amp_gauss=1
width=160
amp_diff=0.008

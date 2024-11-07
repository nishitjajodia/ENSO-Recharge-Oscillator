# -*- coding: utf-8 -*-
"""
Initialize
"""
import numpy as np

#non-dimentionalisong constants
T_nd=7.5 # Kelvin
h_nd=150 # meter
t_nd=2 # months

#Critical Values
mu_c=2/3
w_c=np.sqrt(3/32)
t_c=t_nd*2*np.pi/w_c # time period

#other constants
b0=2.5
gamma=0.75
c=1
r=0.25
alpha=0.125
mu_ann=0.2
#mu_0=0.75
tau=12/2 #dimentionalised
fann=0.02
fran=0.2

#Initial conditions
Tinit=1.125/T_nd # Kelvin
hinit=0/h_nd 
taucor=(1/30)/t_nd
    
    # -*- coding: utf-8 -*-
import numpy as np
from Constants import *
import matplotlib.pyplot as plt
import TimeSchemes as TS

#%%plotting
def RKplots(time,title,nt,dt,T_nd,h_nd,t_nd,E1,E2,mu,mu_0,en,T_purtinit=0,h_purtinit=0,varymu=False,varyE1=False,purturb=False):
    '''
        Parameters
        ----------
        time : 1d array
            time-series.
        title : string
            title of plots.
        nt : int
            number of time-steps.
        dt : float
            size of time-steps.
        T_nd : float
            non-dimentionalising constant fot temperature.
        h_nd : TYPE
            non-dimentionalising constant fot depth.
        t_nd : TYPE
            non-dimentionalising constant fot time.
        E1 : float
           Represents random wind stress forcing added the system.
        E2 : float
           Represents random heating added the system.
        mu : float
            coupling coefficient.
        mu_0 : float
        en : float
           Degree of non-linearity.
        T_purtinit : float, optional
            purturbed Tinit . The default is 0.
        h_purtinit : float, optional
            purturbed hinit. The default is 0.
        varymu : bool, optional
            Decides if mu is varied. The default is False.
        varyE1 : bool, optional
            Decides if E1 is varied. The default is False.
        purturb : bool, optional
            decides if ensemble is created. The default is False.
    
        Returns
        -------
        None.

    ''' 
    # Intialise arrays
    Te_nd= np.zeros(nt)
    hw_nd= np.zeros(nt)
    Te_nd[0]=Tinit
    hw_nd[0]=hinit
    if purturb==True:
        Te_nd[0]=T_purtinit
        hw_nd[0]=h_purtinit
        
    for i in range(1,nt):             
        Te_nd[i], hw_nd[i]= TS.RK4(Te_nd[i-1],hw_nd[i-1],dt,en,E1,E2,mu,mu_0,i,varymu,varyE1)

    Te=Te_nd*T_nd  # dimentionalised
    hw=hw_nd*h_nd  # dimentionalised

    if purturb==False:    
        fig, (ax1, ax2) = plt.subplots(1,2,figsize=(20,4))
        fig.suptitle(title)
        ax1.plot(time,Te,label='$T_E$')
        ax1.plot(time,hw/10,label='$h_w/10$')
        ax1.set_xlabel("Time (Months)")
        ax1.set_ylabel("$T_E$, $h_w/10$")
        ax1.set_title("(a) Time series plots for $T_E$ and $h_w/10$")
        ax1.legend()
        ax2.plot(Te,hw/10,alpha=0.7,color='green',label="Phase Plot")
        ax2.set_xlabel('[$T_E$]')
        ax2.set_ylabel('[$h_w/10]')
        ax2.set_title("(b) $T_E-h_w/10$ Phase plots")
        ax2.scatter(Te[-1],hw[-1]/10,color='red',label='End Point')
        ax2.legend()
        ax1.grid(True,which="Both")
        ax2.grid(True,which="Both")
    
    return Te,hw

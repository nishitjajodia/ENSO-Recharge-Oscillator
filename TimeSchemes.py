# -*- coding: utf-8 -*-
"""
functions
"""
import ENSO as enso
from Constants import *

def RK4(T,h,dt,en,E1,E2,mu,mu_0,i,varymu=False,varyE1=False):    
    '''
    Calculates 
    T_new= Te at next time step 
    and 
    h_new-hw at next time step 
    using RK4

    Parameters
    ----------
    h : float
        Thermocline anomaly.
    T : float
        SST anomaly.
    dt : float
        size of time-step.
    en : float
       Degree of non-linearity.
    E1 : float
       Represents random wind stress forcing added the system.
    E2 : float
       Represents random heating added the system.
    mu : float
        coupling coefficient.
    mu_0 : float
 
    i : int
        index of current time-step.
    varymu : bool, optional
        Decides if mu is varied. The default is False.
    varyE1 : bool, optional
        Decides if E1 is varied. The default is False.

    Returns
    -------
    T_new : float
        SST anomaly at next time step.
    h_new : float
        thermocline depth anomaly at next time step.

    '''
    b=b0*mu
    R=gamma*b -c
    i=i-1
    if varymu==True:
        mu1=mu_0*( 1 + mu_ann*np.cos( (2*np.pi*(dt*i)/tau) - 5*np.pi/6) )
        b1=b0*mu1
        R1=gamma*b1 -c
        
        mu2=mu_0*( 1 + mu_ann*np.cos( (2*np.pi*(dt*(i+0.5))/tau) - 5*np.pi/6) )
        b2=b0*mu2
        R2=gamma*b2 -c
        
        mu4=mu_0*( 1 + mu_ann*np.cos( (2*np.pi*(dt*(i+1))/tau) - 5*np.pi/6) )
        b4=b0*mu4
        R4=gamma*b4 -c
        
    else:
        mu1,mu2,mu4=mu,mu,mu
        R1,R2,R4=R,R,R
        b1,b2,b4=b,b,b
    
    if varyE1==True:
        W1=np.random.uniform(-1,1)
        E11= fann*np.cos(2*np.pi*dt*i/tau) + fran*W1*taucor/dt
        
        W2=np.random.uniform(-1,1)
        E12= fann*np.cos(2*np.pi*dt*(i+0.5)/tau) + fran*W2*taucor/dt
        
        W4=np.random.uniform(-1,1)
        E14= fann*np.cos(2*np.pi*dt*(i+1)/tau) + fran*W4*taucor/dt
    else:
        E11,E12,E14=E1,E1,E1
        
        
    k1=enso.dTE(R1,h,T,en,b1,gamma,E11,E2)
    l1=enso.dhw(r,h,T,alpha,b1,E11)
    
    k2=enso.dTE(R2,h+ (l1*dt)/2,T+ (k1*dt)/2,en,b2,gamma,E12,E2)
    l2=enso.dhw(r,h+ (l1*dt)/2,T+ (k1*dt)/2,alpha,b2,E12)
    
    k3=enso.dTE(R2,h+ (l2*dt)/2,T+ (k2*dt)/2,en,b2,gamma,E12,E2)
    l3=enso.dhw(r,h+ (l2*dt)/2,T+ (k2*dt)/2,alpha,b2,E12)
    
    k4=enso.dTE(R4,(h+ l3*dt),(T+ k3*dt),en,b4,gamma,E14,E2)
    l4=enso.dhw(r,(h+ l3*dt),(T+ k3*dt),alpha,b4,E14)
    
    T_new= T + (dt/6) *(k1 + 2*k2 + 2*k3 + k4)
    h_new= h + (dt/6) *(l1 + 2*l2 + 2*l3 + l4)

    return T_new, h_new
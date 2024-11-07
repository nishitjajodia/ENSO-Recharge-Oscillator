# -*- coding: utf-8 -*-
"""
ENSO
"""
def dTE(R,h,T,en,b,gamma,E1,E2):
    '''
    Calculates rate of change of SST anomaly dTe/dt

    Parameters
    ----------
    R : float
        Collectivelydescribes the Bjerknes positive feedback process
    h : float
        Thermocline anomaly.
    T : float
        SST anomaly.
    en : float
        Degree of non-linearity.
    b : float
        Measure of thermocline slope.
    gamma : float
        Specifies the feedback of the thermocline gradient on the SST gradient..
    E1 : float
        Represents random wind stress forcing added the system.
    E2 : float
        Represents random heating added the system.

    Returns
    -------
    d_dT : float
        Rate of change of SST anomaly dTe/dt.
    '''
    d_dT= R*T + gamma*h - en*((h + b*T)**3)+ gamma*E1 + E2
    return d_dT

def dhw(r,h,T,alpha,b,E1):
    '''
    Rate of change of thermocline depth anomaly dhw/dt

    Parameters
    ----------
    r : float
        represents damping of the upper ocean heat content.
    h : float
        Thermocline anomaly.
    T : float
        SST anomaly
    alpha : float
        relates enhanced easterly wind stress to the recharge of ocean heat content.
    b : float
        Measure of thermocline slope.
E1 : float
    Represents random wind stress forcing added the system.
    
    Returns
    -------
    d_dh : float
        Rate of change of thermocline depth anomaly dhw/dt.

    '''
    d_dh=-r*h -alpha*b*T -alpha*E1
    return d_dh


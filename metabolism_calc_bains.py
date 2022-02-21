#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 19:22:32 2021

@author: sukrit
"""
"""
Purpose of this code is to make the plots to help determine if metabolisms (here methanogenesis, ammoniagenesis) are profitable at runaway conditions (lots of the metabolic product, possibly high temperature)
"""
########################
###Import useful libraries
########################
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import pdb

########################
###Define useful constants, all in CGS (via http://www.astro.wisc.edu/~dolan/constants.html)
########################

#Unit conversions
km2m=1.e3 #1 km in m
km2cm=1.e5 #1 km in cm
cm2km=1.e-5 #1 cm in km
amu2g=1.66054e-24 #1 amu in g
bar2atm=0.9869 #1 bar in atm
Pa2bar=1.e-5 #1 Pascal in bar
bar2Pa=1.e5 #1 bar in Pascal
deg2rad=np.pi/180.
bar2barye=1.e+6 #1 Bar in Barye (the cgs unit of pressure)
barye2bar=1.e-6 #1 Barye in Bar
micron2m=1.e-6 #1 micron in m
micron2cm=1.e-4 #1 micron in cm
metricton2kg=1000. #1 metric ton in kg

#Fundamental constants
c=2.997924e10 #speed of light, cm/s
h=6.6260755e-27 #planck constant, erg/s
k=1.380658e-16 #boltzmann constant, erg/K
sigma=5.67051e-5 #Stefan-Boltzmann constant, erg/(cm^2 K^4 s)
R_earth=6371.*km2m#radius of earth in m
R_sun=69.63e9 #radius of sun in cm
AU=1.496e13#1AU in cm


########################
###Load Amend & Shock 2001 delta-G_f^{T,P_sat}
#"Standard states used in these calculations are as follows: for solids and water, unit activity of the pure compound at any temperature and pressure; for gases, unit fugacity of the pure gas at any temperature and 1 bar; and for aqueous solutes, unit activity of a hypothetical 1 molal solution referenced to in¢nite dilution at any temperature and pressure."
#Units of delta-G are kJ/mol
#Temp units are K (converted)
########################
temp_abscissa=np.array([2.0, 18.0, 25.0, 37.0, 45.0, 55.0, 70.0, 85.0, 100.0])+273.15 #Convert C values in Amend & Shock 2001 to real values. 

deltaG_F_PSat_T_H2O_l=interp.interp1d(temp_abscissa, np.array([-235.64, -236.70, -237.18, -238.04, -238.63, -239.39, -240.57, -241.81,-243.08]), kind='linear') #water, liquid (unit activity) #A&S 2001 Table 4.1

deltaG_F_PSat_T_H2_g=interp.interp1d(temp_abscissa, np.array([2.98, 0.91, 0.0, -1.57, -2.63, -3.96, -5.97, -8.00, -10.05]), kind='linear') #H2, gas #A&S 2001 Table 4.1
deltaG_F_PSat_T_H2_aq=interp.interp1d(temp_abscissa, np.array([18.89, 18.11, 17.72, 16.99, 16.46, 15.76, 14.62, 13.39, 12.08]), kind='linear') #H2, aqueous #A&S 2001 Table 4.1

deltaG_F_PSat_T_N2_g=interp.interp1d(temp_abscissa, np.array([4.38, 1.34, 0.0, -2.31, -3.85, -5.79, -8.72, -11.66, -14.63]), kind='linear') #N2, gas #A&S 2001 Table 5.1
deltaG_F_PSat_T_N2_aq=interp.interp1d(temp_abscissa, np.array([20.15, 18.84, 18.18, 16.98, 16.12, 14.99, 13.18, 11.24, 9.19]), kind='linear') #N2, aqueous #A&S 2001 Table 5.1

deltaG_F_PSat_T_NH3_g=interp.interp1d(temp_abscissa, np.array([-12.06, -15.11, -16.45, -18.77, -20.33, -22.28, -25.23, -28.20, -31.20]), kind='linear') #NH3, gas
deltaG_F_PSat_T_NH3_aq=interp.interp1d(temp_abscissa, np.array([-24.30, -25.96, -26.71, -28.02, -28.92, -30.06, -31.82, -33.63, -35.50]), kind='linear') #NH3, aqueous

deltaG_F_PSat_T_CO2_g=interp.interp1d(temp_abscissa, np.array([-389.48, -392.87, -394.36, -396.93, -398.66, -400.83, -404.10, -407.40, -410.73]), kind='linear') #CO2, gas #A&S 2001 Table 7.1
deltaG_F_PSat_T_CO2_aq=interp.interp1d(temp_abscissa, np.array([-383.51, -385.17, -385.98, -387.44, -388.48, -389.84, -391.99, -394.26, -396.66]), kind='linear') #CO2, aqueous #A&S 2001 Table 7.1

deltaG_F_PSat_T_CH4_g=interp.interp1d(temp_abscissa, np.array([-46.47, -49.42, -50.72, -52.97, -54.47, -56.36, -59.22, -62.10, -65.02]), kind='linear') #CH4, gas #A&S 2001 Table 7.1
deltaG_F_PSat_T_CH4_aq=interp.interp1d(temp_abscissa, np.array([-32.71, -33.87, -34.46, -35.57, -36.38, -37.47, -39.23, -41.12, -43.17]), kind='linear') #CH4, aqueous #A&S 2001 Table 7.1

########################
###Load minimum biomass power production criterion, from Seager et al. 2013, in turn sourced from Tijhuis et al. 1993
########################

def P_me(T):
    E=6.94E1 #kJ mol^-1
    R=8.314E-3 #kJ mol^-1 K^-1
    A=2.2E7 #kJ gram biomass^-1 s^-1
    return A*np.exp(-E/(R*T))

life_threshold_deltaG=-20.0 #kJ mol^-1. If reaction does not generate at least 10 kJ/mol, life not known to be able to use. From Hoehler 2004. -20 kJ/mol is conservative.

########################
###Define functions to calculate delta-G for the various reactions
########################

def deltaG_methanogenesis_modern_Earth(T, conc_CO2, conc_H2, pCH4):
    """
    Equation: CO2 (aq) + 4H2 (aq) --> CH4 (g)+ 2H2O (l)
    Assumes M, M, bar, 1bar total pressure.
    T must be array (even if just 1 entry)
    """
    R=8.314E-3 #kJ mol^-1 K^-1
    
    deltaG_0=deltaG_F_PSat_T_CH4_g(T)+2.0*deltaG_F_PSat_T_H2O_l(T) - (deltaG_F_PSat_T_CO2_aq(T) + 4.0*deltaG_F_PSat_T_H2_aq(T))# Standard free energy
    
    Q=(pCH4*(1.0)**2.0)/(conc_CO2*conc_H2**4.0) #Reaction quotient Q
    
    deltaG=deltaG_0+R*T*np.log(Q)
    
    deltaG[deltaG > life_threshold_deltaG]=np.nan #if at least 10 kJ/mol not generated, bugs not known to live.
    
    return deltaG

def deltaG_methanogenesis_early_Earth(T, pCO2, pH2, pCH4):
    """
    Equation: CO2 (g) + 4H2 (g) --> CH4 (g)+ 2H2O (l)
    Assumes bar, 1bar total pressure, for gases.
    T must be array (even if just 1 entry)
    """
    R=8.314E-3 #kJ mol^-1 K^-1
    
    deltaG_0=deltaG_F_PSat_T_CH4_g(T)+2.0*deltaG_F_PSat_T_H2O_l(T) - (deltaG_F_PSat_T_CO2_g(T) + 4.0*deltaG_F_PSat_T_H2_g(T))# Standard free energy
    
    Q=(pCH4*(1.0)**2.0)/(pCO2*pH2**4.0) #Reaction quotient Q
    
    deltaG=deltaG_0+R*T*np.log(Q)
    
    deltaG[deltaG > life_threshold_deltaG]=np.nan #if at least 10 kJ/mol not generated, bugs not known to live.
    
    return deltaG

def deltaG_ammoniagenesis(T, pN2, pH2, pNH3):
    """
    Equation: N2 (g) + 3H2 (g) --> 2NH3 (g)
    Assumes bar, 1bar total pressure, for gases.
    T must be array (even if just 1 entry)
    """
    R=8.314E-3 #kJ mol^-1 K^-1
    
    deltaG_0=2.0*deltaG_F_PSat_T_NH3_g(T) - (deltaG_F_PSat_T_N2_g(T) + 3.0*deltaG_F_PSat_T_H2_g(T))# Standard free energy
    
    Q=(pNH3**2.0)/(pN2*pH2**3.0) #Reaction quotient Q
    
    deltaG=deltaG_0+R*T*np.log(Q)
    
    deltaG[deltaG > life_threshold_deltaG]=np.nan #if at least 10 kJ/mol not generated, bugs not known to live.

    return deltaG


########################
###Plot methanogenesis
########################
temps=np.linspace(276.0, 373.0,num=50) #Temperature abscissa in K

pCH4=1.6E-2 #bar; peak of runaway

fig, ax=plt.subplots(1, figsize=(8., 6.), sharex=False)
markersizeval=5.
#ax[0].set_title(r'CO, N$_2$-CO$_2$ atm')
ax.plot(temps, P_me(temps)/-deltaG_methanogenesis_modern_Earth(temps, 3.0E-3, 3.0E-9, 1.6E-6), linewidth=2, linestyle='-', color='red', label=r'Modern Earth, pCH$_4=1.6\times10^{-6}$ bar', marker='s', markersize=markersizeval)
ax.plot(temps, P_me(temps)/-deltaG_methanogenesis_modern_Earth(temps, 3.0E-3, 3.0E-9, 2.0E-2), linewidth=2, linestyle='--', color='hotpink', label=r'Modern Earth, pCH$_4=2\times10^{-2}$ bar', marker='s', markersize=markersizeval)
ax.plot(temps, P_me(temps)/-deltaG_methanogenesis_early_Earth(temps, 0.9, 1.0E-3, 6.0E-5), linewidth=2, linestyle='-', color='blue', label=r'CO$_2$-N$_2$ Atm., pCH$_4=6\times10^{-5}$ bar', marker='s', markersize=markersizeval)
ax.plot(temps, P_me(temps)/-deltaG_methanogenesis_early_Earth(temps, 0.9, 1.0E-3, 2.0E-2), linewidth=2, linestyle='--', color='cyan', label=r'CO$_2$-N$_2$ Atm., pCH$_4=2\times10^{-2}$ bar', marker='s', markersize=markersizeval)

ax.axvline(288.0, linestyle='-', color='green', label=r'$T_{surf}=$288K')
ax.axvline(301.0, linestyle='--', color='red', label=r'$T_{surf}=$301K')

ax.set_xscale('linear')
ax.set_yscale('log')
ax.set_ylabel(r'Min. CH$_4$ Prod. Rate (mol g biomass$^{-1}$ s$^{-1}$)', fontsize=14)
ax.set_xlabel(r'Surface Temperature (K)', fontsize=14)  
plt.tick_params(labelsize=14)

ax.legend(loc='best', ncol=1, borderaxespad=0., fontsize=12)   

ax.axhline(1.0E-4, linestyle=':', color='purple')

# ax.set_xlim([3.0E+8, 3.0E+12])
# ax.set_ylim([1.0E-6, 1.0])


plt.savefig('./Plots/plot_methanogenesis.pdf', orientation='portrait',format='pdf')

########################
###Plot "ammoniagenesis"
########################
temps=np.linspace(278.0, 373.0,num=50) #Temperature abscissa in K

pN2=0.1 #bar
pH2=0.9 #bar
pNH3=0.001 #bar; peak of runaway

fig, ax=plt.subplots(1, figsize=(8., 6.), sharex=False)
markersizeval=5.
#ax[0].set_title(r'CO, N$_2$-CO$_2$ atm')
ax.plot(temps, 2.0*P_me(temps)/-deltaG_ammoniagenesis(temps, pN2, pH2, pNH3), linewidth=2, linestyle='-', color='black', marker='s', markersize=markersizeval, label=r'pNH$_3=1\times10^{-3}$ bar') #Factor of 2 accounts for fact that 2 NH3 are produced per N2 consumed. 

ax.set_xscale('linear')
ax.set_yscale('log')
ax.set_ylabel(r'Min. NH$_3$ Prod. Rate (mol g biomass$^{-1}$ s$^{-1}$)', fontsize=14)
ax.set_xlabel(r'Surface Temperature (K)', fontsize=14)  
ax.axvline(288.0, linestyle='-', color='green', label=r'$T_{surf}=$288K')
ax.axvline(299.0, linestyle='--', color='red', label=r'$T_{surf}=$299K')  
ax.legend(loc='best', ncol=1, borderaxespad=0., fontsize=12) 

plt.tick_params(labelsize=14)
ax.axhline(1.0E-4, linestyle=':', color='purple')

# ax.set_xlim([3.0E+8, 3.0E+12])
# ax.set_ylim([1.0E-6, 1.0])


plt.savefig('./Plots/plot_ammoniagenesis.pdf', orientation='portrait',format='pdf')



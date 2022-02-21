"""
Purpose of this code is to read and plot the output of the Hu+2012 photochemical code
"""
########################
###Import useful libraries
########################
import numpy as np
import scipy.optimize
import scipy.integrate
import scipy.interpolate as interp
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
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

hc=1.98644582e-9 #erg nm


########################
###Establish key
########################

#Corrected for Hu 1-indexing vs Python 0-indexing
ind_o=1-1 #O
ind_h=3-1 #H
ind_oh=4-1 #OH
ind_ph3 = 112-1 #PH3
ind_ph3d = 113-1 #PH3D

ind_so2=43-1
ind_so=42-1
ind_ch4=21-1
ind_h2s=45-1
ind_h2=53-1
ind_h2o=7-1

ind_no=12-1
ind_n2o=11-1
ind_nh3=9-1
ind_n2h2=106-1
ind_nh2=10-1

ind_co2=52-1
ind_n2=55-1


ind_s8=79-1
ind_s8a=111-1
ind_ch4o=24-1
ind_c2h2=27-1
ind_c2h4=29-1
ind_ocs=49-1

ind_o1d=56-1
ind_co=20-1
ind_o2=54-1
ind_c2h6=31-1
ind_h2o2=6-1
ind_h2so4=73-1
ind_h2so4a=78-1
ind_ch2o=22-1
ind_o3=2-1
ind_ho2=5-1

def get_conc_std_data(base_file, ind_list):
    """
    Takes:
        -name of concentration_STD file
        -np.array of hu code species, specified by their index as above
    Returns:
        -np.array of surface mixing ratios
        -np.array of column-averaged mixing ratios
    """
    surface_mixing_ratios=np.zeros(np.shape(ind_list))
    column_averaged_mixing_ratios=np.zeros(np.shape(ind_list))

    ########################
    ###Read in base data
    ########################
    base_data=np.genfromtxt(base_file, skip_header=2, skip_footer=0, unpack=False) #Import mapping between numerical ID in code and species name.
    
    z_center_base=base_data[:,0] # Center of altitude bins, km 
    T_z_base=base_data[:,3] # Temperature(z), in K
    P_z_base=base_data[:,4]*Pa2bar*bar2barye # Pressure(z), in Pa converted to Barye
    n_z_s_base=base_data[:,5:] #Number concentrations of the 111 chemical species, in cm**-3, as a function of (altitude, species)

    ###Get molar concentrations
    ###NOTE: May (probably) need to exclude condensed-phase species for molar concentration calculation...probably doesn't matter most of the time, but formally required and mioght matter in some weird edge cases.
    n_z_base=np.sum(n_z_s_base,1) #sum number densities across species. This is a profile for the whole atmosphere.
    n_z_bulkatm_base=P_z_base/(k*T_z_base)
    mc_z_s_base=np.zeros(np.shape(n_z_s_base))
    mr_z_s_base=np.zeros(np.shape(n_z_s_base))

    num_s=np.shape(n_z_s_base)[1]

    for ind2 in range(0, num_s):
        mc_z_s_base[:,ind2]=n_z_s_base[:,ind2]/n_z_base#molar concentration of each species.
        mr_z_s_base[:,ind2]=n_z_s_base[:,ind2]/n_z_bulkatm_base#mixing ratio of each species.


    ########################
    ###Print key parameters
    ######################## 
    for ind3 in range(0, len(ind_list)):
        species_ind=ind_list[ind3]
        surface_mixing_ratios[ind3]=mr_z_s_base[0,species_ind]
        column_averaged_mixing_ratios[ind3]=np.sum(mr_z_s_base[:,species_ind]*n_z_bulkatm_base)/np.sum(n_z_bulkatm_base)
    return surface_mixing_ratios, column_averaged_mixing_ratios

# def return_reaction_rates(concSTD, ChemReac, reaction_label_list):
#     """
#     Takes filenames for a concentrationSTD.dat and ChemicalRate.dat file from the Hu code, and a list of Hu code reaction IDs.
#     From concentrationSTD.dat file extracts z scale. From ChemicalRate file takes reaction rate (cm**-3 s**-1).
#     Returns: z_centers (km), delta_zs (cm), vertically resolved rates for each reaction ID (cm**-3 s**-1), column integrated reaction rate for each reaction ID (cm**-2 s**-1)
#     """
#     ###Get z scale
#     tp=np.genfromtxt(concSTD, skip_header=2, skip_footer=0, unpack=False) #Import simulation results
#     z_centers=tp[:,0] # Center of altitude bins, km 
#     deltazs=(tp[:,2]-tp[:,1]) #width of altitude bins, km    

#     ###Get Chemical Reaction Rates
#     chemlabels=np.genfromtxt(ChemReac, skip_header=1, skip_footer=0, unpack=True, usecols=(0), dtype=None)
#     chemrates=np.genfromtxt(ChemReac, skip_header=1, skip_footer=0, unpack=False, dtype=float)[:,1:]
    
#     ###Extract relevant chemical reaction rates.
    
#     vert_reac_rates={} #Dict which will hold the vertical reaction rates (cm*-3 s**-1)
#     colint_reac_rates={} #Dict which will hold the column-integrated reaction rates (cm**-2 s**-1)
    
#     for target_chemlabel in reaction_label_list: #Step through the chemical reactions we are interested in...
#         for ind in range(0, len(chemlabels)): #Step through each of the chemical reactions in the FILE
#             chemlabel=chemlabels[ind].decode('UTF-8') #decode the indth chemical reaction in the FILE
#             if chemlabel==target_chemlabel: #If they match, this is what we want
#                 vert_reac_rates[target_chemlabel]=chemrates[ind,:] #cm**-3 s**-1
#                 colint_reac_rates[target_chemlabel]=np.sum(chemrates[ind,:]*deltazs*km2cm) #cm**-2 s**-1

#     return z_centers, deltazs, vert_reac_rates, colint_reac_rates

# def get_integrated_uv(file, AU, termination_wavelength):
#     """
#     Takes: 
#         -file given stellar irradiance, normalized to 1 AU equivalent (solar constant)
#         -"Distance" in AU, the way the Hu code does it (so 1.6 stays constant despite moving to M-dwarf)
#         -Termination wavelength, nm
#     """
#     wav, toa=np.genfromtxt(file, skip_header=0, skip_footer=0,usecols=(0,1), unpack=True) #units: nm, W m**-2 nm**-1
#     toa_scaled=toa/AU**2.0 #Correct for planet-star distance. 
#     toa_scaled_units=toa_scaled*1000.0  #change units from W m**-2 nm**-1 to erg s**-1 cm**-2 nm**-1
#     toa_scaled_units_phot=toa_scaled_units*wav/hc #Convert from erg s**-1 cm**-2 nm**-1 to photons s**-1 cm**-2 nm**-1
    
#     toa_func=interp.interp1d(wav, toa_scaled_units_phot, kind='linear')
    
#     integrated_uv=scipy.integrate.quad(toa_func, np.min(wav), termination_wavelength, epsabs=0., epsrel=1.e-3, limit=100000)[0]
    
#     return integrated_uv #photons cm^-2 s^-1. 

def plot_runaway_nh3():
    """
    Plot runaway of NH3 only
    """
    
    ########################    
    ###Process NH3, GJ876, standard (low) deposition
    ########################    
    file_prefix='./scenario_library/GJ876/H2-Full-NH3/NH3_'
    conc_file_list=np.array(['1E8', '3E8','1E9','3E9', '1E10', '2E10', '3E10'])
    file_suffix='/ConcentrationSTD.dat'    
    gj876_nh3_flux_list=np.array([1.0E+8, 3.0E+8, 1.0E+9, 3.0E+9, 1.0E+10, 2.0E+10, 3.0E+10])

    ###   
    #Extract concentration data
    ###
    num_file=len(conc_file_list)
    gj876_nh3_surface_mixing_ratios=np.zeros(np.shape(gj876_nh3_flux_list))
    gj876_nh3_column_averaged_mixing_ratios=np.zeros(np.shape(gj876_nh3_flux_list))

    for ind in range(0, num_file):
        conc_file=file_prefix+conc_file_list[ind]+file_suffix
       
        ###Get NH3 mixing ratios (surface and column-averaged)
        gj876_nh3_surface_mixing_ratios[ind], gj876_nh3_column_averaged_mixing_ratios[ind]=get_conc_std_data(conc_file, np.array([ind_nh3]))
    
    ########################    
    ###Process NH3, Sun, standard (low) deposition
    ########################    
    file_prefix='./scenario_library/Sun/H2-Full-NH3/NH3_'
    conc_file_list=np.array(['1E8', '3E8','1E9','3E9', '1E10', '3E10', '1E11', '3E11'])
    file_suffix='/ConcentrationSTD.dat'    
    sun_nh3_flux_list=np.array([1.0E+8, 3.0E+8, 1.0E+9, 3.0E+9, 1.0E+10, 3.0E+10, 1.0E+11, 3.0E+11])

    ###   
    #Extract concentration data
    ###
    num_file=len(conc_file_list)
    sun_nh3_surface_mixing_ratios=np.zeros(np.shape(sun_nh3_flux_list))
    sun_nh3_column_averaged_mixing_ratios=np.zeros(np.shape(sun_nh3_flux_list))

    for ind in range(0, num_file):
        conc_file=file_prefix+conc_file_list[ind]+file_suffix
       
        ###Get NH3 mixing ratios (surface and column-averaged)
        sun_nh3_surface_mixing_ratios[ind], sun_nh3_column_averaged_mixing_ratios[ind]=get_conc_std_data(conc_file, np.array([ind_nh3]))

    ########################    
    ###Process NH3, GJ876, high deposition
    ########################    
    file_prefix='./scenario_library/GJ876/H2-Full-NH3/NH3_'
    conc_file_list=np.array(['1E8', '3E8','1E9','3E9', '1E10', '2E10', '3E10', '6E10'])
    file_suffix='_highdep/ConcentrationSTD.dat'    
    gj876_nh3_highdep_flux_list=np.array([1.0E+8, 3.0E+8, 1.0E+9, 3.0E+9, 1.0E+10, 2.0E+10, 3.0E+10, 6.0E+10])

    ###   
    #Extract concentration data
    ###
    num_file=len(conc_file_list)
    gj876_nh3_highdep_surface_mixing_ratios=np.zeros(np.shape(gj876_nh3_highdep_flux_list))
    gj876_nh3_highdep_column_averaged_mixing_ratios=np.zeros(np.shape(gj876_nh3_highdep_flux_list))

    for ind in range(0, num_file):
        conc_file=file_prefix+conc_file_list[ind]+file_suffix
       
        ###Get mixing ratios (surface and column-averaged)
        gj876_nh3_highdep_surface_mixing_ratios[ind], gj876_nh3_highdep_column_averaged_mixing_ratios[ind]=get_conc_std_data(conc_file, np.array([ind_nh3]))
        
    ########################    
    ###Process NH3, GJ876, low deposition, Koll TP, 170K stratosphere
    ########################    
    file_prefix='./scenario_library/GJ876/H2-Full-NH3/NH3_'
    conc_file_list=np.array(['1E8', '3E8','1E9','3E9', '6E9', '1E10', '3E10'])
    file_suffix='_KollTP/ConcentrationSTD.dat'    
    gj876_nh3_kolltp_flux_list=np.array([1.0E+8, 3.0E+8, 1.0E+9, 3.0E+9, 6.0E+9,1.0E+10, 3.0E+10])

    ###   
    #Extract concentration data
    ###
    num_file=len(conc_file_list)
    gj876_nh3_kolltp_surface_mixing_ratios=np.zeros(np.shape(gj876_nh3_kolltp_flux_list))
    gj876_nh3_kolltp_column_averaged_mixing_ratios=np.zeros(np.shape(gj876_nh3_kolltp_flux_list))

    for ind in range(0, num_file):
        conc_file=file_prefix+conc_file_list[ind]+file_suffix
       
        ###Get mixing ratios (surface and column-averaged)
        gj876_nh3_kolltp_surface_mixing_ratios[ind], gj876_nh3_kolltp_column_averaged_mixing_ratios[ind]=get_conc_std_data(conc_file, np.array([ind_nh3]))
        
    ########################    
    ###Process NH3, GJ876, low deposition, Koll TP, 210K stratosphere
    ########################    
    file_prefix='./scenario_library/GJ876/H2-Full-NH3/NH3_'
    conc_file_list=np.array(['1E8', '3E8','1E9','3E9', '1E10', '2E10', '4E10', '6E10'])
    file_suffix='_KollTP_hotstrat/ConcentrationSTD.dat'    
    gj876_nh3_kolltp_hotstrat_flux_list=np.array([1.0E+8, 3.0E+8, 1.0E+9, 3.0E+9, 1.0E+10, 2.0E+10, 4.0E+10, 6.0E+10])

    ###   
    #Extract concentration data
    ###
    num_file=len(conc_file_list)
    gj876_nh3_kolltp_hotstrat_surface_mixing_ratios=np.zeros(np.shape(gj876_nh3_kolltp_hotstrat_flux_list))
    gj876_nh3_kolltp_hotstrat_column_averaged_mixing_ratios=np.zeros(np.shape(gj876_nh3_kolltp_hotstrat_flux_list))

    for ind in range(0, num_file):
        conc_file=file_prefix+conc_file_list[ind]+file_suffix
       
        ###Get mixing ratios (surface and column-averaged)
        gj876_nh3_kolltp_hotstrat_surface_mixing_ratios[ind], gj876_nh3_kolltp_hotstrat_column_averaged_mixing_ratios[ind]=get_conc_std_data(conc_file, np.array([ind_nh3]))
        
    ########################
    ###Plot
    ########################
    fig, ax=plt.subplots(1, figsize=(8., 8.), sharex=False)
    markersizeval=5.
    
    ax.plot(gj876_nh3_flux_list, gj876_nh3_column_averaged_mixing_ratios,linewidth=3, linestyle='-', color='red', label='GJ876', marker='s', markersize=markersizeval)
    ax.plot(gj876_nh3_highdep_flux_list, gj876_nh3_highdep_column_averaged_mixing_ratios,linewidth=1.5, linestyle='--', color='red', label='GJ876, High Dep. BCs.', marker='s', markersize=markersizeval)
    ax.plot(gj876_nh3_kolltp_flux_list, gj876_nh3_kolltp_column_averaged_mixing_ratios,linewidth=1.5, linestyle='--', color='magenta', label=r'GJ876, $T_{surf}$=299K', marker='s', markersize=markersizeval)
    ax.plot(gj876_nh3_kolltp_hotstrat_flux_list, gj876_nh3_kolltp_hotstrat_column_averaged_mixing_ratios,linewidth=1.5, linestyle='--', color='hotpink', label=r'GJ876, $T_{surf}$=299K , $T_{strat}$=210K', marker='s', markersize=markersizeval)

    ax.plot(sun_nh3_flux_list, sun_nh3_column_averaged_mixing_ratios,linewidth=1.5, linestyle='-', color='gold', label='Sun', marker='s', markersize=markersizeval)

    ax.axvline(4.0E+10, color='green', linestyle=':', linewidth=4, label=r'Modern Earth Net O$_2$ Production') #taken from Zahnle+2006
    
    ax.axvspan(0.75*1.44E10, 1.25*1.44E10, color='purple', alpha=0.2, label=r'Modern  NH$_3$ Flux to Atmosphere') #derived from estimate of Bouwman+1997 of 54 Tg/year, including 25% uncertainties.

    ax.axvspan(0.75*1.44E10/5, 1.25*1.44E10/2, color='blue', alpha=0.2, label=r'Pre-Industrial NH$_3$ Flux to Atmosphere') #derived from estimate of Bouwman+1997 of 54 Tg/year, including 25% uncertainties.

    abscissa=np.logspace(6,13)
    ax.plot(abscissa, abscissa/(2.5E19*1.0E-7), linestyle='--', color='black', label=r'v$_{dep}$-limited')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel(r'Column-Averaged Mixing Ratio', fontsize=14)
    ax.set_xlabel(r'Surface Flux (cm$^{-2}$ s$^{-1}$)', fontsize=14)  
    ax.legend(loc='upper center', ncol=2, borderaxespad=0., fontsize=12, bbox_to_anchor=(0.5,1.29))   

    ax.set_xlim([1.0E+8, 1.0E+11])
    ax.set_ylim([1.0E-10, 1.0E-2])
    plt.tick_params(labelsize=14)
    fig.tight_layout()
#
#
##    ax[1].axvspan(approximation_nh3_1, 1, facecolor='black', alpha=0.05)       
##    ax[1].axvspan(0.09, 1, facecolor='blue', alpha=0.05)         
    plt.savefig('./Plots/plot_runaway_nh3.pdf', orientation='portrait',format='pdf')


def plot_runaway_solar():
    """
    Plot runaway of NH3 (H2-dominated atmosphere), CH4 (CO2-N2 atmosphere), CO (N2-CO2 atmosphere), O2 (N2-CO2 atmosphere).
    """
    ########################    
    ###Process CO, Sun
    ########################    
    file_prefix='./scenario_library/Sun/N2-CO2-Full-CO/CO_'
    conc_file_list=np.array(['3E8', '1E9','3E9','1E10','3E10','1E11','3E11','5E11','6E11','8E11','1E12','2E12'])
    file_suffix='/ConcentrationSTD.dat'    
    sun_co_flux_list=np.array([3.0E+8, 1.0E+9, 3.0E+9, 1.0E+10, 3.0E+10, 1.0E+11, 3.0E+11, 5.0E+11, 6.0E+11,  8.0E+11, 1.0E+12,2.0E+12])

    ###   
    #Extract concentration data
    ###
    num_file=len(conc_file_list)
    sun_co_surface_mixing_ratios=np.zeros(np.shape(sun_co_flux_list))
    sun_co_column_averaged_mixing_ratios=np.zeros(np.shape(sun_co_flux_list))

    for ind in range(0, num_file):
        conc_file=file_prefix+conc_file_list[ind]+file_suffix
       
        ###Get CO mixing ratios (surface and column-averaged)
        sun_co_surface_mixing_ratios[ind], sun_co_column_averaged_mixing_ratios[ind]=get_conc_std_data(conc_file, np.array([ind_co]))

    ########################    
    ###Process CO, GJ876
    ########################    
    file_prefix='./scenario_library/GJ876/N2-CO2-Full-CO/CO_'
    conc_file_list=np.array(['3E8', '1E9','3E9','6E9', '1E10','2E10', '3E10','1E11', '3E11', '1E12', '2E12'])
    file_suffix='/ConcentrationSTD.dat'    
    gj876_co_flux_list=np.array([3.0E+8, 1.0E+9, 3.0E+9, 6.0E+9, 1.0E+10, 2.0e+10, 3.0E+10, 1.0E+11, 3.0E+11, 1.0E+12, 2.0E+12])

    ###   
    #Extract concentration data
    ###
    num_file=len(conc_file_list)
    gj876_co_surface_mixing_ratios=np.zeros(np.shape(gj876_co_flux_list))
    gj876_co_column_averaged_mixing_ratios=np.zeros(np.shape(gj876_co_flux_list))

    for ind in range(0, num_file):
        conc_file=file_prefix+conc_file_list[ind]+file_suffix
       
        ###Get CO mixing ratios (surface and column-averaged)
        gj876_co_surface_mixing_ratios[ind], gj876_co_column_averaged_mixing_ratios[ind]=get_conc_std_data(conc_file, np.array([ind_co]))
        
    ########################    
    ###Process O2, Sun
    ########################    
    file_prefix='./scenario_library/Sun/N2-CO2-Full-O2/O2_'
    conc_file_list=np.array(['1E8','3E8','1E9','3E9', '6E9', '1E10','2E10', '6E10', '9E10', '4E11', '6E11', '1E12'])
    file_suffix='/ConcentrationSTD.dat'    
    sun_o2_flux_list=np.array([1.0E+8, 3.0E+8, 1.0E+9, 3.0E+9, 6.0E+9, 1.0E+10, 2.0E+10, 6.0E+10, 9.0E+10, 4.0E+11, 6.0E+11, 1.0E+12]) #omitted 3E10 b/c solution doesn't re-converge quite as fast as I'd prefer. 

    ###   
    #Extract concentration data
    ###
    num_file=len(conc_file_list)
    sun_o2_surface_mixing_ratios=np.zeros(np.shape(sun_o2_flux_list))
    sun_o2_column_averaged_mixing_ratios=np.zeros(np.shape(sun_o2_flux_list))

    for ind in range(0, num_file):
        conc_file=file_prefix+conc_file_list[ind]+file_suffix
       
        ###Get O2 mixing ratios (surface and column-averaged)
        sun_o2_surface_mixing_ratios[ind], sun_o2_column_averaged_mixing_ratios[ind]=get_conc_std_data(conc_file, np.array([ind_o2]))

    ########################    
    ###Process O2, GJ876
    ########################    
    file_prefix='./scenario_library/GJ876/N2-CO2-Full-O2/O2_'
    conc_file_list=np.array(['1E8','3E8','1E9','3E9', '1E10','2E10', '3E10', '6E10', '1E11', '3E11', '1E12'])
    file_suffix='/ConcentrationSTD.dat'    
    gj876_o2_flux_list=np.array([1.0E+8, 3.0E+8, 1.0E+9, 3.0E+9, 1.0E+10, 2.0E+10, 3.0E+10, 6.0E+10, 1.0E+11, 3.0E+11, 1.0E+12])

    # file_prefix='./scenario_library/GJ876/H2-Full-O2/O2_'
    # conc_file_list=np.array(['1E8', '1E9','3E9', '1E10','3E10', '1E11', '3E11', '1E12'])
    # file_suffix='/ConcentrationSTD.dat'    
    # gj876_o2_flux_list=np.array([1.0E+8, 1.0E+9, 3.0E+9, 1.0E+10, 3.0E+10, 1.0E+11, 3.0E+11, 1.0E+12])
    ##   
    #Extract concentration data
    ###
    num_file=len(conc_file_list)
    gj876_o2_surface_mixing_ratios=np.zeros(np.shape(gj876_o2_flux_list))
    gj876_o2_column_averaged_mixing_ratios=np.zeros(np.shape(gj876_o2_flux_list))

    for ind in range(0, num_file):
        conc_file=file_prefix+conc_file_list[ind]+file_suffix
       
        ###Get O2 mixing ratios (surface and column-averaged)
        gj876_o2_surface_mixing_ratios[ind], gj876_o2_column_averaged_mixing_ratios[ind]=get_conc_std_data(conc_file, np.array([ind_o2]))
        
    ########################    
    ###Process CH4, Sun
    ########################    
    file_prefix='./scenario_library/Sun/CO2-Full-CH4/CH4_'
    conc_file_list=np.array(['3E8','1E9','3E9', '1E10', '2E10', '3E10', '6E10', '1E11', '3E11'])
    file_suffix='/ConcentrationSTD.dat'    
    sun_ch4_flux_list=np.array([3.0E+8, 1.0E+9, 3.0E+9, 1.0E+10, 2.0E+10, 3.0E+10, 6.0E+10, 1.0E+11, 3.0E+11])

    ###   
    #Extract concentration data
    ###
    num_file=len(conc_file_list)
    sun_ch4_surface_mixing_ratios=np.zeros(np.shape(sun_ch4_flux_list))
    sun_ch4_column_averaged_mixing_ratios=np.zeros(np.shape(sun_ch4_flux_list))

    for ind in range(0, num_file):
        conc_file=file_prefix+conc_file_list[ind]+file_suffix
       
        ###Get CH4 mixing ratios (surface and column-averaged)
        sun_ch4_surface_mixing_ratios[ind], sun_ch4_column_averaged_mixing_ratios[ind]=get_conc_std_data(conc_file, np.array([ind_ch4]))

    ########################    
    ###Process CH4, GJ876
    ########################    
    file_prefix='./scenario_library/GJ876/CO2-Full-CH4/CH4_'
    conc_file_list=np.array(['3E8','1E9','3E9', '1E10', '3E10', '1E11', '3E11'])
    file_suffix='/ConcentrationSTD.dat'    
    gj876_ch4_flux_list=np.array([3.0E+8, 1.0E+9, 3.0E+9, 1.0E+10, 3.0E+10, 1.0E+11, 3.0E+11])

    # conc_file_list=np.array(['3E8','1E9','3E9', '1E10', '3E10', '1E11', '3E11'])
    # file_suffix='_low_O2_CO/ConcentrationSTD.dat'    
    # gj876_ch4_flux_list=np.array([3.0E+8, 1.0E+9, 3.0E+9, 1.0E+10, 3.0E+10, 1.0E+11, 3.0E+11])

    ###   
    #Extract concentration data
    ###
    num_file=len(conc_file_list)
    gj876_ch4_surface_mixing_ratios=np.zeros(np.shape(gj876_ch4_flux_list))
    gj876_ch4_column_averaged_mixing_ratios=np.zeros(np.shape(gj876_ch4_flux_list))

    for ind in range(0, num_file):
        conc_file=file_prefix+conc_file_list[ind]+file_suffix
       
        ###Get CH4 mixing ratios (surface and column-averaged)
        gj876_ch4_surface_mixing_ratios[ind], gj876_ch4_column_averaged_mixing_ratios[ind]=get_conc_std_data(conc_file, np.array([ind_ch4]))
        
    ########################
    ###Plot
    ########################
    fig, ax=plt.subplots(3, figsize=(8., 10.), sharex=False)
    markersizeval=5.
    ax[0].set_title(r'CO, N$_2$-CO$_2$ atm')
    ax[0].plot(sun_co_flux_list, sun_co_column_averaged_mixing_ratios, linewidth=2, linestyle=':', color='pink', label=r'Sun', marker='s', markersize=markersizeval)
    ax[0].plot(gj876_co_flux_list, gj876_co_column_averaged_mixing_ratios, linewidth=2, linestyle='-', color='hotpink', label=r'GJ876', marker='s', markersize=markersizeval)
    # ax[0,0].axvline(get_integrated_uv('./Data/solar00.txt', 1.0, 230.), color='black', linestyle='--', label=r'J$_{H_{2}O}$ (Sun)') #taken from 
    # ax[0,0].axvline(get_integrated_uv('./Data/GJ876_00.txt', 1.0, 230.), color='black', linestyle='-', label=r'J$_{H_{2}O}$ (GJ876)') #taken from 
    

    
    ax[1].set_title(r'O$_2$, N$_2$-CO$_2$ atm')
    ax[1].plot(sun_o2_flux_list, sun_o2_column_averaged_mixing_ratios,linewidth=2, linestyle=':', color='green', label='Sun', marker='s', markersize=markersizeval)
    ax[1].plot(gj876_o2_flux_list, gj876_o2_column_averaged_mixing_ratios,linewidth=2, linestyle='-', color='darkgreen', label='GJ876', marker='s', markersize=markersizeval)
    
    ax[2].set_title(r'CH$_4$, CO$_2$-N$_2$ atm')
    ax[2].plot(sun_ch4_flux_list, sun_ch4_column_averaged_mixing_ratios, linewidth=2, linestyle=':', color='blue', label='Sun', marker='s', markersize=markersizeval)
    ax[2].plot(gj876_ch4_flux_list, gj876_ch4_column_averaged_mixing_ratios, linewidth=2, linestyle='-', color='darkblue', label='GJ876', marker='s', markersize=markersizeval)

    # ax[1,1].set_title(r'NH$_3$, H$_2$-N$_2$ atm.')
    # ax[1,1].plot(sun_nh3_flux_list, sun_nh3_column_averaged_mixing_ratios,linewidth=2, linestyle=':', color='black', label='Sun', marker='s', markersize=markersizeval)
    # ax[1,1].plot(gj876_nh3_flux_list, gj876_nh3_column_averaged_mixing_ratios,linewidth=2, linestyle='-', color='black', label='GJ876', marker='s', markersize=markersizeval)

    abscissa=np.logspace(6,13)
    for ind in range(0,3):
        ax[ind].plot(abscissa, abscissa/(2.5E19*1.0E-7), linestyle='--', color='red', label=r'v$_{dep}$-limited')
        ax[ind].set_xscale('log')
        ax[ind].set_yscale('log')
        ax[ind].set_ylabel(r'Column-Averaged Mixing Ratio')
        ax[ind].set_xlabel(r'Surface Flux')  
        ax[ind].legend(loc='upper left', ncol=1, borderaxespad=0., fontsize=10)   

    ax[0].set_xlim([3.0E+8, 3.0E+12])
    ax[0].set_ylim([1.0E-6, 1.0])
    
    ax[1].set_xlim([3.0E+8, 1.0E+12])
    ax[1].set_ylim([1.0E-8, 1.0])
    
    ax[2].set_xlim([3.0E+8, 1.0E+11])
    ax[2].set_ylim([1.0E-6, 1.0])


 
    plt.subplots_adjust(hspace=0.4)

    plt.savefig('./Plots/plot_runaway_solar_gj876_co_o2_ch4.pdf', orientation='portrait',format='pdf')

    plt.show()    
    

#
# =============================================================================
###############################
###Run
###############################

plot_runaway_nh3()
plot_runaway_solar()


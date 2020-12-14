#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 11:40:16 2019

@author: shungokoyama
"""

import pandas as pd
import xlrd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from matplotlib import ticker, cm
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
import scipy.ndimage
import numpy.ma as ma
import matplotlib.gridspec as gridspec
import h5py

####
#This file is to make figures used in my paper.
###
df_reg = pd.read_csv("/Users/shungokoyama/Programming/result/stat/summary_reg_Tmod_forplot_latest.csv",delimiter=',')
CO2array = df_reg['CO2pressure'][24:48] #it doesnt include index 11 #6.3-3000mbar

df_red=pd.read_csv("/Users/shungokoyama/Programming/result/stat/summary_conv_Tmod_forplot_latest.csv",delimiter=',')


#Main result
#This function creates the figure of H escape flux histories and some specific species mixing ratio over time
#vertically adjusted each other
#plot_Hesc_species("Hesc_12e8_Tmod/Hesc_depfluxes_CO2_500mbar_Tsurf_240_Oesc_1.2e8to2.4e8_depv_0.02_nofixedCO2_noHOCO_Tmod.h5","CO2_12e8_Tmod/CO2_500mbar_Tsurf_240_Oesc_1.2e8to2.4e8_depv_0.02_nofixedCO2_noHOCO_Tmod.h5",['CO','O2'],2.4e8)
def plot_Hesc_species(file_for_Hesc,file_for_species,speciesname_array,Oescflux_after):
    #file = file_for_Hesc
    #h5files = []
    #array_species_overtime = species_mixing_overtime(file_for_species,speciesname,1599)
    fluxvals_mat = []
    h5file = h5py.File(file_for_Hesc,"r")
    time_array = h5file["fluxes/times"].value #get time array here
    time_array = time_array/3.15e7
    
    #create 2 subplots vertically adjusted
    fig, axs = plt.subplots(2,1,sharex=True,dpi=150)
    fig.subplots_adjust(hspace=0)
    
    ##ax0 H escape flux
    fluxvals_mat = h5file["fluxes/fluxvals"].value
    H2O2dep=h5file["depositions/H2O2"].value
    O3dep=h5file["depositions/O3"].value
    HO2dep=h5file["depositions/HO2"].value
    Ototalflux_redox=H2O2dep + O3dep*3 + HO2dep*1.5 + Oescflux_after
    
    axs[0].plot(time_array, fluxvals_mat[0:-1], color='k',label='H escape')
    axs[0].plot(time_array,Ototalflux_redox[0:-1], color='k',linestyle=':',label='net O outflux')
    #axs[0].set_xlabel('time(s)')
    axs[0].set_ylabel(r"Flux$[cm^{-2}s^{-1}]$")
    axs[0].set_ylim((2e8,2e9))
    #axs[0].set_ylim((1e8,5e12))
    axs[0].set_yscale("log")
    #axs[0].tick_params(axis='y')
    axs[0].legend(prop={'size':9.5})
    
    #ax2 species mixing ratio
    axs[1].set_ylabel('Mixing ratio')  # we already handled the x-label with ax1
    axs[1].set_xlabel('Time [year]')
    axs[1].set_yscale("log")
    axs[1].set_xscale("log")
    #axs[1].set_ylim((1e-7,3e-4))
    axs[1].set_ylim((2e-4,1.2e-3))
    for speciesname in speciesname_array:
        array_species_overtime = species_mixing_overtime(file_for_species,speciesname,len(time_array-1))
        axs[1].plot(time_array,array_species_overtime,label=speciesname,color=species_color[speciesname])
        #axs[1].tick_params(axis='y')
        print(array_species_overtime[0]-array_species_overtime[len(time_array)-1])
    
    ##setting for both ax1 and ax2
    axs[1].legend(prop={'size':9.5})
    #axs[1].legend(prop={'size':9.5},loc='upper left')
    axs[1].set_xlim((1e2,1e9))
    #plt.xscale("log")
    #plt.title('H escape flux [H20:9.5pr]')
    plt.show()

###########################################################
##Regulation with partial pressure###
## Main figure of the paper 2 縦長バージョン ############
fig=plt.figure(dpi=300,figsize=(4,3))
outer = gridspec.GridSpec(1, 3, wspace=0.0, hspace=0.1)
for i in range(3):

    inner = gridspec.GridSpecFromSubplotSpec(2, 1,
                    subplot_spec=outer[i], wspace=0.15, hspace=0)
    i=3*i+1
    
    ax2=plt.Subplot(fig, inner[0])
    ax2.plot(df_red['CO2pressure'][i*24:24*(i+1)],np.array(df_red['fO2'][i*24:24*(i+1)],dtype=np.float)*CO2array,label='O2',linewidth=0.8,linestyle='--')
    ax2.plot(df_red['CO2pressure'][i*24:24*(i+1)],np.array(df_red['fCO'][i*24:24*(i+1)],dtype=np.float)*CO2array,label='CO',linewidth=0.8,linestyle='-.')
    ax2.plot(df_red['CO2pressure'][i*24:24*(i+1)],np.array(df_red['fH2'][i*24:24*(i+1)],dtype=np.float)*CO2array,label='H2',linewidth=0.8,linestyle=':')
    #ax2.plot(CO2array,np.array(df_red['fOx'][i*22:22*(i+1)],dtype=np.float)*CO2array,label='pOx',linewidth='0.5',color='gray')
    ax2.plot(CO2array,np.abs(np.array(df_red['fOx'][i*24:24*(i+1)],dtype=np.float)*CO2array),label='|pOx|',marker='.',linewidth='1.0',markersize='1',color='black')
    
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_ylim((5e-6,8e1))
    ax2.set_xlim((5,4e3))
    plt.setp(ax2.get_xticklabels(), visible=False)
    fig.subplots_adjust(hspace=0)

    ax2.text(7,20,'T$_{\mathrm{surf}}$:'+str(i*10 + 200)+'K', size='6')

    ax1 = plt.Subplot(fig, inner[1])
    ax1.plot(df_reg['CO2pressure'][i*24:24*(i+1)],np.array(df_reg['regulationtimescale'][i*24:24*(i+1)],dtype=np.float)/3.1536e7,linewidth='1.0',markersize='1',color='k',label='regulation timescale',marker='.')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim((5,4e3))
    ax1.set_ylim((5e3,2e9))
    if i==1:
        ax1.set_ylabel('Regulation Time [years]',size='5')
        ax2.set_ylabel('Partial Pressure [mbar]',size='5')
    else:
        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.setp(ax1.get_yticklabels(), visible=False)
    ax2.tick_params(labelsize=5)
    ax1.set_xlabel(r'pCO$_2$ [mbar]',size='5')
    #ax1.set_ylabel('Regulation Timescale[years]',size='4')
    #plt.legend()
    #ax1.text(10,10,r'T$_{surf}$:'+str(i*10 + 200)+'K', size='6')
    ax1.tick_params(labelsize=5)
    #ax2.legend(prop={'size': 3})
    #ax1.legend(prop={'size': 3})
    fig.add_subplot(ax1)
    fig.add_subplot(ax2)
    fig.subplots_adjust(hspace=0)
ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size': 6})
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.7),prop={'size': 6})
#######################


df_red=pd.read_csv("/Users/shungokoyama/Programming/result/stat/summary_conv_Tmod_forplot_20191203.csv",delimiter=',')
#df_reg=pd.read_csv("/Users/shungokoyama/Programming/result/stat/summary_conv_Tmod_forplot_20191203.csv",delimiter=',')
CO2array = df_red['CO2pressure'][0:22] #it doesnt include index 11 #6.3-3000mbar

######################
#O2 mixing ratio vs temperature
plt.figure(dpi=150)
O2mixing_vs_temp = np.array(df_red['fO2'][0:-1:22],dtype=np.float)
temp_array = np.array(df_red['Ts'][0:-1:22],dtype=np.float)
plt.plot(temp_array,O2mixing_vs_temp, label=r"O$_2$",marker='.',color='k',linewidth=2)
plt.yscale('log')
plt.ylim((1e-6,1e-2))
plt.ylabel('O$_2$ Mixing ratio',size='15')
plt.xlabel('Surface temperature[K]',size='15')
plt.tick_params(labelsize=12)
plt.text(220,6e-3,'pCO$_2$(initial): 6.3mbar',size='12')
#plt.legend()
######################


#####CO possibility of early Mars#########################
#CO mixing ratio vs CO2 for various surface temperature
plt.figure(dpi=150)
i=0
plt.plot(CO2array,np.array(df_red['fCO'][i*22:22*(i+1)],dtype=np.float),
             label=str(i*10 + 200)+'K',marker='.',color='k',linestyle='--',markersize='5',fillstyle='none')
i=1
plt.plot(CO2array,np.array(df_red['fCO'][i*22:22*(i+1)],dtype=np.float),
             label=str(i*10 + 200)+'K',marker='.',color='k',linestyle='-.',markersize='5',fillstyle='none')
i=2
plt.plot(CO2array,np.array(df_red['fCO'][i*22:22*(i+1)],dtype=np.float),
             label=str(i*10 + 200)+'K',marker='.',color='k',linestyle=':',markersize='5',fillstyle='none')
i=3
plt.plot(CO2array,np.array(df_red['fCO'][i*22:22*(i+1)],dtype=np.float),
             label=str(i*10 + 200)+'K',marker='.',color='k',markersize='5',fillstyle='none')
plt.xlabel("Atmospheric CO2 [mbar]",size=15)
plt.ylabel("CO Mixing ratio",size=15)
plt.legend()
plt.tick_params(labelsize=12)
plt.xscale('log')
plt.yscale('log')
#########################################


### sensitivity of CO depostion velocity ######################################
#pCO2 = 100 mbar
#Ts=200 K
#CO2 fixed
df_COdep = pd.read_csv("/Users/shungokoyama/Programming/result/stat/COdep.csv",delimiter=',',header=None)
COmixing_1G = df_COdep[11] #Vdep = 1e-4~1e-11
plt.figure(dpi=150)
vdep=np.array([1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11])
plt.plot(vdep,COmixing_1G,marker='.',color='k',linestyle='--',markersize='10',fillstyle='none')
plt.xlabel(r"Deposition velocity [cm $s^{-1}$]",size=15)
plt.ylabel("CO Mixing ratio",size=15)
plt.xscale('log')
plt.yscale('log')
plt.tick_params(labelsize=12)
plt.ylim((1e-4,1))
###############################################################################


########regulation timescale vs O escape rate#######
#CO2: 1bar
#Ts: 250K
#This graph shows O escape sensitivity test of early Mars.
plt.figure(dpi=150)
reg_array_240=np.array([8.228414098500866e13,9.035798679480902e13,9.710376957998898e13,9.92E+13,1.2314749526993917e14,1.4324900837921494e14])/3.15e7
reg_array_250=np.array([58240129368216.11,70230076908267.24,75473186683228.2,76018560560614.81,78239750923735.7,77678441738026.81])/3.15e7
Oesc_array=np.array([1e7,5e7,1e8,1.2e8,5e8,1e9])
plt.plot(Oesc_array,reg_array_250,marker='.',label='CO2:1bar , Ts=250K',color='k')
plt.ylim((1e4,1e10))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"O escape flux $[cm^{-2}s^{-1}]$")
plt.ylabel('Regulation timescale [years]')
plt.legend()
#############################################

####regulation vs exospheric temperature(effusion velocity)#############
#Tinf = 240-400K
#CO2:1bar, Ts=250K
plt.figure(dpi=150)
reg_array_Tinf=np.array([7.6e13,4.01E+13,3.73E+13,3.44E+13,3.20395E+13])/3.15e7
Tinf_array=np.array([240,280,320,360,400])
plt.plot(Tinf_array,reg_array_Tinf,marker='.',label='CO2:1bar , Ts=250K',color='k')
plt.ylim((1e4,1e10))
#plt.xscale('log')
plt.yscale('log')
plt.tick_params(labelsize=12)
plt.xlabel("Exospheric Temperature [K]",size=15)
plt.ylabel('Regulation timescale [years]',size=15)
#plt.title("CO2:1bar, Ts:250K")
#plt.legend()
#############################################


#####temperature profile and water vapor profile################
#temperature
plt.figure(dpi=150)
temp_data = pd.read_csv("/Users/shungokoyama/Programming/result/stat/temp_profile_10mbar.csv")
plt.plot(temp_data['T200'],temp_data['alt200']/1e5,label="Ts=200K",color="0.1")
plt.plot(temp_data['T220'],temp_data['alt220']/1e5,label="Ts=220K",color="0.3")
plt.plot(temp_data['T240'],temp_data['alt240']/1e5,label="Ts=240K",color="0.5")
plt.plot(temp_data['T260'],temp_data['alt260']/1e5,label="Ts=260K",color="0.7")
plt.plot(temp_data['T280'],temp_data['alt280']/1e5,label="Ts=280K",color="0.9")
plt.xlabel("Temperature[K]",fontsize=15)
plt.ylabel("Altitude[km]",fontsize=15)
plt.ylim((0,250))
plt.tick_params(labelsize=12)
plt.legend()

#water vapor
plt.figure(dpi=150)
def H2O_abundance(file_2):
    file2 = file_2
    h5file2 = h5py.File(file2,"r") #read a h5 file
    initial_mat2 = h5file2["n_current/n_current_mat"].value
    alt = h5file2["n_current/alt"].value

    return(initial_mat2[20,:]),(initial_mat2[22,:]), alt[1:-1]

H2O200, CO2200, alt200 = H2O_abundance("/Users/shungokoyama/Programming/result/stat/conv_Tmod_Tsurf200_Oesc1.2e8/converged_Tsurf_200_CO2_10mbar_Oesc_1.2e8_depv_0.02_nofixedCO2_noHOCO_Tmod_1billion.h5")
H2O220, CO2220, alt220 = H2O_abundance("/Users/shungokoyama/Programming/result/stat/conv_Tmod_Tsurf220_Oesc1.2e8/converged_Tsurf_220_CO2_10mbar_Oesc_1.2e8_depv_0.02_nofixedCO2_noHOCO_Tmod_1billion.h5")
H2O240, CO2240, alt240 = H2O_abundance("/Users/shungokoyama/Programming/result/stat/conv_Tmod_Tsurf240_Oesc1.2e8/converged_Tsurf_240_CO2_10mbar_Oesc_1.2e8_depv_0.02_nofixedCO2_noHOCO_Tmod_1billion.h5")
H2O260, CO2260, alt260 = H2O_abundance("/Users/shungokoyama/Programming/result/stat/conv_Tmod_Tsurf260_Oesc1.2e8/converged_Tsurf_260_CO2_10mbar_Oesc_1.2e8_depv_0.02_nofixedCO2_noHOCO_Tmod_1billion.h5")
H2O280, CO2280, alt280 = H2O_abundance("/Users/shungokoyama/Programming/result/stat/conv_Tmod_Tsurf280_Oesc1.2e8/converged_Tsurf_280_CO2_10mbar_Oesc_1.2e8_depv_0.02_nofixedCO2_noHOCO_Tmod_1billion.h5")
plt.plot(H2O200,alt200/1e5,label="Ts=200K",color="0.1")
plt.plot(H2O220,alt220/1e5,label="Ts=220K",color="0.3")
plt.plot(H2O240,alt240/1e5,label="Ts=240K",color="0.5")
plt.plot(H2O260,alt260/1e5,label="Ts=260K",color="0.7")
plt.plot(H2O280,alt280/1e5,label="Ts=280K",color="0.9")
plt.xscale('log')
plt.ylabel("Altitude[km]",fontsize=15)
plt.ylim((0,250))
plt.xlabel(r"Number density$ [cm^{-3}]$",fontsize=15)
plt.tick_params(labelsize=12)
plt.legend()




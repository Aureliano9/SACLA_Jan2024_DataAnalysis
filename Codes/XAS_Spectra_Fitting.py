# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 16:06:00 2022

@author: nurekeye

Plotting and reorganizing XAS Spectra data from SACLA
csv files are generated by sorting code
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
from collections import defaultdict
from scipy.optimize import curve_fit
import lmfit as lf



def custom_arctan(x,x0,y0,amp,width_arctan):
    return y0 + amp*(0.5 + np.arctan((x-x0) / width_arctan) / np.pi)

def custom_lorentzian(x,xc,area_lor,width_lor):
    return (2*area_lor / np.pi) * (width_lor / (4*(x - xc)**2 + width_lor**2) )


RunNumber = np.array([1145431])
Runs = RunNumber

Directory = 'C:/Users/nurekeye/Desktop/Projects/SACLA_01-2024/DataAnalysis/Codes/' + str(int(RunNumber)) + '/'
Save_Directory = Directory


# Create empty dictionaries to fill them with Pandas DataFrames, each key corresponds to a run
data = {} #timing tool corrected
datau = {} #uncorrected
for i in Runs:
    data[str(i)] = pd.read_csv(Directory+str(i)+'_0.csv', delimiter=',').replace(' ',np.nan).astype(float)
    datau[str(i)] = pd.read_csv(Directory+str(i)+'_0_raw.csv', delimiter=',').replace(' ',np.nan).astype(float)

fit_res = {}
fit_report = {}
fit_range_u = 13477
fit_range_l = 13458

#%% Laser ON
for k in Runs:
    x_axis_on = data[str(k)].Energy_on[(~np.isnan(data[str(k)].TFY_on)) & (data[str(k)].Energy_on<fit_range_u) & (data[str(k)].TFY_on_std!=0)]
    y_axis_on = data[str(k)].TFY_on[(~np.isnan(data[str(k)].TFY_on)) & (data[str(k)].Energy_on<fit_range_u) & (data[str(k)].TFY_on_std!=0)]
    stddev = data[str(k)].TFY_on_std[(~np.isnan(data[str(k)].TFY_on)) & (data[str(k)].Energy_on<fit_range_u) & (data[str(k)].TFY_on_std!=0)]
    
    gmodel = lf.Model(custom_arctan) + lf.Model(custom_lorentzian)
    params = lf.create_params( x0=dict(value=13472, vary=True, min=13470, max=13480), #Edge in eV
                          y0=dict(value=min(y_axis_on), vary=True, min=-0.1, max=0.2), #Background level in a.u.
                          amp=dict(value=0.22, vary=False, min=0.01, max=2), #Arctangent function amplitude in a.u.
                          width_arctan=dict(value=2.5, vary=True, min=0.8, max=5), #Arctangent function width in ev
                          xc=dict(value=13464, vary=True, min=13460, max=13480), #Lorentzian peak position in eV
                          area_lor=dict(value=0.1, vary=True, min=0.01, max=10), #Lorentzain amplitude in a.u.
                          width_lor=dict(value=3, vary=True, min=0.1, max=5)) #Lorentzian width in eV
    
    result = gmodel.fit(y_axis_on, params=params, x=x_axis_on, 
                        weights = 1/stddev,
                        method = 'nelder')
    fit_res[str(k)] = result #save fit results in internal lmfit format
    fit_report[str(k)] = result.fit_report() #save fit report

    figs, axes = plt.subplots(2,1, sharex=True, height_ratios=[5,1], num=k, figsize=[16,8], dpi=200)
    figs.subplots_adjust(hspace=0)
    axes[0].plot(x_axis_on, y_axis_on, '.', label='data')
    axes[0].plot(x_axis_on, fit_res[str(k)].init_fit, '--', label='initial fit')
    axes[0].plot(x_axis_on, fit_res[str(k)].best_fit, '-', label='best fit')
    comps = fit_res[str(k)].eval_components(x=x_axis_on)
    dely = fit_res[str(k)].eval_uncertainty(sigma=3)
    axes[0].fill_between(x_axis_on, fit_res[str(k)].best_fit-dely, fit_res[str(k)].best_fit+dely,
                            color="#C5C9C7", label=r'3-$\sigma$ band')
    axes[0].legend(loc='lower right')
    axes[0].set_ylabel('Intensity (a.u.)')
    axes[0].text(2*min(x_axis_on)-max(x_axis_on), min(y_axis_on), s=fit_res[str(k)].fit_report(), fontsize=7)

    axes[1].plot(x_axis_on, fit_res[str(k)].residual)
    # axes[1].set_ylabel('Resid')
    axes[1].set_xlim([2*min(x_axis_on)-max(x_axis_on), max(x_axis_on)])
    axes[1].set_xlabel('Energy (eV)')
    
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()
    # plt.savefig(Save_Directory +'/'+str(k) + '_LaserON_fit.png',dpi=300) #saving figures in .png file for 2.9 uJ
    print(fit_res[str(k)].values)
    plt.close()

    x_axis2_on = data[str(k)].Energy_on[(~np.isnan(data[str(k)].TFY_on)) &  (data[str(k)].TFY_on_std!=0)]
    y_axis2_on = data[str(k)].TFY_on[(~np.isnan(data[str(k)].TFY_on)) & (data[str(k)].TFY_on_std!=0)]
    stddev2 = data[str(k)].TFY_on_std[(~np.isnan(data[str(k)].TFY_on)) & (data[str(k)].TFY_on_std!=0)]
    x0 = fit_res[str(k)].values['x0']
    y0 = fit_res[str(k)].values['y0']
    amp = fit_res[str(k)].values['amp']
    width_arctan = fit_res[str(k)].values['width_arctan']
    xc = fit_res[str(k)].values['xc']
    area_lor = fit_res[str(k)].values['area_lor']
    width_lor = fit_res[str(k)].values['width_lor']
    y_arctan = custom_arctan(x_axis2_on, x0, y0, amp, width_arctan)
    y_lor = custom_lorentzian(x_axis2_on, xc, area_lor, width_lor)
    y_fit = y_arctan+y_lor
    y_resid = (y_arctan + y_lor - y_axis2_on) / stddev2
    figs, axes = plt.subplots(2,1, sharex=True, height_ratios=[5,1], num=k*2, figsize=[16,8], dpi=200)
    figs.subplots_adjust(hspace=0)
    axes[0].plot(x_axis2_on, y_axis2_on, '.', label='data')
    axes[0].plot(x_axis2_on, y_arctan + y_lor, '-', label='best fit')
    axes[0].legend(loc='lower right')
    axes[0].set_ylabel('Intensity (a.u.)')
    axes[0].text(2*min(x_axis2_on)-max(x_axis2_on), min(y_axis2_on), s=fit_res[str(k)].fit_report(), fontsize=7)
    axes[1].plot(x_axis2_on, y_resid)
    # axes[1].set_ylabel('Resid')
    axes[1].set_xlim([2*min(x_axis2_on)-max(x_axis2_on), max(x_axis2_on)])
    axes[1].set_xlabel('Energy (eV)')
    plt.savefig(Save_Directory +'/'+str(k) + '_LaserON_fit.png',dpi=300) #saving figures in .png file for 2.9 uJ
    plt.close()
    
#%% Laser OFF
for k in Runs:
    x_axis_off = data[str(k)].Energy_off[(~np.isnan(data[str(k)].TFY_off)) & (data[str(k)].Energy_off<fit_range_u) & (data[str(k)].TFY_off_std!=0)]
    y_axis_off = data[str(k)].TFY_off[(~np.isnan(data[str(k)].TFY_off)) & (data[str(k)].Energy_off<fit_range_u) & (data[str(k)].TFY_off_std!=0)]
    stddev = data[str(k)].TFY_off_std[(~np.isnan(data[str(k)].TFY_off)) & (data[str(k)].Energy_off<fit_range_u) & (data[str(k)].TFY_off_std!=0)]
    
    gmodel = lf.Model(custom_arctan)
    params = lf.create_params( x0=dict(value=13472, vary=True, min=13470, max=13475), #Edge in eV
                          y0=dict(value=min(y_axis_off), vary=True, min=-0.1, max=0.2), #Background level in a.u.
                          amp=dict(value=0.22, vary=False, min=0.01, max=0.2), #Arctangent function amplitude in a.u.
                          width_arctan=dict(value=1, vary=True, min=0.8, max=2)) #Arctangent function width in ev
    
    result = gmodel.fit(y_axis_off, params=params, x=x_axis_off, 
                        weights = 1/stddev,
                        method = 'nelder')
    fit_res[str(k)] = result #save fit results in internal lmfit format
    fit_report[str(k)] = result.fit_report() #save fit report

    figs, axes = plt.subplots(2,1, sharex=True, height_ratios=[5,1], num=k*3, figsize=[16,8], dpi=200)
    figs.subplots_adjust(hspace=0)
    axes[0].plot(x_axis_off, y_axis_off, '.', label='data')
    axes[0].plot(x_axis_off, fit_res[str(k)].init_fit, '--', label='initial fit')
    axes[0].plot(x_axis_off, fit_res[str(k)].best_fit, '-', label='best fit')
    comps = fit_res[str(k)].eval_components(x=x_axis_off)
    dely = fit_res[str(k)].eval_uncertainty(sigma=3)
    axes[0].fill_between(x_axis_off, fit_res[str(k)].best_fit-dely, fit_res[str(k)].best_fit+dely,
                            color="#C5C9C7", label=r'3-$\sigma$ band')
    axes[0].legend(loc='lower right')
    axes[0].set_ylabel('Intensity (a.u.)')
    axes[0].text(2*min(x_axis_off)-max(x_axis_off), min(y_axis_off), s=fit_res[str(k)].fit_report(), fontsize=7)
    axes[1].plot(x_axis_off, fit_res[str(k)].residual)
    # axes[1].set_ylabel('Resid')
    axes[1].set_xlim([2*min(x_axis_off)-max(x_axis_off), max(x_axis_off)])
    axes[1].set_xlabel('Energy (eV)')
    
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()
    # plt.savefig(Save_Directory +'/'+str(k) + '_LaserOFF_fit2.png',dpi=300) #saving figures in .png file for 2.9 uJ
    print(fit_res[str(k)].values)
    plt.close()

         
    x_axis2_off = data[str(k)].Energy_off[(~np.isnan(data[str(k)].TFY_off)) &  (data[str(k)].TFY_off_std!=0)]
    y_axis2_off = data[str(k)].TFY_off[(~np.isnan(data[str(k)].TFY_off)) & (data[str(k)].TFY_off_std!=0)]
    stddev2 = data[str(k)].TFY_off_std[(~np.isnan(data[str(k)].TFY_off)) & (data[str(k)].TFY_off_std!=0)]
    x0 = fit_res[str(k)].values['x0']
    y0 = fit_res[str(k)].values['y0']
    amp = fit_res[str(k)].values['amp']
    width_arctan = fit_res[str(k)].values['width_arctan']
    y_arctan = custom_arctan(x_axis2_off, x0, y0, amp, width_arctan)
    y_fit = y_arctan
    y_resid = (y_arctan - y_axis2_off) / stddev2
    fig, axes = plt.subplots(2,1, sharex=True, height_ratios=[5,1], num=k*4, figsize=[16,8], dpi=200)
    figs.subplots_adjust(hspace=0)
    axes[0].plot(x_axis2_off, y_axis2_off, '.', label='data')
    axes[0].plot(x_axis2_off, y_arctan, '-', label='fit')
    axes[0].legend(loc='lower right')
    axes[0].set_ylabel('Intensity (a.u.)')
    axes[0].text(2*min(x_axis2_off)-max(x_axis2_off), min(y_axis2_off), s=fit_res[str(k)].fit_report(), fontsize=7)
    axes[1].plot(x_axis2_off, y_resid)
    # axes[1].set_ylabel('Resid')
    axes[1].set_xlim([2*min(x_axis2_off)-max(x_axis2_off), max(x_axis2_off)])
    axes[1].set_xlabel('Energy (eV)')
    plt.savefig(Save_Directory +'/'+str(k) + '_LaserOFF_fit.png',dpi=300) #saving figures in .png file for 2.9 uJ
    plt.close()

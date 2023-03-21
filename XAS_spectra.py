#This code is for plotting XAS spectra from SACLA beamtime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
from collections import defaultdict

def spectra_extractor(RunNumber, DataDirectory, SaveFolder, ExtraComment, FemtosecondInPls=6.671, 
                 TimeZero=1375, Attenuation=0):

    f = h5py.File(DataDirectory + str(RunNumber) + '.h5', 'r')
    I1 = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_13_in_volt'][:]
    I0_1 = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_14_in_volt'][:]
    I0_2 = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_15_in_volt'][:]
    I_OpticalLaser = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_4_in_volt'][:]
    OpticalAttenuator = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/eh_2_optical_ND_filter_stage_position'][:]
    MonochromatorAngle = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_1/eh_1_CC1_rotational_stage_position'][:]
    LaserStatus = f['/run_' + str(RunNumber) + '/event_info/bl_3/lh_1/laser_pulse_selector_status'][:]
    TagList_I = f['/run_' + str(RunNumber) + '/event_info/tag_number_list'][:]
    
    ''' Activate this section only for run 1145398. It has 2 optical delays (4373 and 61336), so it has to be sorted over them. 
    OpticalDelay = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/eh_2_optical_delay_stage_position'][:]
    filter_arr = OpticalDelay > 5000
    new_arr =  OpticalDelay[filter_arr] 
    I1 = I1[filter_arr]
    I0_1 = I0_1[filter_arr]
    I0_2 = I0_2[filter_arr]
    I_OpticalLaser = I_OpticalLaser[filter_arr]
    OpticalAttenuator = OpticalAttenuator[filter_arr]
    MonochromatorAngle = MonochromatorAngle[filter_arr]
    LaserStatus = LaserStatus[filter_arr]
    TagList_I = TagList_I[filter_arr]
    '''
    
    MonochromatorAngles = np.array(sorted(list(set(MonochromatorAngle)))) #Find out Angle values

    #Create empty arrays for detector values
    I1_on = []
    I0_1_on = []
    I0_2_on = []
    I1_on_std = []
    I0_1_on_std = []
    I0_2_on_std = []
    I1_I0_1_covariance_on = []
    I1_I0_2_covariance_on = []
    I0_1_I0_2_covariance_on = []
    I1_I0_1_covariance_norm_on = []
    I1_I0_2_covariance_norm_on = []
    I0_1_I0_2_covariance_norm_on = []
    N_values_on = []

    I1_off = []
    I0_1_off = []
    I0_2_off = []
    I1_off_std = []
    I0_1_off_std = []
    I0_2_off_std = []
    I1_I0_1_covariance_off = []
    I1_I0_2_covariance_off = []
    I0_1_I0_2_covariance_off = []
    I1_I0_1_covariance_norm_off = []
    I1_I0_2_covariance_norm_off = []
    I0_1_I0_2_covariance_norm_off = []
    N_values_off = []

    for ds in range(len(MonochromatorAngles)):
        #sum over LaserON shots
        N_values_on_temp = 0
        I1_on_temp = []
        I0_1_on_temp = []
        I0_2_on_temp = []
        for k in range(len(LaserStatus)):
            if ( (MonochromatorAngle[k] == MonochromatorAngles[ds])
                and (LaserStatus[k] == 1)
                and (OpticalAttenuator[k] == Attenuation)
                and (np.isnan(I1[k]) == False)
                and (np.isnan(I0_1[k]) == False)
                and (np.isnan(I0_2[k]) == False) 
                and (0.9>I1[k]>0.01)
                and (0.9>I0_1[k]>0.01)
                and (0.9>I0_2[k]>0.01) ):
                N_values_on_temp+= 1
                I1_on_temp.append(I1[k])
                I0_1_on_temp.append(I0_1[k])
                I0_2_on_temp.append(I0_2[k])
        I1_on.append(np.mean(I1_on_temp))
        I0_1_on.append(np.mean(I0_1_on_temp))
        I0_2_on.append(np.mean(I0_2_on_temp))
        I1_on_std.append(np.std(I1_on_temp))
        I0_1_on_std.append(np.std(I0_1_on_temp))
        I0_2_on_std.append(np.std(I0_2_on_temp))
        I1_I0_1_covariance_on.append(np.cov(I1_on_temp,I0_1_on_temp,bias=True)[1,0])
        I1_I0_2_covariance_on.append(np.cov(I1_on_temp,I0_2_on_temp,bias=True)[1,0])
        I0_1_I0_2_covariance_on.append(np.cov(I0_1_on_temp,I0_2_on_temp,bias=True)[1,0])
        I1_I0_1_covariance_norm_on.append(np.corrcoef(I1_on_temp,I0_1_on_temp)[1,0])
        I1_I0_2_covariance_norm_on.append(np.corrcoef(I1_on_temp,I0_2_on_temp)[1,0])
        I0_1_I0_2_covariance_norm_on.append(np.corrcoef(I0_1_on_temp,I0_2_on_temp)[1,0])
        N_values_on.append(N_values_on_temp)        

        #sum over LaserOFF shots
        N_values_off_temp = 0
        I1_off_temp = []
        I0_1_off_temp = []
        I0_2_off_temp = []
        for l in range(len(LaserStatus)):
            if ( (MonochromatorAngle[l] == MonochromatorAngles[ds])
                and (LaserStatus[l] == 0)
                and (OpticalAttenuator[l] == Attenuation)
                and (np.isnan(I1[l]) == False)
                and (np.isnan(I0_1[l]) == False)
                and (np.isnan(I0_2[l]) == False) 
                and (0.9>I1[l]>0.01)
                and (0.9>I0_1[l]>0.01)
                and (0.9>I0_2[l]>0.01) ):            
                N_values_off_temp+= 1
                I1_off_temp.append(I1[l])
                I0_1_off_temp.append(I0_1[l])
                I0_2_off_temp.append(I0_2[l])
        I1_off.append(np.mean(I1_off_temp))
        I0_1_off.append(np.mean(I0_1_off_temp))
        I0_2_off.append(np.mean(I0_2_off_temp))
        I1_off_std.append(np.std(I1_off_temp))
        I0_1_off_std.append(np.std(I0_1_off_temp))
        I0_2_off_std.append(np.std(I0_2_off_temp))
        I1_I0_1_covariance_off.append(np.cov(I1_off_temp,I0_1_off_temp,bias=True)[1,0])
        I1_I0_2_covariance_off.append(np.cov(I1_off_temp,I0_2_off_temp,bias=True)[1,0])
        I0_1_I0_2_covariance_off.append(np.cov(I0_1_off_temp,I0_2_off_temp,bias=True)[1,0])
        I1_I0_1_covariance_norm_off.append(np.corrcoef(I1_off_temp,I0_1_off_temp)[1,0])
        I1_I0_2_covariance_norm_off.append(np.corrcoef(I1_off_temp,I0_2_off_temp)[1,0])
        I0_1_I0_2_covariance_norm_off.append(np.corrcoef(I0_1_off_temp,I0_2_off_temp)[1,0])
        N_values_off.append(N_values_off_temp)

    #Convert into numpy arrays for usefulness
    I1_on = np.array(I1_on)
    I0_1_on = np.array(I0_1_on)
    I0_2_on = np.array(I0_2_on)
    I1_on_std = np.array(I1_on_std)
    I0_1_on_std = np.array(I0_1_on_std)
    I0_2_on_std = np.array(I0_2_on_std)
    I1_I0_1_covariance_on = np.array(I1_I0_1_covariance_on)
    I1_I0_2_covariance_on = np.array(I1_I0_2_covariance_on)
    I0_1_I0_2_covariance_on = np.array(I0_1_I0_2_covariance_on)
    I1_I0_1_covariance_norm_on = np.array(I1_I0_1_covariance_norm_on)
    I1_I0_2_covariance_norm_on = np.array(I1_I0_2_covariance_norm_on)
    I0_1_I0_2_covariance_norm_on = np.array(I0_1_I0_2_covariance_norm_on)
    N_values_on = np.array(N_values_on)

    I1_off = np.array(I1_off)
    I0_1_off = np.array(I0_1_off)
    I0_2_off = np.array(I0_2_off)
    I1_off_std = np.array(I1_off_std)
    I0_1_off_std = np.array(I0_1_off_std)
    I0_2_off_std = np.array(I0_2_off_std)
    I1_I0_1_covariance_off = np.array(I1_I0_1_covariance_off)
    I1_I0_2_covariance_off = np.array(I1_I0_2_covariance_off)
    I0_1_I0_2_covariance_off = np.array(I0_1_I0_2_covariance_off)
    I1_I0_1_covariance_norm_off = np.array(I1_I0_1_covariance_norm_off)
    I1_I0_2_covariance_norm_off = np.array(I1_I0_2_covariance_norm_off)
    I0_1_I0_2_covariance_norm_off = np.array(I0_1_I0_2_covariance_norm_off)
    N_values_off = np.array(N_values_off)

    #Calculate final TFY and its standard deviation according to error propagation rules
    dqdx = 1 / (I0_1_on + I0_2_on)
    dqdy = - I1_on / (I0_1_on + I0_2_on)**2
    dqdz = - I1_on / (I0_1_on + I0_2_on)**2
    TFY_on = I1_on / (I0_1_on + I0_2_on)
    TFY_on_std = (dqdx*I1_on_std)**2 + (dqdy*I0_1_on_std)**2 + (dqdz*I0_2_on_std)**2 + (2*dqdx*dqdy*I1_I0_1_covariance_on) + (2*dqdx*dqdz*I1_I0_2_covariance_on) + (2*dqdy*dqdz*I0_1_I0_2_covariance_on)
    TFY_on_std = np.sqrt(TFY_on_std)

    dqdx = 1 / (I0_1_off + I0_2_off)
    dqdy = - I1_off / (I0_1_off + I0_2_off)**2
    dqdz = - I1_off / (I0_1_off + I0_2_off)**2
    TFY_off = I1_off / (I0_1_off + I0_2_off)
    TFY_off_std = (dqdx*I1_off_std)**2 + (dqdy*I0_1_off_std)**2 + (dqdz*I0_2_off_std)**2 + (2*dqdx*dqdy*I1_I0_1_covariance_off) + (2*dqdx*dqdz*I1_I0_2_covariance_off) + (2*dqdy*dqdz*I0_1_I0_2_covariance_off)
    TFY_off_std = np.sqrt(TFY_off_std)

    Energy = MonochromatorAngles*(-6.33896967e-03) + 2.68470859e+04

    '''
    create a dictionary and put all the useful arrays into it, then save into a .csv file
    '''
    temp_dict=dict()
    temp_dict['Energy_off'] = Energy
    temp_dict['I1_off'] = I1_off
    temp_dict['I0_1_off'] = I0_1_off
    temp_dict['I0_2_off'] = I0_2_off
    temp_dict['I1_off_std'] = I1_off
    temp_dict['I0_1_off_std'] = I0_1_off
    temp_dict['I0_2_off_std'] = I0_2_off
    temp_dict['I1_I0_1_covariance_off'] = I1_I0_1_covariance_off
    temp_dict['I1_I0_2_covariance_off'] = I1_I0_2_covariance_off
    temp_dict['I0_1_I0_2_covariance_off'] = I0_1_I0_2_covariance_off
    temp_dict['I1_I0_1_covariance_norm_off'] = I1_I0_1_covariance_norm_off
    temp_dict['I1_I0_2_covariance_norm_off'] = I1_I0_2_covariance_norm_off
    temp_dict['I0_1_I0_2_covariance_norm_off'] = I0_1_I0_2_covariance_norm_off
    temp_dict['TFY_off'] = TFY_off
    temp_dict['TFY_off_std'] = TFY_off_std
    temp_dict['N_values_off'] = N_values_off
    
    temp_dict['Energy_on'] = Energy
    temp_dict['I1_on'] = I1_on
    temp_dict['I0_1_on'] = I0_1_on
    temp_dict['I0_2_on'] = I0_2_on
    temp_dict['I1_on_std'] = I1_on
    temp_dict['I0_1_on_std'] = I0_1_on
    temp_dict['I0_2_on_std'] = I0_2_on
    temp_dict['I1_I0_1_covariance_on'] = I1_I0_1_covariance_on
    temp_dict['I1_I0_2_covariance_on'] = I1_I0_2_covariance_on
    temp_dict['I0_1_I0_2_covariance_on'] = I0_1_I0_2_covariance_on
    temp_dict['I1_I0_1_covariance_norm_on'] = I1_I0_1_covariance_norm_on
    temp_dict['I1_I0_2_covariance_norm_on'] = I1_I0_2_covariance_norm_on
    temp_dict['I0_1_I0_2_covariance_norm_on'] = I0_1_I0_2_covariance_norm_on
    temp_dict['TFY_on'] = TFY_on
    temp_dict['TFY_on_std'] = TFY_on_std
    temp_dict['N_values_on'] = N_values_on
    
    #reindex randomly saved columns
    df = pd.DataFrame.from_dict(temp_dict,orient='index').transpose().fillna(' ')
    df = df[['Energy_off','I1_off','I0_1_off','I0_2_off','I1_off_std','I0_1_off_std','I0_2_off_std',
             'I1_I0_1_covariance_off','I1_I0_2_covariance_off','I0_1_I0_2_covariance_off',
             'I1_I0_1_covariance_norm_off','I1_I0_2_covariance_norm_off','I0_1_I0_2_covariance_norm_off',
             'TFY_off','TFY_off_std','N_values_off',
             'Energy_on','I1_on','I0_1_on','I0_2_on','I1_on_std','I0_1_on_std','I0_2_on_std',
             'I1_I0_1_covariance_on','I1_I0_2_covariance_on','I0_1_I0_2_covariance_on',
             'I1_I0_1_covariance_norm_on','I1_I0_2_covariance_norm_on','I0_1_I0_2_covariance_norm_on',
             'TFY_on','TFY_on_std','N_values_on']]
    df.to_csv(str(SaveFolder) + str(RunNumber) + '_' + str(ExtraComment) +'.csv', index=False)
    
    temp_dict2=dict()
    temp_dict2['Energy_off'] = np.array( list( filter( lambda x: x is not None, list( map( ( lambda x,y,z: x if (y<0.5 and z == Attenuation) else None), (MonochromatorAngle*(-6.33896967e-03) + 2.68470859e+04) ,LaserStatus,OpticalAttenuator ) ) ) ) )
    temp_dict2['I1_off'] = np.array( list( filter( lambda x: x is not None, list( map( ( lambda x,y,z: x if (y<0.5 and z == Attenuation) else None), I1,LaserStatus,OpticalAttenuator ) ) ) ) )
    temp_dict2['I0_1_off'] = np.array( list( filter( lambda x: x is not None, list( map( ( lambda x,y,z: x if (y<0.5 and z == Attenuation) else None), I0_1,LaserStatus,OpticalAttenuator ) ) ) ) )
    temp_dict2['I0_2_off'] = np.array( list( filter( lambda x: x is not None, list( map( ( lambda x,y,z: x if (y<0.5 and z == Attenuation) else None), I0_2,LaserStatus,OpticalAttenuator ) ) ) ) )

    temp_dict2['Energy_on'] = np.array( list( filter( lambda x: x is not None, list( map( ( lambda x,y,z: x if (y>0.5 and z == Attenuation) else None), (MonochromatorAngle*(-6.33896967e-03) + 2.68470859e+04) ,LaserStatus,OpticalAttenuator ) ) ) ) )
    temp_dict2['I1_on'] = np.array( list( filter( lambda x: x is not None, list( map( ( lambda x,y,z: x if (y>0.5 and z == Attenuation) else None), I1,LaserStatus,OpticalAttenuator ) ) ) ) )
    temp_dict2['I0_1_on'] = np.array( list( filter( lambda x: x is not None, list( map( ( lambda x,y,z: x if (y>0.5 and z == Attenuation) else None), I0_1,LaserStatus,OpticalAttenuator ) ) ) ) )
    temp_dict2['I0_2_on'] = np.array( list( filter( lambda x: x is not None, list( map( ( lambda x,y,z: x if (y>0.5 and z == Attenuation) else None), I0_2,LaserStatus,OpticalAttenuator ) ) ) ) )

    df = pd.DataFrame.from_dict(temp_dict2,orient='index').transpose().fillna(' ')
    df = df[['Energy_off','I1_off','I0_1_off','I0_2_off',
             'Energy_on','I1_on','I0_1_on','I0_2_on']]
    df.to_csv(str(SaveFolder) + str(RunNumber) + '_' + str(ExtraComment) +'_raw.csv', index=False)
    
    # plt.figure(1)
    # plt.errorbar( Energy, TFY_on, yerr=TFY_on_std )
    # plt.plot( MonochromatorAngles, TFY_on )
    # plt.plot( MonochromatorAngles, (I1_on/(I0_1_on)) )
    # plt.plot( MonochromatorAngles, (I1_on/(I0_2_on)) )
    # plt.xlabel('Energy (pls)')
    # plt.ylabel('Intensity (a.u.)')
    # plt.legend(['Average of 2 detectors', 'Detector 1 divided by 2', 'Detector 2 divided by 2'])


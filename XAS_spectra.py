#This code is for plotting XAS spectra from SACLA beamtime


def spectra_extractor(RunNumber, DataDirectory, SaveFolder, ExtraComment, FemtosecondInPls=6.671, 
                 TimeZero=1375, Attenuation=0):

    f = h5py.File(DataDirectory + 'XAS_' + str(RunNumber) + '.h5', 'r')
    I1 = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_13_in_volt'][:]
    I0_1 = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_14_in_volt'][:]
    I0_2 = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_15_in_volt'][:]
    I_OpticalLaser = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_4_in_volt'][:]
    OpticalAttenuator = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/eh_2_optical_ND_filter_stage_position'][:]
    MonochromatorAngle = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_1/eh_1_CC1_rotational_stage_position'][:]
    LaserStatus = f['/run_' + str(RunNumber) + '/event_info/bl_3/lh_1/laser_pulse_selector_status'][:]
    TagList_I = f['/run_' + str(RunNumber) + '/event_info/tag_number_list'][:]
    
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
                and (OpticalAttenuator[k] == Attenuation)
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
    I1_off = np.array(I1_off)
    I0_1_off = np.array(I0_1_off)
    I0_2_off = np.array(I0_2_off)

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

    # np.savetxt('dataaaa.csv', [p for p in zip((MonochromatorAngles), 
    #                                           (I1_On_ave/(I0_1_On_ave+I0_2_On_ave)), 
    #                                           (I1_Off_ave/(I0_1_Off_ave+I0_2_Off_ave)))], 
    #                                           delimiter=',')

    # plt.figure(1)
    # plt.errorbar( Energy, TFY_on, yerr=TFY_on_std )
    # plt.plot( MonochromatorAngles, TFY_on )
    # plt.plot( MonochromatorAngles, (I1_on/(I0_1_on)) )
    # plt.plot( MonochromatorAngles, (I1_on/(I0_2_on)) )
    # plt.xlabel('Energy (pls)')
    # plt.ylabel('Intensity (a.u.)')
    # plt.legend(['Average of 2 detectors', 'Detector 1 divided by 2', 'Detector 2 divided by 2'])


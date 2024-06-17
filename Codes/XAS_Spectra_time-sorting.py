#This code is for sorting and plotting XAS data from SACLA beamtime

import h5py
import numpy as np
import pandas as pd

def sorting_tool_spectra_time(RunNumber, DataDirectory, DataDirectoryTM, SaveFolder, ExtraComment, FemtosecondInPls=6.671, 
                 TimeZero=1375, BinSize=25, Bin_Threshold=0):
    
    I1 = np.empty(0)
    I0_1 = np.empty(0)
    I0_2 = np.empty(0)
    OpticalDelay = np.empty(0)
    MonochromatorAngle = np.empty(0)
    LaserStatus = np.empty(0)
    TM_data_reduced_pixels = np.empty(0)
        
    for k in RunNumber:
        #Import detector data etc from HDF5 file
        f = h5py.File(DataDirectory + 'XAS_' + str(RunNumber) + '.h5', 'r')
        I1t = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_13_in_volt'][:]
        I0_1t = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_14_in_volt'][:]
        I0_2t = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_15_in_volt'][:]
        I_OpticalLaser = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_4_in_volt'][:]    
        OpticalDelayt = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/eh_2_optical_delay_stage_position'][:]
        MonochromatorAnglet = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_1/eh_1_CC1_rotational_stage_position'][:]
        LaserStatust = f['/run_' + str(RunNumber) + '/event_info/bl_3/lh_1/laser_pulse_selector_status'][:]
        TagList_I = f['/run_' + str(RunNumber) + '/event_info/tag_number_list'][:]
        
        #Import Timing Tool data
        TM_data = np.genfromtxt(DataDirectoryTM + str(RunNumber) + '.csv', delimiter=',',skip_header=2)
        
        '''
        We need to find the tag nubmers of the points which both have a Photodiode measurement
        and a Timing tool corrected measurement. Here we reduce the TM_data table by selecting
        only the rows which have the corresponding tag in the Photodiode data
        '''
        TM_data_reduced_pixelst = [] #For filtered Timing Tool data, values are in pixels
        for fd in range(len(TagList_I)):
            Corresponding_TM_data_index = (np.where(TM_data[:,0]==TagList_I[fd]))[0]
            #Check if the data exists in the TM_data, if not, then the empty array is created and discarded
            DoesDataExist = True
            if not Corresponding_TM_data_index:
                DoesDataExist = False
            #Finds the row with the same tag value.
            if (DoesDataExist == True) and (TagList_I[fd] == TM_data[Corresponding_TM_data_index,0]): 
                TM_data_reduced_pixelst.append(TM_data[Corresponding_TM_data_index])
        
        TM_data_reduced_pixelst = np.array(TM_data_reduced_pixelst)
        TM_data_reduced_pixelst = np.squeeze(TM_data_reduced_pixelst,axis=1)
        
        I1 = np.hstack((I1,I1t))
        I0_1 = np.hstack((I0_1,I0_1t))
        I0_2 = np.hstack((I0_2,I0_2t))
        OpticalDelay = np.hstack((OpticalDelay,OpticalDelayt))
        MonochromatorAngle = np.hstack((MonochromatorAngle,MonochromatorAnglet))
        LaserStatus = np.hstack((LaserStatus,LaserStatust))
        TM_data_reduced_pixels = np.hstack((TM_data_reduced_pixels,TM_data_reduced_pixelst))
        
    MonochromatorAngles_unique = np.array(sorted(list(set(MonochromatorAngle)))) #Find out Angle values

    #This section applies the Timing Tool correction, and adds extra filtering conditions
    I1_off = []
    I1_on = []
    I0_1_off = []
    I0_1_on = []
    I0_2_off = []
    I0_2_on = []
    Delays_on = []
    Delays_off = []
    MonochromatorAngle_on = []
    MonochromatorAngle_off = []
    for sa in range(len(TM_data_reduced_pixels[:,0])):
        '''
        First condition filters out NaN values of timing tool data,
        which is in column 2, which is timing tool edge fit;
        change to column 1 ([sa,1]) to use timing tool edge derivative instead.
        Conditions 2, 3, 4 filter out NaN values of sample and reference photodiodes. 
        Conditions 5 collects the Laser Off status values.
        Conditions 6, 7, 8 select the phodotiode values within indicated range.
        '''
        if ((np.isnan(TM_data_reduced_pixels[sa,2]) == False)
            and (np.isnan(I1[sa]) == False)
            and (np.isnan(I0_1[sa]) == False)
            and (np.isnan(I0_2[sa]) == False)
            and (LaserStatus[sa] == 0)
            and (0.9>I1[sa]>0.01)
            and (0.9>I0_1[sa]>0.01)
            and (0.9>I0_2[sa]>0.01)):          
            I1_off.append(I1[sa])
            I0_1_off.append(I0_1[sa])
            I0_2_off.append(I0_2[sa])
            MonochromatorAngle_off.append(MonochromatorAngle[sa])
            '''
            Here we save the new Timing Tool corrected delays.
            Shift the original delay by expected TimeZero in Pls, then convert
            to femtoseconds. Timing Tool data is shifted by 1000 pixels so that
            the correction is somewhat neutral around 0, and converted to femtoseconds
            The addition and substraction signs are very important, the order shoulnd't
            be changed at the whim. A quote form Katayama et al. (2016):
            'The positive direction in the relative time indicates the earlier 
            irradiation of x-rays with respect to the optical lasers.'
            Also, since the grid is uneven, everythong above 1250 fs is not timing-tool corrected.
            '''
            Delays_off.append( ((OpticalDelay[sa]-TimeZero)*FemtosecondInPls) + 
                              ((1000-TM_data_reduced_pixels[sa,2])*2.6) )

        #Same as above, but Laser On
        elif ((np.isnan(TM_data_reduced_pixels[sa,2]) == False)
            and (np.isnan(I1[sa]) == False)
            and (np.isnan(I0_1[sa]) == False)
            and (np.isnan(I0_2[sa]) == False)
            and (LaserStatus[sa] == 1)
            and (0.9>I1[sa]>0.01)
            and (0.9>I0_1[sa]>0.01)
            and (0.9>I0_2[sa]>0.01)):
            I1_on.append(I1[sa])
            I0_1_on.append(I0_1[sa])
            I0_2_on.append(I0_2[sa])
            MonochromatorAngle_on.append(MonochromatorAngle[sa])
            Delays_on.append( ((OpticalDelay[sa]-TimeZero)*FemtosecondInPls) + 
                              ((1000-TM_data_reduced_pixels[sa,2])*2.6) )

    
    
# =============================================================================
#     #Sort by Delay and convert into numpy array
#     I1_off_time = np.array([p for p in sorted(zip(Delays_off, I1_off))])
#     I1_on_time = np.array([p for p in sorted(zip(Delays_on, I1_on))])
#     I0_1_off_time = np.array([p for p in sorted(zip(Delays_off, I0_1_off))])
#     I0_1_on_time = np.array([p for p in sorted(zip(Delays_on, I0_1_on))])
#     I0_2_off_time = np.array([p for p in sorted(zip(Delays_off, I0_2_off))])
#     I0_2_on_time = np.array([p for p in sorted(zip(Delays_on, I0_2_on))])
# =============================================================================
        
    #Now is the rebinning of the sorted, time-corrected data points
    I1_off_rebinned = []
    I0_1_off_rebinned = []
    I0_2_off_rebinned = []
    I1_off_rebinned_std = []
    I0_1_off_rebinned_std = []
    I0_2_off_rebinned_std = []
    I1_I0_1_covariance_off = []
    I1_I0_2_covariance_off = []
    I0_1_I0_2_covariance_off = []
    I1_I0_1_covariance_norm_off = []
    I1_I0_2_covariance_norm_off = []
    I0_1_I0_2_covariance_norm_off = []
    Delays_off_rebinned = []
    MonochromatorAngle_off_rebinned = []
    N_values_off = []

    I1_on_rebinned = []
    I0_1_on_rebinned = []
    I0_2_on_rebinned = []
    I1_on_rebinned_std = []
    I0_1_on_rebinned_std = []
    I0_2_on_rebinned_std = []
    I1_I0_1_covariance_on = []
    I1_I0_2_covariance_on = []
    I0_1_I0_2_covariance_on = []
    I1_I0_1_covariance_norm_on = []
    I1_I0_2_covariance_norm_on = []
    I0_1_I0_2_covariance_norm_on = []
    Delays_on_rebinned = []
    MonochromatorAngle_on_rebinned = []
    N_values_on = []
    
    Bins_on = list(range(-200, int(max(Delays_on))+BinSize, BinSize))
    Bins_off = Bins_on
    
    for gf in range(len(Bins_off)):
        N_values = 0
        I1_temp = np.empty(0)
        I0_1_temp = np.empty(0)
        I0_2_temp = np.empty(0)
        for hg in range(len(Delays_off)):
            if 0<=Delays_off[hg]-Bins_off[gf]<BinSize:
                N_values+= 1
                I1_temp = np.append(I1_temp,I1_off[hg])
                I0_1_temp = np.append(I0_1_temp,I0_1_off[hg])
                I0_2_temp = np.append(I0_2_temp,I0_2_off[hg])
        if N_values>Bin_Threshold:
            I1_off_rebinned.append(np.mean(I1_temp))
            I0_1_off_rebinned.append(np.mean(I0_1_temp))
            I0_2_off_rebinned.append(np.mean(I0_2_temp))
            I1_off_rebinned_std.append(np.std(I1_temp))
            I0_1_off_rebinned_std.append(np.std(I0_1_temp))
            I0_2_off_rebinned_std.append(np.std(I0_2_temp))
            I1_I0_1_covariance_off.append(np.cov(I1_temp,I0_1_temp,bias=True)[1,0])
            I1_I0_2_covariance_off.append(np.cov(I1_temp,I0_2_temp,bias=True)[1,0])
            I0_1_I0_2_covariance_off.append(np.cov(I0_1_temp,I0_2_temp,bias=True)[1,0])
            I1_I0_1_covariance_norm_off.append(np.corrcoef(I1_temp,I0_1_temp)[1,0])
            I1_I0_2_covariance_norm_off.append(np.corrcoef(I1_temp,I0_2_temp)[1,0])
            I0_1_I0_2_covariance_norm_off.append(np.corrcoef(I0_1_temp,I0_2_temp)[1,0])
            Delays_off_rebinned.append(Bins_off[gf])
            N_values_off.append(N_values)
            
    for gf in range(len(Bins_on)):
        N_values = 0
        I1_temp = []
        I0_1_temp = []
        I0_2_temp = []
        for hg in range(len(Delays_on)):
            if 0<=Delays_on[hg]-Bins_on[gf]<BinSize:
                N_values+= 1
                I1_temp.append(I1_on[hg])
                I0_1_temp.append(I0_1_on[hg])
                I0_2_temp.append(I0_2_on[hg])
        if N_values>Bin_Threshold:
            I1_on_rebinned.append(np.mean(I1_temp))
            I0_1_on_rebinned.append(np.mean(I0_1_temp))
            I0_2_on_rebinned.append(np.mean(I0_2_temp))
            I1_on_rebinned_std.append(np.std(I1_temp))
            I0_1_on_rebinned_std.append(np.std(I0_1_temp))
            I0_2_on_rebinned_std.append(np.std(I0_2_temp))
            I1_I0_1_covariance_on.append(np.cov(I1_temp,I0_1_temp,bias=True)[1,0])
            I1_I0_2_covariance_on.append(np.cov(I1_temp,I0_2_temp,bias=True)[1,0])
            I0_1_I0_2_covariance_on.append(np.cov(I0_1_temp,I0_2_temp,bias=True)[1,0])
            I1_I0_1_covariance_norm_on.append(np.corrcoef(I1_temp,I0_1_temp)[1,0])
            I1_I0_2_covariance_norm_on.append(np.corrcoef(I1_temp,I0_2_temp)[1,0])
            I0_1_I0_2_covariance_norm_on.append(np.corrcoef(I0_1_temp,I0_2_temp)[1,0])
            Delays_on_rebinned.append(Bins_on[gf])
            N_values_on.append(N_values)
    
    I1_off_rebinned = np.array(I1_off_rebinned)
    I0_1_off_rebinned = np.array(I0_1_off_rebinned)
    I0_2_off_rebinned = np.array(I0_2_off_rebinned) 
    I1_off_rebinned_std = np.array(I1_off_rebinned_std)
    I0_1_off_rebinned_std = np.array(I0_1_off_rebinned_std)
    I0_2_off_rebinned_std = np.array(I0_2_off_rebinned_std)
    I1_I0_1_covariance_off = np.array(I1_I0_1_covariance_off)
    I1_I0_2_covariance_off = np.array(I1_I0_2_covariance_off)
    I0_1_I0_2_covariance_off = np.array(I0_1_I0_2_covariance_off)
    I1_I0_1_covariance_norm_off = np.array(I1_I0_1_covariance_norm_off)
    I1_I0_2_covariance_norm_off = np.array(I1_I0_2_covariance_norm_off)
    I0_1_I0_2_covariance_norm_off = np.array(I0_1_I0_2_covariance_norm_off)
    
    I1_on_rebinned = np.array(I1_on_rebinned)
    I0_1_on_rebinned = np.array(I0_1_on_rebinned)
    I0_2_on_rebinned = np.array(I0_2_on_rebinned)
    I1_on_rebinned_std = np.array(I1_on_rebinned_std)
    I0_1_on_rebinned_std = np.array(I0_1_on_rebinned_std)
    I0_2_on_rebinned_std = np.array(I0_2_on_rebinned_std)
    I1_I0_1_covariance_on = np.array(I1_I0_1_covariance_on)
    I1_I0_2_covariance_on = np.array(I1_I0_2_covariance_on)
    I0_1_I0_2_covariance_on = np.array(I0_1_I0_2_covariance_on)
    I1_I0_1_covariance_norm_on = np.array(I1_I0_1_covariance_norm_on)
    I1_I0_2_covariance_norm_on = np.array(I1_I0_2_covariance_norm_on)
    I0_1_I0_2_covariance_norm_on = np.array(I0_1_I0_2_covariance_norm_on)
    
    Delays_off_rebinned = np.array(Delays_off_rebinned)
    Delays_on_rebinned = np.array(Delays_on_rebinned)
    N_values_off = np.array(N_values_off) 
    N_values_on = np.array(N_values_on) 
    
    #Calculate final TFY and its standard deviation according to error propagation rules
    dqdx = 1 / (I0_1_on_rebinned + I0_2_on_rebinned)
    dqdy = - I1_on_rebinned / (I0_1_on_rebinned + I0_2_on_rebinned)**2
    dqdz = - I1_on_rebinned / (I0_1_on_rebinned + I0_2_on_rebinned)**2
    TFY_on = I1_on_rebinned / (I0_1_on_rebinned + I0_2_on_rebinned)
    TFY_on_std = (dqdx*I1_on_rebinned_std)**2 + (dqdy*I0_1_on_rebinned_std)**2 + (dqdz*I0_2_on_rebinned_std)**2 + (2*dqdx*dqdy*I1_I0_1_covariance_on) + (2*dqdx*dqdz*I1_I0_2_covariance_on) + (2*dqdy*dqdz*I0_1_I0_2_covariance_on)
    TFY_on_std = np.sqrt(TFY_on_std)
    
    dqdx = 1 / (I0_1_off_rebinned + I0_2_off_rebinned)
    dqdy = - I1_off_rebinned / (I0_1_off_rebinned + I0_2_off_rebinned)**2
    dqdz = - I1_off_rebinned / (I0_1_off_rebinned + I0_2_off_rebinned)**2
    TFY_off = I1_off_rebinned / (I0_1_off_rebinned + I0_2_off_rebinned)
    TFY_off_std = (dqdx*I1_off_rebinned_std)**2 + (dqdy*I0_1_off_rebinned_std)**2 + (dqdz*I0_2_off_rebinned_std)**2 + (2*dqdx*dqdy*I1_I0_1_covariance_off) + (2*dqdx*dqdz*I1_I0_2_covariance_off) + (2*dqdy*dqdz*I0_1_I0_2_covariance_off)
    TFY_off_std = np.sqrt(TFY_off_std)
    
    
    '''
    create a dictionary and put all the useful arrays into it, then save into a .csv file
    '''
    temp_dict=dict()
    temp_dict['Delays_off'] = Delays_off_rebinned
    temp_dict['I1_off'] = I1_off_rebinned
    temp_dict['I0_1_off'] = I0_1_off_rebinned
    temp_dict['I0_2_off'] = I0_2_off_rebinned
    temp_dict['I1_off_std'] = I1_off_rebinned_std
    temp_dict['I0_1_off_std'] = I0_1_off_rebinned_std
    temp_dict['I0_2_off_std'] = I0_2_off_rebinned_std
    temp_dict['I1_I0_1_covariance_off'] = I1_I0_1_covariance_off
    temp_dict['I1_I0_2_covariance_off'] = I1_I0_2_covariance_off
    temp_dict['I0_1_I0_2_covariance_off'] = I0_1_I0_2_covariance_off
    temp_dict['I1_I0_1_covariance_norm_off'] = I1_I0_1_covariance_norm_off
    temp_dict['I1_I0_2_covariance_norm_off'] = I1_I0_2_covariance_norm_off
    temp_dict['I0_1_I0_2_covariance_norm_off'] = I0_1_I0_2_covariance_norm_off
    temp_dict['TFY_off'] = TFY_off
    temp_dict['TFY_off_std'] = TFY_off_std
    temp_dict['N_values_off'] = N_values_off
    
    temp_dict['Delays_on'] = Delays_on_rebinned
    temp_dict['I1_on'] = I1_on_rebinned
    temp_dict['I0_1_on'] = I0_1_on_rebinned
    temp_dict['I0_2_on'] = I0_2_on_rebinned
    temp_dict['I1_on_std'] = I1_on_rebinned_std
    temp_dict['I0_1_on_std'] = I0_1_on_rebinned_std
    temp_dict['I0_2_on_std'] = I0_2_on_rebinned_std
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
    df = df[['Delays_off','I1_off','I0_1_off','I0_2_off','I1_off_std','I0_1_off_std','I0_2_off_std',
             'I1_I0_1_covariance_off','I1_I0_2_covariance_off','I0_1_I0_2_covariance_off',
             'I1_I0_1_covariance_norm_off','I1_I0_2_covariance_norm_off','I0_1_I0_2_covariance_norm_off',
             'TFY_off','TFY_off_std','N_values_off',
             'Delays_on','I1_on','I0_1_on','I0_2_on','I1_on_std','I0_1_on_std','I0_2_on_std',
             'I1_I0_1_covariance_on','I1_I0_2_covariance_on','I0_1_I0_2_covariance_on',
             'I1_I0_1_covariance_norm_on','I1_I0_2_covariance_norm_on','I0_1_I0_2_covariance_norm_on',
             'TFY_on','TFY_on_std','N_values_on']]
    df.to_csv(str(SaveFolder) + str(RunNumber) + '_' + str(BinSize) + 'fs_' + str(ExtraComment) +'.csv', index=False)
    
    #same as previous step but now without summing the bins
    temp_dict2=dict()
    temp_dict2['Delays_off'] = Delays_off
    temp_dict2['I1_off'] = I1_off
    temp_dict2['I0_1_off'] = I0_1_off
    temp_dict2['I0_2_off'] = I0_2_off
    temp_dict2['Delays_on'] = Delays_on
    temp_dict2['I1_on'] = I1_on
    temp_dict2['I0_1_on'] = I0_1_on
    temp_dict2['I0_2_on'] = I0_2_on

    df2 = pd.DataFrame.from_dict(temp_dict2,orient='index').transpose().fillna(' ')
    df2 = df2[['Delays_off','I1_off','I0_1_off','I0_2_off','Delays_on','I1_on','I0_1_on','I0_2_on']]
    df2.to_csv(str(SaveFolder) + str(RunNumber) + '_' + str(BinSize) + 'fs_' + str(ExtraComment) +'_uncorrected.csv', index=False)


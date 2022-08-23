#This code is for sorting and plotting XAS data from SACLA beamtime, without Timing tool correction

import h5py
import matplotlib.pyplot as plt
import numpy as np

def sorting_tool_simple(RunNumber, DataDirectory, SaveFolder, ExtraComment, FemtosecondInPls=6.671, 
                 TimeZero=1375, Attenuation=0):


    f = h5py.File(DataDirectory + 'XAS_' + str(RunNumber) + '.h5', 'r')
    I1 = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_13_in_volt'][:]
    I0_1 = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_14_in_volt'][:]
    I0_2 = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_15_in_volt'][:]
    OpticalAttenuator = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/eh_2_optical_ND_filter_stage_position'][:]
    OpticalDelay = f['/run_' + str(RunNumber) + '/event_info/bl_3/eh_2/eh_2_optical_delay_stage_position'][:]
    LaserStatus = f['/run_' + str(RunNumber) + '/event_info/bl_3/lh_1/laser_pulse_selector_status'][:]
    
    Delays = np.array(sorted(list(set(OpticalDelay)))) #Find out delay values
    
    #Create empty arrays for detector values
    I1_On_ave = []
    I0_1_On_ave = []
    I0_2_On_ave = []
    I1_Off_ave = []
    I0_1_Off_ave = []
    I0_2_Off_ave = []
    for ds in range(len(Delays)):
        #sum over LaserON shots
        intensity1 = 0
        intensity2 = 0
        intensity3 = 0
        for k in range(len(LaserStatus)):
            if (OpticalDelay[k] == Delays[ds]) and (OpticalAttenuator[k] == Attenuation) and (LaserStatus[k] == 1):
                intensity1 = intensity1 + I1[k]
                intensity2 = intensity2 + I0_1[k]
                intensity3 = intensity3 + I0_2[k]
        I1_On_ave.append(intensity1)
        I0_1_On_ave.append(intensity2)
        I0_2_On_ave.append(intensity3)
        #sum over LaserOFF shots
        intensity01 = 0
        intensity02 = 0
        intensity03 = 0
        for l in range(len(LaserStatus)):
            if (OpticalDelay[l] == Delays[ds]) and (OpticalAttenuator[l] == Attenuation) and (LaserStatus[l] == 0):
                intensity01 = intensity01 + I1[l]
                intensity02 = intensity02 + I0_1[l]
                intensity03 = intensity03 + I0_2[l]
        I1_Off_ave.append(intensity01)
        I0_1_Off_ave.append(intensity02)
        I0_2_Off_ave.append(intensity03)
    
    #Convert into numpy arrays for usefulness
    I1_On_ave = np.array(I1_On_ave)
    I0_1_On_ave = np.array(I0_1_On_ave)
    I0_2_On_ave = np.array(I0_2_On_ave)
    I1_Off_ave = np.array(I1_Off_ave)
    I0_1_Off_ave = np.array(I0_1_Off_ave)
    I0_2_Off_ave = np.array(I0_2_Off_ave)
    
    np.savetxt(str(SaveFolder) + str(RunNumber) + '_Simple_' + str(ExtraComment) +'.csv',
                                             [p for p in zip(((Delays-TimeZero)*FemtosecondInPls), 
                                              I1_Off_ave,
                                              I0_1_Off_ave,
                                              I0_2_Off_ave,
                                              I1_On_ave,
                                              I0_1_On_ave,
                                              I0_2_On_ave)], 
                                              delimiter=',',
                                              header='Delays,I1_off_ave,I0_1_off_ave,I0_2_off_ave,I1_on_ave,I0_1_on_ave,I0_2_on_ave')
    

    '''
    plt.plot( ((Delays-TimeZero)*PlsInFemtosecond), (I1_On_ave/(I0_1_On_ave+I0_2_On_ave)) )
    plt.plot( ((Delays-TimeZero)*PlsInFemtosecond), (I1_On_ave/(I0_1_On_ave*2)) )
    plt.plot( ((Delays-TimeZero)*PlsInFemtosecond), (I1_On_ave/(I0_2_On_ave*2)) )
    plt.xlabel('Time (fs)')
    plt.ylabel('Intensity (a.u.)')
    plt.legend(['Average of 2 detectors', 'Detector 1 divided by 2', 'Detector 2 divided by 2'])
    
    plt.xlim([-600, 800])
    plt.ylim([0,60])
    '''

#This code is for plotting Intensity data from SACLA beamtime to see how the X-Ray power develops over time

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#Select the Run Number and the directory where the file is located
RunNumber = np.array(np.arange(1145439,1145444,1))
DataDirectory = 'C:/Users/nurekeye/Desktop/Projects/SACLA_05-2022/DataAnalysisTry/soap/XAS_Kinetics/'

#Set some initial parameters
I1_total = []
I0_1_total = []
I0_2_total = []
I_OL_total = []
OpticalAttenuator_total = []
LaserStatus_total = []
TagList_total = []
TagInitial = []

for i in RunNumber:
    #Import detector data etc from HDF5 file
    f = h5py.File(DataDirectory + str(i) + '.h5', 'r')
    I1 = f['/run_' + str(i) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_13_in_volt'][:]
    I0_1 = f['/run_' + str(i) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_14_in_volt'][:]
    I0_2 = f['/run_' + str(i) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_15_in_volt'][:]
    I_OL = f['/run_' + str(i) + '/event_info/bl_3/eh_2/photodiode/photodiode_user_4_in_volt'][:]
    OpticalAttenuator = f['/run_' + str(i) + '/event_info/bl_3/eh_2/eh_2_optical_ND_filter_stage_position'][:]
    LaserStatus = f['/run_' + str(i) + '/event_info/bl_3/lh_1/laser_pulse_selector_status'][:]
    TagList = f['/run_' + str(i) + '/event_info/tag_number_list'][:]
    
    I1_total += I1.tolist()
    I0_1_total += I0_1.tolist()
    I0_2_total += I0_2.tolist()
    I_OL_total += I_OL.tolist()
    OpticalAttenuator_total += OpticalAttenuator.tolist()
    LaserStatus_total += LaserStatus.tolist()
    TagList_total += TagList.tolist()
    TagInitial.append(TagList[0])
    
df = pd.DataFrame(list(zip(I1_total,I0_1_total,I0_2_total,I_OL_total,OpticalAttenuator_total,LaserStatus_total,TagList_total)),
                  columns=['I1','I0_1','I0_2','I_OL','OpticalAttenuator','LaserStatus','TagList'])
df.to_csv('IntensityOverTime',index=False)

#plots
#df2 = df.dropna()
#crap = df2['I0_1'].rolling(100).mean()
#for i in range(len(RunNumber)):
#    plt.axvline(x=TagInitial[i],color='k',linestyle='dotted')
#    plt.text(TagInitial[i],0.1,str(RunNumber[i]),rotation=90)

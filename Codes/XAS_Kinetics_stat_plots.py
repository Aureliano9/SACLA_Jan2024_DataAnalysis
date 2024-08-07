import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

RunNumber = np.array([1359438,1359439,1359440,1359441])
# RunNumber = np.array([1359438, 1359439, 1359440, 1359441, 1359442, 1359443, 1359444, 1359445, 1359447, 1359448, 
#                       1359449, 1359450, 1359451, 1359453, 1359454, 1359455]) #400 nm Kinetics
# RunNumber = np.array([1359527, 1359528, 1359529, 1359530, 1359531, 1359532, 1359533, 1359534, 1359535, 1359536, 
#                       1359537, 1359538, 1359539, 1359540, 1359541, 1359542, 1359543, 1359544]) #200 nm Kinetics

Directory = 'C:/Users/nurekeye/Desktop/Projects/SACLA_01-2024/DataAnalysis/2024_05_13/res/laser-off-non-binned/'
# Save_Directory = Directory

# Create empty dictionaries to fill them with Pandas DataFrames, each key corresponds to a run
data = {} #timing tool corrected, rebinned
datau = {} #timing tool corrected, but non-binned
datas = {} #without timing tool correction, binned
Runs = RunNumber

for i in Runs:
    # i=1145442
    data[str(i)] = pd.read_csv(Directory+str(i)+'_25fs_0.csv', delimiter=',').replace(' ',np.nan).astype(float)
    datau[str(i)] = pd.read_csv(Directory+str(i)+'_25fs_0_uncorrected.csv', delimiter=',').replace(' ',np.nan).astype(float)
    datas[str(i)] = pd.read_csv(Directory+str(i)+'_simple_0.csv', delimiter=',').replace(' ',np.nan).astype(float)

j=0
for k in Runs:
    os.mkdir(Directory+str(k)+'/')
    Save_Directory = Directory+str(k)+'/'
    #Detector vs Detector scatter plots, and histograms  
    plt.figure(1, figsize=[10,10])
    plt.subplot(2,2,1)
    plt.plot( datau[str(k)].I1_on, datau[str(k)].I0_1_on,'.')
    plt.plot( datau[str(k)].I1_on, datau[str(k)].I0_2_on,'.')
    plt.plot( datau[str(k)].I0_1_on, datau[str(k)].I0_2_on,'.')
    plt.xlabel('Intensity (a.u.)')
    plt.ylabel('Intensity (a.u.)')
    plt.legend(['I1_on vs I0_1_on', 'I1_on vs I0_2_on', 'I0_1_on vs I0_2_on'])
    plt.title('Laser ON, Run='+str(k))
    plt.subplot(2,2,2)
    plt.hist(datau[str(k)].I1_on,bins=100,histtype='step')
    plt.hist(datau[str(k)].I0_1_on,bins=100,histtype='step')
    plt.hist(datau[str(k)].I0_2_on,bins=100,histtype='step')
    plt.xlabel('Intensity (a.u.)')
    plt.ylabel('Number of points')
    plt.legend(['I1_on', 'I0_1_on', 'I0_2_on'])
    plt.title('Laser ON, Run='+str(k))
    plt.subplot(2,2,3)
    plt.plot( datau[str(k)].I1_off, datau[str(k)].I0_1_off,'.')
    plt.plot( datau[str(k)].I1_off, datau[str(k)].I0_2_off,'.')
    plt.plot( datau[str(k)].I0_1_off, datau[str(k)].I0_2_off,'.')
    plt.xlabel('Intensity (a.u.)')
    plt.ylabel('Intensity (a.u.)')
    plt.legend(['I1_off vs I0_1_off', 'I1_off vs I0_2_off', 'I0_1_off vs I0_2_off'])
    plt.title('Laser OFF, Run='+str(k))
    plt.subplot(2,2,4)
    plt.hist(datau[str(k)].I1_off,bins=100,histtype='step')
    plt.hist(datau[str(k)].I0_1_off,bins=100,histtype='step')
    plt.hist(datau[str(k)].I0_2_off,bins=100,histtype='step')
    plt.xlabel('Intensity (a.u.)')
    plt.ylabel('Number of points')
    plt.legend(['I1_off', 'I0_1_off', 'I0_2_off'])
    plt.title('Laser OFF, Run='+str(k))
    plt.savefig(Save_Directory+str(k)+'det_vs_det_non-binned.png',dpi=300)
    plt.close()
    
    #Separate detector statistic vs Delay, rebinned
    plt.figure(2, figsize=[15,15])
    plt.subplots_adjust(hspace=0.5, wspace=0.2)
    #Delay vs I1_on, I0_1_on and I0_1_on, rebinned
    plt.subplot(3,2,1)
    plt.plot( data[str(k)].Delays_on, data[str(k)].I1_on)
    plt.plot( data[str(k)].Delays_on, data[str(k)].I0_1_on)
    plt.plot( data[str(k)].Delays_on, data[str(k)].I0_2_on)
    plt.xlabel('Delays (fs)')
    plt.ylabel('Intensity (a.u.)')
    plt.legend(['I1_on', 'I0_1_on', 'I0_2_on'])
    plt.title('Separate detector signal, Laser ON, after rebinning')
    #Delay vs I1_off, I0_1_off and I0_1_off, rebinned
    plt.subplot(3,2,2)
    plt.plot( data[str(k)].Delays_off, data[str(k)].I1_off)
    plt.plot( data[str(k)].Delays_off, data[str(k)].I0_1_off)
    plt.plot( data[str(k)].Delays_off, data[str(k)].I0_2_off)
    plt.xlabel('Delays (fs)')
    plt.ylabel('Intensity (a.u.)')
    plt.legend(['I1_off', 'I0_1_off', 'I0_2_off'])
    plt.title('Separate detector signal, Laser OFF, after rebinning')
    #Delay vs standard deviations of I1_on, I0_1_on and I0_1_on, rebinned
    plt.subplot(3,2,3)
    plt.plot( data[str(k)].Delays_on, data[str(k)].I1_on_std)
    plt.plot( data[str(k)].Delays_on, data[str(k)].I0_1_on_std)
    plt.plot( data[str(k)].Delays_on, data[str(k)].I0_2_on_std)
    plt.xlabel('Delays (fs)')
    plt.ylabel('Intensity (a.u.)')
    plt.legend(['I1_on', 'I0_1_on', 'I0_2_on'])
    plt.title('Separate detector standard deviation, Laser ON, after rebinning')
    #Delay vs standard deviations of I1_off, I0_1_off and I0_1_off, rebinned
    plt.subplot(3,2,4)
    plt.plot( data[str(k)].Delays_off, data[str(k)].I1_off_std)
    plt.plot( data[str(k)].Delays_off, data[str(k)].I0_1_off_std)
    plt.plot( data[str(k)].Delays_off, data[str(k)].I0_2_off_std)
    plt.xlabel('Delays (fs)')
    plt.ylabel('Intensity (a.u.)')
    plt.legend(['I1_off', 'I0_1_off', 'I0_2_off'])
    plt.title('Separate detector standard deviation, Laser OFF, after rebinning')
    #Delay vs I1_on, I0_1_on and I0_1_on with error bars, rebinned
    plt.subplot(3,2,5)
    plt.errorbar( data[str(k)].Delays_on, data[str(k)].I1_on, yerr=data[str(k)].I1_on_std, capsize=5)
    plt.errorbar( data[str(k)].Delays_on, data[str(k)].I0_1_on, yerr=data[str(k)].I0_1_on_std, capsize=5)
    plt.errorbar( data[str(k)].Delays_on, data[str(k)].I0_2_on, yerr=data[str(k)].I0_2_on_std, capsize=5)
    plt.xlabel('Delays (fs)')
    plt.ylabel('Intensity (a.u.)')
    plt.legend(['I1_on', 'I0_1_on', 'I0_2_on'])
    plt.title('Separate detector signal with error bars, Laser ON, after rebinning')
    #Delay vs I1_off, I0_1_off and I0_1_off with error bars, rebinned
    plt.subplot(3,2,6)
    plt.errorbar( data[str(k)].Delays_off, data[str(k)].I1_off, yerr=data[str(k)].I1_off_std, capsize=5)
    plt.errorbar( data[str(k)].Delays_off, data[str(k)].I0_1_off, yerr=data[str(k)].I0_1_off_std, capsize=5)
    plt.errorbar( data[str(k)].Delays_off, data[str(k)].I0_2_off, yerr=data[str(k)].I0_2_off_std, capsize=5)
    plt.xlabel('Delays (fs)')
    plt.ylabel('Intensity (a.u.)')
    plt.legend(['I1_off', 'I0_1_off', 'I0_2_off'])
    plt.title('Separate detector signal with error bars, Laser OFF, after rebinning')
    plt.savefig(Save_Directory+str(k)+'det_vs_delay.png',dpi=300)
    plt.close()


    #Delay vs normalized I1_on with error bars, error propagation rules applied. First is standard deviation, second is standard error
    plt.figure(3, figsize=[14,5])
    plt.subplot(1,2,1)
    plt.errorbar( data[str(k)].Delays_on, data[str(k)].TFY_on, yerr=data[str(k)].TFY_on_std, capsize=5)
    plt.errorbar( data[str(k)].Delays_on, data[str(k)].TFY_on, yerr=data[str(k)].TFY_on_std/np.sqrt(data[str(k)].N_values_on), capsize=5)
    plt.xlabel('Delay (fs)')
    plt.ylabel('TFY (a.u.)')
    #plt.xlim([-600,1000])
    plt.legend(['Standard Deviation', 'Standard Error'])
    plt.title('Laser ON, Run='+str(k))
    plt.subplot(1,2,2)
    plt.errorbar( data[str(k)].Delays_off, data[str(k)].TFY_off, yerr=data[str(k)].TFY_off_std, capsize=5)
    plt.errorbar( data[str(k)].Delays_off, data[str(k)].TFY_off, yerr=data[str(k)].TFY_off_std/np.sqrt(data[str(k)].N_values_off), capsize=5)
    plt.xlabel('Delay (fs)')
    plt.ylabel('TFY (a.u.)')
    #plt.xlim([-600,1000])
    plt.legend(['Standard Deviation', 'Standard Error'])
    plt.title('Laser OFF, Run='+str(k))
    plt.savefig(Save_Directory+str(k)+'_kinetics_rebinned.png',dpi=300)
    plt.close()    
    
    #Delay vs normalized I1_on with error bars, error propagation rules applied. First is standard deviation, second is standard error
    plt.figure(4, figsize=[14,5])
    plt.subplot(1,2,1)
    plt.errorbar( data[str(k)].Delays_on, data[str(k)].TFY_on, yerr=data[str(k)].TFY_on_std, capsize=5)
    plt.errorbar( data[str(k)].Delays_on, data[str(k)].TFY_on, yerr=data[str(k)].TFY_on_std/np.sqrt(data[str(k)].N_values_on), capsize=5)
    plt.xlabel('Delay (fs)')
    plt.ylabel('TFY (a.u.)')
    plt.xlim([-600,1300])
    plt.legend(['Standard Deviation', 'Standard Error'])
    plt.title('Laser ON, Run='+str(k))
    plt.subplot(1,2,2)
    plt.errorbar( data[str(k)].Delays_off, data[str(k)].TFY_off, yerr=data[str(k)].TFY_off_std, capsize=5)
    plt.errorbar( data[str(k)].Delays_off, data[str(k)].TFY_off, yerr=data[str(k)].TFY_off_std/np.sqrt(data[str(k)].N_values_off), capsize=5)
    plt.xlabel('Delay (fs)')
    plt.ylabel('TFY (a.u.)')
    plt.xlim([-600,1300])
    plt.legend(['Standard Deviation', 'Standard Error'])
    plt.title('Laser OFF, Run='+str(k))
    plt.savefig(Save_Directory+str(k)+'_kinetics_rebinned_small-scale.png',dpi=300)
    plt.close() 
    
    #Delay vs normalized I1_on with error bars, error propagation rules applied. First is standard deviation, second is standard error
    plt.figure(5, figsize=[14,5])
    plt.subplot(1,2,1)
    plt.errorbar( data[str(k)].Delays_on, data[str(k)].TFY_on-data[str(k)].TFY_off, yerr=np.sqrt(data[str(k)].TFY_on_std**2+data[str(k)].TFY_off_std**2)/np.sqrt(2), capsize=5)
    plt.errorbar( data[str(k)].Delays_on, data[str(k)].TFY_on-data[str(k)].TFY_off, yerr=np.sqrt(data[str(k)].TFY_on_std**2+data[str(k)].TFY_off_std**2)/np.sqrt(2*data[str(k)].N_values_on), capsize=5)
    plt.xlabel('Delay (fs)')
    plt.ylabel('TFY (a.u.)')
    plt.xlim([-600,1300])
    plt.legend(['Standard Deviation', 'Standard Error'])
    plt.title('Laser ON - Laser OFF, Run='+str(k))
    plt.subplot(1,2,2)
    plt.errorbar( data[str(k)].Delays_on, data[str(k)].TFY_on-data[str(k)].TFY_off, yerr=np.sqrt(data[str(k)].TFY_on_std**2+data[str(k)].TFY_off_std**2)/np.sqrt(2), capsize=5)
    plt.errorbar( data[str(k)].Delays_on, data[str(k)].TFY_on-data[str(k)].TFY_off, yerr=np.sqrt(data[str(k)].TFY_on_std**2+data[str(k)].TFY_off_std**2)/np.sqrt(2*data[str(k)].N_values_on), capsize=5)
    plt.xlabel('Delay (fs)')
    plt.ylabel('TFY (a.u.)')
    #plt.xlim([-600,1000])
    plt.legend(['Standard Deviation', 'Standard Error'])
    plt.title('Laser ON - Laser OFF, Run='+str(k))
    plt.savefig(Save_Directory+str(k)+'_kinetics_rebinned_difference.png',dpi=300)
    plt.close() 
    
    #Separate detector statistic vs delay
    plt.figure(6, figsize=[16,12])
    plt.subplot(2,3,1)
    plt.scatter(datau[str(k)].Delays_on, datau[str(k)].I0_2_on,color='#1f77b4',marker='.')
    plt.title('I0_2_ON')
    plt.ylabel('Intensity (a.u.)')
    plt.xlabel('Delay (fs)')
    plt.xlim([-600,1300])
    plt.subplot(2,3,2)
    plt.scatter(datau[str(k)].Delays_on, datau[str(k)].I0_1_on,color='#ff7f0e',marker='.')
    plt.title('I0_1_ON')
    plt.ylabel('Intensity (a.u.)')
    plt.xlabel('Delay (fs)')
    plt.xlim([-600,1300])
    plt.subplot(2,3,3)
    plt.scatter(datau[str(k)].Delays_on, datau[str(k)].I1_on,color='#2ca02c',marker='.')
    plt.title('I1_ON')
    plt.ylabel('Intensity (a.u.)')
    plt.xlabel('Delay (fs)')
    plt.xlim([-600,1300])
    plt.subplot(2,3,4)
    plt.scatter(datau[str(k)].Delays_off, datau[str(k)].I0_2_off,color='#1f77b4',marker='.')
    plt.title('I0_2_OFF')
    plt.ylabel('Intensity (a.u.)')
    plt.xlabel('Delay (fs)')
    plt.xlim([-600,1300])
    plt.subplot(2,3,5)
    plt.scatter(datau[str(k)].Delays_off, datau[str(k)].I0_1_off,color='#ff7f0e',marker='.')
    plt.title('I0_1_OFF')
    plt.ylabel('Intensity (a.u.)')
    plt.xlabel('Delay (fs)')
    plt.xlim([-600,1300])
    plt.subplot(2,3,6)
    plt.scatter(datau[str(k)].Delays_off, datau[str(k)].I1_off,color='#2ca02c',marker='.')
    plt.title('I1_OFF')
    plt.ylabel('Intensity (a.u.)')
    plt.xlabel('Delay (fs)')
    plt.xlim([-600,1300])
    plt.savefig(Save_Directory+str(k)+'_scatter-of-each-detector_vs-delay.png',dpi=300)
    plt.close()
    
    plt.figure(7, figsize=[16,6])
    plt.subplot(1,2,1)
    plt.scatter(datau[str(k)].Delays_on, datau[str(k)].I1_on/(datau[str(k)].I0_1_on+datau[str(k)].I0_2_on),color='k',marker='.')
    plt.title('I1/(I0_1+I0_2)')
    plt.ylabel('Intensity (a.u.)')
    plt.xlabel('Delay (fs)')
    plt.xlim([-600,1300])
    plt.subplot(1,2,2)
    plt.errorbar(data[str(k)].Delays_on, data[str(k)].TFY_on, yerr=data[str(k)].TFY_on_std)
    plt.xlim([-600,1300])
    plt.title('Rebinned')
    plt.ylabel('Intensity (a.u.)')
    plt.xlabel('Delay (fs)')
    plt.savefig(Save_Directory+str(k)+'_normalized-scatter.png',dpi=300)
    plt.close()

    #Correlation between signal and number of shots
    figs, axes = plt.subplots(3,1, num=k, sharex=True)
    axes[0].scatter(data[str(k)].Delays_on, (data[str(k)].TFY_on)/data[str(k)].TFY_on_std)
    axes[0].set_title('Mean/Standard Deviation')
    axes[1].scatter(data[str(k)].Delays_on, np.sqrt(data[str(k)].N_values_on), color='r')
    axes[1].set_title('sqrt(Shots)')
    axes[2].scatter(data[str(k)].Delays_on, ((data[str(k)].TFY_on)/data[str(k)].TFY_on_std) / np.sqrt(data[str(k)].N_values_on), color='g')
    axes[2].set_title('(Mean/Standard Deviation) / (sqrt(Shots))')
    plt.savefig(Save_Directory+str(k)+'_photon-stat.png',dpi=300)
    plt.close()

    plt.figure(8, figsize=[10,5])
    plt.errorbar(data[str(k)].Delays_on, data[str(k)].TFY_on, yerr=data[str(k)].TFY_on_std)
    plt.errorbar(datas[str(k)].Delays_on, datas[str(k)].TFY_on, yerr=datas[str(k)].TFY_on_std)
    plt.xlim([-600,1300])
    plt.title('Comparison of Kinetics with and without Timing Tool correction')
    plt.ylabel('Intensity (a.u.)')
    plt.xlabel('Delay (fs)')
    plt.legend(['with Timing Tool correction', 'without Timing Tool correction'])
    plt.savefig(Save_Directory+str(k)+'_TMA-correctd_vs_TMA-noncorrected.png',dpi=300)
    plt.close()
    
    plt.figure(9, figsize=[10,5])
    plt.plot(data[str(k)].Delays_on, data[str(k)].TFY_on)
    plt.plot(datas[str(k)].Delays_on, datas[str(k)].TFY_on)
    plt.xlim([-600,1300])
    plt.title('Comparison of Kinetics with and without Timing Tool correction')
    plt.ylabel('Intensity (a.u.)')
    plt.xlabel('Delay (fs)')
    plt.legend(['with Timing Tool correction', 'without Timing Tool correction'])
    plt.savefig(Save_Directory+str(k)+'_TMA-correctd_vs_TMA-noncorrected2.png',dpi=300)
    plt.close()
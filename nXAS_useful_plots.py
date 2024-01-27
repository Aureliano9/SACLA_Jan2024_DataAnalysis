import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

Directory ='C:/Users/nurekeye/Desktop/Projects/SACLA_05-2022/DataAnalysis/2023_06_21/'
Save_Directory ='C:/Users/nurekeye/Desktop/Projects/SACLA_01-2024/plot/'

Runs = np.array([1145439])

# Create empty dictionaries to fill them with Pandas DataFrames, each key corresponds to a run
data = {} #timing tool corrected
datau = {} #uncorrected

for i in Runs:
    # i=1145442
    data[str(i)] = pd.read_csv(Directory+str(i)+'_25fs_0.csv', delimiter=',').replace(' ',np.nan).astype(float)
    datau[str(i)] = pd.read_csv(Directory+str(i)+'_25fs_0_uncorrected.csv', delimiter=',').replace(' ',np.nan).astype(float)

j=0
for k in Runs:
    #Detector vs Detector scatter plots, and histograms  
    plt.figure(1, figsize=[10,10])
    plt.subplot(2,2,1)
    plt.plot( datau[str(k)].I1_on, datau[str(k)].I0_1_on,'.')
    plt.plot( datau[str(k)].I1_on, datau[str(k)].I0_2_on,'.')
    plt.plot( datau[str(k)].I0_1_on, datau[str(k)].I0_2_on,'.')
    plt.xlabel('Intensity (a.u.)')
    plt.ylabel('Intensity (a.u.)')
    plt.legend(['I1_on vs I0_1_on', 'I1_on vs I0_1_on', 'I0_1_on vs I0_2_on'])
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
    plt.legend(['I1_off vs I0_1_off', 'I1_off vs I0_1_off', 'I0_1_off vs I0_2_off'])
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
    plt.figure(1, figsize=[15,15])
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
    plt.figure(1, figsize=[14,5])
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
    
    #Separate detector statistic vs delay
    plt.figure(1, figsize=[16,6])
    plt.subplot(1,3,1)
    plt.scatter(datau[str(k)].Delays_on, datau[str(k)].I0_2_on,color='#1f77b4',marker='.')
    plt.title('I0_2')
    plt.ylabel('Intensity (a.u.)')
    plt.xlabel('Delay (fs)')
    plt.xlim([-600,1300])
    plt.subplot(1,3,2)
    plt.scatter(datau[str(k)].Delays_on, datau[str(k)].I0_1_on,color='#ff7f0e',marker='.')
    plt.title('I0_1')
    plt.ylabel('Intensity (a.u.)')
    plt.xlabel('Delay (fs)')
    plt.xlim([-600,1300])
    plt.subplot(1,3,3)
    plt.scatter(datau[str(k)].Delays_on, datau[str(k)].I1_on,color='#2ca02c',marker='.')
    plt.title('I1')
    plt.ylabel('Intensity (a.u.)')
    plt.xlabel('Delay (fs)')
    plt.xlim([-600,1300])
    plt.savefig(Save_Directory+str(k)+'_scatter-of-each-detector_vs-delay.png',dpi=300)
    plt.close()
    
    plt.figure(2, figsize=[16,6])
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




'''
#Delay vs I1_off, I0_1_off and I0_1_off, rebinned
plt.plot( Delays_off_time_rebinned, I1_off_time_rebinned)
plt.plot( Delays_off_time_rebinned, I0_1_off_time_rebinned)
plt.plot( Delays_off_time_rebinned, I0_2_off_time_rebinned)
plt.legend(['I1', 'I0_1', 'I0_2'])

#Delay vs standard deviations of I1_off, I0_1_off and I0_1_off, rebinned
plt.plot( Delays_off_time_rebinned, I1_off_time_rebinned_std)
plt.plot( Delays_off_time_rebinned, I0_1_off_time_rebinned_std)
plt.plot( Delays_off_time_rebinned, I0_2_off_time_rebinned_std)
plt.legend(['I1', 'I0_1', 'I0_2'])

#Delay vs I1_off, I0_1_off and I0_1_off with error bars, rebinned
plt.errorbar( Delays_off_time_rebinned, I1_off_time_rebinned, yerr=I1_off_time_rebinned_std, uplims=True, lolims=True)
plt.errorbar( Delays_off_time_rebinned, I0_1_off_time_rebinned, yerr=I0_1_off_time_rebinned_std, uplims=True, lolims=True)
plt.errorbar( Delays_off_time_rebinned, I0_2_off_time_rebinned, yerr=I0_2_off_time_rebinned_std, uplims=True, lolims=True)
plt.legend(['I1', 'I0_1', 'I0_2'])

#Delay vs I1_on, I0_1_on and I0_1_on, rebinned
plt.plot( Delays_on_time_rebinned, I1_on_time_rebinned)
plt.plot( Delays_on_time_rebinned, I0_1_on_time_rebinned)
plt.plot( Delays_on_time_rebinned, I0_2_on_time_rebinned)
plt.legend(['I1', 'I0_1', 'I0_2'])

#Delay vs standard deviations of I1_on, I0_1_on and I0_1_on, rebinned
plt.plot( Delays_on_time_rebinned, I1_on_time_rebinned_std)
plt.plot( Delays_on_time_rebinned, I0_1_on_time_rebinned_std)
plt.plot( Delays_on_time_rebinned, I0_2_on_time_rebinned_std)
plt.legend(['I1', 'I0_1', 'I0_2'])

#Delay vs I1_on, I0_1_on and I0_1_on with error bars, rebinned
plt.errorbar( Delays_on_time_rebinned, I1_on_time_rebinned, yerr=I1_on_time_rebinned_std, uplims=True, lolims=True)
plt.errorbar( Delays_on_time_rebinned, I0_1_on_time_rebinned, yerr=I0_1_on_time_rebinned_std, uplims=True, lolims=True)
plt.errorbar( Delays_on_time_rebinned, I0_2_on_time_rebinned, yerr=I0_2_on_time_rebinned_std, uplims=True, lolims=True)
plt.legend(['I1', 'I0_1', 'I0_2'])

#Delay vs normalized I1_on, I1_off and their difference, rebinned
plt.plot( Delays_on_time_rebinned, (I1_on_time_rebinned/(I0_1_on_time_rebinned+I0_2_on_time_rebinned)) )
plt.plot( Delays_off_time_rebinned, (I1_off_time_rebinned/(I0_1_off_time_rebinned+I0_2_off_time_rebinned)) )
plt.plot( Delays_on_time_rebinned[0:85], (I1_on_time_rebinned[0:85]/(I0_1_on_time_rebinned[0:85]+I0_2_on_time_rebinned[0:85]))-(I1_off_time_rebinned[0:85]/(I0_1_off_time_rebinned[0:85]+I0_2_off_time_rebinned[0:85])) )

#Delay vs normalized I1_on with error bars, error propagation rules applied, rebinned
I1 = I1_on_time_rebinned/(I0_1_on_time_rebinned+I0_2_on_time_rebinned)
I1_std = np.sqrt( (I1_on_time_rebinned_std/(I0_1_on_time_rebinned+I0_2_on_time_rebinned))**2 + (-I1_on_time_rebinned*I0_1_on_time_rebinned_std/((I0_1_on_time_rebinned+I0_2_on_time_rebinned)**2))**2 + (-I1_on_time_rebinned*I0_2_on_time_rebinned_std/((I0_1_on_time_rebinned+I0_2_on_time_rebinned)**2))**2 )
plt.errorbar( Delays_on_time_rebinned, I1, yerr=I1_std, uplims=True, lolims=True)

#I1_on vs I0_1_on, I0_2_on
plt.plot( I1_on_time_rebinned, I0_1_on_time_rebinned,'.')
plt.plot( I1_on_time_rebinned, I0_2_on_time_rebinned,'.')
'''

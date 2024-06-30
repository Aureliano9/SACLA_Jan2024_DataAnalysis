import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

# RunNumber = np.array([1359403]) #400 nm -3300 pls
# RunNumber = np.array([1359466, 1359467, 1359468, 1359469, 1359470, 1359471, 1359472, 1359473, 1359474, 1359475, 1359476]) #400 nm -501 pls
# RunNumber = np.array([1359404, 1359405, 1359406, 1359407, 1359408, 1359409, 1359410, 1359411, 1359412, 1359413, 1359414, 1359415, 1359416, 1359417]) #400 nm -495 pls
# RunNumber = np.array([1359462, 1359463, 1359464, 1359465]) #400 nm -487 pls
# RunNumber = np.array([1359395, 1359396, 1359397, 1359398, 1359399, 1359400, 1359401, 1359402,1359430, 1359431, 1359432, 1359433, 1359434, 1359435]) #400 nm -480 pls (should be max signal)
# RunNumber = np.array([1359477, 1359478, 1359479, 1359480, 1359481, 1359482, 1359483, 1359484, 1359485, 1359486]) #400 nm -435 pls
# RunNumber = np.array([1359419, 1359420, 1359421, 1359422, 1359424, 1359425, 1359426, 1359427, 1359428]) #400 nm -342 pls

# RunNumber = np.array([1359503, 1359504, 1359505, 1359506, 1359507, 1359508, 1359509, 1359510, 1359511, 1359512, 1359513, 1359514]) #200 nm 0 fs
# RunNumber = np.array([1359495, 1359496, 1359497, 1359498, 1359499, 1359500, 1359501]) #200 nm 100 fs
# RunNumber = np.array([1359516, 1359517, 1359518, 1359519, 1359520, 1359521, 1359522, 1359523, 1359524, 1359525, 1359526]) #200 nm 1000 fs
# RunNumber = np.array([1359556]) #200 nm 100 ps
# RunNumber = np.array([1359554, 1359555]) #200 nm 300 ps
RunNumber = np.array([1359545, 1359546, 1359547, 1359548, 1359549, 1359550, 1359551, 1359552, 1359553]) #200 nm 800 ps

# RunNumber = np.array([1359527])

Directory = 'C:/Users/nurekeye/Desktop/Projects/SACLA_01-2024/DataAnalysis/2024_06_17/Spectra_200nm/'
Save_Directory = Directory
ExtraComment = '200nm_800ps'
# Create empty dictionaries to fill them with Pandas DataFrames, each key corresponds to a run
data = {} #timing tool corrected, rebinned
datau = {} #timing tool corrected, but non-binned
datas = {} #without timing tool correction, binned
Runs = RunNumber

for i in Runs:
    data[str(i)] = pd.read_csv(Directory+str(i)+'_0.csv', delimiter=',').replace(' ',np.nan).astype(float)

Energy_on = []
TFY_on = []
TFY_on_std = []
N_values_on = []
Energy_off = []
TFY_off = []
TFY_off_std = []
N_values_off = []

for k in Runs:
    Energy_on = Energy_on + list(data[str(k)].Energy_on)  
Energy_on = np.unique(np.array(Energy_on))
Energy_on = Energy_on[~np.isnan(Energy_on)]

for k in Runs:
    Energy_off = Energy_off + list(data[str(k)].Energy_off)  
Energy_off = np.unique(np.array(Energy_off))
Energy_off = Energy_off[~np.isnan(Energy_off)]

for j in Energy_on:
    TFY_on_temp = []
    TFY_on_std_temp = []
    N_values_temp = []
    for h in Runs:
        if j in data[str(h)].Energy_on.values:
            ind = list(data[str(h)].Energy_on).index(j)
            TFY_on_temp.append(data[str(h)].TFY_on[ind])
            TFY_on_std_temp.append(data[str(h)].TFY_on_std[ind])
            N_values_temp.append(data[str(h)].N_values_on[ind])
    TFY_on_temp = np.array(TFY_on_temp)
    TFY_on_std_temp = np.array(TFY_on_std_temp/np.sqrt(N_values_temp))
    N_values_temp = np.array(N_values_temp)
    
    TFY_on.append(sum(TFY_on_temp*N_values_temp)/sum(N_values_temp))
    TFY_on_std.append( np.sqrt( (sum(N_values_temp*TFY_on_std_temp**2)/sum(N_values_temp)) + np.std(TFY_on_std_temp)**2) )
    N_values_on.append(sum(N_values_temp))
    
for j in Energy_off:
    TFY_off_temp = []
    TFY_off_std_temp = []
    N_values_temp = []
    for h in Runs:
        if j in data[str(h)].Energy_off.values:
            ind = list(data[str(h)].Energy_off).index(j)
            TFY_off_temp.append(data[str(h)].TFY_off[ind])
            TFY_off_std_temp.append(data[str(h)].TFY_off_std[ind])
            N_values_temp.append(data[str(h)].N_values_off[ind])
    TFY_off_temp = np.array(TFY_off_temp)
    TFY_off_std_temp = np.array(TFY_off_std_temp/np.sqrt(N_values_temp))
    N_values_temp = np.array(N_values_temp)
    
    TFY_off.append(sum(TFY_off_temp*N_values_temp)/sum(N_values_temp))
    TFY_off_std.append( np.sqrt( (sum(N_values_temp*TFY_off_std_temp**2)/sum(N_values_temp)) + np.std(TFY_off_std_temp)**2) )
    N_values_off.append(sum(N_values_temp))    
    
temp_dict=dict()
temp_dict['Energy_on'] = Energy_on
temp_dict['TFY_on'] = TFY_on
temp_dict['TFY_on_std'] = TFY_on_std
temp_dict['N_values_on'] = N_values_on
temp_dict['Energy_off'] = Energy_off
temp_dict['TFY_off'] = TFY_off
temp_dict['TFY_off_std'] = TFY_off_std
temp_dict['N_values_off'] = N_values_off

df = pd.DataFrame.from_dict(temp_dict,orient='index').transpose().fillna(' ')
df = df[['Energy_on','TFY_on','TFY_on_std','N_values_on','Energy_off','TFY_off','TFY_off_std','N_values_off']]
# df.to_csv(str(Save_Directory) + str(ExtraComment) + '_averaged.csv', index=False)

plt.errorbar(Energy_on, TFY_on, TFY_on_std)
plt.errorbar(Energy_off, TFY_off, TFY_off_std)

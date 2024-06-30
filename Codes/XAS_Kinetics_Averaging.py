import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

RunNumber = np.array([1359438, 1359439, 1359440, 1359441]) #400 nm Lorentzian Kinetics
# RunNumber = np.array([1359442, 1359443, 1359444, 1359445]) #400 nm Edge Kinetics
# RunNumber = np.array([1359447, 1359448, 1359449, 1359450]) #400 nm Absorption Max Kinetics
# RunNumber = np.array([1359451, 1359453, 1359454, 1359455]) #400 nm Edge DIff Max Kinetics

# RunNumber = np.array([1359527, 1359528, 1359529, 1359530, 1359539, 1359540, 1359541, 1359542, 1359543, 1359544]) #200 nm Lorentzian Kinetics
# RunNumber = np.array([1359531, 1359532, 1359533, 1359534]) #200 nm Edge Kinetics
# RunNumber = np.array([1359535, 1359536, 1359537, 1359538]) #200 nm Edge Diff Max Kinetics

# RunNumber = np.array([1359527])

Directory = 'C:/Users/nurekeye/Desktop/Projects/SACLA_01-2024/DataAnalysis/2024_06_17/Kinetics_400nm/'
Save_Directory = Directory
ExtraComment = 'Lorentzian_400nm'
# Create empty dictionaries to fill them with Pandas DataFrames, each key corresponds to a run
data = {} #timing tool corrected, rebinned
datau = {} #timing tool corrected, but non-binned
datas = {} #without timing tool correction, binned
Runs = RunNumber

for i in Runs:
    data[str(i)] = pd.read_csv(Directory+str(i)+'_25fs_0.csv', delimiter=',').replace(' ',np.nan).astype(float)

Delays = []
TFY = []
TFY_std = []
N_values = []
for k in Runs:
    Delays = Delays + list(data[str(k)].Delays_on)  
Delays = np.unique(np.array(Delays))
Delays = Delays[~np.isnan(Delays)]

for j in Delays:
    TFY_temp = []
    TFY_std_temp = []
    N_values_temp = []
    for h in Runs:
        if j in data[str(h)].Delays_on.values:
            ind = list(data[str(h)].Delays_on).index(j)
            # TFY_temp.append(data[str(h)].TFY_on[ind])
            # TFY_std_temp.append(data[str(h)].TFY_on_std[ind])
            # N_values_temp.append(data[str(h)].N_values_on[ind])
            
            TFY_temp.append(data[str(h)].TFY_on[ind]-data[str(h)].TFY_off[ind])
            TFY_std_temp.append(np.sqrt(data[str(h)].TFY_on_std[ind]**2+data[str(h)].TFY_off_std[ind]**2))
            N_values_temp.append(data[str(h)].N_values_on[ind]+data[str(h)].N_values_off[ind])
            
    TFY_temp = np.array(TFY_temp)
    TFY_std_temp = np.array(TFY_std_temp/np.sqrt(N_values_temp))
    N_values_temp = np.array(N_values_temp)
    
    TFY.append(sum(TFY_temp*N_values_temp)/sum(N_values_temp))
    TFY_std.append( np.sqrt( (sum(N_values_temp*TFY_std_temp**2)/sum(N_values_temp)) + np.std(TFY_std_temp)**2) )
    N_values.append(sum(N_values_temp))
    
temp_dict=dict()
temp_dict['Delays_on'] = Delays
temp_dict['TFY_on'] = TFY
temp_dict['TFY_on_std'] = TFY_std
temp_dict['N_values_on'] = N_values


df = pd.DataFrame.from_dict(temp_dict,orient='index').transpose().fillna(' ')
df = df[['Delays_on','TFY_on','TFY_on_std','N_values_on']]
# df.to_csv(str(Save_Directory) + str(ExtraComment) + '_averaged.csv', index=False)

plt.errorbar(Delays, TFY, TFY_std)

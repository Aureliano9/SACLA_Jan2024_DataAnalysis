#Input Parameters
import time
start_time = time.time()

from XAS_Kinetics_sorting import sorting_tool
from XAS_Kinetics_sorting import sorting_tool_laser_off_is_nonbinned
from XAS_Kinetics_sorting import sorting_tool_simple
import numpy as np
import os

# Select the Run Number, in a loop if needed. When using range function, last point isn't included, but first point is
RunNumber = np.array([1359439,1359440,1359441])
# RunNumber = np.array([1359438, 1359439, 1359440, 1359441, 1359442, 1359443, 1359444, 1359445, 1359447, 1359448, 
#                       1359449, 1359450, 1359451, 1359453, 1359454, 1359455]) #400 nm Kinetics
# RunNumber = np.array([1359527, 1359528, 1359529, 1359530, 1359531, 1359532, 1359533, 1359534, 1359535, 1359536, 
#                       1359537, 1359538, 1359539, 1359540, 1359541, 1359542, 1359543, 1359544]) #200 nm Kinetics

#Select the directory where the file is located, as well as TM data and saveing folder
DataDirectory = 'C:/Users/nurekeye/Desktop/Projects/SACLA_01-2024/DataAnalysis/2024_05_13/'
DataDirectoryTM = 'C:/Users/nurekeye/Desktop/Projects/SACLA_01-2024/DataAnalysis/2024_05_13/'
SaveFolder = 'C:/Users/nurekeye/Desktop/Projects/SACLA_01-2024/DataAnalysis/2024_05_13/res/'
# os.mkdir(SaveFolder)

#Parameters for sorting and rebinning
FemtosecondInPls = 6.671 #Femtosends in 1 pls
TimeZero = -480 #estimated Time zero in pls
BinSize = 25 #Set the bin size for rebinning in fs
Bin_Threshold = 0 #Threshold for selecting the good bins. If the number of points in that bin is less than Bin_Threshold+1, then that bin is rejected
ExtraComment = '0' #Add comment at the end of the generated file



for k in RunNumber:
    # sorting_tool(k, DataDirectory, DataDirectoryTM, SaveFolder, ExtraComment, FemtosecondInPls, TimeZero, BinSize, Bin_Threshold)
    sorting_tool_laser_off_is_nonbinned(k, DataDirectory, DataDirectoryTM, SaveFolder, ExtraComment, FemtosecondInPls, TimeZero, BinSize, Bin_Threshold)
    sorting_tool_simple(k, DataDirectory, SaveFolder, ExtraComment, FemtosecondInPls, TimeZero)
print('done')
print("--- %s seconds ---" % (time.time() - start_time))

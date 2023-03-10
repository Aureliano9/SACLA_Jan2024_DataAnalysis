#Input Parameters
import time
start_time = time.time()

import XAS_sorting_tool
import XAS_sorting_tool_simple
import numpy as np

#Select the directory where the file is located, as well as TM data
DataDirectory = '/home/nurekeyev/nurekeyev/2022A8049_backup/process_1/Merge/'
DataDirectoryTM = '/home/nurekeyev/nurekeyev/'
SaveFolder = '/home/nurekeyev/nurekeyev/SavedResults/2022_08_24/'

#ExtraComment = 'test' #Add comment at the end of the generated file

FemtosecondInPls = 6.671 #Femtosends in 1 pls
TimeZero = 1375 #estimated Time zero in pls

BinSize = 25 #Set the bin size for rebinning in fs

#Threshold for selecting the good bins. 
#If the number of points in that bin is less than Bin_Threshold+1, then that bin is rejected
Bin_Threshold = 0 

#Select an optical attenuator value. 0 for 5 uJ of laser energy, -2500 for 2.9 uJ
Attenuation = 0
# Select the Run Number, in a loop if needed. WHen using range function, last point isn't included, but first point is
RunNumber = np.array([1145374,1145375,1145376,1145377,1145378,1145381,1145382,1145383,1145384,
		      1145385,1145387,1145388,1145411,1145436,1145437,1145438,1145439,1145440,1145441,1145442,1145443,1145444,1145445,1145446]) #5 uJ runs

RunNumber = np.array([1145374,1145375,1145376,1145377,1145378,1145379,1145380,1145381,1145382,1145383,1145384,1145385,
		      1145387,1145389,1145409,1145410,1145436,1145437,1145438,1145439,1145440,1145441,1145442,1145443,1145444,1145445,1145446]) #2.9 uJ runs

RunNumber = np.array(range(1145436,1145448,1))
#RunNumber = np.array([1145436])
for k in RunNumber:
    ExtraComment = Attenuation
    XAS_sorting_tool.sorting_tool(k, DataDirectory, DataDirectoryTM, SaveFolder, ExtraComment, FemtosecondInPls, TimeZero, Attenuation, BinSize, Bin_Threshold)
    XAS_sorting_tool_simple.sorting_tool_simple(k, DataDirectory, SaveFolder, ExtraComment, FemtosecondInPls, TimeZero, Attenuation)

#Full possible inputs:
#XAS_sorting_tool.sorting_tool(RunNumber,DataDirectory,DataDirectoryTM,SaveFolder,ExtraComment,FemtosecondInPls,TimeZero,Attenuation,BinSize,Bin_Threshold)
#XAS_sorting_tool_simple.sorting_tool_simple(RunNumber, DataDirectory, SaveFolder, ExtraComment, FemtosecondInPls, TimeZero, Attenuation)

print('done')
print("--- %s seconds ---" % (time.time() - start_time))

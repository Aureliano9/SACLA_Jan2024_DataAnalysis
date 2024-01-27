#Input Parameters
import time
start_time = time.time()

import nXAS_sorting_tool
import numpy as np

#Select the directory where the file is located, as well as TM data and saveing folder
DataDirectory = '/home/nurekeyev/nurekeyev/2022A8049_backup/process_1/Merge/'
DataDirectoryTM = '/home/nurekeyev/nurekeyev/'
SaveFolder = '/home/nurekeyev/nurekeyev/SavedResults/2022_08_24/'

#Parameters for sorting and rebinning
FemtosecondInPls = 6.671 #Femtosends in 1 pls
TimeZero = 1375 #estimated Time zero in pls
BinSize = 25 #Set the bin size for rebinning in fs
Bin_Threshold = 0 #Threshold for selecting the good bins. If the number of points in that bin is less than Bin_Threshold+1, then that bin is rejected
ExtraComment = '' #Add comment at the end of the generated file

# Select the Run Number, in a loop if needed. When using range function, last point isn't included, but first point is
RunNumber = np.array([1145436])

for k in RunNumber:
    XAS_sorting_tool.sorting_tool(k, DataDirectory, DataDirectoryTM, SaveFolder, ExtraComment, FemtosecondInPls, TimeZero, BinSize, Bin_Threshold)

print('done')
print("--- %s seconds ---" % (time.time() - start_time))

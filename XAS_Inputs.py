#Input Parameters
import time
start_time = time.time()

import sorting_tool
import numpy as np

#Select the directory where the file is located, as well as TM data
DataDirectory = '/home/nurekeyev/nurekeyev/2022A8049_backup/process_1/Merge/'
DataDirectoryTM = '/home/nurekeyev/nurekeyev/'
SaveFolder = '/home/nurekeyev/nurekeyev/SavedResults/2022_08_23/'

ExtraComment = 'test' #Add comment at the end of the generated file

#Set some initial parameters
FemtosecondInPls = 6.671 #Femtosends in 1 pls
TimeZero = 1375 #estimated Time zero in pls

BinsSize = [10,20,30,40,50] #Set the bin size for rebinning in fs
BinsSize = np.array(BinsSize)
#Threshold for selecting the good bins. 
#If the number of points in that bin is less than Bin_Threshold+1, then that bin is rejected
Bin_Threshold = 0 

#Select an optical attenuator value, 0 for 0%, -2500 for 50%
Attenuation = 0
# Select the Run Number, in a loop if needed
RunNumber = 1145445

for i in range(len(BinsSize)):
	BinSize = BinsSize[i]
	XAS_sorting_tool.sorting_tool(RunNumber,DataDirectory,DataDirectoryTM,SaveFolder,ExtraComment,FemtosecondInPls,TimeZero,Attenuation,BinSize,Bin_Threshold)

print('done')
print("--- %s seconds ---" % (time.time() - start_time))

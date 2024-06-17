#Input Parameters
import time
start_time = time.time()
from XAS_Spectra_sorting import spectra_extractor
import numpy as np


# Select the Run Number, in a loop if needed. When using range function, last point isn't included, but first point is
#RunNumber = np.array([1359395])
# RunNumber = np.array([1359395, 1359396, 1359397, 1359398, 1359399, 1359400, 1359401, 1359402, 1359403, 1359404, 
#                       1359405, 1359406, 1359407, 1359408, 1359409, 1359410, 1359411, 1359412, 1359413, 1359414, 
#                       1359415, 1359416, 1359417, 1359419, 1359420, 1359421, 1359422, 1359424, 1359425, 1359426, 
#                       1359427, 1359428, 1359430, 1359431, 1359432, 1359433, 1359434, 1359435, 1359462, 1359463, 
#                       1359464, 1359465, 1359466, 1359467, 1359468, 1359469, 1359470, 1359471, 1359472, 1359473, 
#                       1359474, 1359475, 1359476, 1359477, 1359478, 1359479, 1359480, 1359481, 1359482, 1359483, 
#                       1359484, 1359485, 1359486]) #400 nm Spectra

RunNumber = np.array([1359495, 1359496, 1359497, 1359498, 1359499, 1359500, 1359501, 1359503, 1359504, 1359505, 
                      1359506, 1359507, 1359508, 1359509, 1359510, 1359511, 1359512, 1359513, 1359514, 1359516, 
                      1359517, 1359518, 1359519, 1359520, 1359521, 1359522, 1359523, 1359524, 1359525, 1359526,
                      1359545, 1359546, 1359547, 1359548, 1359549, 1359550, 1359551, 1359552, 1359553, 1359554, 
                      1359555, 1359556]) #200 nm Spectra

#Select the directory where the file is located
DataDirectory = 'C:/Users/nurekeye/Desktop/Projects/SACLA_01-2024/DataAnalysis/Codes/'
SaveFolder = 'C:/Users/nurekeye/Desktop/Projects/SACLA_01-2024/DataAnalysis/Codes/'


ExtraComment = '0' #Add comment at the end of the generated file

for k in RunNumber:
	spectra_extractor(k, DataDirectory, SaveFolder, ExtraComment)

print('done')
print("--- %s seconds ---" % (time.time() - start_time))


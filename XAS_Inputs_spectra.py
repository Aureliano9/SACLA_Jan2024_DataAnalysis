#Input Parameters
import time
start_time = time.time()
from XAS_spectra import spectra_extractor
import numpy as np



#Select the directory where the file is located
DataDirectory = '/home/nurekeyev/nurekeyev/2022A8049_backup/process_1/Merge/'
SaveFolder = '/home/nurekeyev/nurekeyev/SavedResults/2022_08_24/'

#ExtraComment = 'test' #Add comment at the end of the generated file

FemtosecondInPls = 6.671 #Femtosends in 1 pls
TimeZero = 1375 #estimated Time zero in pls

#Select an optical attenuator value. 0 for 5 uJ of laser energy, -2500 for 2.9 uJ
Attenuation = 0
# Select the Run Number, in a loop if needed. WHen using range function, last point isn't included, but first point is
RunNumber = np.array([1145366,1145367,1145390,1145392,1145394,1145395,1145396,1145398,
		      1145403,1145417,1145418,1145423,1145424,1145425,1145427,1145430,1145431,1145434]) #5 uJ runs
RunNumber = np.array([1145372,1145373,1145391,1145393,1145397,1145398,1145399,1145404,
		      1145405,1145406,1145407,1145408,1145419,1145420,1145421,1145422,1145426,1145432,1145433]) #2.9 uJ runs
#RunNumber = np.array(range(1145436,1145448,1))
#RunNumber = np.array([1145436])
for k in RunNumber:
	ExtraComment = Attenuation
	spectra_extractor(k, DataDirectory, SaveFolder, ExtraComment, FemtosecondInPls, TimeZero, Attenuation)

print('done')
print("--- %s seconds ---" % (time.time() - start_time))

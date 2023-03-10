#Input Parameters
import time
start_time = time.time()

import spectra_extractor
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
from collections import defaultdict


#Select the directory where the file is located
DataDirectory = '/home/nurekeyev/nurekeyev/2022A8049_backup/process_1/Merge/'
SaveFolder = '/home/nurekeyev/nurekeyev/SavedResults/2022_08_24/'

#ExtraComment = 'test' #Add comment at the end of the generated file

FemtosecondInPls = 6.671 #Femtosends in 1 pls
TimeZero = 1375 #estimated Time zero in pls

#Select an optical attenuator value. 0 for 5 uJ of laser energy, -2500 for 2.9 uJ
Attenuation = np.array([0, -2500])
# Select the Run Number, in a loop if needed. WHen using range function, last point isn't included, but first point is
RunNumber = np.array(range(1145436,1145448,1))
#RunNumber = np.array([1145436])
for k in RunNumber:
	for i in Attenuation:
		ExtraComment = i
		spectra_extractor(k, DataDirectory, SaveFolder, ExtraComment, FemtosecondInPls=6.671, TimeZero=1375, Attenuation=i)

print('done')
print("--- %s seconds ---" % (time.time() - start_time))

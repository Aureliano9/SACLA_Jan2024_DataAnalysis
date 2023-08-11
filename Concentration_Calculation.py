# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 16:23:11 2023

@author: nurekeye
"""
import numpy as np
import matplotlib.pyplot as plt


BeamWaistO = 85 #optical beam diameter in um
BeamWaistX = 5 #X-Ray beam diameter in um
PulseEnergyO = 5 #Optical laser pulse energy in uJ
PulseEnergyX = 300 #X-Ray laser pulse energy in uJ
JetThickness = 80 #jet thickness in um
Concentration = 0.0352 #np.linspace(0.001, 0.05, 50) #Concentraion of solute in M
ExtinctionCoef = 10000 #Extinction coefficient at optical laser excitation in L*mol^-1*cm^-1
AbsCoefXMass = 147 #X-Ray Mass Absorption Coefficient in cm^2/g
MolMass = 79.904 * (1.66 * 10**-24) * (6.022 * 10 **23); #Molar Mass in g/mol
PhotonWavelengthO = 201 #Optical laser photon wavelength in nm
PhotonEnergyX = 13.5 #X-Ray laser photon wavelength in keV

ConcFor1OD = 1*10**4 / (ExtinctionCoef*JetThickness) #What should be the concentration so that OD=1
PhotonEnergyO = 6.626 * (10**-17) * 3 / PhotonWavelengthO #Optical laser photon energy in J
PhotonCountO = (10**-6) * PulseEnergyO / PhotonEnergyO #Number of Optical photons per pulse
AbsorptionO = ExtinctionCoef*Concentration*JetThickness*10**-4 #Absorption of optical photons
SoluteCountInVolumeO = (JetThickness*BeamWaistO**2*np.pi/4)*Concentration*6.022*10**8 #Number of solute particles in an excitation volume
ExcitedStateFraction = (PhotonCountO-(PhotonCountO/10**AbsorptionO))/SoluteCountInVolumeO #Fraction of Excited state Solute particles after Optical laser

PhotonCountX = 10**-6*PulseEnergyX/(PhotonEnergyX*10**3*1.602*10**-19) #Number of X-Ray Photons
AbsCrossSectionX = AbsCoefXMass*MolMass/(6.022*10**23) #X-Ray Absorption Cross section in cm^2, otherwise known as sigma in ln(I1/I0)=n*sigma*d
SolutePartsPerVol = Concentration*6.022*10**17 #Number of Solute particles in mm^-3
AbsorptionX = AbsCrossSectionX*SolutePartsPerVol*JetThickness*10**-1 #Absorption of X-Ray photons
XRayPhotonsAbsorbed = PhotonCountX-PhotonCountX/np.exp(1)**AbsorptionX #Number of Absorbed X-Ray photons after Optcail laser
PhotonsOnDetectorX = XRayPhotonsAbsorbed*0.5*10**-3

AbsCrossSectionO = ExtinctionCoef * np.log(10) * 1000 / (6.022*10**23) #Optical Absorption Cross section in cm^2, otherwise known as sigma in ln(I1/I0)=n*sigma*d
PhotonCountO_in_cm2 = (10**8) * 4 * PhotonCountO / (np.pi * BeamWaistO**2) #Number of Optical photons per cm^2
res = AbsCrossSectionO * PhotonCountO_in_cm2
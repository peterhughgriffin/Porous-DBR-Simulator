# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 17:39:58 2019

Example file demonstrating the TMM_Class to simulate the reflectivity of porous DBRs

@author: Peter Griffin
"""
from TMM_Class import DBR
import matplotlib.pyplot as plt
import numpy as np

#%%############################  Setup ############################
#Porosity values at different polarisations
Phi_45 = 0.4
Phi_135 = 0.51

# Set Phi_90 as the average of the two porosity values
Phi_90 = Phi_45 +Phi_135/2

# Number of top layers (With no polarisation effect)
Bot_Layers = 7
# Total number of DBR repeats
Nlayers = 12


#%%############################ Build structures ############################

## Grading parameters
NGrades = 11
Factor = 1
Order = 1/8

# Refractive index
n_Space = 1 # Can change if pores do not contain air
# GaN Template thickness
T_Temp = 3400
# Sapphire substrate refractive index (assumed infinite thickness)
n_Sub = 1.76

path = "Polarisation/"
n_File = "GaN_RefractiveIndex_Barker-o"

Period = 97.3
Phi = 0.37
T_Rat = 0.345

All45 = DBR("All_45",path, Period,T_Rat,Phi_45,Nlayers,T_Temp)
All45.MakeGrades(NGrades,Factor,Order,Phi_45)

All135 = DBR("All_135",path, Period,T_Rat,Phi_135,Nlayers,T_Temp)
All135.MakeGrades(NGrades,Factor,Order,Phi_135)

Bot45 = DBR("Bot_45",path, Period,T_Rat,Phi_45,Bot_Layers,T_Temp)
Bot45.MakeGrades(NGrades,Factor,Order,Phi_45)

Bot135 = DBR("Bot_135",path, Period,T_Rat,Phi_135,Bot_Layers,T_Temp)
Bot135.MakeGrades(NGrades,Factor,Order,Phi_135)

Top90 = DBR("Top_90",path, Period,T_Rat,Phi_90,Nlayers-Bot_Layers,0)
Top90.MakeGrades(NGrades,Factor,Order,Phi_90)

#%%############################ Run sims ############################

# Using basic one DBR structure funbction: Simulate

All45.Simulate(n_Sub,n_Space,n_File)
All135.Simulate(n_Sub,n_Space,n_File)

# Using SimulateParts to simulate a structure where the DBR structure changes
Bot45.SimulateParts(n_Sub,n_Space,n_File,Top90)
Bot135.SimulateParts(n_Sub,n_Space,n_File,Top90)


#%%############################ Plot data ###########################

# Plot raw reflectivity

# Initialise variable to create separate plots
try:
    index=index+1
except:
    index=1

plt.figure(index)


Sim = All45
plt.plot(Sim.Wav,Sim.Rnorm, label = Sim.label)

Sim = All135
plt.plot(Sim.Wav,Sim.Rnorm, label = Sim.label)

Sim = Bot45
plt.plot(Sim.Wav,Sim.Rnorm, label = Sim.label)

Sim = Bot135
plt.plot(Sim.Wav,Sim.Rnorm, label = Sim.label)


plt.ylim([0,1])
plt.xlim([350,650])
plt.legend(prop={'size': 25})

#%% Plot Difference
index=index+1
plt.figure(index)

AllDiff=[]
for a,b in zip(All135.Rnorm,All45.Rnorm):
    AllDiff.append(a-b)

plt.plot(All135.Wav,AllDiff, label = "All")

BotDiff=[]
for a,b in zip(Bot135.Rnorm,Bot45.Rnorm):
    BotDiff.append(a-b)
    
plt.plot(Bot135.Wav,BotDiff, label = Bot_Layers)


plt.ylim([0,1])
plt.xlim([350,650])
plt.legend(prop={'size': 25})

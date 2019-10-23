# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 10:27:06 2019

This code creates a class DBR that represents a reflection experiment on a porous DBR.
The reflection is modelled using the transfer matrix model given structural information about the DBR.
The model is capable of simulating graded interfaces between the layers.

MakeGrades calculates the layer structure to simulate graded interfaces for a layer with given porosity and thickness
Simulate calculates the normal reflectivity of the structure
WriteData produces a csv file of the simulated data.
It uses the tmm package by Steven Byrnes: https://pypi.org/project/tmm/

@author: Peter Griffin
"""
from __future__ import division, print_function, absolute_import
from tmm import coh_tmm
import numpy as np
from numpy import linspace, inf
from scipy.interpolate import interp1d

def porosity_to_n(porosity,GaN_n,air_n):
    """Convert a porosity to a refractive index. using the volume averaging theory"""
    porous_n = np.sqrt((1-porosity)*GaN_n*GaN_n + porosity*air_n*air_n)
    return porous_n

class DBR:
    """Main DBR class"""
    def __init__(self,label,path, Period,T_Rat,Phi,NLayers,T_Temp):
        self.path = path                 # The path that the data will be written to
        self.label = label               # A useful label for this structure
        self.Period = Period             # The period thickness of the DBR
        self.T_Rat = T_Rat               # The thickness ratio of the porous and non-porous layers
        self.T_Por = T_Rat*Period        # The porous layer thickness
        self.T_GaN = Period-self.T_Por   # The non-porous layer thickness
        self.Phi = Phi                   # The overall average porosity of the porous layers
        self.NLayers = NLayers           # The number of DBR pairs (i.e. NLayers= 10 means 10 porous layers and 10 non-porous layers)
        self.T_Temp = T_Temp             # The thickness of the underlying GaN template
        self.nPor=[]                     # A list for storing the refractive indexes of the porous layers, which are graded

        
    def MakeGrades(self,NGrades,Factor,Order,Phi):
        """This function calculates multiple graded layers for each porous layer with an overall average porosity given by Phi"""
        # Check number of Grades is odd
        if NGrades%2==0:
            raise ValueError("NGrades must be odd!")
        
        halfN = int((NGrades+1)/2)
        # Make porosity of graded layers        
        Por_Grad1 = np.fromfunction(lambda i,j: 1+Factor*j**Order,(1,halfN),dtype=int)
        Por_Grad2 = np.fromfunction(lambda i,j: 1+Factor*(halfN-j-2)**Order,(1,halfN-1),dtype=int)
        
        Por_Grad = np.concatenate((Por_Grad1, Por_Grad2), axis=None)
        
        # Adjust the graded n to give correct overall porosity
        a = np.sum(Por_Grad)/(NGrades*self.Phi)
        self.Por_Graded = Por_Grad/a

        # Make thickness of graded layers
        T_Grade = self.T_Por/NGrades
        self.T_Graded = [T_Grade]*NGrades
        
        
    def Simulate(self, n_Sub,n_Space,n_File,n = 2.38):
        """This function simulates the reflectivity of the DBR structure using Rnorm from the TMM package"""
        # n_File is a string that is the filename of one of the included refractive index files for GaN. If it's false then a constant value is used.
        
        # Run the sims
        print(self.Period)
        
        # Using constant value for n_GaN
        if n_File==False:
            # Initialise Wavelength range to model (Resolution is a nm)
            self.Wav=linspace(200,1000,800)
        
            # Make Layer Thickness List
            Repeat = [self.T_GaN] + self.T_Graded
            self.d_list = [inf]+ Repeat*self.NLayers + [self.T_Temp] + [inf]
            
            # initialize lists of y-values to plot
            self.Rnorm=[] 
    
            ############ Calculations ############
            
            n_Por = porosity_to_n(self.Por_Graded,n,n_Space).tolist()
            # Make list of refractive indices
            Repeat = [n] + n_Por
            self.n_list = [1]+ Repeat*self.NLayers +[n] + [n_Sub]
            
            for Lambda in self.Wav:
                # For normal incidence, s and p polarizations are identical.
                # I arbitrarily decided to use 's'.
                self.Rnorm.append(coh_tmm('s',self.n_list, self.d_list, 0, Lambda)['R'])
                
        # Using data file for n_GaN
        else:
            # Load GaN refractive index data
            self.n_File = n_File
            data = np.loadtxt(self.n_File+'.csv', delimiter=',',skiprows=1)
                                        
            [self.Wav, n_GaN] = np.transpose(data) # Transpose data and assign the two columns
        
            self.Wav=self.Wav*1000 # Convert Wavelength to nm (from um)
            
            n_func = interp1d(self.Wav, n_GaN) # This interpolation function allows the look up of n for any wavelength
        
            # Initialise Wavelength range to model (Resolution is a nm)
            self.Wav=linspace(min(self.Wav),min(1000,self.Wav[-1]),num=int((1000-min(self.Wav))))
        
            # Make Layer Thickness List
            Repeat = [self.T_GaN] + self.T_Graded
            self.d_list = [inf]+ Repeat*self.NLayers + [self.T_Temp] + [inf]
            
                # initialize lists of y-values to plot
            self.Rnorm=[] 
        
            ############ Calculations ############
            
            for Lambda in self.Wav:
                n = n_func(Lambda).tolist()
                
                n_Por = porosity_to_n(self.Por_Graded,n,n_Space).tolist()
                self.nPor.append(n_Por)
                # Make list of refractive indices
                Repeat = [n] + n_Por
                self.n_list = [1]+ Repeat*self.NLayers +[n] + [n_Sub]
                
                # For normal incidence, s and p polarizations are identical.
                # I arbitrarily decided to use 's'.
                self.Rnorm.append(coh_tmm('s',self.n_list, self.d_list, 0, Lambda)['R'])
    
    def WriteData(self):
        """This function writes the calculated DBR to a .csv file"""
        
        Header = ',%s %d %% %.2f nm %.3f Ratio\n%d pair DBR\n%.2f nm porous layer\n%.2f nm non-porous layer\n%d %% Porosity\n\n'%(self.label,self.Phi*100, self.Period, self.T_Rat, self.NLayers, self.T_Por, self.T_GaN, self.Phi*100)
        FilePath = '%sTMM_%s_%dPr_%dnm_%d-%d_%dPc.csv'%(self.path,self.label,self.NLayers, self.Period, self.T_Por, self.T_GaN, self.Phi*100)
        Output=np.column_stack((self.Wav, self.Rnorm))
        np.savetxt(FilePath,Output, header= 'TMM Simulation result\n'+Header+'Wavelength,Reflectance\nnm,', comments='',delimiter = ',')
    

    def SimulateParts(self, n_Sub,n_Space,n_File,TopDBR,n = 2.38):
        """This function is an alternative to SImulate for structures made up of two different DBRs. Takes another DBR structure as an additional argument and calculates the reflectivity of the combined structure"""
        # Run the sims
        print(self.Period)
    
        # Using constant value for n_GaN
        if n_File==False:
            # Initialise Wavelength range to model (Resolution is a nm)
            self.Wav=linspace(200,1000,800)
        
            # Make Layer Thickness List
            Repeat = [self.T_GaN] + self.T_Graded
            RepeatTop = [TopDBR.T_GaN] + TopDBR.T_Graded
            self.d_list = [inf]+ RepeatTop*TopDBR.NLayers + Repeat*self.NLayers + [self.T_Temp] + [inf]
            
            # initialize lists of y-values to plot
            self.Rnorm=[] 
    
            ############ Calculations ############
            
            n_Por = porosity_to_n(self.Por_Graded,n,n_Space).tolist()
            n_PorTop = porosity_to_n(TopDBR.Por_Graded,n,n_Space).tolist()
            # Make list of refractive indices
            Repeat = [n] + n_Por
            RepeatTop = [n] + n_PorTop
            self.n_list = [1] + RepeatTop*TopDBR.NLayers + Repeat*self.NLayers +[n] + [n_Sub]
            
            for Lambda in self.Wav:
                # For normal incidence, s and p polarizations are identical.
                # I arbitrarily decided to use 's'.
                self.Rnorm.append(coh_tmm('s',self.n_list, self.d_list, 0, Lambda)['R'])
    
        # Using data file for n_GaN
        else:
            # Load GaN refractive index data
            self.n_File = n_File
            data = np.loadtxt(self.n_File+'.csv', delimiter=',',skiprows=1)
                
            [self.Wav, n_GaN] = np.transpose(data) # Transpose data and assign the two columns
        
            self.Wav=self.Wav*1000 # Convert Wavelength to nm (from um)
            
            n_func = interp1d(self.Wav, n_GaN) # This interpolation function allows the look up of n for any wavelength
        
            # Initialise Wavelength range to model (Resolution is a nm)
            self.Wav=linspace(min(self.Wav),min(1000,self.Wav[-1]),num=int((1000-min(self.Wav))))
        
            # Make Layer Thickness List
            Repeat = [self.T_GaN] + self.T_Graded
            RepeatTop = [TopDBR.T_GaN] + TopDBR.T_Graded
            self.d_list = [inf]+ RepeatTop*TopDBR.NLayers + Repeat*self.NLayers + [self.T_Temp] + [inf]
            
                # initialize lists of y-values to plot
            self.Rnorm=[] 
        
            ############ Calculations ############
            
            for Lambda in self.Wav:
                n = n_func(Lambda).tolist()
                
                n_Por = porosity_to_n(self.Por_Graded,n,n_Space).tolist()
                self.nPor.append(n_Por)
                n_PorTop = porosity_to_n(TopDBR.Por_Graded,n,n_Space).tolist()
                TopDBR.nPor.append(n_PorTop)
                # Make list of refractive indices
                Repeat = [n] + n_Por
                RepeatTop = [n] + n_PorTop
                self.n_list = [1] + RepeatTop*TopDBR.NLayers + Repeat*self.NLayers +[n] + [n_Sub]
                
                # For normal incidence, s and p polarizations are identical.
                # I arbitrarily decided to use 's'.
                self.Rnorm.append(coh_tmm('s',self.n_list, self.d_list, 0, Lambda)['R'])
    


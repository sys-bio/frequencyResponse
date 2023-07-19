# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 11:51:33 2023

@author: hsauro
"""

import roadrunner
import numpy as np
import math
import matplotlib.pyplot as plt

class FreqencyResponse:
    
    def __init__(self, r):
        self.r = r
    
    def getLogSpace (self, startW, numDecades, numPoints):
        d1 = 0
        y = np.zeros (numPoints)
        for i in range (numPoints-1):
            y[i] = i*(numDecades-d1)
            y[i] = y[i]/(numPoints-1)
            y[i] = d1 + y[i]
            y[i] = (10**y[i]) * startW
        y[numPoints-1] = (10**numDecades)*startW
        return y
    
    def getPhase (self, val):
        if (val.real == 0.0) and (val.imag == 0.0):
           return 0.0
        else:
           return math.atan2(val.imag, val.real)


    # Implement instead:
    # This should be able to handle models with conserved cycles
    #   H(jw) = (jw I - Nr dv/dx L )^{-1} Nr dv/dp
    def getFrequencyResponse(self, startFrequency, numberOfDecades, numberOfPoints,
                      parameterId, variableId, useDB=True, useHz = True):
    
        numReactions = self.r.getNumReactions()
        numFloatingSpecies = self.r.getNumFloatingSpecies()
        reactionIds = self.r.getReactionIds()
        speciesIds = self.r.getFloatingSpeciesIds()
        boundarySpeciesIds = self.r.getBoundarySpeciesIds()
        globalParameterIds = self.r.getGlobalParameterIds()
    
        if variableId not in speciesIds:
           raise Exception ("Can't find species:" + variableId + " in the model")     
               
        # Check that parameter name exists
        if parameterId not in boundarySpeciesIds:
           if parameterId not in globalParameterIds:
               raise Exception ("Can't find parameter:" + variableName + " in the model")     
               
        try:
           self.r.steadyState()
        except Exception as e:
            raise Exception ("Can't find the steady state: " + str (e))
               
        # Get the link matrix in case there are conserved moieties,
        # equals the identity matrix if not
        linkMatrix = self.r.getLinkMatrix()
        # Convert Jacobian matrix to complex matrix
        Jac = self.r.getFullJacobian()
        Jac = np.array (Jac)
        Jac = Jac.astype (complex)
        # Compute Jac*linkMatrix
        JacL = np.matmul (Jac, linkMatrix)
    
        # Will need ee for the flux TFs mwhen implemented
        #ee = self.r.getUnscaledElasticityMatrix()
        # Convert stoich matrix to complex matrix
        Nr = self.r.getReducedStoichiometryMatrix()
        Nr = np.array (Nr)
        Nr = Nr.astype (complex)
        
        resultArray = np.zeros ((numberOfPoints, 3))
    
        identMatrix = np.identity(numFloatingSpecies, dtype=complex)
        
        w = self.getLogSpace (startFrequency, numberOfDecades, numberOfPoints)   
        for i, wf in enumerate(w):
               
            # Compute H(jw) = (jw I - Nr dv/dx L )^{-1} Nr dv/dp
            jw_val = complex (0, wf)
            T1 = identMatrix*jw_val
            T2 = T1 - Jac
            inverse = np.linalg.inv(T2)
            T3 = np.matmul(inverse, Nr)
            dvdp = np.zeros(numReactions, dtype=complex)
            for k in range (numReactions): 
                val = self.r.getUnscaledParameterElasticity(reactionIds[k], parameterId);
                dvdp[k] = val+0j
            T4 = np.matmul (T3, dvdp)
            
            # we need to pick out the requested output species
            for  j in range (numFloatingSpecies):
                if (speciesIds[j] == variableId):
                    dw = abs(T4[j])
                    if useDB:
                       dw = 20.0 * math.log10(dw)
                    resultArray[i, 1] = dw
        
                    val = complex (T4[j], 0)
                    phase = self.getPhase(val)
                    resultArray[i, 2] = phase
                    break
    
            if (useHz):
               # Store frequency, convert to Hz by dividing by 2Pi
               resultArray[i,0] = w[i] / (2. * math.pi);
            else:
               # Store frequency, leave as rad/sec
               resultArray[i,0] = w[i]
               
        # This is to prevent discontinous jumps at the 180 degree points
        resultArray[:,2] = np.unwrap(resultArray[:,2])
        # Convert to degrees
        resultArray[:,2] = (180/math.pi)*resultArray[:,2]
        # Store result internally so that it can be used by plot
        self.results = resultArray
        return resultArray

    def plot (self):
    
        plt.subplot(211)
        plt.xscale("log")
        plt.grid(visible=True, which='major', color='k', linestyle='-', alpha=0.3)
        plt.grid(visible=True, which='minor', color='k', linestyle='-', alpha=0.3)
        plt.minorticks_on()
        plt.plot (self.results[:,0], self.results[:,1], 'r', linewidth=1)
        plt.ylabel('Amp')
        
        plt.subplot(212)
        plt.xscale("log")
        plt.grid(visible=True, which='major', color='k', linestyle='-', alpha=0.3)
        plt.grid(visible=True, which='minor', color='k', linestyle='-', alpha=0.3)
        plt.minorticks_on()
        plt.plot (self.results[:,0], self.results[:,2], 'r', linewidth=1)
        plt.ylabel('Phase (Degrees)')
        plt.xlabel ('Frequency')       
        
        plt.show()
    
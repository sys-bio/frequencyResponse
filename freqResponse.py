# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 11:51:33 2023

@author: hsauro
"""

import roadrunner
import numpy as np
import math
import matplotlib.pyplot as plt

class FreqencyResposne:
    
    def __init__(self):
        pass
    
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

  
    # This function is used in roadrunners freq analysis code
    # Note sure if its working quite right
    def getAdjustment(self, z):
        adjustment = 0
        if (z.real >= 0 and z.imag >= 0):
            adjustment = 0
        elif (z.real >= 0 and z.imag < 0):
            adjustment = 360
        elif (z.real < 0 and z.imag >= 0):
            adjustment = 0
        else:
            adjustment = 360
    
        return adjustment


    # Implement instead:
    # This should be able to handle models with conserved cycles
    #   H(jw) = (jw I - Nr dv/dx L )^{-1} Nr dv/dp
    def getFrequencyResponse(self, r, startFrequency, numberOfDecades, numberOfPoints,
                      parameterId, variableId, useDB=True, useHz = True):
           
        try:
           r.steadyState()
        except Exception as e:
            raise Exception ("Can't find the steady state: " + str (e))
            
        # Get the link matrix in case there are conserved moieties,
        # equals the identity matrix if not
        linkMatrix = r.getLinkMatrix()
        # Convert Jacobian matrix to complex matrix
        Jac = r.getFullJacobian()
        Jac = np.array (Jac)
        Jac = Jac.astype (complex)
        # Compute Jac*linkMatrix
        JacL = np.matmul (Jac, linkMatrix)
    
        # Will need ee for the flux TFs mwhen implemented
        #ee = r.getUnscaledElasticityMatrix()
        # Convert stoich matrix to complex matrix
        Nr = r.getReducedStoichiometryMatrix()
        Nr = np.array (Nr)
        Nr = Nr.astype (complex)
    
        numReactions = r.getNumReactions()
        numFloatingSpecies = r.getNumFloatingSpecies()
        reactionIds = r.getReactionIds()
        speciesIds = r.getFloatingSpeciesIds()
        boundarySpeciesIds = r.getBoundarySpeciesIds()
        globalParameterIds = r.getGlobalParameterIds()
    
        if variableId not in speciesIds:
           raise Exception ("Can't find species:" + variableId + " in the model")     
               
        # Check that parameter name exists
        if parameterId not in boundarySpeciesIds:
           if parameterId not in globalParameterIds:
               raise Exception ("Can't find parameter:" + variableName + " in the model")     
               
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
                val = r.getUnscaledParameterElasticity(reactionIds[k], parameterId);
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
                    phase = (180.0 / math.pi) * self.getPhase(val) + self.getAdjustment(val)
                    resultArray[i, 2] = phase
                    break
    
            if (useHz):
               # Store frequency, convert to Hz by dividing by 2Pi
               resultArray[i,0] = w[i] / (2. * math.pi);
            else:
               # Store frequency, leave as rad/sec
               resultArray[i,0] = w[i]
               
        # Store result internally so that it can be used by plot
        self.results = resultArray
        return resultArray

    def plot (self):
    
        plt.subplot(211)
        plt.xscale("log")
        plt.grid(visible=True, which='major', color='k', linestyle='-', alpha=0.3)
        plt.grid(visible=True, which='minor', color='k', linestyle='-', alpha=0.3)
        plt.minorticks_on()
        plt.plot (self.results[:,0], self.results[:,1], 'r')
        plt.ylabel('Amp')
        
        plt.subplot(212)
        plt.xscale("log")
        plt.grid(visible=True, which='major', color='k', linestyle='-', alpha=0.3)
        plt.grid(visible=True, which='minor', color='k', linestyle='-', alpha=0.3)
        plt.minorticks_on()
        plt.plot (self.results[:,0], self.results[:,2], 'r')
        plt.ylabel('Phase')
               
        plt.show()
    
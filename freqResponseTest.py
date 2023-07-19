# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 20:47:13 2023

@author: hsauro
"""

import tellurium as te
import roadrunner
import numpy as np
import matplotlib.pyplot as plt
from freqResponse import *

r = te.loada("""
    $Xo -> S1; k1*Xo
    S1 -> S2; k2*S1
    S2 -> S3; k3*S2
    S3 -> S4; k4*S3
    S4 -> S5; k5*S4
    S5 -> ; k5*S5
    
    k1 = 0.1; k2 = 0.23; k3 = 0.4
    k4 = 0.33; k5 = 0.26
    Xo = 1
""")
 
fr = FreqencyResponse(r)
results = fr.getFrequencyResponse(0.01, 3, 100, 'Xo', 'S5')

fr.plot()
        
        
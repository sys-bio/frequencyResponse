# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 20:47:13 2023

@author: hsauro
"""

import tellurium as te
import roadrunner
import numpy as np
import matplotlib.pyplot as plt
from freqResponse import FrequencyResponse
import unittest

IGNORE_TEST = False
IS_PLOT = True

ROADRUNNER = te.loada("""
    $Xo -> S1; k1*Xo
    S1 -> S2; k2*S1
    S2 -> S3; k3*S2
    S3 -> S4; k4*S3
    S4 -> ; k5*S4
    
    k1 = 0.1; k2 = 0.23; k3 = 0.4
    k4 = 0.33; k5 = 0.26
    Xo = 1
""")

class TestFrequencyResponse(unittest.TestCase):

    def setUp(self):
        self.fr = FrequencyResponse()

    def testGetFrequencyResponse(self):
        if IGNORE_TEST:
            return
        results = self.fr.getFrequencyResponse(ROADRUNNER, 0.01, 3, 100, 'Xo', 'S3')
        if IS_PLOT:
            self.fr.plot()

if __name__ == '__main__':
    unittest.main()
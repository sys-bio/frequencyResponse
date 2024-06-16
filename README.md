# frequencyResponse
Computes the frequency response for a SBML or Antimony biochemmcal reaction model.

This isn't a pip install package.

Instead copy the freqResponse.py to where you will do your work. 

Run the freqResponseTest.py to try it out.

Exmaple of how to use the code:

```
import tellurium as te
import numpy as np
import matplotlib.pyplot as plt
from freqResponse import *

# Two step pathway
r0 = te.loada('''
   $Xo -> S1; k1*Xo
    S1 -> ;   k2*S1
    
    k1 = 5; k2 = 0.01;
    Xo = 1
''')

fr = FreqencyResponse(r0)
results = fr.getSpeciesFrequencyResponse(0.0001, 4, 100, 'k1', 'S1')
fr.plot()
```

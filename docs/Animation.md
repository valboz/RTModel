[Back to **Final assessment**](FinalAssessment.md)

# Animating fits

It is sometimes useful to analyze the fitting process by simple visualization tools. This is possible by the `RTModel.plotmodel` subpackage, which also includes the possibility to show an animation of the fitting process.

Suppose you want to see how the binary-lens fit from the initial condition 190 works for the event `/event001`. The following code will visualize the step chain of the Levenberg-Marquardt fits in the `(s,q)` plane:

```
from plotmodel import *
import matplotlib.pyplot as plt
import os
import numpy as np
import glob

eventname = '/event001'
model = 'LS0190'
plotchain(eventname,model,0,1)
```

The output will look like this

<img src="plotchain.png" width = 400>



[Go to **Archiving and updating models**](Archive.md)

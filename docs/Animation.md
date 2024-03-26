[Back to **Final assessment**](FinalAssessment.md)

# Animating fits

It is sometimes useful to analyze the fitting process by simple visualization tools. This is possible by the `RTModel.plotmodel` subpackage, which also includes the possibility to show an animation of the fitting process.

## Visualizing the step chain in the parameter space

Suppose you want to see how the binary-lens fit from the initial condition 190 works for the event `/event001`. The following code will visualize the step chain of the Levenberg-Marquardt fits in the `(s,q)` plane:

```
from plotmodel import *
import matplotlib.pyplot as plt
import os
import numpy as np
import glob

eventname = '/event001'
model = 'LS0190'
parameter1 = 0
parameter2 = 1
plotchain(eventname,model,parameter1,parameter2)
```

The output will look like this

<img src="plotchain.png" width = 400>

The `plotchain()` function shows the chains of steps taken by the Levenberg-Marquardt fit from the initial condition to the minima found. The blue line is the first fit, while the following fits after the bumping mechanism described in [Fitting](Fitting.md) are shown with different colors. Each minimum is marked by a filled circle. The parameter space is identified by the two values given as `parameter1,parameter2`. The order of parameters is the one specified in [Model categories](ModelCategories.md).

## Animating the fit

The fit process can be animated by the following code

```
filename = 'PreModels/LS0190/LS0190-stepchain0.dat'
plotmodel(eventname = eventname, modelfile = filename, printpars = False, animate = 1,interval = 800)
```

The output is a gif file `ani.gif` generated in the directory `/event001` which looks like this

<img src="ani.gif" width = 900>

[Go to **Archiving and updating models**](Archive.md)

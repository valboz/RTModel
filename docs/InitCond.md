[Back to **Data pre-processing**](DataPreprocessing.md)

# Initial conditions

As explained in [Full modeling run](ModelingRun.md), all steps in the modeling run are performed by external pre-compiled executables. Each of them can be launched separately through the corresponding function available in `RTModel`.

## The `InitCond` module

The second step in the modeling run is the determination of the initial conditions for the following fits. This is performed by a specific external module called `InitCond`. This can be launched by the corresponding function called `InitCond()`:

```
import RTModel
rtm = RTModel.RTModel('/event001')
rtm.Reader()
rtm.InitCond()
```

With this code, we just perform the data pre-processing performed by `Reader` without proceeding to the following steps. In the `/event001` directory you will see the following products appear:
- a new subfolder called `ini`. This contains the file `Reader.ini` file, which contains the current options with which `Reader` has been launched;
- a file named `LCToFit.txt` containing all data points that will be used for modeling after combining all photometry files found in `/Data`;
- a file named `FilterToData.txt` containing the ordered list of names of the datasets used to build `LCToFit.txt`.

After the execution of `Reader`, you may call the `run()` function to complete the modeling run or the `InitCond()` function if you just want to check the results of the next step, which is [Initial conditions](InitCond.md) settings.

## Setting initial conditions

## Options for initial conditions

### The `config_InitCond()` function

The user may specify his/her own options to drive the initial conditions to the desired result by calling the `config_InitCond()` function with the proper options:

```
import RTModel
rtm = RTModel.RTModel('/event001')
rtm.config_InitCond(npeaks = 2, peakthreshold = 10.0, oldmodels = 4, override = None, nostatic = False, onlyorbital = False, usesatellite = 0)
rtm.run()
```

The call to `config_InitCond()` will affect all following executions of the `InitCond` module, whether called through `run()` or `InitCond()`. If you want to change your options, you may call `config_InitCond()` again.

### Description of the options

self.InitCond_npeaks = npeaks # Number of peaks in the observed light curve to be considered for setting initial conditions.
self.InitCond_peakthreshold = peakthreshold # Number of sigmas necessary for a deviation to be identified as a maximum or a minimum.
self.InitCond_oldmodels = oldmodels # Maximum number of old models to include in new run as initial conditions
self.InitCond_override = override # Override peak identification and manually set peak times
self.InitCond_nostatic = nostatic or onlyorbital # No static models will be calculated.
self.InitCond_noparallax = onlyorbital; # Only orbital motion models will be calculated.
self.InitCond_usesatellite = usesatellite; # Satellite to be used for initial conditions. Ground telescopes by default.


[Go to **Fitting**](Fitting.md)

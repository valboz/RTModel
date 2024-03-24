[Back to **Data pre-processing**](DataPreprocessing.md)

# Initial conditions

## The `InitCond` module

The second step in the modeling run is the determination of the initial conditions for the following fits. This task is performed by a specific external module called `InitCond`. This can be launched by the corresponding function called `InitCond()`:

```
import RTModel
rtm = RTModel.RTModel('/event001')
rtm.Reader()
rtm.InitCond()
```

With this code, we first perform the data pre-processing by `Reader` and then we set the initial conditions by `InitCond`. In the `/event001` directory you will see the following products appear:
- in the subfolder `ini/`, the file `InitCond.ini` appears, which contains the current options with which `InitCond` has been launched;
- 

After the execution of `InitCond`, you may call the `run()` function to complete the modeling run or the `launch_fits()` function for the fits you are interested in, as described in [Fitting](Fitting.md).

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

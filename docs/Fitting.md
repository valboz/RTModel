[Back to **Initial conditions**](InitCond.md)

# Fitting

## The `LevMar` module

The third step in the modeling run is fitting of all models from the given initial conditions. This task is performed by a specific external module called `LevMar`. This can be launched by the corresponding function called `LevMar()`:

```
import RTModel
rtm = RTModel.RTModel('/event001')
rtm.Reader()
rtm.InitCond()
rtm.LevMar('PS0000')
```

With this code, we first perform the data pre-processing by `Reader`, we set the initial conditions by `InitCond` and finally launch the fit of the single-lens-single-source model from the first initial condition found in the file `InitCondPS.txt` (see [Initial conditions](InitCond.md)). Initial condition seeds in each `InitCondXX.txt` file are numbered starting from zero (e.g. `PS0011` would be the 12th initial condition).

In the `/event001` directory you will see the following products appear:
- A new subdirectory called `PreModels/` is created with a file `minchi.dat`. This file contains the value of the minimum chi square among all calculated preliminary models contained in the subdirectories within `PreModels/`.
- In this subdirectory there will be a subsubdirectory called `PS0000/`, dedicated to the models resulting from this run of `LevMar`.
- Inside `PS0000/` we will find some files numbered `0.txt`, `1.txt`, ... containing the preliminary models found by this run of `LevMar`.
- There will also be files named `PS0000-stepchain0.txt`, `PS0000-stepchain1.txt` ... containing all steps of the Levenberg-Marquardt fit taken for each of the corresponding models.
- Besides these, there will also be a file `nlc.dat` containing the numbers of models calculated so far and a file `tPS0000.dat` generated at the ecit of the `LevMar` module to mark the successful closure of the module.

After the execution of `LevMar`, you may call the `run()` function to complete the modeling run or continue with other calls to `LevMar()`, depending on your intentions.

## Launching all fits for a specific category

The `LevMar()` function only launches one fit and will be rarely useful to a generic user. More interesting is the `launch_fits()` function, that launches all fits for a specific model category:

```
import RTModel
rtm = RTModel.RTModel('/event001')
rtm.Reader()
rtm.InitCond()
rtm.launch_fits('PS')
```

In this code, the `launch_fits()` will launch all fits of the `'PS'` [model category](ModelCategories.md) using all available processors or those that have been indicated by the users through the `set_processors()` function (see [Modeling Run](ModelingRun.md)). A progress bar will inform the user of the fits completed and those yet to be done.

At the end of the execution of the `launch_fits()` function, the `PreModels/` directory will be populated by all subdirectories corresponding to fits from all initial conditions in `InitCondPS.txt`. The file `minchi.dat` will contain the minimum chi square found so far, with the name of the subdirectory containing the best model.

## Setting initial conditions

In order to set initial conditions for modeling, `InitCond` executes the following steps:

- Spline approximation of each dataset;
- Identification of peaks in each datasets;
- Removal of duplicate peaks;
- Definition of initial conditions for single-lens models by using the main peak(s) and a grid search over the shape parameters;
- Definition of initial conditions for binary-lens models by matching detected peaks to peaks of templates from a libary;
- Addition of old models from previous runs (if any).

The details of each step will be illustrated in a future publication. 

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

Here we describe the options for `InitCond` in detail indicating their default values.

- `peakthreshold = 10.0`: Number of sigmas necessary for a deviation to be identified as a peak in a concave section with respect to a straight line joining the left and right boundaries of the section. A too low value will include noise in the baseline among peaks. A too high vlaue will ignore small anomalies.
- `npeaks = 2`: Number of peaks in the observed light curve to be considered for setting initial conditions. If you choose to use more than 2 peaks you will have many more fits to be run, with greater chances of success but longer computational time.
- `nostatic = False`: If True, static models will not be calculated. This is useful if higher orders are significant and cannot be treated as a simple perturbation of static models. Furthermore, this option is recommended if you have observations from a satellite spaced by a distance of the order of au.
- `onlyorbital = False`: If true, only orbital motion models will be calculated.
- `usesatellite = 0`: Initial conditions are set only considering peaks in the indicated satellite. If zero, ground datasets are used for initial conditions.
- `oldmodels = 4`: If previous runs have been archived, the chosen number of best models from the last previous run are included as initial conditions. This can be useful for refining old models with new data or options.
- `override = None`: If a t-uple is specified here (e.g `(8760.1, 8793.1)`), the elements of the t-uple are taken as peak positions in the data and directly used to define the initial conditions. The whole spline and peak identification procedure is then skipped.

Notice that the options that are not explicitly specified in the call to `config_InitCond()` are always reset to their default values.

### Recording the options

In each modeling run, the options for `InitCond` are stored in the file `InitCond.ini` in the `/ini` subdirectory within the event directory for later reference. If the modeling run is [archived](Archive.md), also the whole `/ini` subdirectory is saved so that the user may check the options used in each modeling run.

[Go to **Model selection**](ModelSelection.md)

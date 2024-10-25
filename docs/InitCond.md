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

With this code, we first perform the data pre-processing by `Reader` and then we set the initial conditions by `InitCond`. In the `/event001` directory you will see that the following products appear:
- In the subfolder `ini/`, the file `InitCond.ini` appears, which contains the current options with which `InitCond` has been launched.
- A new subdirectory called `InitCond/` is created. In this subdirectory there are several text files named `InitCondXX.txt` and some `PreInitCondXX.txt`, with XX replaced by the label of the corresponding [model category](ModelCategories.md).
- Each `InitCondXX.txt` contains the number of initial conditions and the parameters of each initial conditions line by line.
- A file `spline.txt` in the base event folder containing the points of the spline models built by `InitCond` for each dataset.

After the execution of `InitCond`, you may call the `run()` function to complete the modeling run or the `launch_fits()` function for the fits you are interested in, as described in [Fitting](Fitting.md).

## Setting initial conditions

In order to set initial conditions for modeling, `InitCond` executes the following steps:

- Spline approximation of each dataset;
- Identification of peaks in each datasets;
- Removal of duplicate peaks;
- Definition of initial conditions for single-lens models by using the main peak(s) and a grid search over the shape parameters;
- Definition of initial conditions for binary-lens models by matching detected peaks to peaks of templates from a library;
- Addition of old models from previous runs (if any).

The details of each step are illustrated in the [RTModel paper](https://ui.adsabs.harvard.edu/abs/2024A%26A...688A..83B/abstract). 

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

- `peakthreshold = 10.0`: Number of sigmas necessary for a deviation to be identified as a peak in a concave section with respect to a straight line joining the left and right boundaries of the section. A too low value will include noise in the baseline among peaks. A too high value will ignore small anomalies.
- `npeaks = 2`: Number of peaks in the observed light curve to be considered for setting initial conditions. If you choose to use more than 2 peaks you will have many more fits to be run, with greater chances of success but longer computational time.
- `templatelibrary = None`: Alternative user-defined template library file to be used to build initial conditions for binary-lens fits. You may learn more about the customization of template libraries in [Template libraries](TemplateLibraries.md).  
- `override = None`: If a t-uple is specified here (e.g `(8760.1, 8793.1)`), the elements of the t-uple are taken as peak positions in the data and directly used to define the initial conditions. The whole spline and peak identification procedure is then skipped.
- `usesatellite = 0`: Initial conditions are set only considering peaks in the indicated [satellite](Satellite.md). If zero, ground datasets are used for initial conditions.
- `nostatic = False`: If `True`, static models will not be calculated. This is useful if higher orders are significant and cannot be treated as a simple perturbation of static models. Furthermore, this option is recommended if you have observations from a satellite spaced by a distance of the order of au.
- `onlyorbital = False`: If `True`, only orbital motion models will be calculated.
- `modelcategories = ['PS','PX','BS','BO','LS','LX','LO']`: The list of model categories to be fit in this modeling run. Refer to [Model categoreis](ModelCategories.md) for the meaning of the labels.
- `oldmodels = 4`: If previous runs have been [archived](Archive.md), the chosen number of best models from the last previous run are included as initial conditions. This can be useful for refining old models with new data or options.
- `onlyupdate = False`: If `True`, the modeling run will be limited to a quick update of old models found in the last [archived run](Archive.md) as specified by the `oldmodels` option.

All options are also accessible separately as properties of the `RTModel` class. The user may thus modify the value of each option one by one. The names of the properties are the same as the options in `config_InitCond()` with the prefix `InitCond_`, as shown in the example below:

```
rtm.InitCond_peakthreshold = 20.0
rtm.InitCond_nostatic = True
```

An empty call to `config_InitCond()` with no parameters will reset all variables to the default values. Notice that the options that are not explicitly specified in the call to `config_InitCond()` are always reset to their default values. This is also true if you previously used the `recover_options()` function to inherit the options from a previous run (see [Archiving and updating](Archive.md)).

### Recording the options

In each modeling run, the options for `InitCond` are stored in the file `InitCond.ini` in the `/ini` subdirectory within the event directory for later reference. If the modeling run is [archived](Archive.md), also the whole `/ini` subdirectory is saved so that the user may check the options used in each modeling run. The function `recover_options()` can be used to load the options from a previous run.

[Go to **Fitting**](Fitting.md)

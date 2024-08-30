[Back to **Fitting**](Fitting.md)

# Model selection

## The `ModelSelector` module

The fourth step in the modeling run is the selection of best models for a given model category. This task is performed by a specific external module called `ModelSelector`. This module can be launched by the corresponding function called `ModelSelector()`:

```
import RTModel
rtm = RTModel.RTModel('/event001')
rtm.Reader()
rtm.InitCond()
rtm.launch_fits('PS')
rtm.ModelSelector('PS')
```

With this code, we first perform the data pre-processing by `Reader`, we set the initial conditions by `InitCond`, we launch all fits of the single-lens-single-source model with `launch_fits` and then we select the best models within this category with `ModelSelector`.

In the `/event001` directory you will see the following products appear:
- A new subdirectory called `Models/` is created. This will contain the best models for each category.
- One or more files named `PSXXXX-X.txt` containing the details of the selected models. Each model is identified by the label for the model category followed by the number of initial condition and then by the fit number. Any model obtained by user-defined initial conditions will carry the same label specified by the user.
- In addition, in the `/InitCond` subdirectory, some initial conditions files are updated to include more initial conditions obtained by perturbing the best models found in this category. For example, after the single-lens-single-source fits, initial conditions for binary lenses starting from best models found with single lens are added. These are particularly useful to model small anomalies due to planets.

After the execution of `ModelSelector`, you may call the `run()` function to complete the modeling run or continue with other calls to `launch_fits()` and `ModelSelector()`, or going to [final assessment](FinalAssessment.md) with `Finalizer()`, depending on your intentions.

## Model files

Each model file contains:

- The parameters of the model, starting from the non-linear parameters as described in [Model categories](ModelCategories.md), followed by the blend and source fluxes for each dataset, in the order shown in `FilterToData.txt`, and closing with the chi square.
- The 1-sigma error for each parameter as listed in the first line, except for the chi square.
- The covariance matrix for the parameters as used in the fit. Some of them are fit in log scale (see [Model categories](ModelCategories.md)).

## The model selection

The `ModelSelector` module sorts all preliminary model of the chosen category by their chi square. Models with overlapping covariance ellipsoid are discarded as duplicates. Reflections are considered according to the symmetries of the model category.

## Options for model selection

### The `config_ModelSelector()` function

The user may specify his/her own options to drive the initial conditions to the desired result by calling the `config_ModelSelector()` function with the proper options:

```
import RTModel
rtm = RTModel.RTModel('/event001')
rtm.config_ModelSelector(sigmasoverlap = 3.0, sigmachisquare = 1.0, maxmodels = 10)
rtm.run()
```

The call to `config_ModelSelector()` will affect all following executions of the `ModelSelector` module, whether called through `run()` or `ModelSelector()`. If you want to change your options, you may call `config_ModelSelector()` again.

### Description of the options

Here we describe the options for `LevMar` in detail indicating their default values.

- `sigmasoverlap = 3.0`: Maximum number of sigmas to declare overlap.
- `sigmachisquare = 1.0`: Besides the best model, `ModelSelector` retains competing models up to a threshold given by sqrt(2*chisqure), which represents one sigma in the chi square distribution. This threshold can be changed by this option to include less or more competing models in the final selection.
- `maxmodels = 10`: Maximum number of competing models to report

Notice that the options that are not explicitly specified in the call to `config_LevMar()` are always reset to their default values. This is also true if you previously used the `recover_options()` function to inherit the options from a previous run (see [Archiving and updating](Archive.md)).

### Recording the options

In each modeling run, the options for `ModelSelector` are stored in the file `ModelSelector.ini` in the `/ini` subdirectory within the event directory for later reference. If the modeling run is [archived](Archive.md), also the whole `/ini` subdirectory is saved so that the user may check the options used in each modeling run. The function `recover_options()` can be used to load the options from a previous run.

[Go to **Final assessment**](FinalAssessment.md)

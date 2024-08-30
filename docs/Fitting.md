[Back to **Initial conditions**](InitCond.md)

# Fitting

The third step in the modeling run is fitting of all models from the given initial conditions. This task is performed by a specific external module called `LevMar`. For more direct control, the user may either choose to launch all fits of a given category or launch one individual fit from a specified initial condition. We will examine both possibilities in the following.

## The `LevMar` module

The basic fitting module can be launched by the corresponding function called `LevMar()`:

```
import RTModel
rtm = RTModel.RTModel('/event001')
rtm.Reader()
rtm.InitCond()
rtm.LevMar('PS0000')
```

With this code, we first perform the data pre-processing by `Reader`, we set the initial conditions by `InitCond` and finally launch the fit of the single-lens-single-source model from the first initial condition found in the file `InitCondPS.txt` (see [Initial conditions](InitCond.md)). Initial condition seeds in each `InitCondXX.txt` file are numbered starting from zero (e.g. `PS0011` would be the 12th initial condition).

In the `/event001` directory you will see the following products appear:
- A new subdirectory called `PreModels/` is created with a file `minchi.dat`. This file contains the value of the minimum chi square among all  preliminary models contained in the subdirectories within `PreModels/`.
- In the subdirectory `PreModels/` there will be a subsubdirectory called `PS0000/`, dedicated to the models resulting from this run of `LevMar`.
- Inside `PS0000/` we will find some files numbered `0.txt`, `1.txt`, ... containing the preliminary models found by this run of `LevMar`.
- There will also be files named `PS0000-stepchain0.txt`, `PS0000-stepchain1.txt` ... containing all steps of the Levenberg-Marquardt fit taken for each of the corresponding models.
- Besides these, there will also be a file `nlc.dat` containing the numbers of models calculated so far and a file `tPS0000.dat` generated at the exit of the `LevMar` module to mark the successful closure of the module.

After the execution of `LevMar`, you may call the `run()` function to complete the modeling run or continue with other calls to `LevMar()`, depending on your intentions.

## Launching all fits for a specific category

The `LevMar()` function only launches one fit and will be rarely useful to a generic user, except for checking or repeating a specific initial condition. More interesting is the `launch_fits()` function, that launches all fits for a specific model category:

```
import RTModel
rtm = RTModel.RTModel('/event001')
rtm.Reader()
rtm.InitCond()
rtm.launch_fits('PS')
```

In this code, the `launch_fits()` will launch all fits of the `'PS'` [model category](ModelCategories.md) using all available processors or those that have been indicated by the users through the `set_processors()` function (see [Modeling Run](ModelingRun.md)). A progress bar will inform the user of the fits completed and those yet to be done.

At the end of the execution of the `launch_fits()` function, the `PreModels/` directory will be populated by all subdirectories corresponding to fits from all initial conditions in `InitCondPS.txt`. The file `minchi.dat` will contain the minimum chi square found so far, with the name of the subdirectory containing the best model.

The following step would be the [Model selection](ModelSelection.md) within the fitted model category.

## The Levenberg-Marquardt fit

The `LevMar` module executes a number of Levenberg-Marquardt fits from the specified initial condition. Every time a minimum is found, it is filled with a 'bumper'. Any subsequent run hitting the bumper will be bounced off in a different direction in the parameter space. In this way, new minima can be found from the same initial condition.

The details of this fitting strategy are illustrated in the [RTModel paper](https://ui.adsabs.harvard.edu/abs/2024A%26A...688A..83B/abstract). 

## Options for fitting

### The `config_LevMar()` function

The user may specify his/her own options to drive the initial conditions to the desired result by calling the `config_InitCond()` function with the proper options:

```
import RTModel
rtm = RTModel.RTModel('/event001')
rtm.config_LevMar(nfits = 5, timelimit = 600.0, maxsteps = 50, bumperpower = 2.0)
rtm.run()
```

The call to `config_LevMar()` will affect all following executions of the `LevMar` module, whether called through `run()` or `launch_fits()` or `LevMar()`. If you want to change your options, you may call `config_LevMar()` again.

### Description of the options

Here we describe the options for `LevMar` in detail indicating their default values.

- `nfits = 5`: Number of fits executed from the same initial condition.
- `maxsteps = 50`: Maximum number of steps for each fit.
- `timelimit = 600.0`: Maximum time in seconds allowed for the execution of `LevMar`. If the limit is reached, the execution is stopped and the last step is saved as a preliminary model.
- `bumperpower = 2.0`: Size of the bumper in the parameter space expressed in sigmas. The bumper is created with the shape determined by the local covariance matrix. and the size given by this parameter.

Notice that the options that are not explicitly specified in the call to `config_LevMar()` are always reset to their default values.

### Recording the options

In each modeling run, the options for `LevMar` are stored in the file `LevMar.ini` in the `/ini` subdirectory within the event directory for later reference. If the modeling run is [archived](Archive.md), also the whole `/ini` subdirectory is saved so that the user may check the options used in each modeling run.

[Go to **Model selection**](ModelSelection.md)

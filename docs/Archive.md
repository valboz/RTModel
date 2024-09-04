[Back to **Plotting models**](PlotModel.md)

# Archiving and Updating models

## Archiving a run

After the completion of a modeling run, it is possible to save and **archive the results** and run a new modeling run on the same event with different options or data files. This is useful if you have additional data for some event and you want to keep track of previous results. Another possibility is that the first modeling run with the default options was not satisfactory and you want to try different options for the model search but you want to compare old and new results.

```
import RTModel
rtm = RTModel.RTModel()
rtm.set_event('/event001')
rtm.archive_run()
```

In this code the function `archive_run()` moves all the files and subdirectories inside `/event001` to a new subdirectory called `/run-0001`. Only the `/Data` subdirectory is left in the event directory. The event is then ready for a new modeling run.

All consecutive archived runs are stored in directories with increasing number. Optionally, the user may specify a different archive directory by `rtm.archive_run(destination = 'myarchive')`.

## Previously found best models

An important feature of archived runs is that the best models found in the most recent archived run (identified by its progressive number) are included in the set of initial conditions of the new run. In this way, `RTModel` will also have simple **updates of the previous best models** available in the model selection in the new run.

## Re-using options from previous runs

If you want to use the same options as in the previous run, you can use the `recover_options` function. Supposing you have an event with a completed run, the following code

```
import RTModel
rtm = RTModel.RTModel()
rtm.set_event('/event001')
rtm.recover_options()
```

will generate the output
```
Reader ---  ['tau = 0.1', 'binning = 20000', 'otherseasons = 1', 'renormalize = 1', 'thresholdoutliers = 10']
InitCond ---  ['npeaks = 2', 'peakthreshold = 10.0', 'oldmodels = 4', 'usesatellite = 0', 'override = 6757.903 6753.5']
LevMar ---  ['nfits = 5', 'maxsteps = 50', 'timelimit = 600.0', 'bumperpower = 2.0']
ModelSelector ---  ['sigmasoverlap = 3.0', 'sigmachisquare = 1.0', 'maxmodels = 10']
```

This is the list of all options used in the previous run. You will find detailed explanations for each of these options in the sections for advanced users, starting from [Data pre-processing](DataPreprocessing.md). At the moment, we note that by this function you can repeat a run using exactly the same options you used in the previous run. In fact, all options are automatically loaded in your object `rtm` and will be used as they are unless you explicitly change them. 

You may also pick up the options from one of the archived runs by specifying the optional argument `run`
```
rtm.recover_options(run = '/eventoo1/run-0001')
```

[Go to **Satellite datasets**](Satellite.md)

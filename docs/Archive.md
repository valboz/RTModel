[Back to **Plotting models**](PlotModel.md)

# Archiving and Updating models

After the completion of a modeling run, it is possible to save and archive the results and run a new modeling run on the same event with different options or data files. This is useful if you have additional data for some event and you want to keep track of previous results. Another possibility is that the first modeling run with the default options was not satisfactory and you want to try different options for the model search but you want to compare old and new results.

```
import RTModel
rtm = RTModel.RTModel()
rtm.set_event('/event001')
rtm.archive_run()
```

In this code the function `archive_run()` moves all the files and subdirectories inside `/event001` to a new subdirectory called `/run-0001`. Only the `/Data` subdirectory is left in the event directory. The event is then ready for a new modeling run.

All consecutive archived runs are stored in directories with increasing number. Optionally, the user may specify a different archive directory by `rtm.archive_run(destination = 'myarchive')`.

An important feature of archived runs is that the best models found in the most recent archived run (identified by its progressive number) are included in the set of initial conditions of the new run. In this way, `RTModel` will also have simple updates of the previous best models available in the model selection in the new run.

[Go to **Satellite datasets**](Satellite.md)

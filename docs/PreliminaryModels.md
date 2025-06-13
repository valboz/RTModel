[Back to **High-Resolution Imaging**](HighResolutionImaging.md)

# Preliminary Models

By default, `RTModel` cleans up all the preliminary models that have not been selected in order to save disk space. However, these can be useful for careful diagnostics of the modeling process. If you want to see these preliminary models, you should run modeling with the `cleanup = False` option:
```
rtm.run('/event002', cleanup = False)
```

At the end of the modeling run, you will also see the directory

```
PreModels/               # Directory containing all preliminary models calculated by all fits
```

The preliminary models calculated by `RTModel` may occupy 30MB of disk space in nearly 10000 files. Unless you need some specific debugging or some deeper investigation, we suggest to cleanup the event directory from the directory `/PreModels`, thus saving 99% disk space. You may do it after the modeling run by executing the following line:

```
rtm.cleanup_preliminary_models()
```

Besides the preliminary models, in this directory you will also find the 'stepchain' files, which contain the full list of steps in the fit leading to the calculated preliminary models. These can be visualized by [purposed functions](Animation.md#visualizing-the-step-chain-in-the-parameter-space).

[Go to **Animation of the Fit process**](Animation.md)

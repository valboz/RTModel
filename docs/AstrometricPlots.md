[Back to **Astrophotometric fits**](Astrophotometric.md)

# Astrometric Plots

`RTModel` comes with a subpackage that is intended for fast and basic visualization of results with a minimal number of options.

Assuming you have a completed run on some event in its directory `/event001`, we may plot the best models with the following code

```
import RTModel.plotmodel as plm
import glob

event = '/event001'
models = glob.glob(event +'/FinalModels/*')
model = models[0] # let's plot the first of the best models

myplot = plm.plotmodel(eventname = event, modelfile = model)
```

The output will look like this

<img src="figs/fig_astrophot.png" width = 900>


[Go to **High-Resolution Imaging**](HighResolutionImaging.md)

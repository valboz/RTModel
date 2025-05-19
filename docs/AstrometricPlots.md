[Back to **Astrophotometric fits**](Astrophotometric.md)

# Astrometric Plots

The subpackage `RTModel.plotmodel` contains functions for plotting the astrometric centroid trajectory and comparing with the observations.

Assuming you have a completed run on some astrophotometric event in its directory `/astroevent001`, we proceed as for purely photometric events to [plot the light curve](PlotModel.ms) as usual

```
import RTModel.plotmodel as plm
import glob

event = '/astroevent001'
models = glob.glob(event +'/FinalModels/*')
model = models[0] # let's plot the first of the best models

myplot = plm.plotmodel(eventname = event, modelfile = model)
```

The output will look like this

<img src="figs/fig_astrophot.png" width = 900>




[Go to **High-Resolution Imaging**](HighResolutionImaging.md)

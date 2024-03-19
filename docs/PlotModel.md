[Back to **Model categories**](ModelCategories.md)

# Plotting models

`RTModel` comes with a subpackage that is intended for fast and basic visualization of results with a minimal number of options.

Assuming you have a completed run on some event in its directory `/event001`, we may plot the best models with the following code

```
import plotmodel as plm
import glob

event = '/event001'
models = glob.glob(event +'/SelectedModels/*')
model = models[0] # let's plot the first of the best models

myplot = plm.plotmodel(eventname = event, modelfile = model)
```

The output will look like this

<img src="plotmodel_fig1.png" width = 900>

On the left, we have the model light curve with the data points. Residuals are also shown below. the source trajectory and the caustics are shown on the right. The plots are followed by the list of parameters with their errors. For each telescope we also have the blending fraction ($$F_B/F_*$$) and the baseline magnitude. Finally, the chi square for the model is displayed.

Here is a list of options available for the `plotmodel` function:
- `eventname`: Directory of the event prepared along the indications in [Data preparation](DataPreparation.md)
- `modelfile = None`: The file name containing the model we want to plot. The type of model is identified by the first two characters in the filename. E.g. 'LX0000-0.txt' is a file containing a binary lens with parallax (see [Model categories](ModelCategories.md)). The parameters are read from the file. If `modelfile` is not specified, you may plot any kind of models specifying your parameters in input by providing the arguments `model` and `parameters`.
- `model = None`: If modelfile is left blank, you may specify here the model you want to plot following the labels given in [Model categories](ModelCategories.md). E.g. `model = 'LS'`. The parameters of the model should be given through the argument `parameters`
- `parameters = []`: parameters of a user-defined model. The order and meaning of the parameters depends on the model chosen through the argument `model` (see [Model categories](ModelCategories.md)).
- 

The `plotmodel` function returns a `plotmodel` object that can be manipulated to customize the plot further. 







[Go to **Archiving and updating models**](Archive.md)

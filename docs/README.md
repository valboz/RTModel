
# Documentation

In this document we will describe all use cases of the `RTModel` library and provide ready-to-use examples that you can copy/paste to your code. 

## Quick start

Provided you have prepared a folder for your microlensing event (e.g. `/event001`) with all datasets in the subfolder `/event001/Data` as detailed in the section about [sata preparation](DataPreparation.md), in order to analyze the microlensing event with the default settings, just type

```
import RTModel
rtm=RTModel.RTModel()
rtm.run('/event001')
```

Then the modeling starts and you may sit back and watch the progress bars tracking the status of modeling. You will notice that the Single-source-binary-lens fit will take 1-2 hours whereas the other fits are nearly instantaneous. At the end of all fits you will see the final assessment of `RTModel` for the event with a summary of the best chi squares found for each category of models.

The possible models proposed by `RTModel` are collected as separate files in the directory `/event001/SelectedModels`. You may plot them using the `RTModel.plotmodel` package:

```
import plotmodel as plm
import glob
import os

models = glob.glob('/event001/SelectedModels/*')
modelfile = models[0]
model = os.path.basename(modelfile)

plm.plotmodel(eventname = event, modelfile = modelfile, model = model)
```

## Summary

In the following pages, we will describe all functionalities of `RTModel` in detail.

- [Data preparation](DataPreparation.md)

- [Full modeling run](ModelingRun.md)

- [Data pre-processing](DataPreprocessing.md)

- [Initial conditions](InitialConditions.md)

- [Fitting](Fitting.md)

- [Model selection](ModelSelection.md)

- [Final assessment and results](FinalAssessment.md)

- [Plotting models and results](Plotmodel.md)

- [Satellite datasets](Satellites.md)

- [Limb darkening](LimbDarkening.md)

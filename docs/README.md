
# Documentation

In this document we will describe all use cases of the `RTModel` library and provide ready-to-use examples that you can copy/paste to your code. 

## Quick start

Provided you have prepared a folder for your microlensing event (e.g. `/event001`) with all datasets in the subfolder `/event001/Data` as detailed in the section about [data preparation](DataPreparation.md), in order to analyze the microlensing event with the default settings, just type

```
import RTModel
rtm=RTModel.RTModel()
rtm.run('/event001')
```

Then the modeling starts and you may sit back and watch the progress bars tracking the status of modeling. You will notice that all fits are nearly instantaneous except for the Binary-lens-single-source fit that takes 1-2 hours. At the end of all fits, you will see the final assessment of `RTModel` for the event, along with a summary of the best chi squares found for each category of models.

The final model(s) proposed by `RTModel` are collected as separate files in the directory `/event001/SelectedModels`. You may plot them using the `RTModel.plotmodel` package:

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

In this documentation, we describe all functionalities of `RTModel` in detail. A novel user should read the following pages at least.

- [Data preparation](DataPreparation.md)

- [Full modeling run](ModelingRun.md)

- [Model Categories](ModelCategories.md)

- [Plotting models and results](Plotmodel.md)

Additional useful functionalities are discussed in the following pages.

- [Archiving and updating models](ArchivingUpdating.md)

- [Satellite datasets](Satellites.md)

- [Limb darkening](LimbDarkening.md)

Advanced users may attempt a deeper understanding of the modeling steps and optimize `RTModel` by numerous options.

- [Data pre-processing](DataPreprocessing.md)

- [Initial conditions](InitialConditions.md)

- [Fitting](Fitting.md)

- [Model selection](ModelSelection.md)

- [Final assessment and results](FinalAssessment.md)

- [Animating fits](Animation.md)

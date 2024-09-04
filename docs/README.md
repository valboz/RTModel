
# Documentation

In this document we will describe all use cases of the `RTModel` library and provide ready-to-use examples that you can copy/paste to your code. 

## Quick start

Provided you have prepared a folder for your microlensing event (e.g. `/event001`) with all datasets in the subfolder `/event001/Data` as detailed in the section about [data preparation](DataPreparation.md), in order to analyze the microlensing event with the default settings, just type

```
import RTModel
rtm = RTModel.RTModel()
rtm.run('/event001')
```

Then the modeling starts and you may sit back and watch the progress bars tracking the status of modeling. You will notice that all fits are nearly instantaneous except for the Binary-lens-single-source fit that takes 1-2 hours. At the end of all fits, you will see the final assessment of `RTModel` for the event, along with a summary of the best chi squares found for each category of models.

The final model(s) proposed by `RTModel` are collected as separate files in the directory `/event001/FinalModels`. You may plot them using the `RTModel.plotmodel` package:

```
import RTModel.plotmodel as plm
import glob

event = '/event001'
models = glob.glob(event +'/FinalModels/*')
model = models[0]

plm.plotmodel(eventname = event, modelfile = model)
```

## Summary

In this documentation, we describe all functionalities of `RTModel` in detail. A novel user should read the following pages at least.

- [Data preparation](DataPreparation.md)

- [Full modeling run](ModelingRun.md)

- [Model Categories](ModelCategories.md)

- [Plotting models and results](PlotModel.md)

Additional useful functionalities are discussed in the following pages.

- [Archiving and updating models](ArchivingUpdating.md)

- [Satellite datasets](Satellite.md)

- [Limb darkening](LimbDarkening.md)

Advanced users may attempt a deeper understanding of the modeling steps and optimize `RTModel` by numerous options.

- [Data pre-processing](DataPreprocessing.md)

- [Initial conditions](InitCond.md)

- [Fitting](Fitting.md)

- [Model selection](ModelSelection.md)

- [Final assessment and results](FinalAssessment.md)

- [Animating fits](Animation.md)

- [Template library customization](TemplateLibrary.md)

## Success rate

The success rate of `RTModel` has been evaluated on the simulated events created for the [WFIRST data challenge](https://roman.ipac.caltech.edu/docs/street_data_challenge1_results.pdf) by Matthew Penny. Here we report the results from the current version.

### Planetary regime

In this challenge there were 51 binary events with q<0.03 (planetary regime). The results obtained by the current version of  `RTModel` with the default options were the following:
- 42 full successes;
- 6 cases in which an s<1 solution was found instead of the s>1 (or the opposite);
- 1 cases in which a binary solution was preferred;
- 2 case in which the anomaly was not detected at all (too much noise)

### Binary regime

In the stellar binary regime q>0.03 there were 75 events simulated in the data challenge. We note that orbital motion was simulated by assigning the two transverse components, which leads to non-physical trajectories not reproducible by circular orbital motion. Anyway, the results obtained are the following:

- 58 full successes;
- 7 cases in which an s<1 solution was found instead of the s>1;
- 8 cases in which a different binary model was found;
- 2 cases in which the anomaly was too weak or perfectly reproduced by a single-lens with parallax

The success rate is observed to decrease significantly for events with high orbital motion. This should be partly due to the different way this effect is taken into account in the simulation and in our fitting code.

In addition, we note that such events with strong higher order effects would have been better fit with the `nostatic = True` option (see [Initial conditions](InitCond.md)).


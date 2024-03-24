[Back to **Data preparation**](DataPreparation.md)

# Modeling run

## RTModel start-up 

Once data have been prepared according to the format specified in [Data preparation](DataPreparation.md), we can go to Python and write a simple code to perform the full modeling run.

```
import RTModel
rtm = RTModel.RTModel('/event001')
```
Here we have first imported the `RTModel` package, then we have defined an instance to an `RTModel` object. We may optionally specify the name of the event to analyze already in the constructor or leave it blank and indicate it later. The event is identified through the path to the directory prepared before.

The output of the constructor looks like this:
```
*********************
****   RTModel   ****
*********************
Number of processors: 24
```

It is useful to know that by default `RTModel` will use all available processors to run fits in parallel. If you are concerned that your CPU will burn, you may explicitly choose the number of processors you want to use by the following line

```
nprocessors = 8
rtm.set_processors(nprocessors)
```

where `nprocessors` is the number you want to use.

## Launching modeling run(s)

After that, we are ready to launch modeling through the `run()` function. In this example, we are assuming that we just want to use all available processors (see above).

```
import RTModel
rtm = RTModel.RTModel('/event001')
rtm.run()
```

Note that an alternative equivalent writing is the following:
 
```
import RTModel
rtm = RTModel.RTModel()
rtm.run('/event001')
```

The event name can be specified directly in the `run()` function. In principle, more events can be modeled by the same `RTModel` object one after the other:

```
import RTModel
rtm = RTModel.RTModel()
rtm.run('/event001')
rtm.run('/event002')
...
```

Or we can have a list of events to model in a campaign and perform a single loop to model them all.

```
import RTModel
import glob

rtm = RTModel.RTModel()

eventlist = glob.glob(/campaign2024/*)
for event in eventlist:
  rtm.run(event)
```

If you just want to change the name of the event without launching a modeling run, you may use the `set_event()` function

```
rtm.set_event('/event002')
```

A modeling run may take from minutes to hours depending on the number of data points, the complexity of the event, the available CPU. Progress bars are shown to inform the user of the current status of the modeling run. The slowest step is the Binary-lens-single-source fitting, which takes 90% of the time.

## Output of a modeling run

### Assessment on the nature of the event

At the end of the modeling run, a summary is displayed with the best chi square achieved for each model category, the final assessment and the list of proposed models. The same content can be found in a file named `nature.txt` that is created inside the event directory.

### Best models

The best models listed in `nature.txt` are available as text files in the subdirectory `event001/FinalModels`.

Each file contains the list of parameters in the first line, including the background and source fluxes for each telescope and the total chi square. This means that the first line contains `nps + 2 * ntel + 1` values, where `nps` is the number of parameters in the model category (e.g. 7 for binary-lens-single-source) and `ntel` is the number of datasets.

A detailed explanation of parameters for each model category is available in [Model categories](ModelCategories.md)

The second line contains the uncertainty on each of the above listed parameters (except for the chi square). Therefore, this line contains `nps + 2 * ntel` values.

The remaining lines contain the covariance matrix between the model parameters. Therefore, there are `nps` additional lines containing `nps` values each. Note that some parameters are internally fit in logarithmic scale (see [Model categories](ModelCategories.md) for which parameters are fit as ln). The reported covariance matrix is thus calculated on these internal parameters.

### Additional products

Additional partial products of modeling are stored in the event directory. At the end of the modeling run we will find the following files and directories
```
Data/                    # Directory containing the original input data files as described in data preparation
ini/                     # Directory containing text files specifying the options for individual modules
InitCond/                # Directory containing the text files with the initial conditions for fitting
PreModels/               # Directory containing subdirectories with all models calculated by all fits
Models/                  # Directory containing selected models for each category
FinalModels/             # Directory containing the best models as proposed in the final assessment (see above)
LCToFit.txt              # Text file containing the formatted and pre-processed data points
FilterToData.txt         # Table of datasets names
spline.txt               # List of points in the spline approximation used for initial conditions
nature.txt               # Text file containing the final assessment on the event and the list of best models
```

These files will be explained in the following chapters in due time. They can be useful for careful diagnostics of the modeling process. You may vision an event with a completed run among the provided examples: [event001done.zip](/events/event001done.zip).

## Structure of a modeling run

The `run()` function performs a precise sequence of operations on the selected event. The individual operations are performed by external pre-compiled modules called by the function. Here we summarize the individual steps:

1. Data pre-processing
2. Initial conditions setting
3. Fits for a specific category of models
4. Model selection within the category
5. Final assessment on the event

Steps 3 and 4 are repeated for each model category. Details on the model categories and their parameters are given in the [following page](ModelCategories.md)

Each step is described in a dedicated documentation page (see [Summary](README.md)). Each step can be executed separately through a dedicated function. Numerous options are available to change the behavior of each step.

In case the modeling run is interrupted, it can be resumed by executing the `run()` function again on the same event. The function checks whether the expected output of each step is present and in such case it moves to the next step.

[Go to **Model categories**](ModelCategories.md)

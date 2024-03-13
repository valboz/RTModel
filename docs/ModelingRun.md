[Back to **Data preparation**](DataPreparation.md)

# Modeling run

## Constructor and number of processors

Once data have been prepared according to the format specified in [Data preparation](DataPreparation.md), we can go to Python and write a simple code to perform the full modeling run.

```
import RTModel
rtm = RTModel.RTModel('/event001')
```
First we import the `RTModel` package, then we define an instance to an `RTModel` object. We may optionally specify the name of the event to analyze already in the constructor or leave it blank and indicate it later. The event is identified through the path to the directory prepared before.

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

The best models listed in this summary are available as text files in the subdirectory `event001/SelectedModels`.

Each file contains the list of parameters


## Structure of a modeling run

The `run()` function performs a sequence of operations on the selected event and generate some output files containing the final assessment and the proposed best models for the event. The individual operations are performed by external pre-compiled models 


[Go to **Data pre-processing**](DataPreprocessing.md)

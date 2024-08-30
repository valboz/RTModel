[Back to **Archiving and updating models**](Archive.md)

# Satellite datasets

Observations from satellites, if available, are identified through their photometry filename ending with a number, e.g. `Spitzer1.dat`. Each satellite is identified by a different number that is used to select the correct ephemerides table. All other files ending by anything that is not a number is considered as taken from a ground telescope.

In all models without parallax ('PS', 'BS', 'BO', 'LS' (see [Model categories](ModelCategories.md))), satellite datasets are treated in the same way as ground datasets. For models including parallax ('PX', 'LX', 'LO'), the position of the satellite in space is taken into account through its respective ephemerides table. This leads to a different source trajectory with respect to the caustics for the event as seen from the satellite. Therefore, the light curve model for the satellite data can be profoundly different. This difference is extremely useful to fix the parallax parameters and retrieve physical information on the lens distance.

## Ephemerides tables

In order to exploit satellite observations, we need satellite ephemerides covering the observation period. `RTModel` conforms to `VBBinaryLensing` standard, using ephemerides of the satellite in the format given by the [NASA Horizons system](http://ssd.jpl.nasa.gov/horizons.cgi).

In particular, we assume five columns:
- JD
- RA (degrees)
- Dec (degrees)
- Distance from Earth (AU)
- Distance rate change (not really needed but included by default in Horizons).

In order to obtain the table in the correct format in [Horizons](http://ssd.jpl.nasa.gov/horizons.cgi), in Table Settings you should only check flags 1. Astrometric RA & Dec and 20. Observer range & range-rate. Date/time format should be Julian Day calendar, angle format should be decimal degrees, range units should be astronomical units.

An example of a valid satellite ephemerid tables is in [/events/satellite1.txt](/events/satellite1.txt).

The satellite table(s) should be named "satellite*.txt" (with * replaced by a number) and placed in some directory `/satellitedir', which does not need to be inside the event directory. 

As satellite datasets are identified by a number preceding the extension, they are matched to the ephemerides tables with the same number.

In order to let `RTModel` find the ephemerides tables, you should include a call to the function `set_satellite_dir` in your code before the start of the modeling run:

```
import RTModel
rtm = RTModel.RTModel('/OB190033')
rtm.set_satellite_dir('/satellitedir')
rtm.run()
```

## Modeling strategies

By default `RTModel` looks for static models first (without parallax) and then only performs fits with parallax on the best static models found. In case of satellite observations, this strategy may work if the satellite is very close to Earth (Earth or Moon orbiting and L2 are good examples) or for very small parallax values. The small difference between the source trajectories as seen from the Earth and from the satellite can be easily recovered as a higher order refinement.

However, for au-separated satellites, static models are often useless or even driven to completely wrong results if satellite data are forced to follow the same light curve model as for ground telescopes. Then, we should skip static models and start the search with parallax included from the beginning of the search. This is achieved by the `nostatic = True` option in the `config_InitCond()` function:

```
import RTModel
rtm = RTModel.RTModel('/OB190033')
rtm.set_satellite_dir('/satellitedir')
rtm.config_InitCond(nostatic = True)
rtm.run()
```

This is thus the recommended sequence of instructions for modeling an event including satellite data. Further discussion and details are provided in the [Initial conditions](InitCond.md) page.

## Plotting satellite observations

For events including satellite observations, the `plotmodel` call should include the indication of the directory containing the ephemerides tables:

```
plm.plotmodel(eventname = event, modelfile = model, satellitedir = '/satellitedir')
```

The plots will include the light curves for ground datasets and for satellite datasets with different colors. The same colors are used to represent the source trajectories as seen from ground and satellite in the plot on the right:

<img src="plotmodel_fig2.png" width = 900>


[Go to **Limb darkening**](LimbDarkening.md)

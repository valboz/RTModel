[Back to **Summary**](README.md)

# Data preparation

In this page we explain how the data should be prepared for the analysis by `RTModel`. You can find some sample events already prepared in the directory [events](/events). These can be useful to see how an event directory should look like at the beginning or at the end of the modeling run.

## Directory structure

Each microlensing event should have its own dedicated directory somewhere on your PC (e.g. `/event001/`). 

This directory should contain a subdirectory named `/Data`.

The `/Data` directory should contain all photometric time series available for the analysis. Each series (data collected by a single telescope in one filter) corresponds to one file with extension  `.dat`.

## Photometry files

The content of each `.dat` file should be as in the following example:

```
# Mag err HJD-2450000
19.0232 0.012 8370.1223
19.0150 0.011 8370.2421
19.0034 0.011 8370.3697
18.9712 0.010 8370.4911
18.9592 0.011 8370.6114
18.9109 0.009 8370.8234
18.8798 0.009 8371.0092
...

```

The first line is a header specifying the content of the columns. 

`RTModel` accepts both input in magnitudes or fluxes. 

If the header contains the keyword "Mag", then `RTModel` assumes that each line contains magnitude, error and Heliocentric Julian Date - 2450000 for each individual photometric measurement.

If the header DOES NOT contain the keyword "Mag", then `RTModel` assumes that each line contains flux, error and Heliocentric Julian Date - 2450000 for each individual photometric measurement.

## Event coordinates

Event coordinates are contained in a file with extension `.coordinates` (e.g. `event001.coordinates`) placed in the same `/Data` directory along with the photometry files.

The content of the file should be in the form `HH:MM:SS.S +DD:PP:SS.S` for right ascension and declination respectively, e.g. `18:00:23.32 -32:11:09.7`.

## Optional input files

Other optional input files are observations from [satellite](Satellite.md), [limb darkening coefficients](LimbDarkening.md) and [normalizations for datasets](DataPreprocessing.md#forcing-error-bar-normalization), described in the corresponding pages. [Astrophotometric datasetes](Astrophotometric.md) are also discussed in a separate page.

[Go to **Modeling Run**](ModelingRun.md)

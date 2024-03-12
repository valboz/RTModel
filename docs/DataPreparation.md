[Back to **Summary**](README.md)

# Data preparation

## Directory structure

Each microlensing event should have its own dedicated directory (e.g. `/event001/`). 

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

The first line is just a header and is ignored by `RTModel`.

Each data point contains magnitude, error and Heliocentric Julian Date - 2450000.

## Satellite observations

Observations from satellites, if available, are identified through the filename ending with a number, e.g. `Spitzer1.dat`. Each satellite is identified by a different number that is used to select the correct ephemerides table.

All other files ending by anything that is not a number is considered as taken from a ground telescope.

## Event coordinates

Event coordinates are contained in a file with extension `.coordinates` (e.g. `event001.coordinates`) placed in the same `/Data` directory along with the photometry files.

The content of the file should be in the form `HH:MM:SS.S +DD:PP:SS.S` for right ascension and declination respectively, e.g. `18:00:23.32 -32:11:09.7`.

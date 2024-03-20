This directory contains some simulated events with data prepared in the format required by `RTModel` (see [Data preparation](/docs/DataPreparation.md)).

The zipped directory `event001done.zip` contains the same data as `event001.zip` with the full modeling run already completed. It can be used as reference for what you should obtain at the end of a modeling run.

The event in `OB190033.zip` is OGLE-2019-BLG-0033, published in [Herald et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022A%26A...663A.100H/abstract). The zip file contains all the data used for the analysis, including those from the *Spitzer* satellite. It is a good event to start practising with satellite data. The corresponding ephemerides for *Spitzer* are in the file `satellite1.txt` available in this directory. If you run `RTModel` on this event, make sure to select the option `nostatic = True` in `config_InitCond()` (see [Initial conditions](InitCond.md)).

[Back to **Limb darkening**](LimbDarkening.md)

# Data pre-processing

As explained in [Full modeling run](ModelingRun.md), all steps in the modeling run are performed by external pre-compiled executables. Each of them can be launched separately through the corresponding function available in `RTModel`.

## The `Reader` module

The first step in the modeling run is the data pre-processing, which is performed by a specific external module called `Reader`. This can be launched by the corresponding function called `Reader()`:

```
import RTModel
rtm = RTModel.RTModel('/event001')
rtm.Reader()
```

With this code, we just perform the data pre-processing performed by `Reader` without proceeding with the following steps. In the `/event001` directory you will see the following products appear:
- a new subfolder called `ini`. This contains the file `Reader.ini` file, which contains the current options with which `Reader` has been launched;
- a file named `LCToFit.txt` containing all data points that will be used for modeling after combining all photometry files found in `/Data`;
- a file named `FilterToData.txt` containing the ordered list of names of the datasets used to build `LCToFit.txt`.




```
tau # conventional correlation time for consecutive points
binning # maximum number of points left after re-binning
otherseasons # How to use other seasons (0 = Yes, 1 = decrease significance, 2 = remove)
renormalize # Re-normalize error bars if non-zero
thresholdoutliers # Threshold in sigmas for removing outliers
```

[Go to **Initial conditions**](InitCond.md)

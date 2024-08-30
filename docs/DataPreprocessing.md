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

With this code, we just perform the data pre-processing by the `Reader` module without proceeding to the following steps. In the `/event001` directory you will see the following products appear:
- a new subfolder called `ini/`. This contains the file `Reader.ini` file, which contains the current options with which `Reader` has been launched;
- a file named `LCToFit.txt` containing all data points that will be used for modeling after combining all photometry files found in `/Data`;
- a file named `FilterToData.txt` containing the ordered list of names of the datasets used to build `LCToFit.txt`.

After the execution of `Reader`, you may call the `run()` function to complete the modeling run or the `InitCond()` function if you just want to check the results of the next step, which is [Initial conditions](InitCond.md) settings.

## Pre-processing operations

The `Reader` module performs several operations on the original data to improve the chances of a successful modeling run and limit the computational time. All these operations can be fully controlled by the user through specific options. In particular:

- Data are re-binned until the number of data points is contained within the desired amount;
- Error bars are re-normalized based on the assessment of local scatter;
- Outliers are removed;
- Datasets left with less than 4 points are ignored.

For details about the algorithms used in this pre-processing, please refer to the [RTModel paper](https://ui.adsabs.harvard.edu/abs/2024A%26A...688A..83B/abstract).

## Options for pre-processing

### The `config_Reader()` function

The user may specify his/her own options to drive the pre-processing to the desired result by calling the `config_Reader()` function with the proper options:

```
import RTModel
rtm = RTModel.RTModel('/event001')
rtm.config_Reader(binning = 4000, tau = 0.1, otherseasons = 1, renormalize = 1, thresholdoutliers = 10)
rtm.run()
```

The call to `config_Reader()` will affect all following executions of the `Reader` module, whether called through `run()` or `Reader()`. If you want to change your options, you may call `config_Reader()` again.

### Description of the options

Here we describe the options in detail with their default values:

- `binning = 4000`: the maximum number of data points you want to model. If the original datasets total to less than `binning` they are left untouched.
- `tau = 0.1`: The timescale (in days) used for the assessment of local scatter and for re-binning. In a first approximation, `RTModel` considers variations below `tau` as possible scatter.
- `otherseasons = 1`: how to treat seasons other than the season containing the peak value: 0 for including all seasons, 1 to downgrade the significance of other seasons in the re-binning process, 2 to remove all seasons other than the peak season.
- `renormalize = 1`: if non-zero, all datasets are re-normalized based on the scatter assessment.
- `thresholdoutliers = 10`: threshold in sigmas to remove outliers.

Notice that the options that are not explicitly specified in the call to `config_Reader()` are always reset to their default values. This is also true if you previously used the `recover_options()` function to inherit the options from a previous run (see [Archiving and updating](Archive.md)).

### How to switch off pre-processing

If you want to avoid any modifications to the original data, you may switch off all pre-processing by the following call
```
rtm.config_Reader(binning = 1000000, renormalize = 0, thresholdoutliers = 1000000)
```

Re-binning and outliers removal do not intervene if set to very high numbers and renormalization is switched off.

### Recording the options

In each modeling run, the options for `Reader` are stored in the file `Reader.ini` in the `/ini` subdirectory within the event directory for later reference. If the modeling run is [archived](Archive.md), also the whole `/ini` subdirectory is saved so that the user may check the options used in each modeling run. The function `recover_options()` can be used to load the options from a previous run.

## Forcing error-bar normalization

With the `renormalize` option, the user may choose between no re-normalization and error-bar re-normalization based on local scatter. However, the user may wish to use his/her own normalization factors. These can be provided by adding the file `Normalizations.txt` to the `/Data` subdirectory. This file should contain the normalization factors for each dataset one per line:
```
1.23
0.97
```

In this example, there are just two datasets. After switching-off re-normalization with the `renormalize = 0` option, the error bars will be multiplied by these factors in all fits. The order of the datasets is the same that can be found in the `FilterToData.txt` file, generated by the `Reader` module.

[Go to **Initial conditions**](InitCond.md)

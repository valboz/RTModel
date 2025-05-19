[Back to **Final assessment**](FinalAssessment.md)

# Astrophotometric modeling

Space telescopes or adaptive optics facilities may reach astrometric precisions of the order of milliarcsecond, thus being sensitive to the astrometric shift of the images during the microlensing event. Follow-up observations may also measure the source proper motion providing constraints on the microlensing event. If the lens is luminous, it may be resolved as a separate object, enabling a direct measurement of the relative proper motion. All this additional information can be incorporated in `RTModel` via the inclusion of astrophotometric datasets, which trigger a combined astrophometric fit, in which additional physical parameters are retrieved based on the astrometric measurements.

## Astrophotometric datasets

We have already discussed purely [photometric datasets](DataPreparation.md), which should be prepared with thre columns: magnitude, error, HJD. We recall an example here for convenience:

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

If we have Declination and Right Ascension of our microlensing event with the respective uncertainties, we will have four additional columns: Declination, error on declination, Right Ascension, error on right ascension. Therefore, an astrophotometric dataset will look as follows:

```
# Mag err HJD-2450000 Dec errDec RA errRA
19.0232 0.012 8370.1223 1.2534 1.0 0.0165 1.0
19.0150 0.011 8370.2421 1.7510 1.1 -0.5422 1.1
19.0034 0.011 8370.3697 1.1190 1.1 -0.1981 1.1
18.9712 0.010 8370.4911 1.4281 1.0 0.2119 1.0
18.9592 0.011 8370.6114 1.3005 0.9 0.3982 1.0
18.9109 0.009 8370.8234 1.6233 1.0 0.5121 1.0
18.8798 0.009 8371.0092 2.0223 1.2 0.9411 1.1
...

```

As usual, the first line is just a header and is ignored by `RTModel`. Dec and RA indicate the angular displacements in milliarcseconds from a fixed reference point (decided by the observer) in the North and East directions respectively.

In the directory [`/events`](/events) we have an example of an astrophotometric event ([`astroevent001.zip'](/events/astroevent001.zip)) with one purely photometric dataset and an astrophotometric dataset that can be useful to understand how to prepare such datasets.

We caution that for the purpose of [binning](DataPreprocessing.md#pre-processing-operations) the astrometric undertainties are not taken into account. Re-binning affects photometry and astrometry at the same time.

## Astrophotometric models

If `RTModel` finds an astrophotometric dataset, it automatically includes four additional parameters in the fit, which are necessary to model the astrometric centroid trajectory:

| Parameter | Meaning | 
| --- | --- |
| muS_Dec | Source proper motion component in the North direction in milliarcseconds per year | 
| muS_RA | Source proper motion component in the East direction in milliarcseconds per year |
| piS | Geometric parallax of the source in milliarcseconds |
| thetaE | Einstein angle in milliarcseconds |

Furthermore, the microlensing parallax components piN and piE are always included in the model, which means that no static models are fitted. This is because the lens proper motion is obtained from the source proper motion using the information on the relative proper motion hidden in the standard microlensing parameters. In definitive, the [model categories](ModelCategories.md) used by default in astrophotometric fits are `['PX','BO','LX','LO']`, with a total of 10, 14, 13, 16 parameters respectively. The Keplerian fit 'LK' with 18 parameters can be added by the user, if desired.

## Results

The results of an astrophotometric modeling run are displayed in the end with the chi square including both the photometric and the astrometric contributions. The assessment is contained in `nature.txt` with the list of models which are available as text files in the subdirectory `/FinalModels`, similarly to purely [photometric fits](ModelingRun.md#best-models). However, the model files also contain the information relative to astrometric parameters.

In particular, the list of parameters in the first line contains `nps + 4 * ntel + 1` values, where `nps` is the number of parameters in the model category including the 4 astrometric parameters (e.g. 13 for binary-lenses) and `ntel` is the number of datasets. For each dataset, besides the background and source flux, we also have the centroid position at time t0 in Dec and RA, which is obtained by the fitting. Finally, we have the chi square.

The second line contains the uncertainty on each of the above listed parameters (except for the chi square). Therefore, this line contains `nps + 4 * ntel` values.

The remaining lines contain the covariance matrix between the model parameters. Therefore, there are `nps` lines containing `nps` values each. 

[Go to **Astrometric plots**](AstrometricPlots.md)

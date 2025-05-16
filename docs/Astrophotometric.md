[Back to **Final assessment**](FinalAssessment.md)

# Astrophotometric datasets

Space telescopes or adaptive optics facilities may reach astrometric precisions of the order of milliarcsecond, thus being sensitive to the astrometric shift of the images during the microlensing event. Follow-up observations may also measure the source proper motion providing constraints on the microlensing event. If the lens is luminous, it may be resolved as a separate object, enabling a direct measurement of the relative proper motion.

All this additional information can be incorporated in `RTModel` via the inclusion of astrophotometric datasets. We have already discussed purely [photometric datasets](DataPreparation.md), which should be prepared with thre columns: magnitude, error, HJD. We recall an example here for convenience:

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

If we have Declination and Right Ascension of our microlensing event with the respective uncertainties, we will have four additional columns: Declination, error, Right Ascension, error. An astrophotometric dataset will look as follows:

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

As usual, the first line is just a header and is inessential



[Go to **High Resolution Imaging**](HighResolutionImaging.md)

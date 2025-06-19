[Back to **Modeling run**](ModelingRun.md)

# Model categories

In the current version, `RTModel` fits 7 different model categories by default. Each of them has a speficic label with two characters. This table summarizes the main information

| Label | Model | Number of parameters |
| --- | --- | --- |
| PS | Single-lens-single-source | 4 |
| PX | Single-lens-single-source with parallax | 6 |
| BS | Single-lens-binary-source | 7 |
| BO | Single-lens-binary-source with xallarap | 10 |
| LS | Binary-lens-single-source | 7 |
| LX | Binary-lens-single-source with parallax | 9 |
| LO | Binary-lens-single-source with circular orbital motion | 12 |

## Including or excluding model categories

By default, `RTModel` fits all the model categories listed above to the data, providing models for all of them and comparing their chi square to make its final assessment. However, the user may specify which model categories to fit and even include additional models not proposed in the default modeling run by specifying the corresponding option, as described in [Initial conditions](InitCond.md). 

## Additional model categories

Some additional model categories are available but not included in the default modeling run. They can be included by adding their label to the list of model categories as detailed in [Initial conditions](InitCond.md). At the moment, we have only one additional category, but more are on their ways.

| Label | Model | Number of parameters |
| --- | --- | --- |
| LK | Binary-lens-single-source with eccentric orbital motion | 14 |

The following sections provide details about parameters and conventions of all categories listed above.

Note that some parameters are fit in logarithmic scale by `RTModel` for better performance. The only consequence for users is that the covariance matrix contained in model files is written for the ln parameters.

## Single-lens-single-source (PS)

This is the basic microlensing model with a single-lens and a single-source. Note that `RTModel` only considers finite-source models. In most cases, there will be no constraints on the source radius anyway.

The parameters of this model are the following. We also indicate if the parameter is internally fit in logarithmic (ln) scale (see above).

| Number | Parameter | Meaning | ln |
| --- | --- | --- | --- |
| 1 | u0 | Impact parameter normalized to Einstein angle | X |
| 2 | tE | Einstein time in days | X |
| 3 | t0 | Closest approach time in HJD |  |
| 4 | rho | Source radius normalized to Einstein angle | X |

We note that only positive values of `u0` are considered in this model since a distinction only comes when parallax is included.

## Single-lens-single-source with parallax (PX)

This is the same as before with parallax vector included. 

| Number | Parameter | Meaning | ln |
| --- | --- | --- | --- |
| 1 | u0 | Impact parameter normalized to Einstein angle |  |
| 2 | tE | Einstein time in days | X |
| 3 | t0 | Closest approach time in HJD |  |
| 4 | rho | Source radius normalized to Einstein angle | X |
| 5 | piN| Parallax component along North |  |
| 6 | piE | Parallax component along East |  |

Note that the positive and negative impact parameters correspond to different models when we include parallax.
The parallax components are calculated with respect to the reference time t0.

## Single-lens-binary-source (BS)

Here we consider two sources and one lens. The two sources have different fluxes and different finite radii. Here are the parameters:

| Number | Parameter | Meaning | ln |
| --- | --- | --- | --- |
| 1 | tE | Einstein time in days | X |
| 2 | FR | Flux ratio of the secondary to the primary source | X |
| 3 | u01 | Impact parameter of the primary source |  |
| 4 | u02 | Impact parameter of the secondary source |  |
| 5 | t01 | Closest approach time of the primary source |  |
| 6 | t02 | Closest approach time of the secondary source |  |
| 7 | rho1 | Source radius for the primary | X |

Note that both sources may have positive or negative impact parameters.

The source radius is only given for the primary star, while the secondary star has a source radius calculated by the relation rho2 = rho1 * FR^(0.225), coming from approximate stellar mass-luminosity-radius relations for solar-type stars. This exponent can be changed by operating on the mass-radius and mass-luminosity exponents, as explained in [Fitting](Fitting.md).

## Single-lens-binary-source with xallarap (BO)

Here we also include circular orbital motion of the two sources around a common center of mass:

| Number | Parameter | Meaning | ln |
| --- | --- | --- | --- |
| 1 | tE | Einstein time in days | X |
| 2 | FR | Flux ratio of the secondary to the primary source | X |
| 3 | u01 | Impact parameter of the primary source |  |
| 4 | u02 | Impact parameter of the secondary source |  |
| 5 | t01 | Closest approach time of the primary source |  |
| 6 | t02 | Closest approach time of the secondary source |  |
| 7 | rho1 | Source radius for the primary | X |
| 8 | piN| Parallax component along North |  |
| 9 | piE | Parallax component along East |  |
| 10 | gamma1 | Angular velocity parallel to the lens axis |  |
| 11 | gamma2 | Angular velocity perpendicular to the lens axis |  |
| 12 | gammaz | Angular velocity along the line of sight |  |

See [VBMicrolensing Binary Sources](https://github.com/valboz/VBMicrolensing/blob/master/docs/python/BinarySources.md) for a detailed explanation.

The mass ratio between the two component is calculated as qs = FR^(0.25) and the radius of the secondary is rho2 = rho1*FR^(0.225). These exponents can be changed by operating on the mass-radius and mass-luminosity exponents, as explained in [Fitting](Fitting.md).

Finally, note that we may also [turn off the light from the secondary component](Fitting.md), if we want to model a source orbiting a dark object.

## Binary-lens-single-source (LS)

This is the "static" binary-lens model:

| Number | Parameter | Meaning | ln |
| --- | --- | --- | --- |
| 1 | s | Separation between the lenses in Einstein radii | X |
| 2 | q | Mass ratio of the secondary to the primary lens | X |
| 3 | u0 | Impact parameter normalized to Einstein angle |  |
| 4 | alpha | Angle between the source velocity and the vector pointing from the secondary to the primary lens |  |
| 5 | rho | Source radius normalized to Einstein angle | X |
| 6 | tE | Einstein time in days | X |
| 7 | t0 | Closest approach time in HJD to the barycenter |  |

For details about the source trajectory parameterization, please check [VBMicrolensing Light Curves](https://github.com/valboz/VBMicrolensing/blob/master/docs/python/LightCurves.md).

## Binary-lens-single-source with parallax (LX)

Here we also include parallax as we did before for the single-lens case:

| Number | Parameter | Meaning | ln |
| --- | --- | --- | --- |
| 1 | s | Separation between the lenses in Einstein radii | X |
| 2 | q | Mass ratio of the secondary to the primary lens | X |
| 3 | u0 | Impact parameter normalized to Einstein angle |  |
| 4 | alpha | Angle between the source velocity and the vector pointing from the secondary to the primary lens |  |
| 5 | rho | Source radius normalized to Einstein angle | X |
| 6 | tE | Einstein time in days | X |
| 7 | t0 | Closest approach time in HJD to the barycenter |  |
| 8 | piN| Parallax component along North |  |
| 9 | piE | Parallax component along East |  |


## Binary-lens-single-source with circular orbital motion (LO)

Finally, we also explore circular orbital motion for our binary lens:

| Number | Parameter | Meaning | ln |
| --- | --- | --- | --- |
| 1 | s | Separation between the lenses in Einstein radii | X |
| 2 | q | Mass ratio of the secondary to the primary lens | X |
| 3 | u0 | Impact parameter normalized to Einstein angle |  |
| 4 | alpha | Angle between the source velocity and the vector pointing from the secondary to the primary lens |  |
| 5 | rho | Source radius normalized to Einstein angle | X |
| 6 | tE | Einstein time in days | X |
| 7 | t0 | Closest approach time in HJD to the barycenter |  |
| 8 | piN| Parallax component along North |  |
| 9 | piE | Parallax component along East |  |
| 10 | gamma1 | Angular velocity parallel to the lens axis |  |
| 11 | gamma2 | Angular velocity perpendicular to the lens axis |  |
| 12 | gammaz | Angular velocity along the line of sight |  |

The three components of the orbital motion are expressed at time t0 in units of days^-1.

The component gamma1 is equivalent to ds/dt/s. The component gamma2 is dalpha/dt. The third component gammaz is generally poorly constrained, but is required to move along a physical circular orbit.

More details are available at [VBMicrolensing Orbital Motion](https://github.com/valboz/VBMicrolensing/blob/master/docs/python/OrbitalMotion.md).

The subpackage [plotmodel](PlotModel.md) contains a function for translating from the fitting parameters to conventional orbital elements.

## Binary-lens-single-source with eccentric orbital motion (LK)

Eccentric orbital motion is not fit by default, but can be requested by setting the corresponding option in the [Initial Conditions](InitCond.md).

| Number | Parameter | Meaning | ln |
| --- | --- | --- | --- |
| 1 | s | Separation between the lenses in Einstein radii | X |
| 2 | q | Mass ratio of the secondary to the primary lens | X |
| 3 | u0 | Impact parameter normalized to Einstein angle |  |
| 4 | alpha | Angle between the source velocity and the vector pointing from the secondary to the primary lens |  |
| 5 | rho | Source radius normalized to Einstein angle | X |
| 6 | tE | Einstein time in days | X |
| 7 | t0 | Closest approach time in HJD to the barycenter |  |
| 8 | piN| Parallax component along North |  |
| 9 | piE | Parallax component along East |  |
| 10 | gamma1 | Angular velocity parallel to the lens axis |  |
| 11 | gamma2 | Angular velocity perpendicular to the lens axis |  |
| 12 | gammaz | Angular velocity along the line of sight |  |
| 13 | sz_s | Separation along the line of sight over projected separation |  |
| 14 | a_s3d | Semimajor axis over total separation |  |

The two additional parameters are sufficient to completely define an eccentric orbit. More details are available at [VBMicrolensing Orbital Motion](https://github.com/valboz/VBMicrolensing/blob/master/docs/python/OrbitalMotion.md).

The subpackage [plotmodel](PlotModel.md) contains a function for translating from the fitting parameters to conventional orbital elements.

[Go to **Plotting models**](PlotModel.md)


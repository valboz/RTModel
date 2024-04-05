[Back to **Modeling run**](ModelingRun.md)

# Model categories

In the current version, `RTModel` offers 7 different model categories. Each of them has a speficic label with two characters. This table summarizes the main information

| Label | Model | Number of parameters |
| --- | --- | --- |
| PS | Single-lens-single-source | 4 |
| PX | Single-lens-single-source with parallax | 6 |
| BS | Single-lens-binary-source | 7 |
| BO | Single-lens-binary-source with xallarap | 10 |
| LS | Binary-lens-single-source | 7 |
| LX | Binary-lens-single-source with parallax | 9 |
| LO | Binary-lens-single-source with orbital motion | 12 |

The following sections provide details about parameters and conventions.

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

The source radius is only given for the primary star, while the secondary star has a source radius calculated by the relation rho2 = rho1 * FR^(0.225), coming from approximate stellar mass-luminosity-radius relations for solar-type stars. In general, we do not expect that both sources have a detectable finite-size effect, so this relation is enforced with the only purpose to avoid the exploration of grossly unphysical models.

## Single-lens-binary-source with xallarap (BO)

Here we also include circular orbital motion of the two sources around a common center of mass:

| Number | Parameter | Meaning | ln |
| --- | --- | --- | --- |
| 1 | u01 | Impact parameter of the primary source |  |
| 2 | t01 | Closest approach time of the primary source |  |
| 3 | tE | Einstein time in days | X |
| 4 | rho1 | Source radius for the primary | X |
| 5 | xi1 | Xallarap component parallel to the source velocity | |
| 6 | xi2 | Xallarap component orthogonal to the source velocity | |
| 7 | omega | Orbital angular velocity in days^-1 |  |
| 8 | i | Inclination of the orbital plane in radians |  |
| 9 | phi | Phase of the orbit from the passage on the line of nodes | |
| 10 | qs | Mass ratio of the secondary to the primary | X |

The position of the secondary source is calculated from the position of the primary at any time. See [VBBinaryLensing Binary Sources](https://github.com/valboz/VBBinaryLensing/blob/master/docs/BinarySources.md) for a detailed explanation.

Annual parallax is not considered in this model because it can be mimicked by xallarap, as well known and would only induce bad degeneracies in a higher dimensional parameter space. However, satellite observations may distinguish between parallax and xallarap.

Note that the flux ratio is calculated as FR = qs^4 and the radius of the secondary is rho2 = rho1 * qs^0.9. As before, we do not expect that these mass-radius-luminosity relations are strictly needed except for avoiding the exploration of unphysical models.

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

For details about the source trajectory parameterization, please check [VBBinaryLensing Light Curves](https://github.com/valboz/VBBinaryLensing/blob/master/docs/LightCurves.md).

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


## Binary-lens-single-source with orbital motion (LO)

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
| 11 | gamma2 | Angulr velocity perpendicular to the lens axis |  |
| 12 | gammaz | Angular velocity along the line of sight |  |

The three components of the orbital motion are expressed at time t0 in units of days^-1.

The component gamma1 is equivalent to ds/dt/s. The component gamma2 is dalpha/dt. The third component gammaz is generally poorly constrained, but is required to move along a physical circular orbit.

More details are available at [VBBinaryLensing Orbital Motion](https://github.com/valboz/VBBinaryLensing/blob/master/docs/OrbitalMotion.md).

[Go to **Plotting models**](PlotModel.md)


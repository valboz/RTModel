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
| 5 | pi_N| Parallax component along North |  |
| 6 | pi_E | Parallax component along East |  |

Note that the positive and negative impact parameters correspond to different models when we include parallax.
The parallax components are calculated with respect to the reference time t0.

## Single-lens-binary-source (BS)

Here we consider two sources and one lens. The two sources have different fluxes and different finite radii. Here are the parameters:

| Number | Parameter | Meaning | ln |
| --- | --- | --- | --- |
| 1 | tE | Einstein time in days | X |
| 2 | FR | Flux ratio of the secondary over primary source | X |
| 3 | u01 | Impact parameter of the primary source |  |
| 4 | u02 | Impact parameter of the secondary source |  |
| 5 | t01 | Closest approach time of the primary source |  |
| 6 | t02 | Closest approach time of the secondary source |  |
| 7 | rho1 | Source radius for the primary | X |

Note that both sources may have positive or negative impact parameters.

The source radius is only given for the primary star, while the secondary star has a source radius calculated by the relation rho2 = rho1 * FR^(3.6), coming from approximate stellar mass-luminosity-radius relations in the solar range. In general, we do not expect that both sources have a detectable finite-size effect, so this relation is enforced with the only purpose to avoid the exploration of grossly unphysical models.


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
| 10 | qs | Mass ratio of the secondary over the primary | X |

The position of the secondary source is calculated from the position of the primary at any time. See [VBBinaryLensing Binary Sources](https://github.com/valboz/VBBinaryLensing/blob/master/docs/BinarySources.md) for a detailed explanation.

Note that the flux ratio is calculated as FR = qS^4 and the radius of the secondary is rho2 = rho1 * qs^0.9. As before, we do not expect that these mass-radius-luminosity relations are strictly needed except for avoiding the exploration of unphysical models.





[Go to **Data pre-processing**](DataPreprocessing.md)


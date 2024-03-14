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


[Go to **Data pre-processing**](DataPreprocessing.md)


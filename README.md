# RTModel
`RTModel` is a package for modeling and interpretation of microlensing events. It uses photometric time series collected from ground and/or space telescopes to propose one or more possible models among the following:
- Single-lens-single-source microlensing (i.e. Paczynski)
- Single-lens-binary-source microlensing (with or without xallarap)
- Binary-lens-single-source microlensing (including planetary microlensing, parallax and orbital motion)

All models include the finite-size of the source(s).

The modeling strategy is based on a grid search in the parameter space for single-lens models, whereas a **template library** for binary-lens models is used including all possible geometries of the source trajectory with respect to the caustics. In addition to this global search, planets are searched where maximal deviations from a Paczynski model occurs. 

The library is in the form of a standard Python package that launches specific subprocesses for different tasks. Model fitting is executed in **parallel** exploiting available processors in the machine. The full modeling may take from one to three hours depending on the event and on the machine speed. The results of modeling are given in the form of a text **assessment file**; in addition, **final models** are made available with their parameters and covariance matrices.

`RTModel` also includes a subpackage **`RTModel.plotmodel`** that allows an immediate visualization of models and the possibility to review each individual fitting process as an animated gif.

## Attribution

`RTModel` has been created by Valerio Bozza (University of Salerno) as a product of many years of direct experience on microlensing modeling (see [RTModel webpage](http://www.fisica.unisa.it/GravitationAstrophysics/RTModel.htm)). 

Any scientific use of `RTModel` should be acknowledged by citing the paper [V.Bozza, A&A 688 (2024) 83](https://ui.adsabs.harvard.edu/abs/2024A%26A...688A..83B/abstract), describing all the algorithms behind the code.

We are grateful to Greg Olmschenk, who revised the package installation in order to make it as cross-platform as possible.

## Installation

The easiest way to install `RTModel` is through `pip`. 

First clone this repository.

Then go to the repository directory and type

```
pip install .
```

In alternative, you may directly install it from PyPI without cloning this repository:

```
pip install RTModel
```

Currently, `RTModel` works on Linux, Windows and MacOS, requiring Python >= 3.6. 
A C++ compiler compatible with C++17 standard is needed for installation.
`RTModel` also incorporates version 3.7 of [`VBBinaryLensing`](https://github.com/valboz/VBBinaryLensing).

## Documentation
Full [documentation for the use of RTModel](docs/README.md) is available.

In the directory [events](events) we provide some microlensing data on which you may practise with `RTModel`.

A Jupyter notebook for quick start-up is also available in the [jupyter](jupyter) folder.

## License
`RTModel` is freely available to the community under the 
GNU Lesser General Public License Version 3 included in this repository.



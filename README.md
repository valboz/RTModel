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

A paper is in preparation describing all algorithms behind RTModel in detail.
At the moment, any use of this code for scientific publications should be acknowledged by citing this GitHub repository. 

## Installation

The easiest way to install `RTModel` is through `pip`
```
pip install RTModel
```

Currently, `RTModel` works on Linux and Windows, requiring Python >= 3.6. 
A C++ compiler compatible with C++17 standard is needed for installation.
`RTModel` also incorporates version 3.7 of [`VBBinaryLensing`](https://github.com/valboz/VBBinaryLensing).

Example Jupyter notebooks will be included in `examples/`.

## Documentation
Provisional [documentation for the use of RTModel](docs/README.md) is available and under development.

## License
`RTModel` is freely available to the community under the 
GNU Lesser General Public License Version 3 included in this repository.



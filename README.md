# RTModel
RTModel is a package for modeling and interpretation of microlensing events. It uses photometric time series collected from ground and/or space telescopes to propose one or more possible models among the following:
- Single-lens-single-source microlensing (i.e. Paczynski)
- Single-lens-binary-source microlensing (with or without xallarap)
- Binary-lens-single-source microlensing (including planetary microlensing, parallax and orbital motion)
All models include the finite-size of the source(s).

The modeling strategy is based on a grid search in the parameter space for single-lens models while it is based on a template library for binary-lens models including all possible geometries of the source trajectory with respect to the caustics. In addition to this global search, planets are searched where maximal deviations from a Paczynski model occurs.


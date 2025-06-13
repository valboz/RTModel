[Back to **Limb darkening**](LimbDarkening.md)

# Constrained fits

It is possible to impose constraints on the fit parameters or some pre-defined combinations of them that are particularly important in microlensing observations. Let $f(\mathbf{p})$ be a function of the fit parameters. Suppose we want to force $f(\mathbf{p}) = f_0$ with some tolerance $\sigma$. Then we make the following addition to the $\chi^2$ function:

$\tilde\chi^2 =\chi^2 + \left(\frac{f(\mathbf{p})-f_0}{\sigma}\right)^2$

The modified $\tilde\chi^2$ including a gaussian constraint on $f$ is then used in minimization. In `RTModel` it is also possible to consider asymmetric constraints with different tolerances on either side of $f_0$ following the form 

$\tilde\chi^2 =\chi^2 + \left(\frac{f(\mathbf{p})-f_0}{\sigma_l}\right)^2\Theta(f_0-f(\mathbf{p})) + \left(\frac{f(\mathbf{p})-f_0}{\sigma_r}\right)^2\Theta(f(\mathbf{p})-f_0)$

## Implementation

The practical implementation in `RTModel` occurs through the function `set_constraints` as in the following example

```
import RTModel
rtm = RTModel.RTModel('/event001')

myconstraints = [['tE', 23.0, -1.0, +1.0],
                 ['log_rho', -3, -2, +0.1],
                 ['g_1', 0.0, -0.01, 1.e100]]
rtm.set_constraints(myconstraints)

rtm.run()
```

Here we have built a list of constraints on several parameters or functions that will be included in the modified  $\tilde\chi^2$ as above. Each constraint is in the form of a list following the general form `[function_name, f_0, sigma_l, sigma_r]`. The numbers for `f_0, sigma_l, sigma_r` express the asymmetric constraint as explained above. Note that one-sided constraints can be easily obtained if one of the sigmas is set to extremely large values, as in the third constraint in the example above.

The following sections explain the syntax for the `function_name` in detail.

## Names of the constrained functions

Each entry in the list of constraints is a list starting with the name of the parameter or the function that is going to be constrained. We have four possibilities available to the user. 

### Constraining a parameter

A single fit parameter can be just indicated as constrained function. For the names of the parameters, refer to [Model Categories](ModelCategories.md). Note that if the constraint is defined on a parameter that is not relevant for a given model category, the constraint will just be ignored for that specific model category.

### Constraining the logarithm of a parameter

In the second constraint in the example above we see how to constrain the logarithm (base 10) of a given parameter. The name of the parameter should just be preceded by the prefix `'log_'`.

### Constraining blending

The third constraint in the example above shows how to constrain the blending parameter for telescope number 1 (the telescope list starts from 0). We remind that in `RTModel` the blending parameter is defined as the ratio of the background flux to the source flux: $g_1 = F_{background,1}/F_{source,1}$. Note that in the specific case, non-negative blending has been requested.

### Other functions

To constrain other functions of the parameters, you may choose from the following list, to be updated on demand by users:

- 'muangle': The angle $\psi_\mu$ in radians of the relative lens-source proper motion from the North direction taken counterclockwise. This angle affects the parallax parameters according to

$\arctan \frac{\pi_E}{\pi_N} =\psi_\mu$

This direction can be obtained from high-resolution observations which have been able to resolve the lens and the source individually. The constraint on the magnitude of the proper motion becomes a constraint on $t_E$ when combined with constraints on $\theta_E$ coming from the source study and the finite source effect.

- 't*': The source radius crossing time

$t_E \rho_* =t_*$

## Recovering previous constraints

The constraints used in a given modeling run are stored in the file `Constraints.ini` in the subdirectory `/ini` in the event directory.

The function `recover_options()` illustrated in [Archiving and updating models](Archive.md) recovers previously used constraints for a given event. These can be modified by manipulating the `rtm.constraints` object and will be used in the next modeling run. 

[Go to **Data pre-processing**](DataPreprocessing.md)

[Back to **Limb darkening**](LimbDarkening.md)

# Constrained fits

It is possible to impose constraints on the fit parameters or some pre-defined combinations of them that are particularly important in microlensing observations. Let $f(\mathbf{p})$ be a function of the fit parameters. Suppose we want to force $f(\mathbf{p}) = f_0$ with some tolerance $\sigma$. Then we make the following addition to the $\chi^2$ function:

$\tilde\chi^2 =\chi^2 + \left(\frac{f(\mathbf{p})-f_0}{\sigma}\right)$

The modified $\tilde\chi^2$ including a gaussian constraint on $f$ is then used in minimization. In `RTModel` it is also possible to consider asymmetric constraints of the form 

$\tilde\chi^2 =\chi^2 + \left(\frac{f(\mathbf{p})-f_0}{\sigma_l}\right)\Theta(f_0-f(\mathbf{p})) + \left(\frac{f(\mathbf{p})-f_0}{\sigma_r}\right)$\Theta(f(\mathbf{p})-f_0)$


[Go to **Data pre-processing**](DataPreprocessing.md)

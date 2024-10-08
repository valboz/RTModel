[Back to **Limb darkening**](LimbDarkening.md)

# Constrained fits

It is possible to impose constraints on the fit parameters or some pre-defined combinations of them that are particularly important in microlensing observations. Let $f(\mathbf{p})$ be a function of the fit parameters. Suppose we want to force $f(\mathbf{p}) = f_0$ with some tolerance $\sigma$. Then the $\chi^2$ function is modified as

$\chi^2 = \sum\limits_i \left(\frac{y_i-model_i(\mathbf{p})}{\sigma_i\right)^2 + \left(\frac{f(\mathbf{p})-f_0}{\sigma}\right)$

[Go to **Data pre-processing**](DataPreprocessing.md)

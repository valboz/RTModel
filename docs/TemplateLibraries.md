[Back to **Animating fits**](Animation.md)

# Template libraries

As explained before, the module [InitCond](InitCond.md) sets the initial conditions for binary-lens fits by matching the peaks found in the observed datasets to the peaks of the templates in a library. This matching provides the values for $t_0$ and $t_E$, while all remaining parameters are read from the template.

The default library used by `InitCond` is available [here](/RTModel/data/TemplateLibrary.txt). The first line indicates the number of templates, while the following lines contain the information for each template one by one. Each line contains the parameters $s$ $q$ $u_0$ $\alpha$ $\rho$ followed by the time of two peaks $t_1$ $t_2$. If the template has more than two peaks, the template is repeated in the following line with different choices of the times of the peaks.

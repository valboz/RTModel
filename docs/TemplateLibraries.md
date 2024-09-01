[Back to **Animating fits**](Animation.md)

# Template libraries

As explained before, the module [InitCond](InitCond.md) sets the initial conditions for binary-lens fits by matching the peaks found in the observed datasets to the peaks of the templates in a library. This matching provides the values for $t_0$ and $t_E$, while all remaining parameters are read from the template.

The default library used by `InitCond` is available [here](../../RTModel/data/TemplateLibrary.txt).

[Back to **Model selection**](ModelSelection.md)

# Final assessment

## The `Finalizer` module

The final step in the modeling run is the classification of the microlensing event through comparison of models obtained for each category. This task is performed by a specific external module called `Finalizer`. This module can be launched by the corresponding function called `Finalizer()`:

```
import RTModel
rtm = RTModel.RTModel('/event001')
rtm.Reader()
rtm.InitCond()
rtm.launch_fits('PS')
rtm.ModelSelector('PS')
...
rtm.Finalizer()
```

With this code, we first perform all steps as detailed in the previous sections and then we make the final assessment through `Finalizer()`.

In the `/event001` directory you will see the following products appear:
- A new subdirectory called `FinalModels/` is created. This will contain the final proposed models for the microlensing event, which are just copies of those appearing in the directory `Models/` that have passed all critieria.
- A file `nature.txt` containing a summary of the chi square achieved in each model category, a final assessment, and the list of the competing models that have passed all criteria and copied to the directory `FinalModels/`.

The execution of `Finalizer` closes the modeling run. You may then proceed to [plotting](PlotModel.md) the best models.

## Making the final assessment

The `Finalizer` module compares models of different categories according to Wilks' theorem, which states that the level at which a model with p additional parameters is preferred can be assessed by the chi square distribution with p degrees of freedom. More details will be given in our future publication. Models that are not nested one within the other are compared with the softer threshold given by the 1-sigma level in the chi square distribution.

There are no available options for `Finalizer()` since an unsatisfied user may easily vet the models in the directory `Models/` according to his/her own criteria.

[Go to **Animating fits**](Animation.md)

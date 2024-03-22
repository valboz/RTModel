[Back to **Archiving and updating models**](Archive.md)

# Limb Darkening

The main purpose of `RTModel` is to find preliminary models of microlensing events that can be later studied and refined with all additional information. In general, to meet this goal, it is not necessary to consider the source limb darkening and a uniform brightness profile is adopted in all calculations.

However, there are known exceptions in which limb darkening is essential to avoid claiming planet detections from improperly modeled finite-source effects. In fact, in high-magnification single-lens events it may happen that the difference between the uniform and limb-darkened model can be mimicked by a distortion of the central caustic induced by a planet.

Limb darkening can be included in the whole `RTModel` modeling run just by adding one more file named `LimbDarkening.txt` to the `/Data` subdirectory of the event.

The file `LimbDarkening.txt` should contain the values of the limb darkening coefficients for each telescope



[Go to **Data pre-processing**](DataPreprocessing.md)

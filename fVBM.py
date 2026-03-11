import importlib.util
import os

spec = importlib.util.find_spec("VBMicrolensing")

if spec is None or spec.origin is None:
    print("! Error: VBMicrolensing non installed!")
    raise SystemExit(1)

package_dir = os.path.dirname(spec.origin).replace("\\", "/")
print(package_dir)

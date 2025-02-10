import os
import inspect

try:
    import VBMicrolensing
    package_dir = os.path.dirname(inspect.getfile(VBMicrolensing)).replace('\\','/')

 #   with open("vbmp.txt", "w") as f:
 #       f.write(package_dir)
        
    print(f"{package_dir}")

except ImportError:
    print("! Error: VBMicrolensing non installed!")
    raise


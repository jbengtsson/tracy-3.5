import os, sys

sys.path.append(os.getenv('TRACY_LIB')+'/tracy/lib')
import libtracy

tst

# ext_modules = [
#     Pybind11Extension(
#         "python_example",
#         sorted(glob("src/*.cc")),  # Sort source files for reproducibility
#     ),
# ]

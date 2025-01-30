import sys
import importlib

# Import the main Pybind11 module (mypackage.pybind_module)
pybind_module = importlib.import_module("siconos.numerics.pynumerics")


# Dynamically attach attributes from the Pybind11 main module
__all__ = []  # Keep track of exported names
for attr in dir(pybind_module):
    if not attr.startswith("_"):  # Skip private attributes
        globals()[attr] = getattr(pybind_module, attr)
        __all__.append(attr)



# Dynamically handle submodules
submodules = []

for submodule in submodules:
    try:
        # Import the submodule
        imported_submodule = importlib.import_module(f"siconos.numerics.pynumerics.{submodule}")

        # Attach submodule to the `mypackage` namespace
        globals()[submodule] = imported_submodule

        # Add submodule to `__all__`
        __all__.append(submodule)

        # Optional: Make submodules directly importable
        sys.modules[f"siconos.mechanics.{submodule}"] = imported_submodule

    except ModuleNotFoundError:
        # Handle missing submodules gracefully
        pass

# Hide the internal Pybind11 module
del sys.modules["siconos.numerics.pynumerics"]
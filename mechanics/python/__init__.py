import sys
import importlib

# Import the main Pybind11 module (pymechanics)
pybind_module = importlib.import_module("siconos.mechanics.pymechanics")

# Dynamically attach attributes from the Pybind11 main module
__all__ = []  # Track exported names
for attr in dir(pybind_module):
    if not attr.startswith("_"):  # Skip private attributes
        globals()[attr] = getattr(pybind_module, attr)
        __all__.append(attr)

# Recursive function to handle submodules
def attach_submodules(base_module, base_name):
    """
    Recursively attach submodules to the namespace.
    """
    for attr in dir(base_module):
        if not attr.startswith("_") and isinstance(getattr(base_module, attr), type(importlib)):
            # Check if it's a submodule and import it
            submodule_name = f"{base_name}.{attr}"
            try:
                imported_submodule = importlib.import_module(submodule_name)

                # Attach submodule to the current namespace
                globals()[attr] = imported_submodule
                __all__.append(attr)

                # Optional: Make submodules directly importable
                sys.modules[submodule_name.replace(".pymechanics", "")] = imported_submodule

                # Recursively handle deeper submodules
                attach_submodules(imported_submodule, submodule_name)

            except ModuleNotFoundError:
                # Gracefully handle missing submodules
                pass


# Attach the top-level submodules
attach_submodules(pybind_module, "siconos.mechanics.pymechanics")

# Hide the internal Pybind11 module
del sys.modules["siconos.mechanics.pymechanics"]
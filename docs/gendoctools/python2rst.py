"""Tools used during configuration and build process of
Siconos documentation.

Process:

python files --> rst (sphinx-apidoc) --> html (sphinx)

Siconos is a program dedicated to modeling, simulation and control
 of non smooth dynamical systems.

 Copyright 2024 INRIA.

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
"""

from pathlib import Path

# To display a module on the main python api reference, add a doc in this dict.
# Other modules won't appear on the page but will be available in the toc.
modules_docs = {
    "siconos.externals": "API or tools related to external software libraries used by Siconos.",
    "siconos.numerics": (
        "collection of low-level algorithms for solving algebra and "
        "optimization problems."
    ),
    "siconos.geometry": "geometry tools.",
    "siconos.modeling": "high-level API to model the nonsmooth dynamical systems.",
    "siconos.integrators": "classes and tools for time integration of the dynamics "
    "(one-step integrators).",
    "siconos.nonsmooth_formulations": "classes and tools dedicated to the description"
    " of the nonsmooth problem (formulation and solver).",
    "siconos.simulation": "classes and tools dedicated to the description of the simulation"
    " of a NSDS (time discretization, ...).",
    "siconos.control": "control toolbox.",
    "siconos.mechanics": "toolbox for collision detection and joints",
    "siconos.mechanisms": "toolbox for collision detection and joints"
    " (legacy version, won’t be sustained in long term)",
    "siconos.io": "input/output tools (HDF5, VTK, ...).",
}


def build_python_api_main(outputdir):
    """Parse existing rst files generated for python API
    by sphinx-apidoc and collect them into reference/index.rst

    Called by make rst_api

    Parameters
    ----------
    outputdir : Path()
         sphinx directory which contains rst files
         generated for the api by sphinx-apidoc
    """
    outputdir = Path(outputdir)
    index = outputdir / "python_api.rst"
    rst_dir = outputdir / "python"

    with index.open("w") as f:
        f.write(".. _python_api:\n\n")
        title = "Siconos Python API reference"
        f.write(title + "\n")
        f.write("=" * len(title) + "\n\n")
        f.write(
            "This is the documentation of "
            "`Python <https://www.python.org/>`_ interface to Siconos.\n\n"
        )

        for name in modules_docs:
            rst = rst_dir / f"{name}.rst"

            # if rst.exists():
            #    f.write(f"   python/{name}\n")
        f.write("\n")

        f.write("Main packages\n")
        f.write("-------------\n\n")

        # Main packages: those documented in modules_docs
        # that will be listed on the page.
        # Other modules will be available through the toc on the left
        # thanks to hidden toctree below

        first_indent = "   "
        second_indent = 2*first_indent

        f.write(".. grid:: 3\n")
        f.write(f"{first_indent}:gutter: 3\n\n")
        for name in modules_docs:
            rst = rst_dir / f"{name}.rst"
            if rst.exists():
                f.write(f"{first_indent}.. grid-item-card:: {name}\n")
                f.write(f"{second_indent}:link: python/{name}\n")
                f.write(f"{second_indent}:class-card: sd-bg-code\n")
                f.write(f"{second_indent}:link-type: doc\n\n")
                f.write(f"{second_indent}{modules_docs[name]}\n\n")
                # f.write(f"* :doc:`{name} <python/{name}>` : {modules_docs[name]}\n")
        f.write("\n")

        f.write("Access to the complete reference of all modules is available")
        f.write(" in the left menu under 'siconos.package'.\n\n")

        f.write(".. toctree::\n")
        f.write("   :hidden:\n")
        f.write("   :maxdepth: 2\n\n")

        f.write("   python/siconos\n\n")

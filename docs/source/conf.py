# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath("../../"))
sys.path.insert(0, os.path.abspath("../../deepchem_server"))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "DeepChem Server"
copyright = "2024, DeepChem Server Team"
author = "DeepChem Server Team"
release = "1.0.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.githubpages",
]

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# -- Options for autodoc extension -------------------------------------------

# This value selects what content will be inserted into the main body of an autoclass directive.
autoclass_content = "both"

# This value selects how the signature will be displayed for the class defined by autoclass directive.
autodoc_class_signature = "mixed"

# This value controls the docstrings inheritance.
autodoc_inherit_docstrings = True

# This value controls the behavior of sphinx.ext.autodoc when it encounters a member without docstrings.
autodoc_undoc_members = True

# This value controls whether or not to display the source code link in the generated documentation.
autodoc_member_order = "bysource"

linkcheck_report_timeouts_as_broken = False

# Mock imports for packages that are not available during documentation build
autodoc_mock_imports = [
    "deepchem",
    "deepchem.feat",
    "deepchem.data",
    "deepchem.models",
    "deepchem.utils",
    "deepchem.metrics",
    "deepchem.trans",
    "deepchem.hyper",
    "deepchem.dock",
    "deepchem.metalearning",
    "deepchem.rl",
    "rdkit",
    "rdkit.Chem",
    "rdkit.DataStructs",
    "fastapi",
    "fastapi.responses",
    "fastapi.exceptions",
    "uvicorn",
    "pydantic",
    "starlette",
    "starlette.responses",
    "starlette.status",
    "requests",
    "requests_toolbelt",
    "pandas",
    "numpy",
    "scipy",
    "sklearn",
    "joblib",
    "pickle",
    "json",
    "yaml",
    "matplotlib",
    "seaborn",
    "plotly",
    "thefuzz",
]

# -- Options for intersphinx extension ---------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "deepchem": ("https://deepchem.readthedocs.io/en/latest/", None),
}

# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True

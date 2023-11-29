import os
import sys

sys.path.insert(0, os.path.abspath('../../modules/'))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Tull echelle reduction'
copyright = '2023, Daniel Krolikowski'
author = 'Daniel Krolikowski'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.napoleon', 'autoapi.extension']

autoapi_dirs = ['../../modules']
autoapi_type = "python"
autoapi_keep_files = True
autoapi_add_toctree_entry = False

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinxdoc'
html_static_path = ['_static']
html_theme_options = { 'sidebarwidth': '25em' }

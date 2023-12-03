import os
import sys
import sphinx_rtd_theme

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

extensions = ['sphinx.ext.napoleon', 'autoapi.extension', 'sphinx_rtd_theme']

autoapi_dirs = ['../../modules']
autoapi_type = "python"
autoapi_keep_files = True
autoapi_add_toctree_entry = False

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# html_theme = 'sphinxdoc'
# html_static_path = ['_static']
# #html_css_files = ['custom.css']
# html_theme_options = { 'sidebarwidth': '28em' }

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['custom.css']

html_theme_options = {
    'collapse_navigation': False,
}


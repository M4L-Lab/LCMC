import os
import sys

# Add the 'lcmc' folder to the Python path so Sphinx can find your code
sys.path.insert(0, os.path.abspath('../lcmc'))

project = 'LCMC'
copyright = '2023, Rajib K. Musa'
author = 'Rajib K. Khan'
release = '1.0.0'

extensions = ['sphinx.ext.autodoc']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


html_theme = 'alabaster'
html_static_path = ['_static']
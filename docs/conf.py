import os, sys
sys.path.insert(0, os.path.abspath(".."))

project = "ndvi2gif"
author = "Diego García Díaz"
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "myst_nb",            # si prefieres nbsphinx: "nbsphinx"
    "sphinx_copybutton",
]
html_theme = "pydata_sphinx_theme"

html_theme_options = {
    "github_url": "https://github.com/tuusuario/ndvi2gif",  # Cambia por tu repo
    "show_toc_level": 2,
    "navbar_align": "left",
}

# CSS personalizado para notebooks
html_css_files = [
    "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.1.1/css/all.min.css",
]

# Docstrings estilo NumPy
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# Autosummary & Autodoc
autosummary_generate = True
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "inherited-members": False,
    # no mostramos privados en la doc pública:
    "private-members": False,
}

# ¡importantes! para que RTD compile sin deps pesadas/EE
autodoc_mock_imports = [
    "ee", "geemap", "geopandas", "fiona", "rasterio",
    "shapely", "seaborn", "matplotlib", "scipy"
]

# Notebooks: no ejecutar en RTD
nb_execution_mode = "off"        # myst-nb
# Si usas nbsphinx: nbsphinx_execute = 'never'

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", {}),
    "numpy": ("https://numpy.org/doc/stable/", {}),
    "matplotlib": ("https://matplotlib.org/stable/", {}),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", {}),
}
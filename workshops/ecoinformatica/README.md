# Taller Ndvi2Gif — Jornadas de Ecoinformática

**Sevilla, 2 de octubre de 2026, 16:00-17:00 (Taller IV)**

Material del taller. Tres notebooks, todos sobre Doñana, pensados para abrirse en Google
Colab sin instalar nada en tu ordenador.

## Antes de empezar

Hace falta una cuenta de Google Earth Engine con un proyecto de Google Cloud asociado. Si
no la tienes, la guía de esta carpeta lo explica paso a paso:

- [`cuenta-earth-engine.pdf`](cuenta-earth-engine.pdf)

En cada notebook hay una línea que tienes que cambiar por tu proyecto:

```python
PROJECT = 'tu-proyecto-de-earth-engine'
```

En Colab, además, descomenta la primera celda (`!pip install -q ndvi2gif`) y la línea
`ee.Authenticate()` la primera vez.

## Los notebooks

| | notebook | de qué va |
|---|---|---|
| 1 | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Digdgeo/Ndvi2Gif/blob/master/workshops/ecoinformatica/notebooks/01_rois_y_estadisticos.ipynb) [`01_rois_y_estadisticos.ipynb`](notebooks/01_rois_y_estadisticos.ipynb) | Las cinco formas de dar un ROI, el estadístico como decisión, y cuatro décadas de inundación de la marisma con MNDWI |
| 2 | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Digdgeo/Ndvi2Gif/blob/master/workshops/ecoinformatica/notebooks/02_hidroperiodo.ipynb) [`02_hidroperiodo.ipynb`](notebooks/02_hidroperiodo.ipynb) | `HydroperiodAnalyzer`: días de inundación, anomalías entre ciclos y fiabilidad temporal (IRT) |
| 3 | [![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Digdgeo/Ndvi2Gif/blob/master/workshops/ecoinformatica/notebooks/03_fenologia.ipynb) [`03_fenologia.ipynb`](notebooks/03_fenologia.ipynb) | `SpatialPhenologyAnalyzer`: SOS, POS, EOS y LOS a caballo del Guadalquivir — el arrozal de Isla Mayor y el mosaico de cultivos de Lebrija, dos calendarios en la misma imagen |

El primero está hecho a partir del notebook de ejemplos `ndvi2gif extended version` del
repositorio; los otros dos son versiones cortas, en castellano, de tutoriales que están
completos (y en inglés) en el libro.

Lo que no da tiempo a ver en el taller, pero está en el libro con el mismo nivel de
detalle: **luces nocturnas** (DMSP + VIIRS sobre Doñana, y los cambios de alumbrado que
fingen una cosecha), **incendios** (cuarenta años de Landsat sobre la península, y por qué
la serie entre décadas mide el archivo y no el monte), clasificación multisensor, SAR y
calidad de aguas.

## Para seguir después del taller

- **El libro**: https://digdgeo.github.io/Ndvi2Gif/ — tutoriales, referencia de los más de
  40 índices, opciones de ROI, clasificación, SAR, luces nocturnas e incendios.
- **El repositorio**: https://github.com/Digdgeo/Ndvi2Gif
- **Instalación local**: `pip install ndvi2gif` o `conda install -c conda-forge ndvi2gif`
- **El artículo** (JOSS, 2026): https://doi.org/10.21105/joss.10654
- Dudas y problemas: https://github.com/Digdgeo/Ndvi2Gif/issues

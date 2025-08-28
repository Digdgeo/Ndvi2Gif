import os
import pytest
import matplotlib

# Forzar backend no interactivo para CI/headless
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

# Cerrar figuras automáticamente tras cada test
@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")

# Omitir tests marcados como "ee" salvo que se pida explícitamente
def pytest_collection_modifyitems(config, items):
    if os.environ.get("NDVI2GIF_RUN_EE_TESTS") == "1":
        return
    skip_ee = pytest.mark.skip(
        reason="Set NDVI2GIF_RUN_EE_TESTS=1 to run Earth Engine-dependent tests."
    )
    for item in items:
        if "ee" in item.keywords:
            item.add_marker(skip_ee)
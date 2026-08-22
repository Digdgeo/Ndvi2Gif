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


# ---------------------------------------------------------------------
# Earth Engine session
# ---------------------------------------------------------------------
# Building an NdviSeasonality or an S1ARDProcessor already talks to the Earth
# Engine client, so even the tests that make no EE calls of their own need an
# initialized session. A bare ee.Initialize() only works when the stored
# credentials carry a Cloud project; otherwise the project has to come from
# EARTHENGINE_PROJECT or GOOGLE_CLOUD_PROJECT.

# Tests that run without any Earth Engine session at all
NO_EE_TESTS = {
    "test_public_api_imports",
    "test_direct_class_imports",
    "test_timeseries_analyzer_analyze_trend_on_dataframe",
}


def _init_earthengine():
    """Initialize Earth Engine if needed. Returns (ready, reason)."""
    try:
        import ee
    except ImportError as e:
        return False, str(e)

    try:
        ee.Number(1).getInfo()  # already initialized elsewhere
        return True, None
    except Exception:
        pass

    project = os.environ.get("EARTHENGINE_PROJECT") or os.environ.get(
        "GOOGLE_CLOUD_PROJECT"
    )
    attempts = [{"project": project}] if project else []
    attempts.append({})

    reason = None
    for kwargs in attempts:
        try:
            ee.Initialize(**kwargs)
            ee.Number(1).getInfo()
            return True, None
        except Exception as e:
            reason = str(e)
    return False, reason


# Omitir tests marcados como "ee" salvo que se pida explícitamente, y los que
# necesitan cliente de Earth Engine cuando no hay sesión disponible
def pytest_collection_modifyitems(config, items):
    ready, reason = _init_earthengine()

    skip_ee = pytest.mark.skip(
        reason="Set NDVI2GIF_RUN_EE_TESTS=1 to run Earth Engine-dependent tests."
    )
    skip_no_session = pytest.mark.skip(
        reason=(
            f"Earth Engine session not available ({reason}). Authenticate with "
            "`earthengine authenticate`, or set EARTHENGINE_PROJECT to a "
            "registered Cloud project."
        )
    )
    run_ee_tests = os.environ.get("NDVI2GIF_RUN_EE_TESTS") == "1"

    for item in items:
        if "ee" in item.keywords and not run_ee_tests:
            item.add_marker(skip_ee)
        elif not ready and item.name not in NO_EE_TESTS:
            item.add_marker(skip_no_session)

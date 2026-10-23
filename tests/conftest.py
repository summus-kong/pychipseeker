from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


@pytest.fixture
def examples_dir() -> Path:
    return EXAMPLES

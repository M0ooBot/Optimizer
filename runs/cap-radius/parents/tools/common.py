"""Shared helpers for the parent tools."""

import importlib.util
import sys
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
PARENTS = TOOLS.parent
REPO = PARENTS.parents[2]
CONVERTER_SRC = Path("/home/m0obot/repos/Airfoil-Converter/src")
SHAPES = ["flat", "twisted", "swept55", "fin40", "fin5", "wingb"]
VARIANTS = ["sharp", "blunt"]


def load_parent(shape, variant):
    path = PARENTS / shape / variant / "parent.py"
    name = f"parent_{shape}_{variant}"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def converter():
    if str(CONVERTER_SRC) not in sys.path:
        sys.path.insert(0, str(CONVERTER_SRC))
    from airfoil_converter import geometry
    return geometry

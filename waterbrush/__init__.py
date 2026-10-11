"""waterbrush: a deterministic, physics-based watercolor simulator.

Model: Curtis et al. 1997, "Computer-Generated Watercolor" (shallow-water
layer, pigment-deposition layer, capillary layer), composited with
Kubelka-Munk optics. Physical parameters and their sources are listed in
notes/research.md.
"""

__version__ = "0.0.1"

from .paper import Paper
from .pigment import Pigment, BURNT_SIENNA
from .sim import Sim, Params

__all__ = ["Paper", "Pigment", "BURNT_SIENNA", "Sim", "Params"]

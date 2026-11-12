"""Compatibility alias for the neutral :mod:`analysis_data.core` boundary."""

import sys

from analysis_data import core as _core

# Preserve attribute patching and private helper imports used by older clients.
# A re-export wrapper would create a second module namespace and silently make
# patches target the wrong persistence globals.
sys.modules[__name__] = _core

"""Deprecated: use ``import roadrunner_nitro as roadrunner``, which is libroadrunner's Python module on the Rust
engine. This alias will be removed in the next release."""
import warnings

from roadrunner_nitro import *  # noqa: F401,F403
from roadrunner_nitro import RoadRunner  # noqa: F401
from .roadrunner.extended_roadrunner import turbo  # noqa: F401

warnings.warn("tellurium_nitro.nitro is deprecated: use 'import roadrunner_nitro as roadrunner'.",
              DeprecationWarning, stacklevel=2)

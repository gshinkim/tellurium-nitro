"""Drop-in replacement for libroadrunner backed by the Rust engine: ``import tellurium.nitro as roadrunner``.

Everything is libroadrunner's, except ``RoadRunner``: it runs simulate/steadyState/MCA/sensitivities in
roadrunner_rs and falls back to roadrunner for anything the Rust engine can't do exactly the same way.
"""
from roadrunner import *  # noqa: F401,F403
from .roadrunner.extended_roadrunner import ExtendedRoadRunner as RoadRunner, turbo  # noqa: F401

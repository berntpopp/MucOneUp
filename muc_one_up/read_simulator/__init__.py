#!/usr/bin/env python3
"""
MUC1 read simulation package.

This package provides a modular implementation of the read simulation pipeline
for the MucOneUp project, handling the generation of simulated sequencing reads
from MUC1 gene region sequences with VNTR variants.
"""


def __getattr__(name: str):
    """Lazy import to avoid pulling in heavy dependencies at import time."""
    if name == "simulate_reads_pipeline":
        from .pipeline import simulate_reads_pipeline

        return simulate_reads_pipeline
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = ["simulate_reads_pipeline"]

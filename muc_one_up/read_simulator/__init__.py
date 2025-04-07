#!/usr/bin/env python3
"""
MUC1 read simulation package.

This package provides a modular implementation of the read simulation pipeline
for the MucOneUp project, handling the generation of simulated sequencing reads
from MUC1 gene region sequences with VNTR variants.
"""

from .pipeline import simulate_reads_pipeline

__all__ = ["simulate_reads_pipeline"]

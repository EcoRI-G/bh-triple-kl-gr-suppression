"""
Exploratory parameter-space survey of hierarchical black hole triples.

This script performs a small grid of N-body simulations using:

- REBOUND (WHFast integrator)
- REBOUNDx (general relativistic precession)
- Orbit-averaged gravitational-wave decay (Peters 1964)

Purpose:
--------
To identify regions of parameter space where Kozai–Lidov oscillations
lead to high eccentricity and potential merger, and to compare these
outcomes with and without relativistic precession.

Output:
-------
bh_triple_survey_results_gw_fast.csv
"""

from __future__ import annotations

import math
import itertools
from dataclasses import dataclass, asdict
from typing import List

import numpy as np
import pandas as pd

import rebound
import reboundx


# ============================================================
# Constants
# ============================================================

C_AU_PER_YR = 63241.077
G_AU3_YR2_MSUN = 4.0 * math.pi**2


# ============================================================
# Data structures
# ============================================================

@dataclass
class TripleIC:
    """Initial conditions for survey runs."""

    m1: float = 30.0
    m2: float = 20.0
    m3: float = 10.0

    a_in: float = 0.2
    e_in: float = 0.1

    a_out: float = 3.0
    e_out: float = 0.4
    inc_out_deg: float = 80.0


# ============================================================
# Simulation setup
# ============================================================

def build_simulation(ic: TripleIC, use_gr: bool) -> rebound.Simulation:
    """Construct simulation."""

    sim = rebound.Simulation()
    sim.units = ("AU", "yr", "Msun")
    sim.integrator = "whfast"
    sim.dt = 0.01

    sim.add(m=ic.m1)
    sim.add(m=ic.m2, a=ic.a_in, e=ic.e_in)
    sim.add(
        m=ic.m3,
        a=ic.a_out,
        e=ic.e_out,
        inc=math.radians(ic.inc_out_deg),
    )

    sim.move_to_com()

    if use_gr:
        rebx = reboundx.Extras(sim)
        gr = rebx.load_force("gr_full")
        gr.params["c"] = C_AU_PER_YR
        rebx.add_force(gr)

    return sim

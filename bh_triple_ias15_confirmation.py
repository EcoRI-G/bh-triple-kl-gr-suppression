"""
IAS15 confirmation integrations for hierarchical black hole triples.

This script performs high-accuracy integrations of selected representative
systems using:

- REBOUND N-body integrator
- IAS15 adaptive high-order integrator
- REBOUNDx general relativistic precession (gr_full)
- Orbit-averaged gravitational-wave decay (Peters 1964)

Purpose:
--------
To confirm the dynamical behaviour identified in the exploratory survey,
specifically the suppression of Kozai–Lidov-driven mergers by relativistic
precession.

Output:
-------
bh_triple_ias15_confirmation_results.csv

Dependencies:
-------------
numpy, pandas, rebound, reboundx
"""

from __future__ import annotations

import math
from dataclasses import dataclass, asdict
from typing import Dict, List

import numpy as np
import pandas as pd

import rebound
import reboundx


# ============================================================
# Constants
# ============================================================

C_AU_PER_YR = 63241.077        # Speed of light in AU/yr
G_AU3_YR2_MSUN = 4.0 * math.pi**2
EPS = 1e-14


# ============================================================
# Data structures
# ============================================================

@dataclass
class TripleIC:
    """Initial conditions for a hierarchical triple system."""

    label: str

    m1: float = 30.0
    m2: float = 20.0
    m3: float = 10.0

    a_in: float = 0.2
    e_in: float = 0.1

    a_out: float = 8.0
    e_out: float = 0.4
    inc_out_deg: float = 90.0

    Omega_in_deg: float = 0.0
    omega_in_deg: float = 0.0
    f_in_deg: float = 0.0

    Omega_out_deg: float = 0.0
    omega_out_deg: float = 0.0
    f_out_deg: float = 0.0


@dataclass
class RunResult:
    """Container for simulation output."""

    label: str
    physics: str
    gw_on: int

    initial_stable: int
    stability_ratio: float
    stability_ratio_crit: float

    merged: int
    merger_time_yr: float
    max_e_in: float
    min_r_in_au: float
    final_a_in_au: float
    final_e_in: float
    t_end_yr: float
    status: str

    m1: float
    m2: float
    m3: float
    a_in: float
    e_in: float
    a_out: float
    e_out: float
    inc_out_deg: float


# ============================================================
# Helper functions
# ============================================================

def deg(x: float) -> float:
    """Convert degrees to radians."""
    return math.radians(x)


def inner_binary_state(sim: rebound.Simulation) -> Dict[str, float]:
    """Return instantaneous separation and orbital elements."""
    p1, p2 = sim.particles[0], sim.particles[1]

    r_vec = np.array([p2.x - p1.x, p2.y - p1.y, p2.z - p1.z])
    v_vec = np.array([p2.vx - p1.vx, p2.vy - p1.vy, p2.vz - p1.vz])

    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)

    mu = G_AU3_YR2_MSUN * (p1.m + p2.m)
    energy = 0.5 * v * v - mu / r

    a = -mu / (2.0 * energy) if abs(energy) > EPS else np.inf

    h_vec = np.cross(r_vec, v_vec)
    e_vec = np.cross(v_vec, h_vec) / mu - r_vec / r
    e = np.linalg.norm(e_vec)

    return {"r": r, "a": a, "e": e}


def peters_decay_step(sim: rebound.Simulation, dt: float) -> None:
    """Apply orbit-averaged GW decay to the inner binary."""
    p1, p2 = sim.particles[0], sim.particles[1]
    state = inner_binary_state(sim)

    a, e = state["a"], state["e"]

    if not np.isfinite(a) or not np.isfinite(e):
        return

    factor = (G_AU3_YR2_MSUN**3 * p1.m * p2.m * (p1.m + p2.m)) / (C_AU_PER_YR**5)

    da = -(64.0 / 5.0) * factor / (a**3 * (1 - e**2)**3.5) * dt
    de = -(304.0 / 15.0) * factor * e / (a**4 * (1 - e**2)**2.5) * dt

    a_new = max(a + da, 1e-8)
    e_new = min(max(e + de, 0.0), 0.999999)

    sim.particles[1].a = a_new
    sim.particles[1].e = e_new


# ============================================================
# Simulation
# ============================================================

def build_simulation(ic: TripleIC, use_gr: bool) -> rebound.Simulation:
    """Construct simulation with optional relativistic precession."""

    sim = rebound.Simulation()
    sim.units = ("AU", "yr", "Msun")
    sim.integrator = "ias15"

    sim.add(m=ic.m1)
    sim.add(m=ic.m2, a=ic.a_in, e=ic.e_in)
    sim.add(
        m=ic.m3,
        a=ic.a_out,
        e=ic.e_out,
        inc=deg(ic.inc_out_deg),
    )

    sim.move_to_com()

    if use_gr:
        rebx = reboundx.Extras(sim)
        gr = rebx.load_force("gr_full")
        gr.params["c"] = C_AU_PER_YR
        rebx.add_force(gr)

    return sim

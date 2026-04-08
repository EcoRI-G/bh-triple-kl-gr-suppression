# Relativistic Suppression of Kozai–Lidov-Driven Mergers in Hierarchical Black Hole Triples

This repository contains the simulation scripts and data products used in the study:

**“Relativistic Suppression of Kozai–Lidov-Driven Mergers in Hierarchical Black Hole Triples”**

---

## Overview

This project investigates the dynamical evolution of hierarchical black hole triple systems using direct N-body simulations. The primary focus is on the competition between:

- Kozai–Lidov (KL) oscillations  
- Relativistic apsidal precession  

The simulations demonstrate that systems which undergo mergers under Newtonian dynamics can remain non-merging when post-Newtonian effects are included.

---

## Repository Contents

### Simulation Scripts

- `bh_triple_survey_gw_fast.py`  
  Exploratory parameter-space survey using the WHFast integrator

- `bh_triple_ias15_confirmation.py`  
  High-accuracy IAS15 integrations for representative systems

---

### Data Files

- `bh_triple_survey_results_gw_fast.csv`  
  Output of the exploratory survey

- `bh_triple_ias15_confirmation_results.csv`  
  Output of the IAS15 confirmation runs

---

## Requirements

Install dependencies using:

```bash
pip install numpy pandas matplotlib rebound reboundx

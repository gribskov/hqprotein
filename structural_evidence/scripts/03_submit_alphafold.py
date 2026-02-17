#!/usr/bin/env python3
#
# 03_submit_alphafold.py - Generate AlphaFold input, submit jobs, collect results
# Alex Baena, Feb 2026
#
# Not needed yet. S. cerevisiae already has pre-computed AlphaFold
# structures from the AlphaFold Database (6,055 CIF files).
#
# This becomes necessary when running predictions on new protein
# sequences (e.g. corrupted models, novel organisms). Will handle
# input generation, job submission (SLURM/Boltz/ColabFold), and
# collecting CIF output for 04_parse_cif.py.
#

# tric-compbio-toolkit

A collection of Python tools developed during my PhD for the computational
analysis of Trimeric Intracellular Cation (TRIC) channels — spanning structure
modeling, molecular dynamics (MD) analysis, water/hydrogen-bond analysis,
sequence analysis, and database processing.

These scripts automate large-scale, repetitive computational-biology workflows:
building and evaluating thousands of homology models, quantifying ion and water
permeation through membrane channels, profiling pore geometry over time, and
running comparative sequence analysis across taxa.

> **Note on paths:** these were written for a specific research environment, so
> file paths and system-specific parameters are defined at the top of each
> script. Set them to your own paths before running. The value here is the
> analysis logic, which is reusable across MD and structural-bioinformatics
> projects.

## Contents

### structure-modeling/
Automation around MODELLER for large-scale homology modeling.
- `every_required_folder.py` — generates per-target folders and auto-writes MODELLER `align2d`/`AutoModel` scripts for batch modeling
- `evaluate.py` — parses model logs to collect scoring (e.g. SOAP/DOPE) across many models
- `execute.py` — batch execution driver for the modeling pipeline

### md-analysis/
Analysis of GROMACS MD trajectories.
- `permeation.py` — tracks ion/water permeation events through the channel pore
- `ion_permeation_calculation.py` — generates per-ion index files and extracts center-of-mass trajectories for K⁺/Cl⁻ permeation
- `temporal_radius_profile.py` — computes pore-radius profiles over time (temporal HOLE-style heatmaps)
- `helix_kink_angle.py` — measures transmembrane helix kink/crossing geometry
- `apo_halo.py` — comparative analysis between apo and bound (holo) states

### water-analysis/
Hydration and hydrogen-bond analysis inside the pore.
- `water_analysis.py` — water occupancy/hydration profiling along the permeation pathway
- `average_water_molecule.py` — average water positions/counts within the channel
- `RING.py` — residue interaction / hydrogen-bond network analysis

### sequence-analysis/
Comparative sequence analysis across TRIC homologs.
- `cluster_specific_pairwise_alignment.py` — cluster-wise pairwise alignment across taxonomic groups

### database/
- `database_analysis.py` — processing and analysis supporting the dbTRIC database

## Requirements
See `requirements.txt`. Core dependencies: NumPy, pandas, Matplotlib, MDAnalysis, Biopython.
External tools used by some scripts: GROMACS, MODELLER, EMBOSS.

## Related work
Companion to [dbTRIC](https://github.com/gigity-gigity/dbTRIC), a full-stack web
database for TRIC channels.

## Author
**Prashant Upadhyay** — PhD candidate, Biological Sciences & Bioengineering,
IIT Kanpur · pupadhyay21@iitk.ac.in


###################################################################################################################################
############################ calculate the ion permeation event for the potassium and chloride ion through the channel ############
###################################################################################################################################
#
# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-monomer ION permeation (K+ / Cl-) through a TRIC channel (MDAnalysis only).
Same physical definition as the corrected water script, so the ion and water
analyses are directly comparable:

  * pore axis = pore-lining residues per subunit (set PROTOMER_SELS), NOT the
    whole-subunit centroid (which is several A off the pathway);
  * NARROW cylinder R_PORE ~ 5-6 A;
  * CENTRE-crossing test: an event fires only if the ion is inside the narrow
    cylinder while |dz| < MID_BAND (~3 A) -> single-file transit, not a lipid/
    edge crossing;
  * time-continuous z (no per-frame PBC fold flicker);
  * event-producing ions highlighted in the plot, over all ions drawn faint,
    so figure and .dat show the SAME ions.

Ions are few, so this is a single pass over ALL ions (the two-pass candidate
filter used for water is a water-scale optimisation and is not needed here).

Validation for a homotrimer: the three pores should give ROUGHLY EQUAL counts.

Requires: MDAnalysis, numpy, matplotlib.
"""

import os
import numpy as np
import MDAnalysis as mda
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG ============================
GRO = "/media/supremeleader/Pantera/simulation/lipid_simulation/5WTR/charmm-gui-7977600322/gromacs/5wtr/center.gro"
XTC = "/media/supremeleader/Pantera/simulation/lipid_simulation/5WTR/charmm-gui-7977600322/gromacs/5wtr/center.xtc"
OUTDIR = "/media/supremeleader/Pantera/simulation/lipid_simulation/5WTR/charmm-gui-7977600322/gromacs/5wtr/ion_permeation"
SYSTEM_NAME = "5WTR"

# CHARMM naming. GROMACS/AMBER ports: "resname K" / "resname CL", etc.
SPECIES = {
    "Potassium": {"sel": "resname POT", "out": "potassium"},
    "Chloride":  {"sel": "resname CLA", "out": "chloride"},
}

PROTEIN_SEL = "protein"
# Use the SAME pore-lining residues you set for the water script:
#   PROTOMER_SELS = ["segid PROA and resid 30 63 66 129 132",
#                    "segid PROB and resid 30 63 66 129 132",
#                    "segid PROC and resid 30 63 66 129 132"]
# None -> falls back to subunit CA centroid (approximate; prints a warning).
PROTOMER_SELS = None
N_PORES = 3

Z_GATE    = 15.0                 # A : |dz| > Z_GATE == committed to a reservoir
R_PORE    = 6.0                  # A : NARROW pore cylinder radius (TUNE ~5-6)
MID_BAND  = 3.0                  # A : must be inside the cylinder while |dz| < MID_BAND to count
LAST_NS   = 500.0                # analyse the final LAST_NS; set None for the whole trajectory
SAMPLE_PS = 10.0                 # ions are few -> keep fine sampling
YLIM      = (-30, 30)
# ================================================================


def whole_mean(coord, L):
    ref = np.median(coord)
    return (coord - L * np.round((coord - ref) / L)).mean()


def build_pore_groups(u):
    prot = u.select_atoms(f"({PROTEIN_SEL}) and name CA")
    if len(prot) == 0:
        raise RuntimeError("No protein CA atoms found - check PROTEIN_SEL.")
    if PROTOMER_SELS:
        groups = [g for g in (u.select_atoms(s) for s in PROTOMER_SELS) if len(g)]
        if len(groups) != N_PORES:
            raise RuntimeError(f"PROTOMER_SELS gave {len(groups)} non-empty groups, "
                               f"expected {N_PORES}.")
        return prot, groups
    print("WARNING: PROTOMER_SELS is None -> subunit CA centroids used as pore axes "
          "(approximate). Set pore-lining residues for trustworthy per-monomer counts.")
    seen, uniq = set(), []
    for seg in prot.segments:
        if len(seg.atoms.select_atoms("name CA")) and seg.segid not in seen:
            seen.add(seg.segid); uniq.append(seg.segid)
    if len(uniq) == N_PORES:
        return prot, [prot.select_atoms(f"segid {sid}") for sid in uniq]
    idx = np.array_split(np.arange(prot.n_atoms), N_PORES)
    return prot, [prot[chunk] for chunk in idx]


def frame_slice(u):
    n = len(u.trajectory)
    t0 = u.trajectory[0].time
    dt = (u.trajectory[1].time - t0) if n > 1 else SAMPLE_PS
    stride = max(1, int(round(SAMPLE_PS / dt))) if dt > 0 else 1
    if LAST_NS is None:
        return 0, stride, dt
    t_end = t0 + dt * (n - 1)
    return max(0, int(np.ceil((t_end - LAST_NS * 1000.0 - t0) / dt))), stride, dt


def pore_axes(groups, box):
    Lx, Ly = box[0], box[1]
    return np.array([[whole_mean(g.positions[:, 0], Lx),
                      whole_mean(g.positions[:, 1], Ly)] for g in groups])


def xy_dist(pos, axes, box):
    Lx, Ly = box[0], box[1]
    dx = pos[:, 0][:, None] - axes[:, 0][None, :]
    dy = pos[:, 1][:, None] - axes[:, 1][None, :]
    dx -= Lx * np.round(dx / Lx)
    dy -= Ly * np.round(dy / Ly)
    return np.sqrt(dx * dx + dy * dy)


def detect_events(dz_ser, r_ser, pore_ser, times):
    """reservoir -> pore CENTRE (inside narrow cylinder, |dz|<MID_BAND) -> opposite reservoir."""
    events = []
    state = 0
    passed = False
    pore_at_mid = -1
    for dz, r, p, t in zip(dz_ser, r_ser, pore_ser, times):
        if (r < R_PORE) and (abs(dz) < MID_BAND):
            passed = True
            pore_at_mid = int(p)
        side = 1 if dz > Z_GATE else (-1 if dz < -Z_GATE else 0)
        if side != 0 and side != state:
            if state != 0 and passed:
                events.append((t, "down" if state == 1 else "up", pore_at_mid))
            passed = False
            state = side
    return events


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    u = mda.Universe(GRO, XTC)
    prot, groups = build_pore_groups(u)
    ions = {name: u.select_atoms(s["sel"]) for name, s in SPECIES.items()}
    for name, ag in ions.items():
        print(f"{name}: {ag.n_atoms} ions")
    start_frame, stride, dt = frame_slice(u)
    npore = len(groups)
    print(f"pore groups: {npore} ({[g.n_atoms for g in groups]} atoms) | "
          f"R_PORE={R_PORE} A, MID_BAND=+/-{MID_BAND} A | from frame {start_frame}, "
          f"every {stride*dt:.0f} ps")

    # ------- single pass: continuous z + per-pore distance for every ion -------
    times, cz_chk = [], []
    dz_all = {n: [] for n in SPECIES}
    r_all = {n: [] for n in SPECIES}
    pore_all = {n: [] for n in SPECIES}
    prev_z = {n: None for n in SPECIES}
    for ts in u.trajectory[start_frame::stride]:
        box = ts.dimensions; Lz = box[2]
        cz = whole_mean(prot.positions[:, 2], Lz); cz_chk.append(cz)
        axes = pore_axes(groups, box)
        for name, ag in ions.items():
            pos = ag.positions
            zraw = pos[:, 2]
            z = zraw - Lz * np.round((zraw - (cz if prev_z[name] is None else prev_z[name])) / Lz)
            prev_z[name] = z
            rmat = xy_dist(pos, axes, box)
            dz_all[name].append(z - cz)
            r_all[name].append(rmat.min(axis=1))
            pore_all[name].append(rmat.argmin(axis=1))
        times.append(ts.time / 1000.0)
    times = np.asarray(times)
    cz_chk = np.asarray(cz_chk)
    print(f"channel centre cz: mean={cz_chk.mean():.2f} A  std={cz_chk.std():.2f} A")

    combined = {}
    for name, s in SPECIES.items():
        dz = np.asarray(dz_all[name]); r = np.asarray(r_all[name]); pr = np.asarray(pore_all[name])
        resids = ions[name].resids
        nion = dz.shape[1]
        per_pore = {p: {"down": 0, "up": 0} for p in range(npore)}
        events, event_cols = [], set()
        for i in range(nion):
            ev = detect_events(dz[:, i], r[:, i], pr[:, i], times)
            if ev:
                event_cols.add(i)
            for (t, d, p) in ev:
                events.append((resids[i], t, d, p))
                if p in per_pore:
                    per_pore[p][d] += 1
        events.sort(key=lambda e: e[1])
        ndown = sum(1 for e in events if e[2] == "down")
        nup = sum(1 for e in events if e[2] == "up")
        combined[name] = (per_pore, len(events), ndown, nup)

        # ---- plot: all ions faint, permeating ions highlighted ----
        plt.rcParams.update({"font.family": "sans-serif", "font.size": 12})
        fig, ax = plt.subplots(figsize=(10, 6))
        for i in range(nion):
            if i not in event_cols:
                ax.plot(times, dz[:, i], lw=0.3, color="0.82", rasterized=True)
        evc = sorted(event_cols)
        cmap = plt.get_cmap("turbo")
        for k, i in enumerate(evc):
            ax.plot(times, dz[:, i], lw=0.9, color=cmap(k / max(1, len(evc) - 1)), rasterized=True)
        ax.set_xlabel("Time (ns)"); ax.set_ylabel("Channel Z-coordinate (\u00c5)")
        ax.set_title(f"{name} permeation in {SYSTEM_NAME} ({len(evc)} permeating ions)")
        ax.set_xlim(times[0], times[-1])
        if YLIM: ax.set_ylim(*YLIM)
        ax.axhline(Z_GATE, ls="--", lw=0.6, color="0.4")
        ax.axhline(-Z_GATE, ls="--", lw=0.6, color="0.4")
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(OUTDIR, f"{s['out']}_{SYSTEM_NAME.lower()}.{ext}"),
                        dpi=200, bbox_inches="tight")
        plt.close(fig)

        # ---- events file ----
        with open(os.path.join(OUTDIR, f"{s['out']}_permeation_events.dat"), "w") as f:
            f.write(f"# {name} permeation in {SYSTEM_NAME}; "
                    f"window={times[0]:.1f}-{times[-1]:.1f} ns; z_gate=+/-{Z_GATE}; "
                    f"R_pore={R_PORE}; mid_band=+/-{MID_BAND}; {npore} pores\n")
            f.write(f"# total={len(events)} down={ndown} up={nup} net_down={ndown-nup}\n")
            f.write("# resid  time_ns  direction  pore_index\n")
            for rid, t, d, p in events:
                f.write(f"{rid} {t:.3f} {d} {p}\n")

        print(f"\n{name} per-monomer permeation:")
        for p in range(npore):
            d = per_pore[p]["down"]; up = per_pore[p]["up"]
            print(f"  pore {p}: total {d+up}  (down {d}, up {up}, net_down {d-up})")
        print(f"  ALL   : total {len(events)}  (down {ndown}, up {nup}, net_down {ndown-nup})")
        tot = [per_pore[p]['down'] + per_pore[p]['up'] for p in range(npore)]
        if max(tot) > 0 and min(tot) < 0.5 * (sum(tot) / npore):
            print("  NOTE: uneven across a homotrimer -> check PROTOMER_SELS / R_PORE.")

    with open(os.path.join(OUTDIR, "permeation_summary.dat"), "w") as f:
        f.write(f"# {SYSTEM_NAME} window={times[0]:.1f}-{times[-1]:.1f} ns  z_gate=+/-{Z_GATE}  "
                f"R_pore={R_PORE}  mid_band=+/-{MID_BAND}  cz_std={cz_chk.std():.2f}\n")
        f.write("# species  pore  down  up  total  net_down\n")
        for name, (per_pore, tot, nd, nu) in combined.items():
            for p in range(npore):
                d = per_pore[p]["down"]; up = per_pore[p]["up"]
                f.write(f"{name} {p} {d} {up} {d+up} {d-up}\n")
            f.write(f"{name} ALL {nd} {nu} {tot} {nd-nu}\n")
    print(f"\nAll outputs in: {OUTDIR}")


if __name__ == "__main__":
    main()
#
# ###################################################################################################################################
# ############################ calculate the ion permeation event for the potassium and chloride ion through the channel ############
# ###################################################################################################################################
#
# ###################################################################################################################################
# ############################ calculate the water permeation through the channel ###################################################
# ###################################################################################################################################
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-monomer WATER permeation through a TRIC channel (MDAnalysis only, no GROMACS).
Last 500 ns by default.

FIX vs the fat-cylinder version that gave 99.995% candidates and an impossible
homotrimer asymmetry (21 / 246 / 50):

  * pore axis on the REAL pathway   : set PROTOMER_SELS to each subunit's
    pore-lining residues (not the whole-subunit centroid, which is several A off);
  * NARROW cylinder                 : R_PORE ~ 5-6 A instead of 12;
  * CENTRE-crossing requirement     : an event fires only if the water is inside
    the narrow cylinder while |dz| < MID_BAND (~3 A), i.e. as it crosses the pore
    centre - "somewhere in the membrane slab" no longer counts;
  * the plot now draws the EVENT-producing waters highlighted, so the figure and
    the .dat describe the SAME molecules.

Validation for a homotrimer: the three pores should give ROUGHLY EQUAL counts.
A >~2x spread means the axis/radius is still wrong.

Requires: MDAnalysis, numpy, matplotlib.
"""

import os
import numpy as np
import MDAnalysis as mda
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG ============================
GRO = "/media/supremeleader/Pantera/simulation/lipid_simulation/5WTR/charmm-gui-7977600322/gromacs/5wtr/center.gro"
XTC = "/media/supremeleader/Pantera/simulation/lipid_simulation/5WTR/charmm-gui-7977600322/gromacs/5wtr/center.xtc"
OUTDIR = "/media/supremeleader/Pantera/simulation/lipid_simulation/5WTR/charmm-gui-7977600322/gromacs/5wtr/water_permeation"
SYSTEM_NAME = "5WTR"

WATER_SEL = "resname TIP3 and name OH2"     # CHARMM TIP3 oxygen; GROMACS SOL -> "resname SOL and name OW"

PROTEIN_SEL = "protein"
# STRONGLY RECOMMENDED: one selection per subunit giving that subunit's PORE-LINING
# residues (the ones your HOLE profile passes through). Example placeholder below -
# replace resid lists with your actual constriction/lining residues per chain:
#   PROTOMER_SELS = ["segid PROA and resid 30 63 66 129 132",
#                    "segid PROB and resid 30 63 66 129 132",
#                    "segid PROC and resid 30 63 66 129 132"]
# If left None it falls back to the subunit CA centroid (approximate; prints a warning).
PROTOMER_SELS = None
N_PORES = 3

Z_GATE    = 15.0                 # A : |dz| > Z_GATE == committed to a reservoir
R_PORE    = 6.0                  # A : NARROW pore cylinder radius (was 12; TUNE ~5-6)
MID_BAND  = 3.0                  # A : must be inside the cylinder while |dz| < MID_BAND to count
Z_BAND    = 22.0                 # A : Pass-1 z-slab half-width (> Z_GATE)
LAST_NS   = 500.0
SAMPLE_PS = 20.0
YLIM      = (-30, 30)
CAND_WARN_FRAC = 0.20            # warn if more than this fraction of waters become candidates
SAMPLE_BG = 300                  # faint background traces for context in the plot
# ================================================================


def whole_mean(coord, L):
    ref = np.median(coord)
    return (coord - L * np.round((coord - ref) / L)).mean()


def build_pore_groups(u):
    prot = u.select_atoms(f"({PROTEIN_SEL}) and name CA")
    if len(prot) == 0:
        raise RuntimeError("No protein CA atoms found - check PROTEIN_SEL.")
    if PROTOMER_SELS:
        groups = [u.select_atoms(s) for s in PROTOMER_SELS]
        groups = [g for g in groups if len(g)]
        if len(groups) != N_PORES:
            raise RuntimeError(f"PROTOMER_SELS gave {len(groups)} non-empty groups, "
                               f"expected {N_PORES}. Check your selections.")
        return prot, groups
    print("WARNING: PROTOMER_SELS is None -> using whole-subunit CA centroids as pore "
          "axes. This is APPROXIMATE and off the true pathway by several A. For "
          "trustworthy per-monomer counts, set PROTOMER_SELS to pore-lining residues.")
    seen, uniq = set(), []
    for seg in prot.segments:
        if len(seg.atoms.select_atoms("name CA")) and seg.segid not in seen:
            seen.add(seg.segid); uniq.append(seg.segid)
    if len(uniq) == N_PORES:
        return prot, [prot.select_atoms(f"segid {sid}") for sid in uniq]
    idx = np.array_split(np.arange(prot.n_atoms), N_PORES)
    return prot, [prot[chunk] for chunk in idx]


def frame_slice(u):
    n = len(u.trajectory)
    t0 = u.trajectory[0].time
    dt = (u.trajectory[1].time - t0) if n > 1 else SAMPLE_PS
    stride = max(1, int(round(SAMPLE_PS / dt))) if dt > 0 else 1
    t_end = t0 + dt * (n - 1)
    start_frame = max(0, int(np.ceil((t_end - LAST_NS * 1000.0 - t0) / dt)))
    return start_frame, stride, dt


def pore_axes(groups, box):
    Lx, Ly = box[0], box[1]
    return np.array([[whole_mean(g.positions[:, 0], Lx),
                      whole_mean(g.positions[:, 1], Ly)] for g in groups])


def xy_dist(pos, axes, box):
    Lx, Ly = box[0], box[1]
    dx = pos[:, 0][:, None] - axes[:, 0][None, :]
    dy = pos[:, 1][:, None] - axes[:, 1][None, :]
    dx -= Lx * np.round(dx / Lx)
    dy -= Ly * np.round(dy / Ly)
    return np.sqrt(dx * dx + dy * dy)


def detect_events(dz_ser, r_ser, pore_ser, times):
    """reservoir -> pore CENTRE (inside narrow cylinder, |dz|<MID_BAND) -> opposite reservoir."""
    events = []
    state = 0
    passed = False
    pore_at_mid = -1
    for dz, r, p, t in zip(dz_ser, r_ser, pore_ser, times):
        if (r < R_PORE) and (abs(dz) < MID_BAND):
            passed = True
            pore_at_mid = int(p)
        side = 1 if dz > Z_GATE else (-1 if dz < -Z_GATE else 0)
        if side != 0 and side != state:
            if state != 0 and passed:
                events.append((t, "down" if state == 1 else "up", pore_at_mid))
            passed = False
            state = side
    return events


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    u = mda.Universe(GRO, XTC)
    prot, groups = build_pore_groups(u)
    waterO = u.select_atoms(WATER_SEL)
    if len(waterO) == 0:
        raise RuntimeError("No water oxygens found - check WATER_SEL.")
    start_frame, stride, dt = frame_slice(u)
    print(f"protein CA: {prot.n_atoms} | pore groups: {len(groups)} "
          f"({[g.n_atoms for g in groups]} atoms) | water O: {waterO.n_atoms}")
    print(f"window: last {LAST_NS:.0f} ns from frame {start_frame}, every {stride*dt:.0f} ps")
    print(f"pore def: R_PORE={R_PORE} A, MID_BAND=+/-{MID_BAND} A (centre-crossing test)")

    # ------------------------- PASS 1: candidate filter -------------------------
    cand_ix = set()
    cz_chk = []
    for ts in u.trajectory[start_frame::stride]:
        box = ts.dimensions; Lz = box[2]
        cz = whole_mean(prot.positions[:, 2], Lz); cz_chk.append(cz)
        axes = pore_axes(groups, box)
        wz = waterO.positions[:, 2]
        dz = wz - cz
        dz -= Lz * np.round(dz / Lz)
        near = np.abs(dz) < Z_BAND
        if not near.any():
            continue
        si = np.nonzero(near)[0]
        rmin = xy_dist(waterO.positions[si], axes, box).min(axis=1)
        interior = (rmin < R_PORE) & (np.abs(dz[si]) < MID_BAND)   # narrow + centre
        if interior.any():
            cand_ix.update(int(x) for x in waterO.ix[si[interior]])
    cz_chk = np.asarray(cz_chk)
    frac = len(cand_ix) / waterO.n_atoms
    print(f"channel centre cz: mean={cz_chk.mean():.2f} A  std={cz_chk.std():.2f} A")
    print(f"candidate waters: {len(cand_ix)} of {waterO.n_atoms}  ({100*frac:.3f}%)")
    if frac > CAND_WARN_FRAC:
        print(f"WARNING: {100*frac:.1f}% of waters are candidates -> cylinder still too "
              f"loose or axis off the pathway. Lower R_PORE / fix PROTOMER_SELS.")
    if not cand_ix:
        print("No water reached a pore centre -> 0 events.")
        return

    cand = u.atoms[sorted(cand_ix)]
    ncand = cand.n_atoms

    # ------------------------- PASS 2: track candidates -------------------------
    times, dz_all, r_all, pore_all = [], [], [], []
    prev_z = None
    for ts in u.trajectory[start_frame::stride]:
        box = ts.dimensions; Lz = box[2]
        cz = whole_mean(prot.positions[:, 2], Lz)
        axes = pore_axes(groups, box)
        pos = cand.positions
        zraw = pos[:, 2]
        z = zraw - Lz * np.round((zraw - (cz if prev_z is None else prev_z)) / Lz)
        prev_z = z
        rmat = xy_dist(pos, axes, box)
        dz_all.append(z - cz)
        r_all.append(rmat.min(axis=1))
        pore_all.append(rmat.argmin(axis=1))
        times.append(ts.time / 1000.0)
    times = np.asarray(times)
    dz_all = np.asarray(dz_all); r_all = np.asarray(r_all); pore_all = np.asarray(pore_all)

    # ------------------------------- counting -------------------------------
    npore = len(groups)
    resids = cand.resids
    per_pore = {p: {"down": 0, "up": 0} for p in range(npore)}
    all_events, event_cols = [], set()
    for i in range(ncand):
        ev = detect_events(dz_all[:, i], r_all[:, i], pore_all[:, i], times)
        if ev:
            event_cols.add(i)
        for (t, d, p) in ev:
            all_events.append((resids[i], t, d, p))
            if p in per_pore:
                per_pore[p][d] += 1
    all_events.sort(key=lambda e: e[1])
    ndown = sum(1 for e in all_events if e[2] == "down")
    nup = sum(1 for e in all_events if e[2] == "up")

    # -------- plot: event waters highlighted over a faint background sample --------
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 12})
    fig, ax = plt.subplots(figsize=(10, 6))
    rng = np.random.default_rng(0)
    bg = rng.choice(ncand, size=min(SAMPLE_BG, ncand), replace=False)
    for i in bg:
        ax.plot(times, dz_all[:, i], lw=0.3, color="0.82", rasterized=True)
    evc = sorted(event_cols)
    cmap = plt.get_cmap("turbo")
    for k, i in enumerate(evc):
        ax.plot(times, dz_all[:, i], lw=0.8,
                color=cmap(k / max(1, len(evc) - 1)), rasterized=True)
    ax.set_xlabel("Time (ns)"); ax.set_ylabel("Channel Z-coordinate (\u00c5)")
    ax.set_title(f"Water permeation in {SYSTEM_NAME} ({len(evc)} permeating waters highlighted)")
    ax.set_xlim(times[0], times[-1])
    if YLIM: ax.set_ylim(*YLIM)
    ax.axhline(Z_GATE, ls="--", lw=0.6, color="0.4")
    ax.axhline(-Z_GATE, ls="--", lw=0.6, color="0.4")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUTDIR, f"water_{SYSTEM_NAME.lower()}.{ext}"),
                    dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ------------------------------- outputs -------------------------------
    with open(os.path.join(OUTDIR, "water_permeation_events.dat"), "w") as f:
        f.write(f"# Water permeation in {SYSTEM_NAME}; last {LAST_NS:.0f} ns; "
                f"z_gate=+/-{Z_GATE} A; R_pore={R_PORE} A; mid_band=+/-{MID_BAND} A; "
                f"{npore} pores; candidates={ncand}\n")
        f.write(f"# total={len(all_events)} down={ndown} up={nup} net_down={ndown-nup}\n")
        f.write("# resid  time_ns  direction  pore_index\n")
        for rid, t, d, p in all_events:
            f.write(f"{rid} {t:.3f} {d} {p}\n")

    with open(os.path.join(OUTDIR, "water_permeation_summary.dat"), "w") as f:
        f.write(f"# {SYSTEM_NAME} last {LAST_NS:.0f} ns  z_gate=+/-{Z_GATE}  "
                f"R_pore={R_PORE}  mid_band=+/-{MID_BAND}  "
                f"window={times[0]:.1f}-{times[-1]:.1f} ns\n")
        f.write("# pore  down  up  total  net_down\n")
        for p in range(npore):
            d = per_pore[p]["down"]; up = per_pore[p]["up"]
            f.write(f"{p} {d} {up} {d+up} {d-up}\n")
        f.write(f"ALL {ndown} {nup} {len(all_events)} {ndown-nup}\n")

    print(f"\nper-monomer water permeation (last {LAST_NS:.0f} ns):")
    for p in range(npore):
        d = per_pore[p]["down"]; up = per_pore[p]["up"]
        print(f"  pore {p}: total {d+up}  (down {d}, up {up}, net_down {d-up})")
    print(f"  ALL   : total {len(all_events)}  (down {ndown}, up {nup}, net_down {ndown-nup})")
    tot = [per_pore[p]['down'] + per_pore[p]['up'] for p in range(npore)]
    if max(tot) > 0 and min(tot) < 0.5 * (sum(tot) / npore):
        print("NOTE: per-pore counts are uneven for a homotrimer -> axis/radius likely "
              "still imperfect; verify PROTOMER_SELS are the true pore-lining residues.")
    print(f"\nAll outputs in: {OUTDIR}")


if __name__ == "__main__":
    main()


# ###################################################################################################################################
# ############################ calculate the water permeation through the channel ###################################################
# ###################################################################################################################################

###################################################################################################################################
############################ contact analysis for few of the water molecules passing through the channel ##########################
###################################################################################################################################

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PERMEATION PATHWAY of each transiting water: the ORDERED, de-duplicated sequence of
residue interactions it forms as it moves along z through the pore, with the water's
z-coordinate at each step. CSV output.

Rules (as requested):
  * one entry per UNIQUE (residue, interaction-subtype), at its FIRST occurrence;
    the same interaction re-forming is NOT logged again;
  * a DIFFERENT subtype with the same residue IS a new entry
    (e.g. water donor to Y29 -OH, then later acceptor / electrostatic with Y29);
  * subtypes: hbond_wdonor, hbond_wacceptor, elec, polar_contact;
  * VDW / nonpolar is ignored;
  * a weak polar_contact is suppressed for a residue that also forms an H-bond or
    electrostatic contact during the same transit;
  * z_rel_A = water-O z relative to the channel centre (sign = which side of the pore).

Thornton criteria: HB D..A<=3.9, H..A<=2.5, angle>=90 (McDonald & Thornton);
ELEC <=4.0 A (Barlow & Thornton). GRO is sufficient (geometric only).

Requires: MDAnalysis, numpy.
"""

import os, csv
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance, calc_angles

# ============================ CONFIG ============================
GRO    = "/media/supremeleader/Pantera/simulation/lipid_simulation/5WTR/charmm-gui-7977600322/gromacs/5wtr/center.gro"
XTC    = "/media/supremeleader/Pantera/simulation/lipid_simulation/5WTR/charmm-gui-7977600322/gromacs/5wtr/center.xtc"
EVENTS = "/media/supremeleader/Pantera/simulation/lipid_simulation/5WTR/charmm-gui-7977600322/gromacs/5wtr/water_permeation/water_permeation_events.dat"
OUTDIR = "/media/supremeleader/Pantera/simulation/lipid_simulation/5WTR/charmm-gui-7977600322/gromacs/5wtr/water_permeation/contacts"

WATER_SEL   = "resname TIP3"
PROTEIN_SEL = "protein"

Z_TRANSIT   = 18.0
WIN_NS      = 2.0
CONTACT_PS  = 10.0

HB_DA_MAX   = 3.9
HB_HA_MAX   = 2.5
HB_ANG_MIN  = 90.0
SB_DIST     = 4.0
POLAR_CUT   = 3.5

SUPPRESS_POLAR_IF_BONDED = True

BACKBONE = {"N", "CA", "C", "O", "OT1", "OT2", "OXT", "HN", "H", "HA", "HA1", "HA2"}
STRONG = {"hbond_wdonor", "hbond_wacceptor", "elec"}
RANK = {"hbond_wdonor": 3, "hbond_wacceptor": 3, "elec": 2, "polar_contact": 1}

# ---- side-chain functional group by (resname, atom); backbone handled separately ----
SC_GROUP = {
    ("SER", "OG"): "hydroxyl", ("THR", "OG1"): "hydroxyl", ("TYR", "OH"): "phenol-hydroxyl",
    ("ASP", "OD1"): "carboxylate", ("ASP", "OD2"): "carboxylate",
    ("GLU", "OE1"): "carboxylate", ("GLU", "OE2"): "carboxylate",
    ("ASN", "OD1"): "amide-carbonyl", ("ASN", "ND2"): "amide-amino",
    ("GLN", "OE1"): "amide-carbonyl", ("GLN", "NE2"): "amide-amino",
    ("LYS", "NZ"): "amine", ("ARG", "NE"): "guanidinium", ("ARG", "NH1"): "guanidinium",
    ("ARG", "NH2"): "guanidinium", ("TRP", "NE1"): "indole-NH",
    ("CYS", "SG"): "thiol", ("MET", "SD"): "thioether",
}
HIS_N = {"ND1", "NE2"}


def func_group(resname, atom):
    if atom in BACKBONE:
        if atom in ("O", "OT1", "OT2", "OXT"):
            return "backbone-carbonyl"
        if atom in ("N", "HN", "H"):
            return "backbone-amide"
        return "backbone"
    if resname in ("HIS", "HSD", "HSE", "HSP", "HIP") and atom in HIS_N:
        return "imidazole-N"
    return SC_GROUP.get((resname, atom), "sidechain-polar")
# ================================================================


def whole_mean(coord, L):
    ref = np.median(coord)
    return (coord - L * np.round((coord - ref) / L)).mean()


def build_protein_donors(u, prot):
    u.trajectory[0]
    Hs = prot.select_atoms("name H*")
    heavy = prot.select_atoms("name N* O* S*")
    if not (len(Hs) and len(heavy)):
        return np.array([], int), np.array([], int)
    pairs, dists = capped_distance(Hs.positions, heavy.positions,
                                   max_cutoff=1.2, box=u.dimensions, return_distances=True)
    best = {}
    for (hi, di), dd in zip(pairs, dists):
        if hi not in best or dd < best[hi][1]:
            best[hi] = (di, dd)
    return (np.array([Hs[h].ix for h in best], int),
            np.array([heavy[best[h][0]].ix for h in best], int))


def frame_interactions(wat, wO, wH, box, acceptors, dH_H, dH_D, u, cations, anions, prot_polar):
    """strongest interaction per residue this frame -> {key:(subtype,patom,watom,dist,ang)}"""
    Opos = wO.positions; Hpos = wH.positions
    cand = {}
    def offer(typ, patom, watom, dist, ha, ang):
        k = (patom.segid, int(patom.resid), patom.resname)
        cur = cand.get(k)
        if cur is None or RANK[typ] > RANK[cur[0]] or (RANK[typ] == RANK[cur[0]] and dist < cur[3]):
            cand[k] = (typ, patom, watom, dist, ha, ang)

    if acceptors.n_atoms:
        pa_, da_ = capped_distance(Opos, acceptors.positions, max_cutoff=HB_DA_MAX,
                                   box=box, return_distances=True)
        for (_, aj), dDA in zip(pa_, da_):
            A = acceptors[aj]
            for h in range(len(wH)):
                ha = float(np.linalg.norm(Hpos[h] - A.position))
                if ha > HB_HA_MAX:
                    continue
                ang = np.degrees(calc_angles(Opos[0], Hpos[h], A.position, box=box))
                if ang >= HB_ANG_MIN:
                    offer("hbond_wdonor", A, wH[h].name, dDA, ha, ang)
                    break
    if len(dH_D):
        dpos = u.atoms[dH_D].positions
        pd_, dd2 = capped_distance(Opos, dpos, max_cutoff=HB_DA_MAX, box=box, return_distances=True)
        for (_, dj), dDA in zip(pd_, dd2):
            D = u.atoms[dH_D[dj]]; H = u.atoms[dH_H[dj]]
            ha = float(np.linalg.norm(H.position - Opos[0]))
            if ha > HB_HA_MAX:
                continue
            ang = np.degrees(calc_angles(D.position, H.position, Opos[0], box=box))
            if ang >= HB_ANG_MIN:
                offer("hbond_wacceptor", D, wO[0].name, dDA, ha, ang)
    if cations.n_atoms:
        pc, dc = capped_distance(Opos, cations.positions, max_cutoff=SB_DIST, box=box, return_distances=True)
        for (_, cj), r in zip(pc, dc):
            offer("elec", cations[cj], wO[0].name, r, None, None)
    if anions.n_atoms and len(wH):
        pan, dan = capped_distance(Hpos, anions.positions, max_cutoff=SB_DIST, box=box, return_distances=True)
        for (hi_, aj), r in zip(pan, dan):
            offer("elec", anions[aj], wH[hi_].name, r, None, None)
    pp, dp = capped_distance(Opos, prot_polar.positions, max_cutoff=POLAR_CUT, box=box, return_distances=True)
    for (_, pj), r in zip(pp, dp):
        offer("polar_contact", prot_polar[pj], wO[0].name, r, None, None)
    return cand


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    u = mda.Universe(GRO, XTC)
    prot = u.select_atoms(PROTEIN_SEL)
    prot_polar = prot.select_atoms("name N* O* S*")
    acceptors = prot.select_atoms("name O* OD1 OD2 OE1 OE2 OG OG1 OH "
                                  "or (resname HIS HSD HSE HSP and name ND1 NE2)")
    cations = prot.select_atoms("(resname LYS and name NZ) or "
                                "(resname ARG and name NE NH1 NH2) or "
                                "(resname HSP HIP and name ND1 NE2)")
    anions = prot.select_atoms("(resname ASP and name OD1 OD2) or "
                               "(resname GLU and name OE1 OE2) or (name OT1 OT2 OXT)")
    dH_H, dH_D = build_protein_donors(u, prot)

    t0 = u.trajectory[0].time
    dt = (u.trajectory[1].time - t0) if len(u.trajectory) > 1 else CONTACT_PS
    stride = max(1, int(round(CONTACT_PS / dt)))
    half = int(round(WIN_NS * 1000.0 / dt))

    events = []
    with open(EVENTS) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.split()
            events.append((int(c[0]), float(c[1]), c[2], int(c[3])))
    print(f"{len(events)} events")

    out = os.path.join(OUTDIR, "water_pathway.csv")
    fcsv = open(out, "w", newline="")
    w = csv.writer(fcsv)
    w.writerow(["water", "order", "time_ns", "z_rel_A", "type",
                "segid", "resid", "resname", "atom", "func_group", "moiety",
                "water_atom", "dist_A", "HA_dist_A", "angle_deg", "direction", "pore"])

    for ev, (resid, t_cross, direction, pore) in enumerate(events):
        wat = u.select_atoms(f"({WATER_SEL}) and resid {resid}")
        if wat.n_atoms == 0:
            continue
        wO = wat.select_atoms("name OH2 OW O")[:1]
        wH = wat.select_atoms("name H*")
        c_frame = int(round((t_cross * 1000.0 - t0) / dt))
        lo, hi = max(0, c_frame - half), min(len(u.trajectory) - 1, c_frame + half)

        # ---- collect the transit timeline (time-ordered) ----
        timeline = []   # (time_ns, z_rel, {key:(subtype,patom,watom,dist,ang)})
        for fr in range(lo, hi + 1, stride):
            ts = u.trajectory[fr]
            box = ts.dimensions; Lz = box[2]
            cz = whole_mean(prot.positions[:, 2], Lz)
            zc = wO.positions[0, 2]
            dz = zc - cz - Lz * round((zc - cz) / Lz)
            if abs(dz) > Z_TRANSIT:
                continue
            fi = frame_interactions(wat, wO, wH, box, acceptors, dH_H, dH_D, u,
                                    cations, anions, prot_polar)
            timeline.append((ts.time / 1000.0, dz, fi))

        # ---- which residues ever form a strong bond (to suppress weak polar) ----
        strong_res = set()
        for _, _, fi in timeline:
            for (seg, rid, rn), (typ, *_ ) in fi.items():
                if typ in STRONG:
                    strong_res.add((seg, rid, rn))

        # ---- build ordered pathway: first occurrence of each (residue, subtype) ----
        seen = set()
        order = 0
        for t, dz, fi in timeline:
            # within a frame, order new interactions by distance (closest first)
            newly = []
            for reskey, (typ, patom, watom, dist, ha, ang) in fi.items():
                if SUPPRESS_POLAR_IF_BONDED and typ == "polar_contact" and reskey in strong_res:
                    continue
                key = (reskey, typ)
                if key in seen:
                    continue
                newly.append((dist, reskey, typ, patom, watom, ha, ang))
            for dist, reskey, typ, patom, watom, ha, ang in sorted(newly, key=lambda x: x[0]):
                seen.add((reskey, typ))
                order += 1
                moiety = "BB" if patom.name in BACKBONE else "SC"
                seg, rid, rn = reskey
                fg = func_group(rn, patom.name)
                w.writerow([f"TIP3_{resid}", order, f"{t:.3f}", f"{dz:.2f}", typ,
                            seg, rid, rn, patom.name, fg, moiety, watom,
                            f"{dist:.2f}", (f"{ha:.2f}" if ha is not None else ""),
                            (f"{ang:.1f}" if ang is not None else ""),
                            direction, pore])

        if (ev + 1) % 25 == 0:
            print(f"  {ev+1}/{len(events)} events done")

    fcsv.close()
    print(f"\nWritten: {out}")


if __name__ == "__main__":
    main()

###################################################################################################################################
############################ contact analysis for few of the water molecules passing through the channel ##########################
###################################################################################################################################

################################## average contact taken by the water for permeation ####################
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CONSENSUS (average) permeation path from water_pathway.csv.

For each direction (down / up) it answers: as a water crosses the channel, which
residues does it MOST OFTEN interact with, in what order along z, and by which
interaction? The path is not a single water's route but the statistical average
over all permeating waters.

Method
  * bin the pore axis (z_rel_A) into slabs of BIN_WIDTH A;
  * per direction, per slab, occupancy(residue) =
        (# distinct waters of that direction that interact with the residue in the slab)
        / (# distinct waters of that direction present in the slab);
  * the average path = the top residue(s) per slab, read entry -> exit
        (down: z high -> low ; up: z low -> high);
  * dominant interaction type + mean distance/angle reported per residue/slab;
  * up-vs-down compared slab by slab.

resid is per-protomer numbering, so aggregating by resid pools the 3 equivalent pores.

Requires: pandas, numpy, matplotlib.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG ============================
CSV       = "/media/supremeleader/Pantera/simulation/lipid_simulation/5WTR/charmm-gui-7977600322/gromacs/5wtr/water_permeation/contacts/water_pathway.csv"   # <- your file
OUTDIR    = os.path.join(os.path.dirname(CSV), "consensus_path")

BIN_WIDTH = 2.0     # A : z-slab thickness
Z_MIN, Z_MAX = None, None   # None -> derived from the data (min/max z_rel_A of the contacts).
                            # Set numbers only if you want to force a fixed window.
Z_PAD = 0.0         # A : optional extra padding beyond the observed min/max
TOP_N     = 3       # residues listed per slab in the consensus table
MIN_OCC   = 0.15    # residue must reach this occupancy in a slab to enter the path string
HEAT_OCC  = 0.20    # residue enters the heatmap if it reaches this occupancy in any slab
DROP_POLAR = False  # True -> ignore polar_contact, keep only hbond/elec
# ================================================================


def resid_label(row):
    return f"{row['resname']}{int(row['resid'])}"


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    df = pd.read_csv(CSV)
    if DROP_POLAR:
        df = df[df.type != "polar_contact"].copy()
    df["res"] = df.resname.astype(str) + df.resid.astype(int).astype(str)

    # --- z-range: data-driven unless explicitly overridden (portable across proteins) ---
    zlo = Z_MIN if Z_MIN is not None else df.z_rel_A.min() - Z_PAD
    zhi = Z_MAX if Z_MAX is not None else df.z_rel_A.max() + Z_PAD
    # snap to whole BIN_WIDTH multiples so slabs are clean and symmetric-friendly
    zlo = np.floor(zlo / BIN_WIDTH) * BIN_WIDTH
    zhi = np.ceil(zhi / BIN_WIDTH) * BIN_WIDTH
    print(f"z-range used: {zlo:.1f} to {zhi:.1f} A  (from data: "
          f"{df.z_rel_A.min():.1f} .. {df.z_rel_A.max():.1f})")

    edges = np.arange(zlo, zhi + BIN_WIDTH, BIN_WIDTH)
    centers = (edges[:-1] + edges[1:]) / 2.0
    df["zbin"] = pd.cut(df.z_rel_A, edges, labels=centers, include_lowest=True)
    df = df.dropna(subset=["zbin"])
    df["zbin"] = df["zbin"].astype(float)

    consensus_rows = []
    path_strings = {}
    occ_tables = {}   # direction -> DataFrame (index res, columns zbin, values occupancy)

    for d in ("down", "up"):
        sub = df[df.direction == d]
        if sub.empty:
            continue
        present = sub.groupby("zbin")["water"].nunique()            # denominator per slab
        # distinct waters per (slab, residue)
        n_res = (sub.groupby(["zbin", "res"])["water"].nunique()
                    .rename("n").reset_index())
        n_res["occ"] = n_res.apply(lambda r: r["n"] / present[r["zbin"]], axis=1)
        # dominant type + mean geometry per (slab, residue) -- explicit, version-robust
        grp = sub.groupby(["zbin", "res"], observed=True)
        mode = lambda s: s.value_counts().idxmax()
        agg = pd.DataFrame({
            "dom_type":  grp["type"].agg(mode),
            "mean_dist": grp["dist_A"].mean(),
            "mean_ang":  grp["angle_deg"].mean(),
            "func":      grp["func_group"].agg(mode),
            "moiety":    grp["moiety"].agg(mode),
        }).reset_index()
        m = n_res.merge(agg, on=["zbin", "res"])

        # occupancy matrix for the heatmap
        occ_tables[d] = m.pivot_table(index="res", columns="zbin", values="occ", fill_value=0.0)

        # entry->exit slab order
        slab_order = sorted(centers, reverse=(d == "down"))

        # consensus table (top N per slab) + path string (top-1 above MIN_OCC)
        path_bits = []
        for z in slab_order:
            g = m[m.zbin == z].sort_values("occ", ascending=False)
            if g.empty:
                continue
            for rank, (_, r) in enumerate(g.head(TOP_N).iterrows(), 1):
                consensus_rows.append({
                    "direction": d, "z_bin": z, "rank": rank,
                    "residue": r["res"], "func_group": r["func"], "moiety": r["moiety"],
                    "occupancy": round(r["occ"], 3), "n_waters": int(r["n"]),
                    "n_present": int(present[z]),
                    "dominant_type": r["dom_type"],
                    "mean_dist_A": round(r["mean_dist"], 2),
                    "mean_angle": (round(r["mean_ang"], 1) if pd.notna(r["mean_ang"]) else ""),
                })
            top = g.iloc[0]
            if top["occ"] >= MIN_OCC:
                path_bits.append(f"z={z:+.0f} {top['res']}[{top['dom_type']},"
                                 f"{top['occ']*100:.0f}%]")
        path_strings[d] = "  ->  ".join(path_bits)

    cons = pd.DataFrame(consensus_rows)
    cons.to_csv(os.path.join(OUTDIR, "consensus_pathway.csv"), index=False)

    # ---- up vs down comparison: top residue per slab ----
    comp_rows = []
    for z in sorted(centers, reverse=True):
        row = {"z_bin": z}
        for d in ("down", "up"):
            g = cons[(cons.direction == d) & (cons.z_bin == z) & (cons["rank"] == 1)]
            if len(g):
                row[f"{d}_res"] = g.iloc[0]["residue"]
                row[f"{d}_type"] = g.iloc[0]["dominant_type"]
                row[f"{d}_occ"] = g.iloc[0]["occupancy"]
            else:
                row[f"{d}_res"] = ""; row[f"{d}_type"] = ""; row[f"{d}_occ"] = np.nan
        row["same_residue"] = (row.get("down_res") == row.get("up_res")
                               and row.get("down_res") != "")
        comp_rows.append(row)
    comp = pd.DataFrame(comp_rows)
    comp.to_csv(os.path.join(OUTDIR, "up_vs_down_comparison.csv"), index=False)

    # readable path strings
    with open(os.path.join(OUTDIR, "consensus_path.txt"), "w") as f:
        for d in ("down", "up"):
            f.write(f"# {d.upper()} permeation - average path (entry -> exit)\n")
            f.write(path_strings.get(d, "(none)") + "\n\n")
        f.write("# top residue per z-slab, down vs up (same_residue flag)\n")
        f.write(comp.to_string(index=False) + "\n")

    # ---- heatmap: occupancy(residue, z) for down, up, and difference ----
    if occ_tables:
        allz = list(centers)
        res_keep = set()
        for d, tbl in occ_tables.items():
            res_keep |= set(tbl.index[(tbl.max(axis=1) >= HEAT_OCC)])
        res_keep = sorted(res_keep)
        # order residues by occupancy-weighted mean z (across both directions)
        def mean_z(res):
            num = den = 0.0
            for tbl in occ_tables.values():
                if res in tbl.index:
                    for z in tbl.columns:
                        num += z * tbl.loc[res, z]; den += tbl.loc[res, z]
            return num / den if den else 0.0
        res_keep.sort(key=mean_z)

        def mat(d):
            tbl = occ_tables.get(d)
            M = np.zeros((len(res_keep), len(allz)))
            if tbl is not None:
                for i, r in enumerate(res_keep):
                    for j, z in enumerate(allz):
                        if r in tbl.index and z in tbl.columns:
                            M[i, j] = tbl.loc[r, z]
            return M
        Md, Mu = mat("down"), mat("up")

        fig, axs = plt.subplots(1, 3, figsize=(15, max(4, 0.32 * len(res_keep))),
                                sharey=True)
        for ax, M, ttl, cm, vlim in [
            (axs[0], Md, "DOWN occupancy", "viridis", (0, 1)),
            (axs[1], Mu, "UP occupancy", "viridis", (0, 1)),
            (axs[2], Md - Mu, "DOWN - UP", "coolwarm", (-0.6, 0.6))]:
            im = ax.imshow(M, aspect="auto", cmap=cm, vmin=vlim[0], vmax=vlim[1],
                           extent=[allz[0], allz[-1], len(res_keep) - 0.5, -0.5])
            ax.set_title(ttl); ax.set_xlabel("z relative to channel centre (\u00c5)")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        axs[0].set_yticks(range(len(res_keep)))
        axs[0].set_yticklabels(res_keep, fontsize=8)
        fig.suptitle("Consensus water-permeation path (residue occupancy along z)")
        fig.tight_layout()
        fig.savefig(os.path.join(OUTDIR, "consensus_path.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)

    print("DOWN path:\n ", path_strings.get("down", "(none)"))
    print("\nUP path:\n ", path_strings.get("up", "(none)"))
    print(f"\nOutputs in: {OUTDIR}")


if __name__ == "__main__":
    main()
################################## average contact taken by the water for permeation ####################


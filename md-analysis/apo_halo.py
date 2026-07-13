#######################################################################
#################################### RMSD #############################
#######################################################################
# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Trimer-averaged C-alpha RMSD (Angstrom): apo (6IYX) vs holo (6IZF).
#
# Quantity ("trimer average RMSD"): each chain fitted INDEPENDENTLY (rot+trans
# least-squares) onto the production-start reference, RMSD per chain, then /3:
#         RMSD_trimer_avg = (R_chainA + R_chainB + R_chainC) / 3
# This is the divide-by-3 analogue of running `gmx rms ... -fit rot+trans` per chain.
#
# Outputs (all in OUTDIR), RMSD in Angstrom, sampled every 0.5 ns, x-ticks every 100 ns:
#   COMPARISON (apo vs holo), raw trace + running average overlaid (xmgrace style):
#     - trimer_avg_rmsd_TM_apo_vs_holo.(png/pdf)        TM region (H1-H7)
#     - trimer_avg_rmsd_protein_apo_vs_holo.(png/pdf)   full protein
#   PER-STRUCTURE DECOMPOSITION (3 chains + trimer average), full protein:
#     - chain_rmsd_decomposition_6IYX.(png/pdf)
#     - chain_rmsd_decomposition_6IZF.(png/pdf)
#   plus .dat tables (time_ns, trimer_avg, chainA, chainB, chainC).
#
# GROMACS equivalence: same rot+trans least-squares (QCP == Kabsch); gmx mass-weights
# by default but that is IDENTICAL for C-alpha (equal masses cancel). gmx reports nm,
# this reports Angstrom (x10).
#
# Requires: MDAnalysis, numpy, matplotlib
# """
#
# import os
# import numpy as np
# import MDAnalysis as mda
# from MDAnalysis.analysis.rms import rmsd as mda_rmsd
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt
# from matplotlib.ticker import MultipleLocator
#
# # ============================ CONFIG ============================
# OUTDIR = "/media/supremeleader/Pantera/simulation/apo_halo"
#
# SYSTEMS = {
#     "6IYX \u2013 apo (no DAG)": {
#         "tag":  "6IYX",
#         "top":  "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/center.gro",
#         "traj": "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/center.xtc",
#         "ref":  None,          # None -> frame 0 of 'top'; else path to production-start .gro
#         "color": "#0072B2",    # blue
#     },
#     "6IZF \u2013 holo (DAG)": {
#         "tag":  "6IZF",
#         "top":  "/media/supremeleader/Pantera/simulation/lipid_simulation/6IZF/charmm-gui-7977903380/gromacs/6izf/center.gro",
#         "traj": "/media/supremeleader/Pantera/simulation/lipid_simulation/6IZF/charmm-gui-7977903380/gromacs/6izf/center.xtc",
#         "ref":  None,
#         "color": "#D55E00",    # vermillion
#     },
# }
#
# helices = {
#     "H1": (23,  39),  "H2": (50,  73),
#     "H3": (84,  98),  "H4": (112, 138), "H5": (144, 170),
#     "H6": (185, 202), "H7": (209, 229),
# }
#
# N_CHAINS            = 3
# DT_TARGET_PS        = 250.0   # sample every 0.25 ns
# TICK_NS             = 100.0   # major x-tick every 100 ns
# MINOR_TICK_NS       = 20.0
# START_AT_ZERO       = True
# USE_MASS_WEIGHTING  = False    # False == gmx for C-alpha-only
#
# # ---- xmgrace-style running average ----
# RUNNING_AVG_NS      = 5.0      # window length of the moving average, in ns
# RUNNING_AVG_CENTERED = False   # False = trailing (xmgrace default behaviour)
#
# # raw per-frame trace is the prominent layer; running average is a thin overlay on top
# RAW_LW     = 0.7
# RAW_ALPHA  = 0.60
# AVG_LW     = 1.2
#
# CHAIN_COLORS = ["#377EB8", "#4DAF4A", "#984EA3"]   # chA, chB, chC
# AVG_COLOR    = "#000000"
#
# AA = ("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL "
#       "HSD HSE HSP HID HIE HIP HISD HISE HISP CYX CYM ASH GLH LYN ARN").split()
# # ================================================================
#
# helix_resids = sorted({r for (lo, hi) in helices.values() for r in range(lo, hi + 1)})
# EXPECTED_TM_PER_CHAIN = len(helix_resids)
#
#
# def chain_blocks(u):
#     aa_sel = "resname " + " ".join(AA)
#     res = u.select_atoms(f"name CA and ({aa_sel})").residues
#     n = len(res)
#     if n == 0:
#         raise ValueError("No protein C-alpha found - check residue names.")
#     if n % N_CHAINS != 0:
#         raise ValueError(f"{n} protein C-alpha not divisible by {N_CHAINS} chains.")
#     per = n // N_CHAINS
#     return [res[k * per:(k + 1) * per] for k in range(N_CHAINS)], res[0].resid
#
#
# def chain_ca_indices(u, which):
#     blocks, base_first = chain_blocks(u)
#     out = []
#     for block in blocks:
#         offset = block[0].resid - base_first
#         sel = [r for r in block if (r.resid - offset) in helix_resids] if which == "tm" else list(block)
#         out.append(np.array([r.atoms.select_atoms("name CA").indices[0] for r in sel], dtype=int))
#     return out
#
#
# def running_average(y, win_pts, centered=False):
#     """Simple moving average of length win_pts (xmgrace 'Running averages').
#     Returns (idx, smoothed) where idx indexes the time array."""
#     win_pts = max(1, int(win_pts))
#     if win_pts == 1 or len(y) < win_pts:
#         return np.arange(len(y)), np.asarray(y, float)
#     kernel = np.ones(win_pts) / win_pts
#     sm = np.convolve(y, kernel, mode="valid")          # len = N - win + 1
#     if centered:
#         start = (win_pts - 1) // 2
#     else:
#         start = win_pts - 1                            # trailing: align to last point of window
#     idx = np.arange(start, start + len(sm))
#     return idx, sm
#
#
# def analyse(label, cfg, which):
#     u    = mda.Universe(cfg["top"], cfg["traj"])
#     refU = mda.Universe(cfg["ref"]) if cfg["ref"] else mda.Universe(cfg["top"])
#     chain_idx = chain_ca_indices(u, which)
#     counts = [len(c) for c in chain_idx]
#     tag = "TM" if which == "tm" else "protein"
#     print(f"  [{cfg['tag']}] {tag}: C-alpha/chain = {counts}", end="")
#     if which == "tm":
#         print(f" (expected {EXPECTED_TM_PER_CHAIN})", end="")
#         if any(c != EXPECTED_TM_PER_CHAIN for c in counts):
#             print("  !! WARNING numbering/range mismatch", end="")
#     print()
#
#     ag_chains  = [u.atoms[c] for c in chain_idx]
#     ref_chains = [refU.atoms[c].positions.copy() for c in chain_idx]
#     w_chains   = [a.masses if USE_MASS_WEIGHTING else None for a in ag_chains]
#
#     n_frames = len(u.trajectory)
#     dt = (u.trajectory[1].time - u.trajectory[0].time) if n_frames >= 2 else DT_TARGET_PS
#     step = max(1, int(round(DT_TARGET_PS / dt))) if dt > 0 else 1
#
#     times, per_chain = [], [[] for _ in range(N_CHAINS)]
#     for ts in u.trajectory[::step]:
#         times.append(ts.time)
#         for k in range(N_CHAINS):
#             per_chain[k].append(mda_rmsd(ag_chains[k].positions, ref_chains[k],
#                                          weights=w_chains[k], center=True, superposition=True))
#     times = np.asarray(times)
#     t_ns = (times - times[0]) / 1000.0 if START_AT_ZERO else times / 1000.0
#     per_chain = np.asarray(per_chain)
#     avg = per_chain.mean(axis=0)        # divide-by-3
#
#     os.makedirs(OUTDIR, exist_ok=True)
#     np.savetxt(os.path.join(OUTDIR, f"{tag}_rmsd_{cfg['tag']}.dat"),
#                np.column_stack([t_ns, avg, per_chain[0], per_chain[1], per_chain[2]]),
#                fmt="%.4f",
#                header=f"time_ns  trimer_avg_A  chainA_A  chainB_A  chainC_A  # {tag} C-alpha RMSD")
#
#     dt_ns = (t_ns[1] - t_ns[0]) if len(t_ns) > 1 else DT_TARGET_PS / 1000.0
#     return {"label": label, "tag": cfg["tag"], "color": cfg["color"],
#             "t": t_ns, "avg": avg, "per_chain": per_chain, "dt_ns": dt_ns}
#
#
# def _style():
#     plt.rcParams.update({
#         "font.family": "sans-serif", "font.size": 13,
#         "axes.linewidth": 1.2, "xtick.direction": "in", "ytick.direction": "in",
#         "xtick.major.size": 6, "ytick.major.size": 6,
#         "xtick.minor.size": 3, "ytick.minor.size": 3,
#         "xtick.top": True, "ytick.right": True, "legend.frameon": False,
#     })
#
#
# def _finish_axes(ax, tmax, ymax):
#     ax.set_xlabel("Time (ns)")
#     ax.set_ylabel("RMSD (\u00c5)")
#     ax.set_xlim(0, tmax)
#     ax.set_ylim(0, ymax * 1.05)
#     ax.margins(0)
#     ax.xaxis.set_major_locator(MultipleLocator(TICK_NS))
#     ax.xaxis.set_minor_locator(MultipleLocator(MINOR_TICK_NS))
#
#
# def make_compare_figure(results, title, outname):
#     _style()
#     fig, ax = plt.subplots(figsize=(7.2, 4.3))
#     tmax = max(r["t"][-1] for r in results)
#     ymax = max(r["avg"].max() for r in results)
#     win = max(1, round(RUNNING_AVG_NS / results[0]["dt_ns"]))
#     for r in results:
#         ax.plot(r["t"], r["avg"], lw=RAW_LW, color=r["color"], alpha=RAW_ALPHA,
#                 zorder=1)                                                  # raw per-frame RMSD (prominent)
#         idx, sm = running_average(r["avg"], win, RUNNING_AVG_CENTERED)     # running average
#         ax.plot(r["t"][idx], sm, lw=AVG_LW, color=r["color"], alpha=1.0,
#                 zorder=3, solid_capstyle="round", label=r["label"])        # thin overlay
#     _finish_axes(ax, tmax, ymax)
#     ax.set_title(title)
#     ax.legend(loc="lower right")
#     fig.tight_layout(pad=0.4)
#     for ext in ("png", "pdf"):
#         fig.savefig(os.path.join(OUTDIR, f"{outname}.{ext}"), dpi=300, bbox_inches="tight")
#     plt.close(fig)
#     print(f"  saved {outname}  (running-avg window = {win} pts = {win*results[0]['dt_ns']:.1f} ns)")
#
#
# def make_decomposition_figure(r, title, outname):
#     _style()
#     fig, ax = plt.subplots(figsize=(7.2, 4.3))
#     tmax = r["t"][-1]
#     ymax = max(r["per_chain"].max(), r["avg"].max())
#     for k in range(N_CHAINS):
#         ax.plot(r["t"], r["per_chain"][k], lw=1.0, color=CHAIN_COLORS[k],
#                 alpha=0.85, label=f"Chain {chr(65+k)}")
#     ax.plot(r["t"], r["avg"], lw=2.0, color=AVG_COLOR, label="Trimer average")
#     _finish_axes(ax, tmax, ymax)
#     ax.set_title(title)
#     ax.legend(loc="lower right", ncol=2)
#     fig.tight_layout(pad=0.4)
#     for ext in ("png", "pdf"):
#         fig.savefig(os.path.join(OUTDIR, f"{outname}.{ext}"), dpi=300, bbox_inches="tight")
#     plt.close(fig)
#     print(f"  saved {outname}")
#
#
# def main():
#     os.makedirs(OUTDIR, exist_ok=True)
#
#     print("== TM region (comparison) ==")
#     tm = [analyse(lbl, cfg, "tm") for lbl, cfg in SYSTEMS.items()]
#     make_compare_figure(tm,
#         "Trimer-Averaged TM-Region C$_\\alpha$ RMSD: apo (6IYX) vs holo (6IZF)",
#         "trimer_avg_rmsd_TM_apo_vs_holo")
#
#     print("== full protein (comparison) ==")
#     fp = [analyse(lbl, cfg, "all") for lbl, cfg in SYSTEMS.items()]
#     make_compare_figure(fp,
#         "Trimer-Averaged Whole-Protein C$_\\alpha$ RMSD: apo (6IYX) vs holo (6IZF)",
#         "trimer_avg_rmsd_protein_apo_vs_holo")
#
#     print("== per-structure chain decomposition: TM region ==")
#     for r in tm:
#         make_decomposition_figure(r,
#             f"Per-Chain & Trimer-Average C$_\\alpha$ RMSD \u2013 {r['tag']} (TM region)",
#             f"chain_rmsd_decomposition_TM_{r['tag']}")
#
#     print("== per-structure chain decomposition: full protein ==")
#     for r in fp:
#         make_decomposition_figure(r,
#             f"Per-Chain & Trimer-Average C$_\\alpha$ RMSD \u2013 {r['tag']} (full protein)",
#             f"chain_rmsd_decomposition_protein_{r['tag']}")
#
#     print(f"\nAll outputs in: {OUTDIR}")
#
#
# if __name__ == "__main__":
#     main()

#######################################################################
#################################### RMSF #############################
#######################################################################
# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Per-residue C-alpha RMSF (Angstrom): apo (6IYX) vs holo (6IZF), with TM helices shaded.
#
# Quantity: each chain is least-squares fitted (rot+trans) to the production-start
# reference, then RMSF (std. dev. of position around the trajectory average) is computed
# per residue. Trimer average = (chainA + chainB + chainC)/3 per residue.
# This reproduces:  printf "<chain>" | gmx rmsf -f center.xtc -s ref.gro -res   (default -fit).
#
# Figures (OUTDIR), RMSF in Angstrom, sampled every 0.25 ns, TM helices shaded:
#   1) rmsf_protein_apo_vs_holo      : trimer-average RMSF, apo vs holo (no smoothing)
#   2) rmsf_chains_6IYX              : chain A/B/C + trimer average (apo)
#   3) rmsf_chains_6IZF              : chain A/B/C + trimer average (holo)
# plus rmsf_<TAG>.dat tables (resid, trimer_avg, chainA, chainB, chainC).
#
# GROMACS note: gmx -fit uses the selected (RMSF) group for the superposition and
# mass-weights it; for C-alpha that equals geometric weighting. gmx reports nm (x10).
#
# Requires: MDAnalysis, numpy, matplotlib
# """
#
# import os
# import numpy as np
# import MDAnalysis as mda
# from MDAnalysis.analysis.align import rotation_matrix
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt
# from matplotlib.ticker import MultipleLocator
#
# # ============================ CONFIG ============================
# OUTDIR = "/media/supremeleader/Pantera/simulation/apo_halo"
#
# SYSTEMS = {
#     "6IYX \u2013 apo (no DAG)": {
#         "tag":  "6IYX",
#         "top":  "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/center.gro",
#         "traj": "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/center.xtc",
#         "ref":  None,          # None -> frame 0 of 'top'; else production-start .gro path
#         "color": "#0072B2",
#     },
#     "6IZF \u2013 holo (DAG)": {
#         "tag":  "6IZF",
#         "top":  "/media/supremeleader/Pantera/simulation/lipid_simulation/6IZF/charmm-gui-7977903380/gromacs/6izf/center.gro",
#         "traj": "/media/supremeleader/Pantera/simulation/lipid_simulation/6IZF/charmm-gui-7977903380/gromacs/6izf/center.xtc",
#         "ref":  None,
#         "color": "#D55E00",
#     },
# }
#
# # TM helix residue ranges (protomer numbering) -> shaded as transparent bands.
# helices = {
#     "H1": (23,  39),  "H2": (50,  73),
#     "H3": (84,  98),  "H4": (112, 138), "H5": (144, 170),
#     "H6": (185, 202), "H7": (209, 229),
# }
#
# N_CHAINS           = 3
# DT_TARGET_PS       = 250.0    # sample every 0.25 ns
# USE_MASS_WEIGHTING = False    # False == gmx for C-alpha-only
# RES_TICK           = 20       # major x-tick spacing (residues)
# RES_MINOR          = 10
# SHOW_HELIX_LABELS  = True
# TM_BAND_COLOR      = "0.55"   # grey
# TM_BAND_ALPHA      = 0.18
#
# CHAIN_COLORS = ["#377EB8", "#4DAF4A", "#984EA3"]   # chA, chB, chC
# AVG_COLOR    = "#000000"
#
# AA = ("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL "
#       "HSD HSE HSP HID HIE HIP HISD HISE HISP CYX CYM ASH GLH LYN ARN").split()
# # ================================================================
#
#
# def chain_ca_indices(u):
#     """Split protein C-alpha into N_CHAINS equal consecutive residue blocks (homotrimer);
#     return list of global atom-index arrays (one Cα per residue, residue order)."""
#     aa_sel = "resname " + " ".join(AA)
#     res = u.select_atoms(f"name CA and ({aa_sel})").residues
#     n = len(res)
#     if n == 0 or n % N_CHAINS != 0:
#         raise ValueError(f"{n} protein C-alpha not divisible by {N_CHAINS} chains.")
#     per = n // N_CHAINS
#     return [np.array([r.atoms.select_atoms("name CA").indices[0]
#                       for r in res[k*per:(k+1)*per]], dtype=int) for k in range(N_CHAINS)]
#
#
# def compute_rmsf(label, cfg):
#     u    = mda.Universe(cfg["top"], cfg["traj"])
#     refU = mda.Universe(cfg["ref"]) if cfg["ref"] else mda.Universe(cfg["top"])
#     chain_idx = chain_ca_indices(u)
#     ag      = [u.atoms[c] for c in chain_idx]
#     # reference (centred) per chain
#     ref_c   = []
#     for c in chain_idx:
#         x = refU.atoms[c].positions.astype(float)
#         ref_c.append(x - x.mean(axis=0))
#     w = [ag[k].masses if USE_MASS_WEIGHTING else None for k in range(N_CHAINS)]
#
#     natoms = [len(c) for c in chain_idx]
#     S1  = [np.zeros((n, 3)) for n in natoms]   # sum of aligned coords
#     S2  = [np.zeros(n)      for n in natoms]   # sum of |aligned|^2
#     nfr = 0
#
#     n_frames = len(u.trajectory)
#     dt = (u.trajectory[1].time - u.trajectory[0].time) if n_frames >= 2 else DT_TARGET_PS
#     step = max(1, int(round(DT_TARGET_PS / dt))) if dt > 0 else 1
#
#     for ts in u.trajectory[::step]:
#         for k in range(N_CHAINS):
#             x = ag[k].positions.astype(float)
#             xc = x - x.mean(axis=0)                       # remove translation
#             R, _ = rotation_matrix(xc, ref_c[k], weights=w[k])
#             xa = xc @ R.T                                 # remove rotation
#             S1[k] += xa
#             S2[k] += (xa * xa).sum(axis=1)
#         nfr += 1
#
#     rmsf = []
#     for k in range(N_CHAINS):
#         mean_pos = S1[k] / nfr
#         var = S2[k] / nfr - (mean_pos * mean_pos).sum(axis=1)   # <|x|^2> - |<x>|^2
#         rmsf.append(np.sqrt(np.clip(var, 0, None)))            # Angstrom
#     rmsf = np.array(rmsf)                                       # (3, nres)
#     avg = rmsf.mean(axis=0)                                     # divide-by-3
#
#     resids = u.atoms[chain_idx[0]].resids                       # protomer numbering (x-axis)
#     print(f"  [{cfg['tag']}] {nfr} frames (stride {step}, ~{step*dt/1000:.3f} ns), "
#           f"{rmsf.shape[1]} residues, <RMSF>={avg.mean():.3f} A")
#
#     os.makedirs(OUTDIR, exist_ok=True)
#     np.savetxt(os.path.join(OUTDIR, f"rmsf_{cfg['tag']}.dat"),
#                np.column_stack([resids, avg, rmsf[0], rmsf[1], rmsf[2]]),
#                fmt=["%d", "%.4f", "%.4f", "%.4f", "%.4f"],
#                header="resid  trimer_avg_A  chainA_A  chainB_A  chainC_A  # per-residue C-alpha RMSF")
#     return {"label": label, "tag": cfg["tag"], "color": cfg["color"],
#             "resids": resids, "rmsf": rmsf, "avg": avg}
#
#
# def _style():
#     plt.rcParams.update({
#         "font.family": "sans-serif", "font.size": 13,
#         "axes.linewidth": 1.2, "xtick.direction": "in", "ytick.direction": "in",
#         "xtick.major.size": 6, "ytick.major.size": 6,
#         "xtick.minor.size": 3, "ytick.minor.size": 3,
#         "xtick.top": True, "ytick.right": True, "legend.frameon": False,
#     })
#
#
# def _add_tm_bands(ax, ymax):
#     for name, (lo, hi) in helices.items():
#         ax.axvspan(lo, hi, color=TM_BAND_COLOR, alpha=TM_BAND_ALPHA, lw=0, zorder=0)
#         if SHOW_HELIX_LABELS:
#             ax.text((lo + hi) / 2.0, ymax * 1.005, name, ha="center", va="bottom",
#                     fontsize=9, color="0.35")
#
#
# def _finish(ax, resids, ymax):
#     ax.set_xlabel("Residue number")
#     ax.set_ylabel("RMSF (\u00c5)")
#     ax.set_xlim(resids.min(), resids.max())
#     ax.set_ylim(0, ymax * 1.05)
#     ax.margins(0)
#     ax.xaxis.set_major_locator(MultipleLocator(RES_TICK))
#     ax.xaxis.set_minor_locator(MultipleLocator(RES_MINOR))
#
#
# def _save(fig, outname):
#     fig.tight_layout(pad=0.4)
#     for ext in ("png", "pdf"):
#         fig.savefig(os.path.join(OUTDIR, f"{outname}.{ext}"), dpi=300, bbox_inches="tight")
#     plt.close(fig)
#     print(f"  saved {outname}")
#
#
# def make_compare_figure(results, outname):
#     _style()
#     fig, ax = plt.subplots(figsize=(8.0, 4.3))
#     ymax = max(r["avg"].max() for r in results)
#     _add_tm_bands(ax, ymax)
#     for r in results:
#         ax.plot(r["resids"], r["avg"], lw=1.4, color=r["color"], label=r["label"], zorder=3)
#     _finish(ax, results[0]["resids"], ymax)
#     ax.set_title("Per-Residue C$_\\alpha$ RMSF (trimer average): apo (6IYX) vs holo (6IZF)")
#     ax.legend(loc="upper right")
#     _save(fig, outname)
#
#
# def make_decomposition_figure(r, outname):
#     _style()
#     fig, ax = plt.subplots(figsize=(8.0, 4.3))
#     ymax = max(r["rmsf"].max(), r["avg"].max())
#     _add_tm_bands(ax, ymax)
#     for k in range(N_CHAINS):
#         ax.plot(r["resids"], r["rmsf"][k], lw=0.9, color=CHAIN_COLORS[k],
#                 alpha=0.85, label=f"Chain {chr(65+k)}", zorder=2)
#     ax.plot(r["resids"], r["avg"], lw=1.8, color=AVG_COLOR, label="Trimer average", zorder=3)
#     _finish(ax, r["resids"], ymax)
#     ax.set_title(f"Per-Residue C$_\\alpha$ RMSF \u2013 {r['tag']} (full protein)")
#     ax.legend(loc="upper right", ncol=2)
#     _save(fig, outname)
#
#
# def main():
#     os.makedirs(OUTDIR, exist_ok=True)
#     print("== computing RMSF ==")
#     res = [compute_rmsf(lbl, cfg) for lbl, cfg in SYSTEMS.items()]
#     print("== figures ==")
#     make_compare_figure(res, "rmsf_protein_apo_vs_holo")
#     for r in res:
#         make_decomposition_figure(r, f"rmsf_chains_{r['tag']}")
#     print(f"\nAll outputs in: {OUTDIR}")
#
#
# if __name__ == "__main__":
#     main()

#######################################################################
#################################### DSSP TM2 #########################
#######################################################################
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DSSP secondary-structure CARPET (residue x time), reproducing gmx do_dssp + xpm2ps.

Uses mdtraj's full 8-state DSSP (compute_dssp simplified=False) -> same SS categories
as mkdssp/xpm2ps (A-Helix H, 3-Helix G, 5-Helix I, B-Sheet E, B-Bridge B, Turn T,
Bend S, Coil), drawn as a colored residue-vs-time matrix.

For each system and each helix (TM2, TM5) it draws the 3 chains as stacked carpets.

Equivalent to:
  export DSSP=/usr/bin/mkdssp
  gmx do_dssp -f center.xtc -s md.gro -n idx.ndx -o helix.xpm -dt 1000 ; gmx xpm2ps ...

Requires: mdtraj, numpy, matplotlib     (pip install mdtraj)
"""

import os
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator

# ============================ CONFIG ============================
OUTDIR = "/media/supremeleader/Pantera/simulation/apo_halo"

SYSTEMS = {
    "6IYX": {
        "top":  "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/center.gro",
        "traj": "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/center.xtc",
    },
    "6IZF": {
        "top":  "/media/supremeleader/Pantera/simulation/lipid_simulation/6IZF/charmm-gui-7977903380/gromacs/6izf/center.gro",
        "traj": "/media/supremeleader/Pantera/simulation/lipid_simulation/6IZF/charmm-gui-7977903380/gromacs/6izf/center.xtc",
    },
}

TM_HELICES = {"TM2": (50, 73), "TM5": (144, 170)}
N_CHAINS     = 3
DT_TARGET_PS = 1000.0    # 1 ns, matching your gmx -dt 1000 (set 250 for 0.25 ns)

# DSSP 8-state -> integer code -> color/label (xpm2ps-style palette)
SS_ORDER  = ["H", "G", "I", "E", "B", "T", "S", "C"]
SS_COLOR  = {"H": "#2166AC", "G": "#7F7F7F", "I": "#9E5DB0", "E": "#D62728",
             "B": "#000000", "T": "#F4D03F", "S": "#2CA02C", "C": "#FFFFFF"}
SS_LABEL  = {"H": "A-Helix", "G": "3-Helix", "I": "5-Helix", "E": "B-Sheet",
             "B": "B-Bridge", "T": "Turn", "S": "Bend", "C": "Coil"}
CODE = {s: i for i, s in enumerate(SS_ORDER)}
# ================================================================


def protein_residue_blocks(traj):
    """Residues that have a CA atom, split into N_CHAINS equal consecutive blocks."""
    top = traj.topology
    ca_res = sorted({top.atom(i).residue.index for i in top.select("name CA")})
    n = len(ca_res)
    if n == 0 or n % N_CHAINS != 0:
        raise ValueError(f"{n} protein residues not divisible by {N_CHAINS}.")
    per = n // N_CHAINS
    blocks = [ca_res[k * per:(k + 1) * per] for k in range(N_CHAINS)]
    base_first = top.residue(blocks[0][0]).resSeq
    return blocks, base_first


def char_to_code(ss_row):
    out = np.full(ss_row.shape, CODE["C"], dtype=int)
    for s, c in CODE.items():
        out[ss_row == s] = c
    out[ss_row == " "] = CODE["C"]          # mdtraj loop
    return out


def load_and_dssp(cfg):
    f0 = md.load_frame(cfg["traj"], 0, top=cfg["top"])
    f1 = md.load_frame(cfg["traj"], 1, top=cfg["top"])
    dt = float(f1.time[0] - f0.time[0])
    stride = max(1, int(round(DT_TARGET_PS / dt))) if dt > 0 else 1
    traj = md.load(cfg["traj"], top=cfg["top"], stride=stride)
    ss = md.compute_dssp(traj, simplified=False)            # (n_frames, n_residues), 8-state
    codes = char_to_code(ss)
    t_ns = traj.time / 1000.0
    blocks, base_first = protein_residue_blocks(traj)
    print(f"    frames={traj.n_frames} (stride {stride}, ~{stride*dt/1000:.3f} ns)")
    return traj, codes, t_ns, blocks, base_first


def helix_columns(traj, block, base_first, lo, hi):
    """Residue indices + protomer numbers of one chain's residues within [lo,hi]."""
    top = traj.topology
    offset = top.residue(block[0]).resSeq - base_first
    cols, nums = [], []
    for ridx in block:
        prot = top.residue(ridx).resSeq - offset
        if lo <= prot <= hi:
            cols.append(ridx); nums.append(prot)
    return np.array(cols), np.array(nums)


def make_carpet(tag, helix, lohi, codes, t_ns, blocks, base_first, traj):
    lo, hi = lohi
    cmap = ListedColormap([SS_COLOR[s] for s in SS_ORDER])
    norm = BoundaryNorm(np.arange(-0.5, len(SS_ORDER) + 0.5), cmap.N)

    fig, axes = plt.subplots(N_CHAINS, 1, figsize=(8.2, 5.4), sharex=True)
    present = set()
    for k in range(N_CHAINS):
        cols, nums = helix_columns(traj, blocks[k], base_first, lo, hi)
        mat = codes[:, cols].T                       # (n_helix_res, n_frames)
        present.update(np.unique(mat).tolist())
        ax = axes[k]
        ax.imshow(mat, aspect="auto", origin="lower", interpolation="nearest",
                  cmap=cmap, norm=norm,
                  extent=[t_ns[0], t_ns[-1], nums.min() - 0.5, nums.max() + 0.5])
        ax.set_ylabel(f"Chain {chr(65+k)}\nresidue")
        ax.yaxis.set_major_locator(MultipleLocator(5))
    axes[-1].set_xlabel("Time (ns)")
    axes[0].set_title(f"Secondary structure \u2013 {tag}, {helix} (res {lo}\u2013{hi})")

    handles = [Patch(facecolor=SS_COLOR[s], edgecolor="0.5", label=SS_LABEL[s])
               for s in SS_ORDER if CODE[s] in present]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    base = os.path.join(OUTDIR, f"dssp_carpet_{tag}_{helix}")
    for ext in ("png", "pdf"):
        fig.savefig(base + "." + ext, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"    saved dssp_carpet_{tag}_{helix}")


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    for tag, cfg in SYSTEMS.items():
        print(f"== {tag} ==")
        traj, codes, t_ns, blocks, base_first = load_and_dssp(cfg)
        for helix, lohi in TM_HELICES.items():
            make_carpet(tag, helix, lohi, codes, t_ns, blocks, base_first, traj)
    print(f"\nAll outputs in: {OUTDIR}")


if __name__ == "__main__":
    main()
#######################################################################
#################################### DSSP TM2 #########################
#######################################################################


#######################################################################
################## SPACE WITH LATERAL FENESTRALS ######################
#######################################################################

#######################################################################
################## SPACE WITH LATERAL FENESTRALS ######################
#######################################################################


#######################################################################
######### DAG OVERALL INTERACTION/ HEAD & TAIL INTERACTION  ###########
#######################################################################

#######################################################################
######### DAG OVERALL INTERACTION/ HEAD & TAIL INTERACTION  ###########
#######################################################################



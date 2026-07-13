#################################################################################
############# code for temporal radius profile using claude  ####################
#################################################################################
#!/usr/bin/env python3
"""
Temporal radius profile of a TRIC channel (trimer) using HOLE via MDAnalysis.

What this script does
---------------------
1. For each of the 3 monomer trajectories, analyses ONLY the last `LAST_NS` ns,
   computing one HOLE pore profile per ns (i.e. every `frames_per_ns` frames).
2. Excludes protein hydrogens from the pore calculation.
3. Runs the three chains (split into frame-blocks) in PARALLEL across CPU cores
   so the HOLE step is much faster than the original serial pipeline.
4. Reads results straight from MDAnalysis (`ha.results.profiles`) instead of
   moving / re-parsing hole*.out text files.
5. Produces:
       - 3 per-monomer temporal radius-profile heatmaps (Chain1/2/3)
       - 1 trimer-AVERAGE temporal radius-profile heatmap
   ...and caches every profile to .npz so you can re-plot without re-running HOLE.

Why HOLE (and not something "faster")
-------------------------------------
HOLE is the gold standard and gives proper, publication-grade radii, so it is
kept. The slowness in the old code was (a) calling HOLE serially over every
frame and (b) the file-juggling + per-segment matplotlib plotting. Both are
fixed here. If you ever want a genuinely different engine, CHAP is the modern
alternative, but it is harder to install and changes the methodology; for
consistency with your previous figures, parallel HOLE is the right choice.
"""

import os
import shutil
import tempfile
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import matplotlib
matplotlib.use("Agg")  # no GUI needed; we save PNG/PDF
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
from matplotlib.ticker import MultipleLocator

import MDAnalysis as mda
from MDAnalysis.analysis import hole2

warnings.filterwarnings("ignore")  # silence MDAnalysis PDB-writer warnings


# ============================= USER SETTINGS =============================
SYSTEM = "6IZF"  # used in titles and output filenames

base_dir   = "/media/supremeleader/Pantera/simulation/lipid_simulation/6IZF/charmm-gui-7977903380/gromacs/6izf"
output_dir = "/media/supremeleader/Pantera/simulation/analysis_2025/temporal_radius_profile_claude"

# 3 chains — same folder, numbered files
chains = {
    "Chain1": {"gro": os.path.join(base_dir, "chain1.gro"),
               "xtc": os.path.join(base_dir, "chain1.xtc")},
    "Chain2": {"gro": os.path.join(base_dir, "chain2.gro"),
               "xtc": os.path.join(base_dir, "chain2.xtc")},
    "Chain3": {"gro": os.path.join(base_dir, "chain3.gro"),
               "xtc": os.path.join(base_dir, "chain3.xtc")},
}

HOLE_EXECUTABLE = "/home/supremeleader/softwares/executable/hole2/exe/hole"

# --- selection: protein heavy atoms only (no hydrogens) -----------------
# Standard GROMACS/AMBER/CHARMM hydrogen names start with 'H'. If your force
# field uses number-leading hydrogen names (e.g. 1HB), change to:
#   "protein and not (name H* or name [0-9]H*)"
SELECT = "protein and not name H*"

# --- analysis window / sampling -----------------------------------------
LAST_NS       = 500     # analyse only the last 500 ns of each trajectory
FRAMES_PER_NS = None    # None -> auto-detect from trajectory dt.
                        # Set an int to force it (your old value was 100,
                        # i.e. 10 ps/frame). One HOLE profile is taken per ns.

# --- HOLE parameters (match your previous run) --------------------------
CVECT       = [0, 0, 1]  # channel axis (membrane normal)
CPOINT      = None       # Same as your trusted old run: let HOLE locate the pore.
                         # (Forcing a start point — including "center_of_geometry"
                         # — made it trace only the upper vestibule for this
                         # monomer, so we don't.)
                         # FALLBACK: if a chain ever traces only part of the pore,
                         # give an explicit lumen coordinate in Å instead, e.g.
                         # CPOINT = (65.968, 67.737, 43.463) from your old script.
END_RADIUS  = 10.0       # stop expanding the sphere beyond this radius (Å)
SAMPLE      = 0.20       # HOLE sample-plane spacing (Å); 0.20 = HOLE/old default
RANDOM_SEED = None       # Same as your trusted old run. A FIXED seed was what
                         # broke the profile (top-only, all dark blue), so we let
                         # HOLE seed itself. This no longer crashes because HOLE's
                         # scratch files now go to a short LOCAL path (SCRATCH_ROOT),
                         # which is the real cause of the original "gseed cannot
                         # open temporary file" error.
                         # If — and only if — you still hit a gseed crash, set this
                         # to an int AND set an explicit CPOINT coordinate above.

# --- common grid used for the heatmaps and the trimer average -----------
DZ            = 0.25     # channel-coordinate bin width (Å)
Z_LIMIT_LOW  = -30.0    # clamp the grid to a sane range (Å)
Z_LIMIT_HIGH =  30.0

# --- parallelism ---------------------------------------------------------
N_WORKERS        = os.cpu_count() or 4
BLOCKS_PER_CHAIN = max(1, N_WORKERS // 3)  # split each chain into N blocks
REUSE_CACHED     = True   # if a chain's .npz already exists, skip HOLE & reuse

# HOLE's scratch files MUST go on a short, LOCAL path (native ext4 etc.).
# Do NOT point this at the Pantera external drive — HOLE's Fortran temp-file
# opens fail there ("machine_dep.f ... No such file or directory").
SCRATCH_ROOT     = tempfile.gettempdir()   # usually /tmp ; or "/dev/shm" for RAM
# ========================================================================


# ---------------------------- HOLE worker -------------------------------
def _hole_worker(job):
    """Run HOLE on one frame-block of one chain. Returns {abs_frame: (z, r)}.

    HOLE writes its Fortran scratch files into its *working directory*. On a long
    path or an external/NTFS/exFAT mount those temp-file opens fail
    ("machine_dep.f ... No such file or directory"), so we chdir into a short,
    local scratch directory for the duration of the run and put tmpdir there too.
    """
    name, gro, xtc, start, stop, step, params = job

    orig_cwd = os.getcwd()
    workdir = tempfile.mkdtemp(prefix=f"hole_{name}_", dir=params["scratch_root"])
    try:
        os.chdir(workdir)                 # HOLE scratch files now land here (local)

        u = mda.Universe(gro, xtc)        # gro/xtc are absolute paths, unaffected
        ha = hole2.HoleAnalysis(
            u,
            select=params["select"],
            cpoint=params["cpoint"],
            cvect=params["cvect"],
            executable=params["executable"],
            end_radius=params["end_radius"],
            sample=params["sample"],
            tmpdir=".",                   # = workdir, the short local scratch dir
        )
        ha.run(start=start, stop=stop, step=step, random_seed=params["random_seed"])

        # results.profiles: {frame_index: recarray(rxn_coord, radius, cen_line_D)}
        profiles = ha.results.profiles if hasattr(ha, "results") else ha.profiles

        out = {}
        for fidx, prof in profiles.items():
            names = prof.dtype.names      # ('rxn_coord', 'radius', 'cen_line_D')
            # Channel coordinate = cen_line_D (HOLE column 2). Your old code used
            # line[2] for the y-axis ("Pore_coordinate"). cen_line_D runs the full
            # length of the pore; rxn_coord/cenxyz.cvec (column 0) collapses into a
            # narrow band for this monomer, which is what produced the broken plot.
            z = np.asarray(prof[names[2]], dtype=float)   # cen_line_D
            r = np.asarray(prof[names[1]], dtype=float)   # radius
            out[int(fidx)] = (z, r)

        try:
            ha.delete_temporary_files()
        except Exception:
            pass
    finally:
        os.chdir(orig_cwd)                # restore cwd (worker is reused for jobs)
        shutil.rmtree(workdir, ignore_errors=True)
    return name, out


# ------------------------- helper functions -----------------------------
def plan_chain(gro, xtc):
    """Work out frames_per_ns, the analysis window, and the analysed frames."""
    u = mda.Universe(gro, xtc)
    n_frames = u.trajectory.n_frames
    dt_ps = float(u.trajectory.dt)  # ps between saved frames
    if FRAMES_PER_NS is not None:
        fpns = int(FRAMES_PER_NS)
    else:
        fpns = max(1, int(round(1000.0 / dt_ps))) if dt_ps > 0 else 100

    step = fpns                       # 1 profile per ns
    window = LAST_NS * fpns           # frames spanning the last LAST_NS ns
    # Anchor the start on a whole-ns boundary so the time axis reads cleanly
    # (e.g. 500.0 ns, not 500.01). For a 0..N trajectory the last LAST_NS ns is
    # the closed interval [T_total - LAST_NS, T_total], i.e. LAST_NS+1 profiles.
    start = max(0, n_frames - 1 - window)
    stop = n_frames
    analysed = list(range(start, stop, step))
    return dict(n_frames=n_frames, dt_ps=dt_ps, fpns=fpns,
                start=start, stop=stop, step=step, analysed=analysed)


def build_matrix(raw, analysed, z_grid):
    """Interpolate every frame's (z, radius) onto z_grid -> 2D array (frames x z)."""
    R = np.full((len(analysed), len(z_grid)), np.nan)
    for i, fidx in enumerate(analysed):
        if fidx not in raw:
            continue
        z, r = raw[fidx]
        if z.size < 2:
            continue
        order = np.argsort(z)
        z, r = z[order], r[order]
        ri = np.interp(z_grid, z, r)
        ri[(z_grid < z.min()) | (z_grid > z.max())] = np.nan  # don't extrapolate
        R[i] = ri
    return R


def create_custom_colormap():
    colors = [
        (0.6, 0, 0),      # darkest red
        (1, 0, 0),        # red
        (1, 0.5, 0.5),    # lighter red
        (1, 0.75, 0),     # dark yellow
        (0, 0.6, 0),      # dark green
        (0.4, 0.7, 0.8),  # light blue
        (0, 0, 0.6),      # dark blue
    ]
    return LinearSegmentedColormap.from_list("custom_cmap", colors)


def plot_profile(time_ns, z_grid, R, title, save_base):
    """Heatmap: x = time (ns), y = channel coordinate (Å), colour = HOLE radius."""
    cmap = create_custom_colormap()
    cmap.set_bad("white")  # NaN (pore not sampled) -> white
    boundaries = [0, 0.6, 1, 1.3, 1.7, 2.2, 3.2, 5]
    norm = BoundaryNorm(boundaries, ncolors=cmap.N, clip=True)

    fig, ax = plt.subplots(figsize=(16, 4))
    mesh = ax.pcolormesh(time_ns, z_grid, R.T, cmap=cmap, norm=norm, shading="auto")

    cbar = fig.colorbar(mesh, ticks=boundaries)
    cbar.ax.set_yticklabels([str(b) for b in boundaries])
    cbar.set_label("HOLE radius (Å)")

    ax.set_ylim(-25, 25)
    ax.yaxis.set_major_locator(MultipleLocator(3))
    ax.set_xlim(time_ns.min(), time_ns.max())
    ax.xaxis.set_major_locator(MultipleLocator(25))   # tick every 25 ns
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Channel coordinate (Å)")
    ax.set_title(title)
    plt.tight_layout()

    fig.savefig(save_base + ".png", dpi=300)
    fig.savefig(save_base + ".pdf", format="pdf", dpi=300)
    plt.close(fig)
    print(f"  saved {os.path.basename(save_base)}.png / .pdf")


# --------------------------------- main ---------------------------------
def main():
    os.makedirs(output_dir, exist_ok=True)
    scratch_root = tempfile.mkdtemp(prefix="tric_hole_", dir=SCRATCH_ROOT)

    params = dict(select=SELECT, cpoint=CPOINT, cvect=CVECT,
                  executable=HOLE_EXECUTABLE, end_radius=END_RADIUS,
                  sample=SAMPLE, scratch_root=scratch_root, random_seed=RANDOM_SEED)

    # ---- 1. plan each chain and build the parallel job list ----
    meta = {}
    jobs = []
    print("Planning analysis windows:")
    for name, p in chains.items():
        m = plan_chain(p["gro"], p["xtc"])
        meta[name] = m
        total_ns = m["n_frames"] * m["dt_ps"] / 1000.0
        print(f"  {name}: {m['n_frames']} frames, dt={m['dt_ps']:.3f} ps "
              f"(~{total_ns:.1f} ns total) | frames/ns={m['fpns']} | "
              f"analysing {len(m['analysed'])} frames "
              f"(frames {m['start']}..{m['stop']} step {m['step']})")

        cache = os.path.join(output_dir, f"{SYSTEM}_{name}_profile.npz")
        if REUSE_CACHED and os.path.exists(cache):
            print(f"    -> cache found, will reuse {os.path.basename(cache)}")
            continue

        # split the analysed-frame list into contiguous blocks for parallelism
        for block in np.array_split(np.array(m["analysed"]), BLOCKS_PER_CHAIN):
            if block.size == 0:
                continue
            b_start, b_stop = int(block[0]), int(block[-1]) + 1
            jobs.append((name, p["gro"], p["xtc"], b_start, b_stop, m["step"], params))

    # ---- 2. run HOLE in parallel ----
    raw = {name: {} for name in chains}
    if jobs:
        print(f"\nRunning HOLE on {len(jobs)} blocks with {N_WORKERS} workers ...")
        with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
            for name, out in ex.map(_hole_worker, jobs):
                raw[name].update(out)
        for name in raw:
            if raw[name]:
                print(f"  {name}: {len(raw[name])} frames analysed")

    shutil.rmtree(scratch_root, ignore_errors=True)

    # ---- 3. assemble per-chain matrices (run cache OR fresh) on a common grid ----
    # First pass: gather all sampled z so every chain shares one grid.
    all_z = []
    loaded_cache = {}
    for name in chains:
        cache = os.path.join(output_dir, f"{SYSTEM}_{name}_profile.npz")
        if not raw[name] and REUSE_CACHED and os.path.exists(cache):
            d = np.load(cache)
            loaded_cache[name] = (d["time"], d["z_grid"], d["R"])
            all_z.append(d["z_grid"])
        else:
            for z, _ in raw[name].values():
                if z.size:
                    all_z.append(z)
    if not all_z:
        raise RuntimeError("No HOLE profiles were produced. Check paths/executable.")
    all_z = np.concatenate([np.atleast_1d(a) for a in all_z])

    z_lo = max(Z_LIMIT_LOW, float(np.floor(all_z.min())))
    z_hi = min(Z_LIMIT_HIGH, float(np.ceil(all_z.max())))
    z_grid = np.arange(z_lo, z_hi + DZ, DZ)

    # Second pass: build/realign each chain's matrix on z_grid and plot it.
    chain_R = {}
    for name in chains:
        if name in loaded_cache:
            t_cached, zg_cached, R_cached = loaded_cache[name]
            # re-interpolate cached R onto the (possibly new) common z_grid
            R = np.vstack([
                np.where((z_grid >= zg_cached.min()) & (z_grid <= zg_cached.max()),
                         np.interp(z_grid, zg_cached, row), np.nan)
                for row in R_cached
            ])
            time_ns = t_cached
        else:
            analysed = meta[name]["analysed"]
            n_ok = sum(1 for f in analysed if f in raw[name])
            if n_ok == 0:
                print(f"  WARNING: no profiles for {name}; skipping.")
                continue
            R = build_matrix(raw[name], analysed, z_grid)
            dt_ps = meta[name]["dt_ps"]
            # Real simulation time of each analysed frame (e.g. 500–1000 ns for
            # the last 500 ns of a 1 µs run), not a 0-based index.
            time_ns = np.array(analysed, dtype=float) * dt_ps / 1000.0
            np.savez(os.path.join(output_dir, f"{SYSTEM}_{name}_profile.npz"),
                     time=time_ns, z_grid=z_grid, R=R)

        chain_R[name] = (time_ns, R)
        plot_profile(time_ns, z_grid, R,
                     f"Temporal radius profile of {SYSTEM} {name}",
                     os.path.join(output_dir, f"{SYSTEM}_{name}_temporal_radius"))

    # ---- 4. trimer average (interpolated onto the shared grid) ----
    if len(chain_R) >= 2:
        n_min = min(R.shape[0] for _, R in chain_R.values())
        stack = np.stack([R[:n_min] for _, R in chain_R.values()], axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)  # all-NaN slices
            R_avg = np.nanmean(stack, axis=0)
        time_avg = next(iter(chain_R.values()))[0][:n_min]  # real MD time (ns)

        np.savez(os.path.join(output_dir, f"{SYSTEM}_trimer_average_profile.npz"),
                 time=time_avg, z_grid=z_grid, R=R_avg)
        plot_profile(time_avg, z_grid, R_avg,
                     f"Trimer-averaged temporal radius profile of {SYSTEM}",
                     os.path.join(output_dir, f"{SYSTEM}_trimer_average_temporal_radius"))
    else:
        print("Not enough chains succeeded to compute a trimer average.")

    print("\nDone. Outputs written to:", output_dir)


if __name__ == "__main__":
    main()

################################################################################
############# code for temporal radius profile using claude ####################
################################################################################


###################################################################
# ############# code for temporal radius profile ####################
# ###################################################################
# ################################## instruction for using the code ########################
# ##read the comments present at lines:
#
# ##########################################################################################
# import matplotlib
# import MDAnalysis as mda
# from MDAnalysis.analysis import hole2
# import numpy as np
# import shutil
# import glob
# import os
# import matplotlib.pyplot as plt
# from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
# from matplotlib.cm import ScalarMappable
#
# import warnings
# # suppress some MDAnalysis warnings when writing PDB files
# warnings.filterwarnings('ignore')
# # File paths
# ############################ provide here the gro and xtc file chain number for which you want to plot the pore profile ###############
# xtc_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/sodium/6iyz/chain3.xtc"
# gro_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/sodium/6iyz/chain3.gro"
# #######################################################################################################################################
# hole_executable = "/home/supremeleader/softwares/executable/hole2/exe/hole"
#
# # cpoint_coords = (65.968,67.737,43.463)     ### if you want to provide the centroid point in your analysis add here
#
# # Load the molecular dynamics trajectory
# u = mda.Universe(gro_file, xtc_file)
#
# # Run HOLE analysis with custom center point
# ha = hole2.HoleAnalysis(
#     u,
#     select='protein',
#     cpoint=None,  # Setting custom center point coordinates
#     executable=hole_executable,
#     cvect=[0,0,1],
#     end_radius= 10
# )
#
# ha.run(step=100)
# # ha.create_vmd_surface(filename='6iyz.vmd')
# # # ha.plot3D()
# #
#
# ############## move output file from a particular location to another #############################################
# def move_specific_txt_files(source_directory, destination_directory):
#     # Ensure the destination directory exists
#     os.makedirs(destination_directory, exist_ok=True)
#
#     # Find all files that start with 'hole' and end with '.txt' in the source directory
#     pattern = os.path.join(source_directory, "hole*.out")
#     files_to_move = glob.glob(pattern)
#
#     for file_path in files_to_move:
#         # Construct the destination path
#         destination_file = os.path.join(destination_directory, os.path.basename(file_path))
#
#         # Move the file
#         shutil.move(file_path, destination_file)
#         # print(f"Moved: {file_path} to {destination_file}")
#
# ############ provide the source and the destination where you want to copy and paste the file ##########################################
# move_specific_txt_files(
#     source_directory="/home/supremeleader/PycharmProjects/Programming",destination_directory="/media/supremeleader/Pantera/simulation/analysis_2024/pore_profile/output")
# ###############################################################################3
#
# # ha.create_vmd_surface(filename='holeanalysis.vmd')
# # ha.delete_temporary_files()
#
# ################### combining data from all the output file into one #################################################
#
# start_string = "cenxyz.cvec      radius  cen_line_D sum{s/(area point sourc"
# output_file_directory = "/media/supremeleader/Pantera/simulation/analysis_2024/pore_profile/output"
# def read_file_between_markers(file_path, start_marker):
#     with open(file_path, 'r') as file:
#         is_between_markers = False
#         result_lines = []
#         for line in file:
#             if start_marker in line:
#                 is_between_markers = True
#                 continue
#
#             if is_between_markers:
#                 if line.strip() == "" or line.startswith("Minimum"):
#                     break  # Stop reading when a blank line or summary line is found
#                 result_lines.append(line)
#     return result_lines
#
# list_of_output_files= os.listdir(output_file_directory)
# sorted_output_files=sorted(list_of_output_files)
#
# ##############################Path for single hole output file ######################
# ############################## provide the name of the outputfile according to the chain number #######################
# output_file_path = "/media/supremeleader/Pantera/simulation/analysis_2024/pore_profile/structures/6iyz_chain3.txt"
# #######################################################################################################################
# ##############################################################################################
# with open(output_file_path, 'w') as output_file:
#     for file in sorted_output_files:
#         output_file.write(f"Frame: {file[4:-4]}\nRadius,Pore_coordinate\n")
#         result_lines = read_file_between_markers(os.path.join(output_file_directory, file), start_string)
#         for line in result_lines:
#
#             line = line.strip().split()
#             if len(line) >= 3:
#                 output_file.write(f"{line[1]},{line[2]}\n")
#             else:
#                 output_file.write(f"Line does not have enough elements: {line}\n")
# ###############################################################################################################
#                     #### make publication quality figure for the monomers not for the trimer ###
#
# def read_data_from_file(file_path):
#     frame_data = {}
#     current_frame = None
#
#     with open(file_path, 'r') as file:
#         for line in file:
#             line = line.strip()
#
#             # Skip header lines like "Radius,Pore_coordinate"
#             if line.startswith("Radius") or line == "":
#                 continue
#
#             if line.startswith("Frame:"):
#                 # Get frame number from the line (after "Frame: ")
#                 current_frame = int(line.split()[1])
#                 frame_data[current_frame] = {"radius": [], "pore_coordinate": []}
#
#             elif line and "," in line:
#                 try:
#                     # Split the line into radius and pore coordinate
#                     radius, pore_coordinate = map(float, line.split(","))
#                     frame_data[current_frame]["radius"].append(radius)
#                     frame_data[current_frame]["pore_coordinate"].append(pore_coordinate)
#                 except ValueError:
#                     # Handle any potential conversion errors gracefully
#                     print(f"Skipping invalid line: {line}")
#
#     return frame_data
#
#
# def create_custom_colormap():
#     # Define a custom colormap that maps radius to specified colors
#     colors = [
#         (0.6, 0, 0),  # Darkest Red (for radius 0 to 0.6)
#         (1, 0, 0),  # Red (for radius 0.6 to 1)
#         (1, 0.5, 0.5),  # Lighter Red (for radius 1 to 1.3)
#         (1, 0.75, 0),  # Dark Yellow (for radius 1.3 to 1.7)
#         (0, 0.6, 0),  # Dark Green (for radius 1.7 to 2.1)
#         (0.4, 0.7, 0.8) ,  # Dark Blue (for radius 2.1 to 3.5)
#         (0, 0, 0.6)  # Black (for radius 3.5 to 4.5)
#     ]
#     return LinearSegmentedColormap.from_list("custom_cmap", colors)
#
#
# def plot_data(frame_data):
#     fig, ax = plt.subplots(figsize=(16, 4))  # Increase the figure size for better visibility
#
#     # Create custom colormap
#     cmap = create_custom_colormap()
#
#     # Define the boundaries and normalization for color mapping
#     boundaries = [0, 0.6, 1, 1.3, 1.7, 2.2, 3.2, 5]  # Updated boundaries for colormap with new black range
#     norm = BoundaryNorm(boundaries, ncolors=cmap.N, clip=True)
#
#     # Prepare a list to hold all segments to be colored
#     segments = []
#
#     # Loop over each frame and prepare data for plotting
#     for frame in sorted(frame_data.keys()):
#         radii = frame_data[frame]["radius"]
#         pore_coordinates = frame_data[frame]["pore_coordinate"]
#
#         for i in range(1, len(pore_coordinates)):
#             x_values = [frame, frame]  # Constant x value for the line segments
#             y_values = [pore_coordinates[i - 1], pore_coordinates[i]]
#
#             color_value = (radii[i - 1] + radii[i]) / 2  # Average radius for color
#             segments.append((x_values, y_values, color_value))
#
#     # Plot each segment with gradient line
#     for x_values, y_values, color_value in segments:
#         ax.plot(x_values, y_values, color=cmap(norm(color_value)),
#                 linewidth=2)  # Increase the line width for better visibility
#
#     # Create a ScalarMappable and add colorbar for radius gradient
#     sm = ScalarMappable(norm=norm, cmap=cmap)
#     sm.set_array([])
#     cbar = plt.colorbar(sm, ticks=boundaries)
#     cbar.ax.set_yticklabels(['0', '0.6', '1', '1.3', '1.7', '2.2', '3.2', '5'])
#     cbar.set_label('Radius')
#
#     # Set y-axis limits to show only the specified pore coordinate range
#     ax.set_ylim(-25, 25)
#
#     # Set y-axis tick spacing
#     ax.yaxis.set_major_locator(plt.MultipleLocator(3))
#
#     # Set x-axis limits to ensure the plot starts at 0 and ends at the last frame number
#     ax.set_xlim(min(frame_data.keys()), max(frame_data.keys()))
#
#     # Set axis labels and title
#     ax.set_xlabel('Frame Number (1 Frame = 1000 ps)')
#     ax.set_ylabel('Channel Coordinate (Å)')
#     ax.set_title('Temporal radius profile of 6IYZ chain3') # change the chain name for figure
#
#     # Adjust layout for publication quality
#     plt.tight_layout()
#     ####################################### save the output file at this path ####################################
#     ################## change the name of the pdf or the png file here in the below address ######################
#     save_path = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/sodium/6iyz/6iyz_hole_chain3.png"
#     plt.savefig(save_path, dpi=300)
#     save_path = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/sodium/6iyz/6iyz_hole_chain3.pdf"
#     plt.savefig(save_path, format='pdf', dpi=300)
#     ##############################################################################################################
#     ###############################################################################
#     # plt.show()
#
# # ha.create_vmd_surface(filename='6iyz_chain1.vmd')  # change the chain name for figure
# ha.delete_temporary_files()
#
# ############# Example usage: Single file containing all the hole output ##################
# ####### change the name of the output file here ##########################################
# file_path = "/media/supremeleader/Pantera/simulation/analysis_2024/pore_profile/structures/6iyz_chain3.txt"
#
# ##########################################################################################
# frame_data = read_data_from_file(file_path)
# #################################################################
# plot_data(frame_data)
#
# ############################### fig modification ############################
#
# # ################################## calculating average temporal radius profile ######################
# # file1="/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/temporal_radius_profile/6iyx_hole_chain1.txt"
# # file2="/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/temporal_radius_profile/6iyx_hole_chain2.txt"
# # file3="/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/temporal_radius_profile/6iyx_hole_chain3.txt"
# #
# # import os
# # import pandas as pd
# # import matplotlib.pyplot as plt
# # from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
# # from matplotlib.cm import ScalarMappable
# #
# # def read_data_from_files(file_list):
# #     all_data = []
# #     for file_path in file_list:
# #         current_frame = None
# #         with open(file_path, 'r') as file:
# #             for line in file:
# #                 line = line.strip()
# #
# #                 if line.startswith("Radius") or line == "":
# #                     continue
# #
# #                 if line.startswith("Frame:"):
# #                     current_frame = int(line.split()[1])
# #
# #                 elif line and "," in line:
# #                     try:
# #                         radius, pore_coordinate = map(float, line.split(","))
# #                         all_data.append({"frame": current_frame, "radius": radius, "pore_coordinate": pore_coordinate})
# #                     except ValueError:
# #                         print(f"Skipping invalid line: {line}")
# #
# #     return pd.DataFrame(all_data)
# #
# #
# # def calculate_average(data):
# #     return data.groupby(['frame', 'pore_coordinate']).radius.mean().reset_index()
# #
# #
# # def save_to_file(data, output_path):
# #     with open(output_path, 'w') as file:
# #         for frame in sorted(data['frame'].unique()):
# #             file.write(f"Frame: {frame:03d}\n")
# #             file.write("Radius,Pore_coordinate\n")
# #             frame_data = data[data['frame'] == frame]
# #             for _, row in frame_data.iterrows():
# #                 file.write(f"{row['radius']:.5f},{row['pore_coordinate']:.5f}\n")
# #             file.write("\n")
# #
# #
# # # Main function
# # def main():
# #     file_list = [
# #         file1,file2,file3
# #     ]
# #     data = read_data_from_files(file_list)
# #     averaged_data = calculate_average(data)
# #     save_to_file(averaged_data, "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/temporal_radius_profile/6iyx_average.txt")
# #     print("Averaged data saved to average.txt")
# #
# # if __name__ == "__main__":
# #     main()
# #
# #
# #
# # def read_data_from_file(file_path):
# #     frame_data = {}
# #     current_frame = None
# #
# #     with open(file_path, 'r') as file:
# #         for line in file:
# #             line = line.strip()
# #
# #             # Skip header lines like "Radius,Pore_coordinate"
# #             if line.startswith("Radius") or line == "":
# #                 continue
# #
# #             if line.startswith("Frame:"):
# #                 # Get frame number from the line (after "Frame: ")
# #                 current_frame = int(line.split()[1])
# #                 frame_data[current_frame] = {"radius": [], "pore_coordinate": []}
# #
# #             elif line and "," in line:
# #                 try:
# #                     # Split the line into radius and pore coordinate
# #                     radius, pore_coordinate = map(float, line.split(","))
# #                     frame_data[current_frame]["radius"].append(radius)
# #                     frame_data[current_frame]["pore_coordinate"].append(pore_coordinate)
# #                 except ValueError:
# #                     # Handle any potential conversion errors gracefully
# #                     print(f"Skipping invalid line: {line}")
# #
# #     return frame_data
# #
# #
# # def create_custom_colormap():
# #     # Define a custom colormap that maps radius to specified colors
# #     colors = [
# #         (0.6, 0, 0),  # Darkest Red (for radius 0 to 0.6)
# #         (1, 0, 0),  # Red (for radius 0.6 to 1)
# #         (1, 0.5, 0.5),  # Lighter Red (for radius 1 to 1.3)
# #         (1, 0.75, 0),  # Dark Yellow (for radius 1.3 to 1.7)
# #         (0, 0.6, 0),  # Dark Green (for radius 1.7 to 2.1)
# #         (0.4, 0.7, 0.8) ,  # Dark Blue (for radius 2.1 to 3.5)
# #         (0, 0, 0.6)  # Black (for radius 3.5 to 4.5)
# #     ]
# #     return LinearSegmentedColormap.from_list("custom_cmap", colors)
# #
# #
# # def plot_data(frame_data):
# #     fig, ax = plt.subplots(figsize=(16, 4))  # Increase the figure size for better visibility
# #
# #     # Create custom colormap
# #     cmap = create_custom_colormap()
# #
# #     # Define the boundaries and normalization for color mapping
# #     boundaries = [0, 0.6, 1, 1.3, 1.7, 2.2, 3.2, 5]  # Updated boundaries for colormap with new black range
# #     norm = BoundaryNorm(boundaries, ncolors=cmap.N, clip=True)
# #
# #     # Prepare a list to hold all segments to be colored
# #     segments = []
# #
# #     # Loop over each frame and prepare data for plotting
# #     for frame in sorted(frame_data.keys()):
# #         radii = frame_data[frame]["radius"]
# #         pore_coordinates = frame_data[frame]["pore_coordinate"]
# #
# #         for i in range(1, len(pore_coordinates)):
# #             x_values = [frame, frame]  # Constant x value for the line segments
# #             y_values = [pore_coordinates[i - 1], pore_coordinates[i]]
# #
# #             color_value = (radii[i - 1] + radii[i]) / 2  # Average radius for color
# #             segments.append((x_values, y_values, color_value))
# #
# #     # Plot each segment with gradient line
# #     for x_values, y_values, color_value in segments:
# #         ax.plot(x_values, y_values, color=cmap(norm(color_value)),
# #                 linewidth=2)  # Increase the line width for better visibility
# #
# #     # Create a ScalarMappable and add colorbar for radius gradient
# #     sm = ScalarMappable(norm=norm, cmap=cmap)
# #     sm.set_array([])
# #     cbar = plt.colorbar(sm, ticks=boundaries)
# #     cbar.ax.set_yticklabels(['0', '0.6', '1', '1.3', '1.7', '2.2', '3.2', '5'])
# #     cbar.set_label('Radius')
# #
# #     # Set y-axis limits to show only the specified pore coordinate range
# #     ax.set_ylim(-25, 25)
# #
# #     # Set y-axis tick spacing
# #     ax.yaxis.set_major_locator(plt.MultipleLocator(3))
# #
# #     # Set x-axis limits to ensure the plot starts at 0 and ends at the last frame number
# #     ax.set_xlim(min(frame_data.keys()), max(frame_data.keys()))
# #
# #     # Set axis labels and title
# #     ax.set_xlabel('Frame Number (1 Frame = 1000 ps)')
# #     ax.set_ylabel('Channel Coordinate (Å)')
# #     ax.set_title('Average temporal radius profile of 6IYX trimer')
# #
# #     # Adjust layout for publication quality
# #     plt.tight_layout()
# #     ####################################### save the output file at this path #####
# #     save_path = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/temporal_radius_profile/6iyx_average.png"
# #     plt.savefig(save_path, dpi=300)
# #     save_path = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/temporal_radius_profile/6iyx_average.pdf"
# #     plt.savefig(save_path, format='pdf', dpi=300)
# #     ###############################################################################
# #     plt.show()
# #
# # ############# Example usage: Single file containing all the hole output ##################
# # file_path = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/temporal_radius_profile/6iyx_average.txt"
# # frame_data = read_data_from_file(file_path)
# # #################################################################
# # plot_data(frame_data)
# #
# # #####################################################################################################
# #
# #
# # ################################ calaculate average ######################################
# # # Gather and print keys for inspection
# # # gathered = ha.gather()
# # # print("Keys in gathered profiles:", gathered.keys())
# # # print("Number of rxn_coords:", len(gathered['rxn_coord']))
# # #
# # # flat = ha.gather(flat=True)
# # # print("Number of rxn_coords in flat gathered:", len(flat['rxn_coord']))
# # #
# # # # Bin the radii and compute mean values
# # # radii, edges = ha.bin_radii(bins=100, range=None)
# # # means, edges = ha.histogram_radii(bins=100, range=None, aggregator=np.mean)
# # #
# # # # Calculate midpoints for plotting
# # # midpoints = 0.5 * (edges[1:] + edges[:-1])
# # #
# # # # Plot the results
# # # plt.plot(midpoints, means)
# # # plt.ylabel(r"Mean HOLE radius $R$ ($\AA$)")
# # # plt.xlabel(r"Pore coordinate $\zeta$ ($\AA$)")
# # # plt.title("Pore Profile Analysis")
# # # plt.show()

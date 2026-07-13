######################################################################################################################
############################# Claude verified helix-helix distance and angle ##################################
######################################################################################################################
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf_backend
import MDAnalysis as mda

# ===========================================================================
# === Setup ===
# ===========================================================================

base_dir   = "/media/supremeleader/Pantera/simulation/lipid_simulation/6IZF/charmm-gui-7977903380/gromacs/6izf"
output_dir = "/media/supremeleader/Pantera/simulation/apo_halo"
os.makedirs(output_dir, exist_ok=True)

# 3 chains — same folder, numbered files
chains = {
    "Chain1": {
        "gro": os.path.join(base_dir, "chain1.gro"),
        "xtc": os.path.join(base_dir, "chain1.xtc"),
    },
    "Chain2": {
        "gro": os.path.join(base_dir, "chain2.gro"),
        "xtc": os.path.join(base_dir, "chain2.xtc"),
    },
    "Chain3": {
        "gro": os.path.join(base_dir, "chain3.gro"),
        "xtc": os.path.join(base_dir, "chain3.xtc"),
    },
}

# Trajectory parameters (from diagnostic run)
TIMESTEP_PS     = 10.0        # ps per frame
CUTOFF_NS       = 500.0       # analyse last 500 ns
STRIDE          = 100         # read every 100th frame → 1 ns sampling
CUTOFF_PS       = (1000.0 - CUTOFF_NS) * 1000.0   # = 500,000 ps

# Sub-directories
for sub in ["plots_angle", "plots_distance", "chain_csv"]:
    os.makedirs(os.path.join(output_dir, sub), exist_ok=True)

# ===========================================================================
# === Helix definitions (6IYX/6IYZ) ===
# ===========================================================================
#
helices = {
    "H1":  (23,  39),  "H2a": (50,  62),  "H2b": (67,  73),
    "H3":  (84,  98),  "H4":  (112, 138), "H5a": (144, 155),
    "H5b": (159, 170), "H6":  (185, 202), "H7":  (209, 229),
}

########### 5WUC/5WUE ##########
# helices = {
#     "H1": (7, 31), "H2a": (34, 47), "H2b": (51, 57), "H3": (67, 86),
#     "H4": (94, 118), "H5a": (120, 133), "H5b": (137, 143),
#     "H6": (152, 170), "H7": (175, 197)
# }
#################################

pairs = [
    (h1, h2)
    for i, h1 in enumerate(helices)
    for h2 in list(helices)[i + 1:]
]
print(f"Helix pairs to analyse: {len(pairs)}")


# ===========================================================================
# === Core geometry functions (all use SVD axis — single consistent method) ===
# ===========================================================================

def compute_svd_axis(coords):
    """Best-fit helix axis via SVD. Returns (centroid, unit_axis)."""
    centroid = coords.mean(axis=0)
    _, _, vh = np.linalg.svd(coords - centroid, full_matrices=False)
    return centroid, vh[0]


def get_endpoints(coords, centroid, axis):
    """Project all atoms; use min/max projections as helix endpoints."""
    proj  = np.dot(coords - centroid, axis)
    return centroid + proj.min() * axis, centroid + proj.max() * axis


def angle_between(v1, v2):
    """Angle in [0°, 90°] between two axes (sign-ambiguity corrected)."""
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    ang = np.degrees(np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0)))
    ang = ang % 360.0
    if ang > 180.0: ang = 360.0 - ang
    if ang > 90.0:  ang = 180.0 - ang
    return ang


def segment_distance(b1, e1, b2, e2):
    """Minimum distance between two finite line segments."""
    u, v   = e1 - b1, e2 - b2
    w0     = b1 - b2
    a, b_, c = np.dot(u, u), np.dot(u, v), np.dot(v, v)
    d, e_  = np.dot(u, w0), np.dot(v, w0)
    denom  = a * c - b_ * b_
    if denom != 0:
        s = (b_ * e_ - c  * d) / denom
        t = (a  * e_ - b_ * d) / denom
    else:
        s, t = 0.0, (e_ / c if c > 1e-8 else 0.0)
    s, t = np.clip(s, 0, 1), np.clip(t, 0, 1)
    return np.linalg.norm((b1 + s * u) - (b2 + t * v))


# ===========================================================================
# === Process each chain ===
# ===========================================================================

# all_results[chain_name][(h1,h2)] = {"frames": arr, "angles": arr, "distances": arr}
all_results = {}

for chain_name, files in chains.items():
    print(f"\n{'='*60}")
    print(f"Processing {chain_name} ...")
    u = mda.Universe(files["gro"], files["xtc"])

    # --- Select frames in last 500 ns using time stamp ---
    selected_indices = [
        ts.frame for ts in u.trajectory
        if ts.time >= CUTOFF_PS
    ]
    # Apply stride within the selected window
    selected_indices = selected_indices[::STRIDE]
    print(f"  Total trajectory frames : {len(u.trajectory)}")
    print(f"  Frames in last 500 ns   : {len(selected_indices)} (after stride={STRIDE})")

    # --- Pre-compute axes for selected frames ---
    axis_data    = {name: [] for name in helices}
    helix_coords = {name: [] for name in helices}
    frame_times  = []

    for idx in selected_indices:
        ts = u.trajectory[idx]
        frame_times.append(ts.time / 1000.0)   # store in ns
        for name, (start, end) in helices.items():
            coords = u.select_atoms(f"name CA and resid {start}-{end}").positions.copy()
            centroid, axis = compute_svd_axis(coords)
            axis_data[name].append((centroid, axis))
            helix_coords[name].append(coords)

    frame_times = np.array(frame_times)

    # --- Save axis files ---
    for name, data in axis_data.items():
        out_path = os.path.join(output_dir, "chain_csv", f"{chain_name}_{name}_Axis.txt")
        with open(out_path, 'w') as f:
            f.write("Time_ns\tCentroid_X\tCentroid_Y\tCentroid_Z\tAxis_X\tAxis_Y\tAxis_Z\n")
            for t, (cen, ax) in zip(frame_times, data):
                f.write(f"{t:.3f}\t{cen[0]:.3f}\t{cen[1]:.3f}\t{cen[2]:.3f}\t"
                        f"{ax[0]:.5f}\t{ax[1]:.5f}\t{ax[2]:.5f}\n")

    # --- Compute angles & distances for every pair ---
    chain_results = {}
    for h1, h2 in pairs:
        angles, distances = [], []
        for i in range(len(frame_times)):
            c1, ax1 = axis_data[h1][i]
            c2, ax2 = axis_data[h2][i]

            angles.append(angle_between(ax1, ax2))

            co1, co2 = helix_coords[h1][i], helix_coords[h2][i]
            b1, e1   = get_endpoints(co1, c1, ax1)
            b2, e2   = get_endpoints(co2, c2, ax2)
            distances.append(segment_distance(b1, e1, b2, e2))

        angles    = np.array(angles)
        distances = np.array(distances)

        # Save per-chain CSV
        csv_path = os.path.join(output_dir, "chain_csv",
                                f"{chain_name}_{h1}-{h2}_angle_distance.csv")
        with open(csv_path, 'w') as f:
            f.write("Time_ns,Crossing_Angle_Degrees,Distance_Ang\n")
            for t, ang, dist in zip(frame_times, angles, distances):
                f.write(f"{t:.3f},{ang:.2f},{dist:.3f}\n")

        chain_results[(h1, h2)] = {
            "frame_times": frame_times,
            "angles":      angles,
            "distances":   distances,
        }

    all_results[chain_name] = chain_results
    print(f"  {chain_name} done — {len(pairs)} pairs saved to chain_csv/")


# ===========================================================================
# === Compute cross-chain average for every pair ===
# ===========================================================================

print("\nComputing cross-chain averages...")
avg_results = {}   # avg_results[(h1,h2)] = {"frame_times", "angles", "distances"}

# Use Chain1 frame_times as reference (all chains have identical time axes)
ref_times = all_results["Chain1"][pairs[0]]["frame_times"]

for h1, h2 in pairs:
    ang_stack  = np.stack([all_results[cn][(h1, h2)]["angles"]    for cn in chains], axis=0)
    dist_stack = np.stack([all_results[cn][(h1, h2)]["distances"] for cn in chains], axis=0)
    avg_ang    = ang_stack.mean(axis=0)
    avg_dist   = dist_stack.mean(axis=0)

    avg_results[(h1, h2)] = {
        "frame_times": ref_times,
        "angles":      avg_ang,
        "distances":   avg_dist,
    }

    # Save average CSV
    csv_path = os.path.join(output_dir, f"AVG_{h1}-{h2}_angle_distance.csv")
    with open(csv_path, 'w') as f:
        f.write("Time_ns,Avg_Crossing_Angle_Degrees,Avg_Distance_Ang\n")
        for t, ang, dist in zip(ref_times, avg_ang, avg_dist):
            f.write(f"{t:.3f},{ang:.2f},{dist:.3f}\n")

print("  Average CSVs saved.")


# ===========================================================================
# === Plotting ===
# ===========================================================================

CHAIN_COLORS = {
    "Chain1": "#2166ac",   # blue
    "Chain2": "#4dac26",   # green
    "Chain3": "#d6604d",   # red-orange
}
AVG_COLOR  = "black"
AVG_LW     = 2.0
CHAIN_LW   = 0.9
CHAIN_ALPHA = 0.65


def plot_pair(ax, h1, h2, quantity):
    """
    Draw one panel: all 3 chains (thin, coloured) + average (thick black)
    + dashed line at grand mean of the average.
    quantity: 'angles' or 'distances'
    """
    is_angle = quantity == "angles"
    ylabel   = "Crossing Angle (°)" if is_angle else "Inter-helix Distance (Å)"

    for chain_name, chain_res in all_results.items():
        res    = chain_res[(h1, h2)]
        times  = res["frame_times"]
        values = res[quantity]
        ax.plot(times, values,
                color=CHAIN_COLORS[chain_name], linewidth=CHAIN_LW,
                alpha=CHAIN_ALPHA, label=chain_name)

    avg_res    = avg_results[(h1, h2)]
    avg_values = avg_res[quantity]
    avg_times  = avg_res["frame_times"]
    grand_mean = avg_values.mean()

    ax.plot(avg_times, avg_values,
            color=AVG_COLOR, linewidth=AVG_LW, label="Trimer Avg")
    ax.axhline(grand_mean, color=AVG_COLOR, linewidth=1.2, linestyle="--",
               label=f"Mean = {grand_mean:.2f}")

    ax.set_xlabel("Time (ns)", fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_title(f"{h1}–{h2}", fontsize=9, fontweight="bold")
    ax.legend(fontsize=6, ncol=2)
    ax.tick_params(labelsize=7)
    ax.set_xlim(avg_times[0], avg_times[-1])


print("\nGenerating plots...")

# ---------- individual PNGs ----------
for h1, h2 in pairs:
    pair_label = f"{h1}-{h2}"

    # Angle PNG
    fig, ax = plt.subplots(figsize=(8, 3.8))
    plot_pair(ax, h1, h2, "angles")
    ax.set_title(f"Helix Crossing Angle  |  {pair_label}  (last 500 ns)",
                 fontsize=10, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "plots_angle", f"{pair_label}_angle.png"),
                dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Distance PNG
    fig, ax = plt.subplots(figsize=(8, 3.8))
    plot_pair(ax, h1, h2, "distances")
    ax.set_title(f"Helix Axis Distance  |  {pair_label}  (last 500 ns)",
                 fontsize=10, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "plots_distance", f"{pair_label}_distance.png"),
                dpi=150, bbox_inches="tight")
    plt.close(fig)

print(f"  {len(pairs)} angle PNGs  → plots_angle/")
print(f"  {len(pairs)} distance PNGs → plots_distance/")

# ---------- combined PDFs (4 panels per page) ----------
def save_combined_pdf(pairs_list, quantity, pdf_path, suptitle):
    n_per_page = 4
    with pdf_backend.PdfPages(pdf_path) as pdf:
        for page_start in range(0, len(pairs_list), n_per_page):
            page_pairs = pairs_list[page_start: page_start + n_per_page]
            n          = len(page_pairs)
            fig, axes  = plt.subplots(2, 2, figsize=(14, 8))
            axes_flat  = axes.flatten()
            for idx, (h1, h2) in enumerate(page_pairs):
                plot_pair(axes_flat[idx], h1, h2, quantity)
            for idx in range(n, n_per_page):
                axes_flat[idx].set_visible(False)
            fig.suptitle(suptitle, fontsize=12, fontweight="bold", y=1.01)
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    print(f"  Saved: {pdf_path}")


save_combined_pdf(
    pairs, "angles",
    os.path.join(output_dir, "ALL_PAIRS_angles.pdf"),
    "Helix Crossing Angles — All Pairs (last 500 ns, 3 chains + average)",
)
save_combined_pdf(
    pairs, "distances",
    os.path.join(output_dir, "ALL_PAIRS_distances.pdf"),
    "Helix Axis Distances — All Pairs (last 500 ns, 3 chains + average)",
)

# ===========================================================================
print("\n" + "="*60)
print("ALL DONE!")
print(f"  Per-chain CSVs (angle+distance) → {output_dir}/chain_csv/")
print(f"  Average CSVs                    → {output_dir}/AVG_*.csv")
print(f"  Per-pair angle PNGs             → {output_dir}/plots_angle/")
print(f"  Per-pair distance PNGs          → {output_dir}/plots_distance/")
print(f"  Combined angle PDF              → {output_dir}/ALL_PAIRS_angles.pdf")
print(f"  Combined distance PDF           → {output_dir}/ALL_PAIRS_distances.pdf")
print("="*60)

######################################################################################################################
############################# Claude verified helix-helix distance and angle #########################################
######################################################################################################################


# import MDAnalysis as mda
# import numpy as np
# import matplotlib.pyplot as plt
#
# gro_file = "/media/supremeleader/Pantera/simulation/analysis_2024/helix_angle/chain3_tm5.gro"
# xtc_file = "/media/supremeleader/Pantera/simulation/analysis_2024/helix_angle/chain3_tm5.xtc"
#
# u = mda.Universe(gro_file, xtc_file)
# helix_part1={}
#
# for ts in u.trajectory[::100]:
#     part1 = u.select_atoms('name CA and resnum 144-170')
#     helix_part1[ts.frame] = []
#     for atom in part1:
#         name = atom.name
#         resid = atom.resid
#         resname = atom.resname
#         x, y, z = atom.position
#         helix_part1[ts.frame].append([resid,resname,x ,y, z])
#
# def calculate_kink_angle(helix_part1):
#     kink_angles = {}
#
#     for timestamp, residues in helix_part1.items():
#         vectors = []
#
#         # Extract coordinates and calculate direction vectors for segments
#         for i in range(1, len(residues) - 1):
#             # Get Cα coordinates for three consecutive residues
#             coord1 = np.array(residues[i - 1][2:5])
#             coord2 = np.array(residues[i][2:5])
#             coord3 = np.array(residues[i + 1][2:5])
#
#             # Calculate vectors between consecutive Cα atoms
#             vec1 = coord2 - coord1
#             vec2 = coord3 - coord2
#
#             # Normalize the bisector vector of vec1 and vec2
#             bisector = (vec1 + vec2) / np.linalg.norm(vec1 + vec2)
#             vectors.append(bisector)
#
#         # Calculate angles between consecutive vectors to find the kink
#         angles = []
#         for j in range(1, len(vectors)):
#             # Compute the angle between two consecutive bisectors
#             dot_product = np.dot(vectors[j - 1], vectors[j])
#             angle = np.arccos(dot_product) * (180 / np.pi)  # Convert to degrees
#             angles.append(angle)
#
#         # Identify the smallest angle as the kink angle
#         kink_angle = min(angles) if angles else None
#         kink_angles[timestamp] = kink_angle
#
#     return kink_angles
#
#
# # Calculate kink angles for the dictionary
# kink_angles = calculate_kink_angle(helix_part1)
#
# # Print kink angles for each timestamp
# # for timestamp, angle in kink_angles.items():
# #     print(f"{int(timestamp/100)}:{angle:.2f}")
#
#  # Output the angles
# with open("/media/supremeleader/Pantera/simulation/analysis_2024/helix_angle/helix5_chain3.txt", "w") as file:
#     for ts, angle in kink_angles.items():
#         file.write(f"{int(ts/100)}: {angle:.2f} \n")
#
# ###################################### figure #######################################
# # Define the path to your file
# file_path = "/media/supremeleader/Pantera/simulation/analysis_2024/helix_angle/helix5_chain3.txt"
#
# # Initialize lists for time and angle
# time_ns = []
# kink_angle_degrees = []
#
# # Read the data from the file
# with open(file_path, 'r') as file:
#     for line in file:
#         # Each line should be in the format "time: angle"
#         time, angle = line.split(":")
#         time_ns.append(int(time.strip()))
#         kink_angle_degrees.append(float(angle.strip()))
#
# # Plotting
# plt.figure(figsize=(12, 5))
# plt.plot(time_ns, kink_angle_degrees, color='b', linestyle='-', linewidth=1.5, markersize=6)
# plt.xlabel("Time (ns)")
# plt.ylabel("Kink Angle (degrees)")
# plt.title("TM Helix-V (Chain-III) Kink Angle over Time")
# plt.xlim([min(time_ns), max(time_ns)])
# ## plt.grid(True)
# plt.savefig("/media/supremeleader/Pantera/simulation/analysis_2024/helix_angle/chain3_helix5.png")
# plt.savefig("/media/supremeleader/Pantera/simulation/analysis_2024/helix_angle/chain3_helix5.pdf")
## plt.show()



########################################################################################
###################### average angle for trimer and the figure #########################
# import numpy as np
# import matplotlib.pyplot as plt
#
# # File paths
# file1 = "/media/supremeleader/Pantera/simulation/analysis_2024/helix_angle/helix5_chain1.txt"
# file2 = "/media/supremeleader/Pantera/simulation/analysis_2024/helix_angle/helix5_chain2.txt"
# file3 = "/media/supremeleader/Pantera/simulation/analysis_2024/helix_angle/helix5_chain3.txt"
# output_file = "/media/supremeleader/Pantera/simulation/analysis_2024/helix_angle/trimer_h5_average.txt"
#
# # Read data from each file
# def read_angles(file_path):
#     angles = {}
#     with open(file_path, 'r') as f:
#         for line in f:
#             timestamp, angle = line.strip().split(':')
#             angles[int(timestamp)] = float(angle)
#     return angles
#
# # Load angles from each file
# angles1 = read_angles(file1)
# angles2 = read_angles(file2)
# angles3 = read_angles(file3)
#
# # Calculate the average of the angles for each timestamp
# timestamps = sorted(angles1.keys())
# averages = {ts: np.mean([angles1[ts], angles2[ts], angles3[ts]]) for ts in timestamps}
#
# # Save the averaged angles to the output file
# with open(output_file, 'w') as f:
#     for ts, avg_angle in averages.items():
#         f.write(f"{ts}: {avg_angle:.2f}\n")
#
# # Print the contents of the new file
# with open(output_file, 'r') as f:
#     print(f.read())
#
# ########### figure ############
#
# # Define the path to your file
# file_path = "/media/supremeleader/Pantera/simulation/analysis_2024/helix_angle/trimer_h5_average.txt"
#
# # Initialize lists for time and angle
# time_ns = []
# kink_angle_degrees = []
#
# # Read the data from the file
# with open(file_path, 'r') as file:
#     for line in file:
#         # Each line should be in the format "time: angle"
#         time, angle = line.split(":")
#         time_ns.append(int(time.strip()))
#         kink_angle_degrees.append(float(angle.strip()))
#
# # Plotting
# plt.figure(figsize=(12, 5))
# plt.plot(time_ns, kink_angle_degrees, color='b', linestyle='-', linewidth=1.5, markersize=6)
# plt.xlabel("Time (ns)")
# plt.ylabel("Kink Angle (degrees)")
# plt.title(" Average TM Helix-V Kink Angle over Time for 6IYX Trimer")
# plt.xlim([min(time_ns), max(time_ns)])
# ## plt.grid(True)
# plt.savefig("/media/supremeleader/Pantera/simulation/analysis_2024/helix_angle/trimer_helix5.png")
# plt.savefig("/media/supremeleader/Pantera/simulation/analysis_2024/helix_angle/trimer_helix5.pdf")
# # plt.show()
#################################################################################################
################################### least square distance ############################################
# print(helix_part2)
#
# ####     """Calculate the centroid of a list of 3D points.""" ######
# def calculate_centroid(points):
#     return np.mean(points, axis=0)
#
# ####  """Calculate the least-squares axis for a set of points."""  ####
# def calculate_least_squares_axis(points):
#     # Step 1: Calculate the centroid
#     centroid = calculate_centroid(points)
#     # Step 2: Center the points around the origin
#     centered_points = points - centroid
#     # Step 3: Construct the covariance matrix
#     covariance_matrix = np.cov(centered_points, rowvar=False)
#     # Step 4: Find the eigenvalues and eigenvectors
#     eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)
#     # Step 5: Select the eigenvector with the largest eigenvalue
#     principal_axis = eigenvectors[:, np.argmax(eigenvalues)]
#     return principal_axis
#
# ### """Calculate the angle in degrees between two axes (direction vectors).""" ###
# def calculate_angle_between_axes(axis1, axis2):
#     # Normalize the vectors
#     unit_axis1 = axis1 / np.linalg.norm(axis1)
#     unit_axis2 = axis2 / np.linalg.norm(axis2)
#     # Calculate the cosine of the angle between the vectors
#     cos_theta = np.dot(unit_axis1, unit_axis2)
#     # Ensure cosine is within the valid range for arccos due to floating point errors
#     cos_theta = np.clip(cos_theta, -1.0, 1.0)
#     # Calculate the angle in radians and convert to degrees
#     angle_radians = np.arccos(cos_theta)
#     angle_degrees = np.degrees(angle_radians)
#     return angle_degrees
#
# # Calculate angles for each timestamp
# angles = {}
# for time_stamp in helix_part1:
#     if time_stamp in helix_part2:  # Ensure both parts have the same timestamp
#         # Convert the lists to numpy arrays
#         points1 = np.array(helix_part1[time_stamp])
#         points2 = np.array(helix_part2[time_stamp])
#
#         # Calculate least-squares axis for each part
#         axis1 = calculate_least_squares_axis(points1)
#         axis2 = calculate_least_squares_axis(points2)
#
#         # Calculate the angle between the two axes
#         angle = calculate_angle_between_axes(axis1, axis2)
#         angles[time_stamp] = angle
#
# # Output the angles
# with open("/media/supremeleader/Pantera/simulation/analysis_2024/helix_angle/helix2_chain1.txt", "w") as file:
#     for ts, angle in angles.items():
#         file.write(f"{int(ts/100)}: {angle:.2f} \n")
#
# ###################################### figure #######################################
# # Define the path to your file
# file_path = "/media/supremeleader/Pantera/simulation/analysis_2024/helix_angle/helix2_chain1.txt"
#
# # Initialize lists for time and angle
# time_ns = []
# kink_angle_degrees = []
#
# # Read the data from the file
# with open(file_path, 'r') as file:
#     for line in file:
#         # Each line should be in the format "time: angle"
#         time, angle = line.split(":")
#         time_ns.append(int(time.strip()))
#         kink_angle_degrees.append(float(angle.strip()))
#
# # Plotting
# plt.figure(figsize=(10, 6))
# plt.plot(time_ns, kink_angle_degrees, color='b', linestyle='-', linewidth=1.5, markersize=6)
# plt.xlabel("Time (ns)")
# plt.ylabel("Kink Angle (degrees)")
# plt.title("TM Helix-II (Chain-I) Kink Angle over Time")
# plt.tight_layout()
# plt.xlim([min(time_ns), max(time_ns)])
# plt.grid(True)
# # plt.show()
# plt.savefig("/media/supremeleader/Pantera/simulation/analysis_2024/helix_angle/helix2_chain1.png")
# plt.savefig("/media/supremeleader/Pantera/simulation/analysis_2024/helix_angle/helix2_chain1.pdf")

#####################################################################################################
# def calculate_kink_angles_with_details(helix_part1):
#     all_kink_angles = {}
#
#     for timestamp, residues in helix_part1.items():
#         vectors = []
#         angle_details = []
#
#         # Extract coordinates and calculate direction vectors for segments
#         for i in range(1, len(residues) - 1):
#             # Get Cα coordinates for three consecutive residues
#             res1, res2, res3 = residues[i - 1], residues[i], residues[i + 1]
#             coord1, coord2, coord3 = np.array(res1[2:5]), np.array(res2[2:5]), np.array(res3[2:5])
#
#             # Calculate vectors between consecutive Cα atoms
#             vec1 = coord2 - coord1
#             vec2 = coord3 - coord2
#
#             # Normalize the bisector vector of vec1 and vec2
#             bisector = (vec1 + vec2) / np.linalg.norm(vec1 + vec2)
#             vectors.append((bisector, res1[1], res2[1], res3[1]))  # Store bisector and residue names
#
#         # Calculate angles between consecutive vectors to find the kink
#         angles = []
#         for j in range(1, len(vectors)):
#             # Compute the angle between two consecutive bisectors
#             vec1, res1_name, res2_name, res3_name = vectors[j - 1]
#             vec2, res4_name, res5_name, res6_name = vectors[j]
#             dot_product = np.dot(vec1, vec2)
#             angle = np.arccos(dot_product) * (180 / np.pi)  # Convert to degrees
#
#             # Store angle and residue details
#             angle_info = {
#                 "angle": angle,
#                 "residues": [(res1_name, res2_name, res3_name), (res4_name, res5_name, res6_name)]
#             }
#             angles.append(angle_info)
#
#         # Identify the minimum angle as the kink angle and print all angles
#         kink_angle = min(angles, key=lambda x: x["angle"]) if angles else None
#         all_kink_angles[timestamp] = {"kink_angle": kink_angle, "angles": angles}
#
#         # Print results
#         print(f"Timestamp {timestamp}:")
#         for idx, angle_info in enumerate(angles):
#             angle_value = angle_info["angle"]
#             res_group1, res_group2 = angle_info["residues"]
#             print(f"  Angle {idx + 1}: {angle_value:.2f}° between residues {res_group1} and {res_group2}")
#
#         if kink_angle:
#             print(f"  Kink Angle: {kink_angle['angle']:.2f}° (Smallest angle)")
#         print()
#
#     return all_kink_angles
#
#
# # Calculate and print kink angles with details for each timestamp
# kink_angles_with_details = calculate_kink_angles_with_details(helix_part1)



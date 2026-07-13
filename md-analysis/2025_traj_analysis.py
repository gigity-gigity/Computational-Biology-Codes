
############################################################################################################################################
######################### Use this code for calculation of inter-helix distance and crossing angle #########################################
############################################################################################################################################
################ this is the final code using most appropriate mathematical principle for calculation of distance and angle ################
############################################################################################################################################
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda

# === Setup ===
gro_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/chain1_new.gro"
xtc_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/chain1_new.xtc"
output_dir = "/media/supremeleader/Pantera/simulation/analysis_2025/comparison/new_results"
os.makedirs(output_dir, exist_ok=True)

# === Helix definitions ===
########## 6IYX/6IYZ #########
helices = {
    "H1": (23, 39), "H2a": (50, 62), "H2b": (67, 73), "H3": (84, 98),
    "H4": (112, 138), "H5a": (144, 155), "H5b": (159, 170),
    "H6": (185, 202), "H7": (209, 229)
}
########### 5WUC/5WUE ##########
# helices = {
#     "H1": (7, 31), "H2a": (34, 47), "H2b": (51, 57), "H3": (67, 86),
#     "H4": (94, 118), "H5a": (120, 133), "H5b": (137, 143),
#     "H6": (152, 170), "H7": (175, 197)
# }
#################################

# === Load universe ===
u = mda.Universe(gro_file, xtc_file)

# === Step 1: Calculate helix axis data ===
axis_data = {name: [] for name in helices}
helix_coords = {name: [] for name in helices}
frames = []

for ts in u.trajectory[::100]:
    frames.append(ts.frame)
    for name, (start, end) in helices.items():
        sel = u.select_atoms(f"name CA and resid {start}-{end}")
        coords = sel.positions
        centroid = coords.mean(axis=0)
        centered = coords - centroid
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
        axis = vh[0]
        axis_data[name].append((centroid, axis))
        helix_coords[name].append(coords.copy())

# === Save helix axis data to individual files ===
for name, data in axis_data.items():
    out_path = os.path.join(output_dir, f"{name}_Axis.txt")
    with open(out_path, 'w') as f:
        f.write("Frame\tCentroid_X\tCentroid_Y\tCentroid_Z\tAxis_X\tAxis_Y\tAxis_Z\n")
        for frame, (centroid, axis) in zip(frames, data):
            f.write(f"{frame}\t{centroid[0]:.3f}\t{centroid[1]:.3f}\t{centroid[2]:.3f}\t"
                    f"{axis[0]:.5f}\t{axis[1]:.5f}\t{axis[2]:.5f}\n")

# === Utilities for angle ===
def angle_between_vectors(v1, v2):
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    dot = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
    return np.degrees(np.arccos(dot))

def quadrant_angle(angle):
    angle = angle % 360.0
    if angle > 180.0:
        angle = 360.0 - angle
    if angle > 90.0:
        angle = 180.0 - angle
    return angle

# === New method for distance ===
def unit_vector(v):
    norm = np.linalg.norm(v)
    return v / norm if norm > 1e-8 else np.zeros_like(v)

def compute_principal_axis(coords):
    if len(coords) < 5:
        raise ValueError("Helix too short to compute principal axis.")
    p = coords[1:-2] + coords[3:] - 2 * coords[:-3]
    p_mean = np.mean(p, axis=0)
    centered = p - p_mean
    inertia_tensor = np.dot(centered.T, centered)
    eigvals, eigvecs = np.linalg.eigh(inertia_tensor)
    axis = eigvecs[:, np.argmin(eigvals)]
    return unit_vector(axis)

def get_endpoints(coords, axis):
    center = coords.mean(axis=0)
    b = center + np.dot(coords[0] - center, axis) * axis
    e = center + np.dot(coords[-1] - center, axis) * axis
    return b, e

def closest_distance_between_axes(b1, e1, b2, e2):
    u = e1 - b1
    v = e2 - b2
    w0 = b1 - b2
    a, b_, c = np.dot(u, u), np.dot(u, v), np.dot(v, v)
    d, e_ = np.dot(u, w0), np.dot(v, w0)
    denom = a * c - b_ * b_
    s = (b_ * e_ - c * d) / denom if denom != 0 else 0
    t = (a * e_ - b_ * d) / denom if denom != 0 else 0
    s = np.clip(s, 0, 1)
    t = np.clip(t, 0, 1)
    p1 = b1 + s * u
    p2 = b2 + t * v
    dist = np.linalg.norm(p1 - p2)
    return dist

# === Step 2 & 3: Calculate angle and distance for each helix pair ===
pairs = [(h1, h2) for i, h1 in enumerate(helices) for h2 in list(helices)[i+1:]]

for h1, h2 in pairs:
    out_combined_file = os.path.join(output_dir, f"{h1}-{h2}_angle_distance.csv")
    with open(out_combined_file, 'w') as f:
        f.write("Frame,Crossing_Angle_Degrees,Distance\n")
        for i, frame in enumerate(frames):
            c1, axis1 = axis_data[h1][i]
            c2, axis2 = axis_data[h2][i]
            angle = angle_between_vectors(axis1, axis2)
            angle = quadrant_angle(angle)
            coords1 = helix_coords[h1][i]
            coords2 = helix_coords[h2][i]
            pa1 = compute_principal_axis(coords1)
            pa2 = compute_principal_axis(coords2)
            b1, e1 = get_endpoints(coords1, pa1)
            b2, e2 = get_endpoints(coords2, pa2)
            dist = closest_distance_between_axes(b1, e1, b2, e2)
            f.write(f"{frame},{angle:.2f},{dist:.3f}\n")

print("All helix analysis (axis, angle, new distance) completed.")

############################################################################################################################################
############################################################################################################################################
######################################## Code for figure and comparing two structures #################################################
############################################################################################################################################

import os
import pandas as pd
import matplotlib.pyplot as plt

folder_6iyx = "/media/supremeleader/Pantera/simulation/analysis_2025/comparison/5wuc"
folder_6iyz = "/media/supremeleader/Pantera/simulation/analysis_2025/comparison/5wue"
output_dir = "/media/supremeleader/Pantera/simulation/analysis_2025/comparison"

# List all files from 6iyx folder
files = sorted([f for f in os.listdir(folder_6iyx) if f.endswith("_angle_distance.csv")])

for file in files:
    path_6iyx = os.path.join(folder_6iyx, file)
    path_6iyz = os.path.join(folder_6iyz, file)

    if not os.path.exists(path_6iyz):
        print(f"Skipping {file}: not found in 6iyz folder")
        continue

    df_6iyx = pd.read_csv(path_6iyx)
    df_6iyz = pd.read_csv(path_6iyz)

    if not df_6iyx['Frame'].equals(df_6iyz['Frame']):
        print(f"Frame mismatch in {file}")
        continue

    frames = df_6iyx['Frame']  # Use actual frame numbers directly

    plt.figure(figsize=(12, 5))

    plt.subplot(1, 2, 1)
    plt.plot(frames, df_6iyx['Distance'], label='5WUC (Closed)', color='blue')
    plt.plot(frames, df_6iyz['Distance'], label='5WUE (Open)', color='red')
    plt.xlabel('Frame')
    plt.ylabel('Distance (Å)')
    plt.title(f'{file[:-19]} Distance')
    plt.grid()
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(frames, df_6iyx['Crossing_Angle_Degrees'], label='5WUC (Closed)', color='blue')
    plt.plot(frames, df_6iyz['Crossing_Angle_Degrees'], label='5WUE (Open)', color='red')
    plt.xlabel('Frame')
    plt.ylabel('Angle (degrees)')
    plt.title(f'{file[:-19]} Crossing Angle')
    plt.grid()
    plt.legend()

    plt.tight_layout()
    out_plot = os.path.join(output_dir, file.replace("_angle_distance.csv", "_comparison.png"))
    plt.savefig(out_plot, dpi=300)
    plt.close()

######################################## Code for figure and comparing two structures #################################################
############################################################################################################################################
############################################################################################################################################


# import MDAnalysis as mda
# import numpy as np
# import matplotlib.pyplot as plt
#
# gro_file2 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1.gro"
# xtc_file2 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1.xtc"
# output_path = "/media/supremeleader/Pantera/simulation/analysis_2025/6iyx/helix_axis/helix1.txt"
#
# x1 = mda.Universe(gro_file2, xtc_file2)
#
# def calculate_helix_axis_and_write(universe, output_file, selection='name CA and resnum 83-100', step=100):
#     with open(output_file, 'w') as f:
#         for ts in universe.trajectory[::step]:
#             atoms = universe.select_atoms(selection)
#             coords = atoms.positions
#
#             # Center the coordinates
#             coords_centered = coords - coords.mean(axis=0)
#
#             # Perform PCA using SVD
#             _, _, vh = np.linalg.svd(coords_centered)
#             axis_vector = vh[0]  # First principal component (unit vector)
#
#             # Write to file
#             f.write(f"{ts.frame}: {axis_vector[0]:.3f},{axis_vector[1]:.3f},{axis_vector[2]:.3f}\n")

# calculate_helix_axis_and_write(x1, output_path)

# (23-39) (51-74) (83-100) (113-135) (144-170) (184-202) (208-232)

############### distance ##########


# import numpy as np
# import matplotlib.pyplot as plt
# # File paths
# file1 = "/media/supremeleader/Pantera/simulation/analysis_2025/6iyx/helix_axis/helix1.txt"
# file2 = "/media/supremeleader/Pantera/simulation/analysis_2025/6iyx/helix_axis/helix3.txt"
# output_file = "/media/supremeleader/Pantera/simulation/analysis_2025/6iyx/helix_axis/helix1-3_new.txt"
#
# # Function to read the vector file into a dictionary
# def read_vector_file(filepath):
#     data = {}
#     with open(filepath, 'r') as f:
#         for line in f:
#             if ":" in line:
#                 timestamp_part, vector_part = line.strip().split(":")
#                 timestamp = int(timestamp_part.strip())
#                 vector = np.array([float(x) for x in vector_part.strip().split(",")])
#                 data[timestamp] = vector
#     return data
#
# # Read vectors from both files
# vectors1 = read_vector_file(file1)
# vectors2 = read_vector_file(file2)
#
# # Find common timestamps
# timestamps = sorted(set(vectors1.keys()) & set(vectors2.keys()))
#
# # Calculate distances and write to output
# with open(output_file, 'w') as out:
#     out.write("# Timestamp : Distance\n")
#     for t in timestamps:
#         v1 = vectors1[t]
#         v2 = vectors2[t]
#         distance = np.linalg.norm(v1 - v2)
#         out.write(f"{t}: {distance:.3f}\n")




# # Read helix axis data
# def read_axis_file(filepath):
#     data = {}
#     with open(filepath, 'r') as f:
#         for line in f:
#             time_part, vector_part = line.strip().split(':')
#             timestamp = int(time_part.strip())
#             vector = np.array([float(x) for x in vector_part.strip().split(',')])
#             data[timestamp] = vector
#     return data
#
# helix1 = read_axis_file(file1)
# helix2 = read_axis_file(file2)
#
# # Calculate distance at each timestamp
# timestamps = sorted(set(helix1.keys()) & set(helix2.keys()))
# distances = []
#
# with open(output_file, 'w') as out:
#     out.write("# Timestamp : Distance (nm)\n")
#     for t in timestamps:
#         v1 = helix1[t]
#         v2 = helix2[t]
#         dist = np.linalg.norm(v1 - v2)  # Units: nanometers (nm)
#         out.write(f"{t}: {dist:.3f} nm\n")
#         distances.append((t, dist))
#
# # Plot distance vs. time
# times, dists = zip(*distances)
# plt.plot(times, dists, marker='o', color='teal')
# plt.xlabel("Time")
# plt.ylabel("Distance between Helix Axes (nm)")
# plt.title("Helix Axis Distance between H1 vs H3")
# plt.grid(True)
# plt.tight_layout()
# plt.show()
# plt.savefig("/media/supremeleader/Pantera/simulation/analysis_2025/6iyz/helix_axis/distances/H6_vs_H7")


########## angle ###########

#
# import numpy as np
# import matplotlib.pyplot as plt
#
# # File paths
# file1 = "/media/supremeleader/Pantera/simulation/analysis_2025/6iyz/helix_axis/helix6_axis.txt"
# file2 = "/media/supremeleader/Pantera/simulation/analysis_2025/6iyz/helix_axis/helix7_axis.txt"
# angle_output_file = "/media/supremeleader/Pantera/simulation/analysis_2025/6iyz/helix_axis/helix-6_helix-7_angle.txt"
#
# # Function to read axis vectors
# def read_axis_file(filepath):
#     data = {}
#     with open(filepath, 'r') as f:
#         for line in f:
#             time_part, vector_part = line.strip().split(':')
#             timestamp = int(time_part.strip())
#             vector = np.array([float(x) for x in vector_part.strip().split(',')])
#             data[timestamp] = vector
#     return data
#
# helix1 = read_axis_file(file1)
# helix2 = read_axis_file(file2)
#
# # Calculate angle at each timestamp
# timestamps = sorted(set(helix1.keys()) & set(helix2.keys()))
# angles = []
#
# with open(angle_output_file, 'w') as out:
#     for t in timestamps:
#         v1 = helix1[t]
#         v2 = helix2[t]
#         # Normalize the vectors
#         v1_norm = v1 / np.linalg.norm(v1)
#         v2_norm = v2 / np.linalg.norm(v2)
#         # Dot product and angle
#         dot = np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0)  # ensure within arccos domain
#         angle_rad = np.arccos(dot)
#         angle_deg = np.degrees(angle_rad)
#         out.write(f"{t}: {angle_deg:.3f}\n")
#         angles.append((t, angle_deg))
#
# # Plotting
# times, angle_values = zip(*angles)
# plt.plot(times, angle_values, marker='o', color='darkorange')
# plt.xlabel("Time")
# plt.ylabel("Angle Between Helix Axes (degrees)")
# plt.title("Angle Between Helix Axes H6 vs H7")
# plt.grid(True)
# plt.tight_layout()
# plt.savefig("/media/supremeleader/Pantera/simulation/analysis_2025/6iyz/helix_axis/angle_fig/H6_vs_H7")



###############################################################################################################



# import MDAnalysis as mda
# import numpy as np
# import matplotlib.pyplot as plt
#
# # File paths
# gro_file2 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1.gro"
# xtc_file2 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1.xtc"
# output_path = "/media/supremeleader/Pantera/simulation/analysis_2025/6iyx/helix_axis/helix3.txt"
#
# # Load trajectory
# x1 = mda.Universe(gro_file2, xtc_file2)
#
# def fit_helix_axis_least_squares(coords):
#     """
#     Given a set of coordinates (Nx3), fit a line using least squares
#     and return the unit direction vector of the best-fit axis.
#     """
#     mean = coords.mean(axis=0)
#     coords_centered = coords - mean
#
#     # Compute covariance matrix
#     cov_matrix = np.dot(coords_centered.T, coords_centered)
#
#     # Get eigenvectors and eigenvalues
#     eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
#
#     # Axis is eigenvector corresponding to largest eigenvalue
#     principal_axis = eigenvectors[:, np.argmax(eigenvalues)]
#
#     # Normalize the axis vector (direction cosines)
#     return principal_axis / np.linalg.norm(principal_axis)
#
# def calculate_helix_axis_and_write(universe, output_file, selection='name CA and resnum 83-100', step=100):
#     with open(output_file, 'w') as f:
#         for ts in universe.trajectory[::step]:
#             atoms = universe.select_atoms(selection)
#             coords = atoms.positions
#
#             # Fit helix axis using least squares method
#             axis_vector = fit_helix_axis_least_squares(coords)
#
#             # Write time and axis (direction cosines)
#             f.write(f"{ts.frame}: {axis_vector[0]:.3f},{axis_vector[1]:.3f},{axis_vector[2]:.3f}\n")
#
# # Run
# calculate_helix_axis_and_write(x1, output_path)

# File paths
import MDAnalysis as mda
# import numpy as np
# import matplotlib.pyplot as plt
# file1 = "/media/supremeleader/Pantera/simulation/analysis_2025/6iyx/helix_axis/helix1.txt"
# file2 = "/media/supremeleader/Pantera/simulation/analysis_2025/6iyx/helix_axis/helix3.txt"
# output_file = "/media/supremeleader/Pantera/simulation/analysis_2025/6iyx/helix_axis/helix1-3.txt"
#
# # Read helix axis data
# def read_axis_file(filepath):
#     data = {}
#     with open(filepath, 'r') as f:
#         for line in f:
#             time_part, vector_part = line.strip().split(':')
#             timestamp = int(time_part.strip())
#             vector = np.array([float(x) for x in vector_part.strip().split(',')])
#             data[timestamp] = vector
#     return data
#
# helix1 = read_axis_file(file1)
# helix2 = read_axis_file(file2)
#
# # Calculate distance at each timestamp
# timestamps = sorted(set(helix1.keys()) & set(helix2.keys()))
# distances = []
#
# with open(output_file, 'w') as out:
#     out.write("# Timestamp : Distance (nm)\n")
#     for t in timestamps:
#         v1 = helix1[t]
#         v2 = helix2[t]
#         dist = np.linalg.norm(v1 - v2)  # Units: nanometers (nm)
#         out.write(f"{t}: {dist:.3f} nm\n")
#         distances.append((t, dist))
#
# # Plot distance vs. time
# times, dists = zip(*distances)
# plt.plot(times, dists, marker='o', color='teal')
# plt.xlabel("Time")
# plt.ylabel("Distance between Helix Axes (nm)")
# plt.title("Helix Axis Distance between H1 vs H3")
# plt.grid(True)
# plt.tight_layout()
# plt.show()

###########################################################################################################


import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
# gro_file2 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1.gro"
# xtc_file2 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1.xtc"
# output_path = "/media/supremeleader/Pantera/simulation/analysis_2025/6iyx/helix_axis/helix3.txt"
#
# def calculate_helix_axis_and_write(universe, output_file, selection='name CA and resnum 83-100', step=100):
#     with open(output_file, 'w') as f:
#         for ts in universe.trajectory[::step]:
#             atoms = universe.select_atoms(selection)
#             coords = atoms.positions
#
#             # Compute centroid of the C-alpha atoms
#             centroid = coords.mean(axis=0)
#
#             # Use least squares line fitting (SVD on centered coords)
#             coords_centered = coords - centroid
#             _, _, vh = np.linalg.svd(coords_centered)
#             axis_vector = vh[0]  # Best-fit line direction (unit vector)
#
#             # Normalize to ensure it's a unit vector
#             axis_vector /= np.linalg.norm(axis_vector)
#
#             # Write direction cosines (unit vector) to file
#             f.write(f"{ts.frame}: {axis_vector[0]:.3f},{axis_vector[1]:.3f},{axis_vector[2]:.3f}\n")
#
# x1 = mda.Universe(gro_file2, xtc_file2)
# calculate_helix_axis_and_write(x1, output_path, selection='name CA and resnum 83-100', step=100)

# file1 = "/media/supremeleader/Pantera/simulation/analysis_2025/6iyx/helix_axis/helix1_axis.txt"
# file2 = "/media/supremeleader/Pantera/simulation/analysis_2025/6iyx/helix_axis/helix3_axis.txt"
# output_file = "/media/supremeleader/Pantera/simulation/analysis_2025/6iyx/helix_axis/helix1-3.txt"
# # Read helix axis data
# def read_axis_file(filepath):
#     data = {}
#     with open(filepath, 'r') as f:
#         for line in f:
#             time_part, vector_part = line.strip().split(':')
#             timestamp = int(time_part.strip())
#             vector = np.array([float(x) for x in vector_part.strip().split(',')])
#             data[timestamp] = vector
#     return data
#
# helix1 = read_axis_file(file1)
# helix2 = read_axis_file(file2)
#
# # Calculate distance at each timestamp
# timestamps = sorted(set(helix1.keys()) & set(helix2.keys()))
# distances = []
#
# with open(output_file, 'w') as out:
#     out.write("# Timestamp : Distance (nm)\n")
#     for t in timestamps:
#         centroid1 = helix1[t]  # 3D position of helix1 centroid at timestamp t
#         centroid2 = helix2[t]  # 3D position of helix2 centroid at timestamp t
#         dist = np.linalg.norm(centroid1 - centroid2)  # Euclidean distance in nanometers
#         out.write(f"{t}: {dist:.3f} nm\n")
#         distances.append((t, dist))
#
#
# # Plot distance vs. time
# times, dists = zip(*distances)
# plt.plot(times, dists, color='teal')
# plt.xlabel("Time")
# plt.ylabel("Distance between Helix Axes (nm)")
# plt.title("Helix Axis Distance between H1 vs H3")
# plt.grid(True)
# plt.tight_layout()
# plt.show()


############################## distance final code ##################

# import MDAnalysis as mda
# import numpy as np
#
# # ----- USER INPUTS -----
# gro_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/chain1.gro"
# xtc_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/chain1.xtc"
# output_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/distance_fig/helix6-7_distance.txt"
#
# helix1_range = (152,170)   # Residue numbers for Helix 1
# helix2_range = (176,197)  # Residue numbers for Helix 2
# # (7,31) (34,57) (67,88) (94,114) (121,143) (152,170) (176,197)
# step = 100  # Analyze every 100th frame
#
# # ----- HELPER FUNCTIONS -----
# def extract_ca_coords(atomgroup):
#     return atomgroup.positions.copy()
#
# def fit_helix_axis(coords):
#     centroid = np.mean(coords, axis=0)
#     centered = coords - centroid
#     _, _, vh = np.linalg.svd(centered)
#     direction = vh[0]
#     return centroid, direction
#
# def shortest_distance_between_lines(p1, v1, p2, v2):
#     cross = np.cross(v1, v2)
#     norm_cross = np.linalg.norm(cross)
#     if norm_cross == 0:
#         return np.linalg.norm(np.cross((p2 - p1), v1)) / np.linalg.norm(v1)
#     distance = np.abs(np.dot((p2 - p1), cross)) / norm_cross
#     return distance
#
# # ----- MAIN ANALYSIS -----
# u = mda.Universe(gro_file, xtc_file)
# resid_offset = u.atoms.residues[0].resid - u.residues[0].resid  # correct indexing
#
# helix1_sel = f"name CA and resid {helix1_range[0] + resid_offset}:{helix1_range[1] + resid_offset}"
# helix2_sel = f"name CA and resid {helix2_range[0] + resid_offset}:{helix2_range[1] + resid_offset}"
#
# helix1_atoms = u.select_atoms(helix1_sel)
# helix2_atoms = u.select_atoms(helix2_sel)
#
# distances = []
#
# with open(output_file, "w") as out:
#     out.write("Frame Distance(Å)\n")
#     for ts in u.trajectory[::step]:
#         coords1 = extract_ca_coords(helix1_atoms)
#         coords2 = extract_ca_coords(helix2_atoms)
#
#         if coords1.shape[0] < 3 or coords2.shape[0] < 3:
#             print(f"Skipping frame {ts.frame} due to insufficient Cα atoms.")
#             continue
#
#         p1, v1 = fit_helix_axis(coords1)
#         p2, v2 = fit_helix_axis(coords2)
#
#         distance = shortest_distance_between_lines(p1, v1, p2, v2)
#         distances.append(distance)
#
#         out.write(f"{ts.frame} {distance:.3f}\n")
#
# # print(f"Distance calculation complete. Results saved to {output_file}")
# import matplotlib.pyplot as plt
#
# # Path to your output file
# output_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/distance_fig/helix6-7_distance.txt"
#
# # Lists to store data
# frames = []
# distances = []
#
# # Read the file and parse data
# with open(output_file, "r") as f:
#     next(f)  # Skip header
#     for line in f:
#         parts = line.strip().split()
#         if len(parts) == 2:
#             frame, distance = int(parts[0]), float(parts[1])
#             frames.append(frame)
#             distances.append(distance)
#
# # Plot
# plt.figure(figsize=(8, 5))
# plt.plot(frames, distances, marker='o', linestyle='-', color='teal')
# plt.title("Shortest Distance Between H6 vs H7", fontsize=14)
# plt.xlabel("Frame", fontsize=12)
# plt.ylabel("Distance (Å)", fontsize=12)
# plt.grid(True)
# plt.tight_layout()
# # plt.show()
# plt.savefig("/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/distance_fig/helix6-7_distance.png")

##################### angle new #############################
# import MDAnalysis as mda
# import numpy as np
# import matplotlib.pyplot as plt
#
# # ----- USER INPUTS -----
# gro_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/chain1.gro"
# xtc_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/chain1.xtc"
# output_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/distance_fig/helix6-7_angle"
#
# helix1_range = (152,170)   # Helix 3
# helix2_range = (176,197)  # Helix 7
# step = 100  # Analyze every 100th frame
#
# # # (7,31) (34,57) (67,88) (94,114) (121,143) (152,170) (176,197)
#
# # ----- HELPER FUNCTIONS -----
# def extract_ca_coords(atomgroup):
#     return atomgroup.positions.copy()
#
# def fit_helix_axis(coords):
#     centroid = np.mean(coords, axis=0)
#     centered = coords - centroid
#     _, _, vh = np.linalg.svd(centered)
#     direction = vh[0]
#     return direction
#
# def angle_between_vectors(v1, v2):
#     dot = np.dot(v1, v2)
#     angle_rad = np.arccos(np.clip(dot / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1.0, 1.0))
#     return np.degrees(angle_rad)
#
# # ----- MAIN ANALYSIS -----
# u = mda.Universe(gro_file, xtc_file)
# resid_offset = u.atoms.residues[0].resid - u.residues[0].resid
#
# helix1_sel = f"name CA and resid {helix1_range[0] + resid_offset}:{helix1_range[1] + resid_offset}"
# helix2_sel = f"name CA and resid {helix2_range[0] + resid_offset}:{helix2_range[1] + resid_offset}"
#
# helix1_atoms = u.select_atoms(helix1_sel)
# helix2_atoms = u.select_atoms(helix2_sel)
#
# frames, angles = [], []
#
# with open(output_file, "w") as out:
#     out.write("Frame Angle(deg)\n")
#     for ts in u.trajectory[::step]:
#         coords1 = extract_ca_coords(helix1_atoms)
#         coords2 = extract_ca_coords(helix2_atoms)
#
#         if coords1.shape[0] < 3 or coords2.shape[0] < 3:
#             print(f"Skipping frame {ts.frame} due to insufficient Cα atoms.")
#             continue
#
#         v1 = fit_helix_axis(coords1)
#         v2 = fit_helix_axis(coords2)
#
#         angle = angle_between_vectors(v1, v2)
#
#         frames.append(ts.frame)
#         angles.append(angle)
#
#         out.write(f"{ts.frame} {angle:.2f}\n")
#
# print(f"Angle calculation complete. Results saved to {output_file}")
#
# # ----- PLOTTING -----
# plt.figure(figsize=(8, 5))
# plt.plot(frames, angles, marker='o', linestyle='--', color='darkred')
# plt.title("Helix 6 vs Helix 7: Angle Over Time", fontsize=14)
# plt.xlabel("Frame", fontsize=12)
# plt.ylabel("Angle (°)", fontsize=12)
# plt.grid(True)
# plt.tight_layout()
# # plt.show()
# plt.savefig("/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/distance_fig/helix6-7_angle_fig.png")







#############################################################
# import matplotlib.pyplot as plt
# file1="/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/distance_angle_fig/helix1-6_5wuc_dis.txt"
# file2= "/media/supremeleader/Pantera/simulation/analysis_2025/6iyx/helix_axis/distances_between_helix/helix1-6_6iyx_dis.txt"
# # Initialize data holders
# frames1, distances1 = [], []
# frames2, distances2 = [], []
# # Read file1
# with open(file1, "r") as f:
#     next(f)  # skip header
#     for line in f:
#         parts = line.strip().split()
#         if len(parts) == 2:
#             frames1.append(int(parts[0]))
#             distances1.append(float(parts[1]))
# # Read file2
# with open(file2, "r") as f:
#     next(f)  # skip header
#     for line in f:
#         parts = line.strip().split()
#         if len(parts) == 2:
#             frames2.append(int(parts[0]))
#             distances2.append(float(parts[1]))
# # Plot
# plt.figure(figsize=(9, 6))
# plt.plot(frames1, distances1, label="5WUC(Prokaryotes)", color='blue', linestyle='-')
# plt.plot(frames2, distances2, label="6IYX(Eukaryotes)", color='red', linestyle='-')
# plt.title("Angle Between TM1 and TM6")
# plt.xlabel("Frame")
# plt.ylabel("Distance (Å)")
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# # Save the figure
# plt.savefig("/media/supremeleader/Pantera/simulation/analysis_2025/comparison/TM_1-6.png")
# # plt.show()  # Uncomment to display

######################################################################################################
######################## correct distance calculation and angle calculation accounting for the handness of the helix ####################

# import MDAnalysis as mda
# import numpy as np
# import matplotlib.pyplot as plt
#
# def compute_principal_axis(coords):
#     """Compute the principal axis (eigenvector with smallest eigenvalue) from Cα coordinates."""
#     n = len(coords)
#     if n < 4:
#         raise ValueError("Not enough atoms to compute virtual p vectors")
#
#     # Construct p vectors as defined in Chothia et al.
#     p = coords[1:-2] + coords[3:] - 2 * coords[:-3]
#     p_mean = np.mean(p, axis=0)
#     p_centered = p - p_mean
#
#     inertia_tensor = np.dot(p_centered.T, p_centered)
#     eigvals, eigvecs = np.linalg.eigh(inertia_tensor)
#     principal_axis = eigvecs[:, np.argmin(eigvals)]
#     return principal_axis / np.linalg.norm(principal_axis)
#
# def get_endpoints(coords, axis):
#     """Return start and end points of a helix projected along its principal axis."""
#     center = coords.mean(axis=0)
#     b = center + np.dot(coords[0] - center, axis) * axis
#     e = center + np.dot(coords[-1] - center, axis) * axis
#     return b, e
#
# def closest_distance_between_axes(b1, e1, b2, e2):
#     """Compute closest points between two helix axes (lines) and their shortest distance."""
#     u = e1 - b1
#     v = e2 - b2
#     w0 = b1 - b2
#     a = np.dot(u, u)
#     b = np.dot(u, v)
#     c = np.dot(v, v)
#     d = np.dot(u, w0)
#     e = np.dot(v, w0)
#
#     denom = a * c - b * b
#     if denom == 0:
#         s, t = 0, 0
#     else:
#         s = (b * e - c * d) / denom
#         t = (a * e - b * d) / denom
#         s = np.clip(s, 0, 1)
#         t = np.clip(t, 0, 1)
#
#     p1 = b1 + s * u
#     p2 = b2 + t * v
#     return np.linalg.norm(p1 - p2), p1, p2
#
# def crossing_angle(b1, p1, p2, b2):
#     """Compute the Chothia-style crossing angle as the dihedral angle between (b1, p1, p2, b2)."""
#     def unit_vector(v):
#         return v / np.linalg.norm(v)
#
#     b1p1 = b1 - p1
#     p1p2 = p2 - p1
#     b2p2 = b2 - p2
#
#     n1 = np.cross(b1p1, p1p2)
#     n2 = np.cross(p1p2, b2p2)
#
#     n1_u = unit_vector(n1)
#     n2_u = unit_vector(n2)
#     p1p2_u = unit_vector(p1p2)
#
#     x = np.dot(n1_u, n2_u)
#     y = np.dot(np.cross(n1_u, n2_u), p1p2_u)
#
#     angle_rad = np.arctan2(y, x)
#     return np.degrees(angle_rad)
#
# # --------- Main Execution ---------
# top = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1.gro"
# traj = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1.xtc"
#
# u = mda.Universe(top, traj)
# helix1 = u.select_atoms("name CA and resid 23-39")
# helix2 = u.select_atoms("name CA and resid 50-74")
#
# distances = []
# angles = []
# times = []
#
# for ts in u.trajectory[::100]:
#     coords1 = helix1.positions
#     coords2 = helix2.positions
#
#     axis1 = compute_principal_axis(coords1)
#     axis2 = compute_principal_axis(coords2)
#
#     b1, e1 = get_endpoints(coords1, axis1)
#     b2, e2 = get_endpoints(coords2, axis2)
#
#     dist, p1, p2 = closest_distance_between_axes(b1, e1, b2, e2)
#     angle = crossing_angle(b1, p1, p2, b2)
#
#     distances.append(dist)
#     angles.append(angle)
#     times.append(ts.time)
#
# # --------- Plotting ---------
# plt.figure(figsize=(12, 5))
#
# plt.subplot(1, 2, 1)
# plt.plot(times, distances, label='Shortest Distance')
# plt.xlabel('Time (ps)')
# plt.ylabel('Distance (Å)')
# plt.title('Helix–Helix Distance')
# plt.grid()
#
# plt.subplot(1, 2, 2)
# plt.plot(times, angles, label='Crossing Angle', color='orange')
# plt.xlabel('Time (ps)')
# plt.ylabel('Angle (degrees)')
# plt.title('Helix–Helix Crossing Angle')
# plt.grid()
#
# plt.tight_layout()
# plt.show()
#
# import pandas as pd
# df = pd.DataFrame({'Time (ps)': times, 'Distance (Å)': distances, 'Crossing Angle (°)': angles})
# df.to_csv("/media/supremeleader/Pantera/simulation/analysis_2025/comparison/helix_helix_distance_angle.csv", index=False)

######################################### new code for positive angle value no fluctuation of distance and angle ########################
# import MDAnalysis as mda
# import numpy as np
# import matplotlib.pyplot as plt
# import pandas as pd
#
# def compute_principal_axis(coords):
#     """Compute the principal axis (eigenvector with smallest eigenvalue) from Cα coordinates."""
#     n = len(coords)
#     if n < 4:
#         raise ValueError("Not enough atoms to compute virtual p vectors")
#
#     # Construct p vectors as defined in Chothia et al.
#     p = coords[1:-2] + coords[3:] - 2 * coords[:-3]
#     p_mean = np.mean(p, axis=0)
#     p_centered = p - p_mean
#
#     inertia_tensor = np.dot(p_centered.T, p_centered)
#     eigvals, eigvecs = np.linalg.eigh(inertia_tensor)
#     principal_axis = eigvecs[:, np.argmin(eigvals)]
#     return principal_axis / np.linalg.norm(principal_axis)
#
# def get_endpoints(coords, axis):
#     """Return start and end points of a helix projected along its principal axis."""
#     center = coords.mean(axis=0)
#     b = center + np.dot(coords[0] - center, axis) * axis
#     e = center + np.dot(coords[-1] - center, axis) * axis
#     return b, e
#
# def closest_distance_between_axes(b1, e1, b2, e2):
#     """Compute closest points between two helix axes (lines) and their shortest distance."""
#     u = e1 - b1
#     v = e2 - b2
#     w0 = b1 - b2
#     a = np.dot(u, u)
#     b = np.dot(u, v)
#     c = np.dot(v, v)
#     d = np.dot(u, w0)
#     e = np.dot(v, w0)
#
#     denom = a * c - b * b
#     if denom == 0:
#         s, t = 0, 0
#     else:
#         s = (b * e - c * d) / denom
#         t = (a * e - b * d) / denom
#         s = np.clip(s, 0, 1)
#         t = np.clip(t, 0, 1)
#
#     p1 = b1 + s * u
#     p2 = b2 + t * v
#     return np.linalg.norm(p1 - p2), p1, p2
#
# def crossing_angle(b1, p1, p2, b2):
#     """Compute the crossing angle as the dihedral angle between (b1, p1, p2, b2)."""
#     def unit_vector(v):
#         return v / np.linalg.norm(v)
#
#     b1p1 = b1 - p1
#     p1p2 = p2 - p1
#     b2p2 = b2 - p2
#
#     n1 = np.cross(b1p1, p1p2)
#     n2 = np.cross(p1p2, b2p2)
#
#     n1_u = unit_vector(n1)
#     n2_u = unit_vector(n2)
#     p1p2_u = unit_vector(p1p2)
#
#     x = np.dot(n1_u, n2_u)
#     y = np.dot(np.cross(n1_u, n2_u), p1p2_u)
#
#     angle_rad = np.arctan2(y, x)
#     return np.degrees(angle_rad)
#
# # --------- Main Execution ---------
# top = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1.gro"
# traj = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1.xtc"
#
# u = mda.Universe(top, traj)
# helix1 = u.select_atoms("name CA and resid 23-39")
# helix2 = u.select_atoms("name CA and resid 50-74")
#
# distances = []
# angles = []
# times = []
#
# for ts in u.trajectory[::100]:
#     coords1 = helix1.positions
#     coords2 = helix2.positions
#
#     axis1 = compute_principal_axis(coords1)
#     axis2 = compute_principal_axis(coords2)
#
#     b1, e1 = get_endpoints(coords1, axis1)
#     b2, e2 = get_endpoints(coords2, axis2)
#
#     dist, p1, p2 = closest_distance_between_axes(b1, e1, b2, e2)
#
#     # Only take the magnitude of the crossing angle
#     angle = abs(crossing_angle(b1, p1, p2, b2))
#
#     distances.append(dist)
#     angles.append(angle)
#     times.append(ts.time)
#
# # --------- Plotting ---------
# plt.figure(figsize=(12, 5))
#
# plt.subplot(1, 2, 1)
# plt.plot(times, distances, label='Shortest Distance')
# plt.xlabel('Time (ps)')
# plt.ylabel('Distance (Å)')
# plt.title('Helix–Helix Distance')
# plt.grid()
#
# plt.subplot(1, 2, 2)
# plt.plot(times, angles, label='Crossing Angle', color='orange')
# plt.xlabel('Time (ps)')
# plt.ylabel('Angle (degrees)')
# plt.title('Helix–Helix Crossing Angle (|°|)')
# plt.grid()
#
# plt.tight_layout()
# plt.show()
#
# # --------- Save CSV ---------
# df = pd.DataFrame({'Time (ps)': times, 'Distance (Å)': distances, 'Crossing Angle (°)': angles})
# df.to_csv("/media/supremeleader/Pantera/simulation/analysis_2025/comparison/helix_helix_distance_angle.csv", index=False)

###################################################################
################################ do not touch this code #############################################
# ##################### code comparing 2 trajectory values of eukaryotes and prokaryotes ###############
# import MDAnalysis as mda
# import numpy as np
# import matplotlib.pyplot as plt
# import pandas as pd
#
# def compute_principal_axis(coords):
#     if len(coords) < 4:
#         raise ValueError("Not enough atoms to compute principal axis.")
#     p = coords[1:-2] + coords[3:] - 2 * coords[:-3]
#     p_mean = np.mean(p, axis=0)
#     p_centered = p - p_mean
#     inertia_tensor = np.dot(p_centered.T, p_centered)
#     eigvals, eigvecs = np.linalg.eigh(inertia_tensor)
#     axis = eigvecs[:, np.argmin(eigvals)]
#     return axis / np.linalg.norm(axis)
#
# def get_endpoints(coords, axis):
#     center = coords.mean(axis=0)
#     b = center + np.dot(coords[0] - center, axis) * axis
#     e = center + np.dot(coords[-1] - center, axis) * axis
#     return b, e
#
# def closest_distance_between_axes(b1, e1, b2, e2):
#     u = e1 - b1
#     v = e2 - b2
#     w0 = b1 - b2
#     a, b, c = np.dot(u, u), np.dot(u, v), np.dot(v, v)
#     d, e = np.dot(u, w0), np.dot(v, w0)
#     denom = a * c - b * b
#     s = (b * e - c * d) / denom if denom != 0 else 0
#     t = (a * e - b * d) / denom if denom != 0 else 0
#     s = np.clip(s, 0, 1)
#     t = np.clip(t, 0, 1)
#     p1 = b1 + s * u
#     p2 = b2 + t * v
#     return np.linalg.norm(p1 - p2), p1, p2
#
# def crossing_angle(b1, p1, p2, b2):
#     def unit(v): return v / np.linalg.norm(v)
#     n1 = np.cross(b1 - p1, p2 - p1)
#     n2 = np.cross(p2 - p1, b2 - p2)
#     n1_u, n2_u, h = unit(n1), unit(n2), unit(p2 - p1)
#     x = np.dot(n1_u, n2_u)
#     y = np.dot(np.cross(n1_u, n2_u), h)
#     return np.degrees(np.arctan2(y, x))
#
# def analyze_trajectory(top, traj, helix1_sel, helix2_sel, label, step=100):
#     u = mda.Universe(top, traj)
#     helix1 = u.select_atoms(helix1_sel)
#     helix2 = u.select_atoms(helix2_sel)
#
#     times, distances, angles = [], [], []
#
#     for ts in u.trajectory[::step]:
#         coords1, coords2 = helix1.positions, helix2.positions
#         axis1, axis2 = compute_principal_axis(coords1), compute_principal_axis(coords2)
#         b1, e1 = get_endpoints(coords1, axis1)
#         b2, e2 = get_endpoints(coords2, axis2)
#         dist, p1, p2 = closest_distance_between_axes(b1, e1, b2, e2)
#         angle = abs(crossing_angle(b1, p1, p2, b2))  # ignore sign
#         times.append(ts.time)
#         distances.append(dist)
#         angles.append(angle)
#
#     df = pd.DataFrame({'Time (ps)': times, 'Distance (Å)': distances, 'Crossing Angle (°)': angles})
#     df.to_csv(f"/media/supremeleader/Pantera/simulation/analysis_2025/comparison/TM1-TM2b_{label}.csv", index=False)
#     return times, distances, angles
#
# # --------- Analyze Both Trajectories ---------
# t1, d1, a1 = analyze_trajectory(
#     top="/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1.gro",
#     traj="/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1.xtc",
#     helix1_sel="name CA and resid 23-39",
#     helix2_sel="name CA and resid 63-73",
#     label="6iyx"
# )
# # (23-39) (50-73) (84-98) (112-138) (143-170) (185-202) (209-229)
# # (50-62, 63-73) (144-155, 156-170)
#
# t2, d2, a2 = analyze_trajectory(
#     top="/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/chain1.gro",
#     traj="/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/chain1.xtc",
#     helix1_sel="name CA and resid 8-29",
#     helix2_sel="name CA and resid 48-57",
#     label="5wuc"
# )
# # (8-29) (35-57) (67-85) (95-114) (121-143) (152-170) (176-197)
# # (35-47, 48-57) (121-133, 134-143)
# # --------- Plotting Comparison ---------
# plt.figure(figsize=(12, 5))
#
# plt.subplot(1, 2, 1)
# plt.plot(t1, d1, label='6iyx (Euk)', color='blue')
# plt.plot(t2, d2, label='5wuc (Pro)', color='red')
# plt.ylabel("Distance (Å)")
# plt.title("TM1–TM2b Distance Comparison")
# plt.legend()
# plt.grid()
#
# plt.subplot(1, 2, 2)
# plt.plot(t1, a1, label='6iyx', color='blue')
# plt.plot(t2, a2, label='5wuc', color='red')
# plt.xlabel("Time (ps)")
# plt.ylabel("Angle (°)")
# plt.title("TM1–TM2b Angle Comparison")
# plt.legend()
# plt.grid()
#
# plt.tight_layout()
# plt.savefig("/media/supremeleader/Pantera/simulation/analysis_2025/comparison/TM1-TM2b.png", dpi=300)
# plt.show()
################################# perfect code ################################################
################################working perfect for most of the helices #########################
############################################## do not touch this code #########################

############# just this code might work ############
########################### working fine and giving zero angles when not able to normalize the zero vector ##########
##################################### dont touch this code as well ###################################################
# (23-39) (50-73) (84-98) (112-138) (143-170) (185-202) (209-229)
# (8-29) (35-57) (67-85) (95-114) (121-143) (152-170) (176-197)

# import MDAnalysis as mda
# import numpy as np
# import matplotlib.pyplot as plt
# import pandas as pd
#
# # ---------- Vector utilities ----------
# def unit_vector(v):
#     norm = np.linalg.norm(v)
#     return v / norm if norm > 1e-8 else np.zeros_like(v)
#
# def compute_principal_axis(coords):
#     """Compute principal axis of a helix using the smallest eigenvector of the virtual p segment cloud."""
#     if len(coords) < 5:
#         raise ValueError("Helix too short to compute principal axis.")
#     # Virtual points
#     p = coords[1:-2] + coords[3:] - 2 * coords[:-3]
#     p_mean = np.mean(p, axis=0)
#     centered = p - p_mean
#     inertia_tensor = np.dot(centered.T, centered)
#     eigvals, eigvecs = np.linalg.eigh(inertia_tensor)
#     axis = eigvecs[:, np.argmin(eigvals)]
#     return unit_vector(axis)
#
# def get_endpoints(coords, axis):
#     """Project first and last CA atoms onto the principal axis to get helix axis endpoints."""
#     center = coords.mean(axis=0)
#     b = center + np.dot(coords[0] - center, axis) * axis
#     e = center + np.dot(coords[-1] - center, axis) * axis
#     return b, e
#
# def closest_distance_between_axes(b1, e1, b2, e2):
#     """Find shortest distance and closest points between two helix axes."""
#     u = e1 - b1
#     v = e2 - b2
#     w0 = b1 - b2
#     a, b, c = np.dot(u, u), np.dot(u, v), np.dot(v, v)
#     d, e = np.dot(u, w0), np.dot(v, w0)
#     denom = a * c - b * b
#     s = (b * e - c * d) / denom if denom != 0 else 0
#     t = (a * e - b * d) / denom if denom != 0 else 0
#     s = np.clip(s, 0, 1)
#     t = np.clip(t, 0, 1)
#     p1 = b1 + s * u
#     p2 = b2 + t * v
#     dist = np.linalg.norm(p1 - p2)
#     return dist, p1, p2
#
# def crossing_angle(b1, p1, p2, b2):
#     """Chothia-style crossing angle: dihedral between b1–p1–p2–b2."""
#     n1 = np.cross(b1 - p1, p2 - p1)
#     n2 = np.cross(p2 - p1, b2 - p2)
#     h = unit_vector(p2 - p1)
#     n1_u = unit_vector(n1)
#     n2_u = unit_vector(n2)
#     x = np.dot(n1_u, n2_u)
#     y = np.dot(np.cross(n1_u, n2_u), h)
#     angle_rad = np.arctan2(y, x)
#     return np.degrees(angle_rad)
#
# # ---------- Analysis wrapper ----------
# def analyze_trajectory(top, traj, helix1_sel, helix2_sel, label, step=100):
#     u = mda.Universe(top, traj)
#     helix1 = u.select_atoms(helix1_sel)
#     helix2 = u.select_atoms(helix2_sel)
#
#     times, distances, angles = [], [], []
#
#     for ts in u.trajectory[::step]:
#         try:
#             coords1 = helix1.positions.copy()
#             coords2 = helix2.positions.copy()
#
#             axis1 = compute_principal_axis(coords1)
#             axis2 = compute_principal_axis(coords2)
#
#             b1, e1 = get_endpoints(coords1, axis1)
#             b2, e2 = get_endpoints(coords2, axis2)
#
#             dist, p1, p2 = closest_distance_between_axes(b1, e1, b2, e2)
#             angle = abs(crossing_angle(b1, p1, p2, b2))  # Use abs for unsigned angle
#
#             times.append(ts.time)
#             distances.append(dist)
#             angles.append(angle)
#         except Exception as e:
#             print(f"[Frame {ts.frame}] Skipped: {e}")
#
#     df = pd.DataFrame({
#         "Time (ps)": times,
#         "Distance (Å)": distances,
#         "Crossing Angle (°)": angles
#     })
#
#     out_csv = f"/media/supremeleader/Pantera/simulation/analysis_2025/comparison/TM1-TM6_{label}.csv"
#     df.to_csv(out_csv, index=False)
#     print(f"Saved: {out_csv}")
#     return times, distances, angles
#
# # ---------- Analyze both systems ----------
# t1, d1, a1 = analyze_trajectory(
#     top="/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1.gro",
#     traj="/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1.xtc",
#     helix1_sel="name CA and resid 23-39",   # TM3 (6iyx)
#     helix2_sel="name CA and resid 185-202",  # TM4 (6iyx)
#     label="6iyx"
# )
# # (23-39) (50-73) (84-98) (112-138) (143-170) (185-202) (209-229)
# t2, d2, a2 = analyze_trajectory(
#     top="/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/chain1.gro",
#     traj="/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/chain1.xtc",
#     helix1_sel="name CA and resid 8-29",    # TM3 (5wuc)
#     helix2_sel="name CA and resid 152-170",   # TM4 (5wuc)
#     label="5wuc"
# )
# # (8-29) (35-57) (67-85) (95-114) (121-143) (152-170) (176-197)
# # ---------- Plotting ----------
# plt.figure(figsize=(12, 5))
#
# # Distance plot
# plt.subplot(1, 2, 1)
# plt.plot(t1, d1, label='6IYX (Euk)', color='blue')
# plt.plot(t2, d2, label='5WUC (Pro)', color='red')
# plt.xlabel("Time (ps)")
# plt.ylabel("Distance (Å)")
# plt.title("TM1–TM6 Distance")
# plt.legend()
# plt.grid()
#
# # Angle plot
# plt.subplot(1, 2, 2)
# plt.plot(t1, a1, label='6IYX (Euk)', color='blue')
# plt.plot(t2, a2, label='5WUC (Pro)', color='red')
# plt.xlabel("Time (ps)")
# plt.ylabel("Crossing Angle (°)")
# plt.title("TM1–TM6 Angle")
# plt.legend()
# plt.grid()
#
# plt.tight_layout()
# plt.savefig("/media/supremeleader/Pantera/simulation/analysis_2025/comparison/TM1-TM6.png", dpi=300)
# plt.show()
###############################################################################################################
########################### working fine and giving zero angles when not able to normalize the zero vector ##########
##################################### dont touch this code as well ###################################################

######################## this section is working very well ##################
###############################################################################
############################################################################
########### helix axis ################################

# import MDAnalysis as mda
# import numpy as np
#
# # Load your structure and trajectory files
# gro_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1.gro"
# xtc_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1.xtc"
# u = mda.Universe(gro_file, xtc_file)
#
# # Select the C-alpha atoms in the helix region
# part1 = u.select_atoms('name CA and resnum 156-170')
# # (23-39) (50-73) (84-98) (112-138) (143-170) (185-202) (209-229)
# # # (50-62, 63-73) (144-155, 156-170)
# # (8-29) (35-57) (67-85) (95-114) (121-143) (152-170) (176-197)
# # # (35-47, 48-57) (121-133, 134-143)
# # Output file to save axis data
# output_file = "/media/supremeleader/Pantera/simulation/analysis_2025/comparison/new_results/h5b_Axis.txt"
# with open(output_file, 'w') as f:
#     f.write("Frame\tCentroid_X\tCentroid_Y\tCentroid_Z\tAxis_X\tAxis_Y\tAxis_Z\n")
#
#     for ts in u.trajectory[::100]:
#         coords = part1.positions
#
#         # Calculate centroid
#         centroid = coords.mean(axis=0)
#
#         # Center coordinates
#         centered = coords - centroid
#
#         # SVD for principal axis (1st principal component)
#         _, _, vh = np.linalg.svd(centered, full_matrices=False)
#         principal_axis = vh[0]  # Direction of maximal spread
#
#         # Write data to file
#         f.write(f"{ts.frame}\t"
#                 f"{centroid[0]:.3f}\t{centroid[1]:.3f}\t{centroid[2]:.3f}\t"
#                 f"{principal_axis[0]:.5f}\t{principal_axis[1]:.5f}\t{principal_axis[2]:.5f}\n")
#
# print(f"Helix axis data written to {output_file}")
#
# ################## crossing angle ############
#
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
#
# def load_helix_axis_file(filepath):
#     """Load helix axis file and extract frame number and axis vector."""
#     df = pd.read_csv(filepath, sep='\t')
#     axis_vectors = df[['Axis_X', 'Axis_Y', 'Axis_Z']].values
#     frames = df['Frame'].values
#     return frames, axis_vectors
#
# def angle_between_vectors(v1, v2):
#     """Compute angle between two vectors in degrees using dot product."""
#     unit_v1 = v1 / np.linalg.norm(v1)
#     unit_v2 = v2 / np.linalg.norm(v2)
#     dot_product = np.dot(unit_v1, unit_v2)
#     angle_rad = np.arccos(np.clip(dot_product, -1.0, 1.0))  # Numerical stability
#     return np.degrees(angle_rad)
#
# def quadrant_angle(angle):
#     """Mimic Chimera's quadrant_angle logic."""
#     while angle < 0.0:
#         angle += 360.0
#     while angle > 360.0:
#         angle -= 360.0
#     if angle > 180.0:
#         angle = 360.0 - angle
#     if angle > 90.0:
#         angle = 180.0 - angle
#     return angle
#
# def plot_crossing_angles(frames, angles, output_image):
#     """Plot and save the crossing angle vs frame."""
#     plt.figure(figsize=(8, 5))
#     plt.plot(frames, angles, marker='o', linestyle='-', color='darkblue', label='Crossing Angle')
#     plt.xlabel("Frame")
#     plt.ylabel("Crossing Angle (degrees)")
#     plt.title("Helix–Helix Crossing Angle vs Frame")
#     plt.grid(True)
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig(output_image, dpi=300)
#     plt.close()
#
# # === File paths ===
# helix1_file = "/media/supremeleader/Pantera/simulation/analysis_2025/comparison/new_results/h5a_Axis.txt"
# helix2_file = "/media/supremeleader/Pantera/simulation/analysis_2025/comparison/new_results/h5b_Axis.txt"
# output_file = "/media/supremeleader/Pantera/simulation/analysis_2025/comparison/new_results/h5a-h5b_crossing_angle.txt"
# output_plot = "/media/supremeleader/Pantera/simulation/analysis_2025/comparison/new_results/h5a-h5b_crossing_angle.png"
#
# # === Load data ===
# frames1, axes1 = load_helix_axis_file(helix1_file)
# frames2, axes2 = load_helix_axis_file(helix2_file)
#
# # === Sanity check ===
# if not np.array_equal(frames1, frames2):
#     raise ValueError("Frame numbers do not match between helix axis files.")
#
# # === Calculate crossing angles, write to file and store for plot ===
# frames = []
# angles = []
#
# with open(output_file, 'w') as f:
#     f.write("Frame\tCrossing_Angle_Degrees\n")
#     for frame, axis1, axis2 in zip(frames1, axes1, axes2):
#         raw_angle = angle_between_vectors(axis1, axis2)
#         adjusted_angle = quadrant_angle(raw_angle)
#         f.write(f"{frame}\t{adjusted_angle:.2f}\n")
#         frames.append(frame)
#         angles.append(adjusted_angle)
#
# print(f"Crossing angles written to: {output_file}")
#
# # === Plotting ===
# plot_crossing_angles(frames, angles, output_plot)
# print(f"Crossing angle plot saved to: {output_plot}")
#
# ############## interhelical distance ####################
#
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
#
# # === Geometry utilities ===
# def cross_product(a, b):
#     return np.cross(a, b)
#
# def point_to_line_distance(p, a, d):
#     """Distance from point p to line (a + td)."""
#     return np.linalg.norm(np.cross(p - a, d) / np.linalg.norm(d))
#
# def shortest_distance_between_lines(p1, d1, p2, d2):
#     """Compute shortest distance between two skew lines."""
#     d1_u = d1 / np.linalg.norm(d1)
#     d2_u = d2 / np.linalg.norm(d2)
#     cross = np.cross(d1_u, d2_u)
#     cross_norm = np.linalg.norm(cross)
#     if np.isclose(cross_norm, 0.0):
#         return point_to_line_distance(p2, p1, d1_u)
#     diff = p2 - p1
#     return np.abs(np.dot(diff, cross)) / cross_norm
#
# # === Load helix axis data ===
# def load_axis_file(filepath):
#     df = pd.read_csv(filepath, sep='\t')
#     centroids = df[['Centroid_X', 'Centroid_Y', 'Centroid_Z']].values
#     axes = df[['Axis_X', 'Axis_Y', 'Axis_Z']].values
#     frames = df['Frame'].values
#     return frames, centroids, axes
#
# # === File paths ===
# helix1_file = "/media/supremeleader/Pantera/simulation/analysis_2025/comparison/new_results/h1_Axis.txt"
# helix2_file = "/media/supremeleader/Pantera/simulation/analysis_2025/comparison/new_results/h2_Axis.txt"
# output_txt = "/media/supremeleader/Pantera/simulation/analysis_2025/comparison/new_results/h1-h2_distance.txt"
# output_plot = "/media/supremeleader/Pantera/simulation/analysis_2025/comparison/new_results/h1-h2_distance.png"
#
# # === Load data ===
# frames1, centroids1, axes1 = load_axis_file(helix1_file)
# frames2, centroids2, axes2 = load_axis_file(helix2_file)
#
# # === Sanity check ===
# if not np.array_equal(frames1, frames2):
#     raise ValueError("Frame numbers do not match between helix axis files.")
#
# # === Compute distances ===
# frames = []
# distances = []
#
# with open(output_txt, 'w') as f:
#     f.write("Frame\tDistance\n")
#     for frame, c1, a1, c2, a2 in zip(frames1, centroids1, axes1, centroids2, axes2):
#         d = shortest_distance_between_lines(c1, a1, c2, a2)
#         frames.append(frame)
#         distances.append(d)
#         f.write(f"{frame}\t{d:.3f}\n")
#
# print(f"Interhelical distances saved to: {output_txt}")
#
# # === Plot ===
# plt.figure(figsize=(8, 5))
# plt.plot(frames, distances, color='darkred', marker='o', linestyle='-', label='Interhelical Distance')
# plt.xlabel("Frame")
# plt.ylabel("Distance (Å)")
# plt.title("Interhelical Distance vs Frame")
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.savefig(output_plot, dpi=300)
# plt.close()
#
# print(f"Plot saved to: {output_plot}")

######################## this section is working very well ##################
###############################################################################
############################################################################

################### combined code of all helix axis/ crossing angle / crossing distance#############################
######################## this section is working very well ##################
###############################################################################
############################################################################
# import os
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import MDAnalysis as mda
#
# # === Setup ===
# gro_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1_new.gro"
# xtc_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1_new.xtc"
# output_dir = "/media/supremeleader/Pantera/simulation/analysis_2025/comparison/new_results"
# os.makedirs(output_dir, exist_ok=True)
#
# # === Helix definitions ===
# helices = {
#     "H1": (23, 39), "H2a": (50, 62), "H2b": (63, 73), "H3": (84, 98),
#     "H4": (112, 138), "H5a": (144, 155), "H5b": (156, 170),
#     "H6": (185, 202), "H7": (209, 229)
# }
#
# # === Load universe ===
# u = mda.Universe(gro_file, xtc_file)
#
# # === Step 1: Calculate helix axis data ===
# axis_data = {name: [] for name in helices}
# frames = []
#
# for ts in u.trajectory[::100]:
#     frames.append(ts.frame)
#     for name, (start, end) in helices.items():
#         sel = u.select_atoms(f"name CA and resid {start}-{end}")
#         coords = sel.positions
#         centroid = coords.mean(axis=0)
#         centered = coords - centroid
#         _, _, vh = np.linalg.svd(centered, full_matrices=False)
#         axis = vh[0]
#         axis_data[name].append((centroid, axis))
#
# # === Save helix axis data to individual files ===
# for name, data in axis_data.items():
#     out_path = os.path.join(output_dir, f"{name}_Axis.txt")
#     with open(out_path, 'w') as f:
#         f.write("Frame\tCentroid_X\tCentroid_Y\tCentroid_Z\tAxis_X\tAxis_Y\tAxis_Z\n")
#         for frame, (centroid, axis) in zip(frames, data):
#             f.write(f"{frame}\t{centroid[0]:.3f}\t{centroid[1]:.3f}\t{centroid[2]:.3f}\t"
#                     f"{axis[0]:.5f}\t{axis[1]:.5f}\t{axis[2]:.5f}\n")
#
# # === Utilities for angle and distance ===
# def angle_between_vectors(v1, v2):
#     v1_u = v1 / np.linalg.norm(v1)
#     v2_u = v2 / np.linalg.norm(v2)
#     dot = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
#     return np.degrees(np.arccos(dot))
#
# def quadrant_angle(angle):
#     angle = angle % 360.0
#     if angle > 180.0:
#         angle = 360.0 - angle
#     if angle > 90.0:
#         angle = 180.0 - angle
#     return angle
#
# def shortest_distance_between_lines(p1, d1, p2, d2):
#     d1_u = d1 / np.linalg.norm(d1)
#     d2_u = d2 / np.linalg.norm(d2)
#     cross = np.cross(d1_u, d2_u)
#     cross_norm = np.linalg.norm(cross)
#     if np.isclose(cross_norm, 0.0):
#         return np.linalg.norm(np.cross(p2 - p1, d1_u)) / np.linalg.norm(d1_u)
#     diff = p2 - p1
#     return np.abs(np.dot(diff, cross)) / cross_norm
#
# # === Step 2 & 3: Calculate angle and distance for each helix pair ===
# pairs = [(h1, h2) for i, h1 in enumerate(helices) for h2 in list(helices)[i+1:]]
#
# for h1, h2 in pairs:
#     out_combined_file = os.path.join(output_dir, f"{h1}-{h2}_angle_distance.csv")
#
#     with open(out_combined_file, 'w') as f:
#         f.write("Frame,Crossing_Angle_Degrees,Distance\n")
#         for frame, (cent1, axis1), (cent2, axis2) in zip(frames, axis_data[h1], axis_data[h2]):
#             angle = angle_between_vectors(axis1, axis2)
#             angle = quadrant_angle(angle)
#             dist = shortest_distance_between_lines(cent1, axis1, cent2, axis2)
#             f.write(f"{frame},{angle:.2f},{dist:.3f}\n")
#
# print("All helix analysis (axis, angle, distance) completed.")




################################## visualization ############################################

# import numpy as np
# import matplotlib.pyplot as plt
# import MDAnalysis as mda
#
# # File paths
# pdb1 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/112ns-5a.pdb"
# pdb2 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/112ns-5b.pdb"
#
# # Load Cα atoms
# u1 = mda.Universe(pdb1)
# u2 = mda.Universe(pdb2)
# ca1 = u1.select_atoms("name CA").positions
# ca2 = u2.select_atoms("name CA").positions
#
# # === Utilities ===
#
# def unit_vector(v):
#     norm = np.linalg.norm(v)
#     return v / norm if norm > 1e-8 else np.zeros_like(v)
#
# def compute_axis_svd(coords):
#     centroid = coords.mean(axis=0)
#     centered = coords - centroid
#     _, _, vh = np.linalg.svd(centered, full_matrices=False)
#     axis = vh[0]
#     return centroid, unit_vector(axis)
#
# def get_endpoints(coords, axis):
#     center = coords.mean(axis=0)
#     b = center + np.dot(coords[0] - center, axis) * axis
#     e = center + np.dot(coords[-1] - center, axis) * axis
#     return b, e
#
# def closest_distance_between_axes(b1, e1, b2, e2):
#     u = e1 - b1
#     v = e2 - b2
#     w0 = b1 - b2
#     a, b_, c = np.dot(u, u), np.dot(u, v), np.dot(v, v)
#     d, e_ = np.dot(u, w0), np.dot(v, w0)
#     denom = a * c - b_ * b_
#     s = (b_ * e_ - c * d) / denom if denom != 0 else 0
#     t = (a * e_ - b_ * d) / denom if denom != 0 else 0
#     s = np.clip(s, 0, 1)
#     t = np.clip(t, 0, 1)
#     p1 = b1 + s * u
#     p2 = b2 + t * v
#     dist = np.linalg.norm(p1 - p2)
#     return dist, p1, p2
#
# def angle_between_vectors(v1, v2):
#     v1_u = unit_vector(v1)
#     v2_u = unit_vector(v2)
#     dot = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
#     angle_deg = np.degrees(np.arccos(dot))
#     return 180 - angle_deg if angle_deg > 90 else angle_deg
#
# # === Compute Axis and Geometry ===
#
# centroid1, axis1 = compute_axis_svd(ca1)
# centroid2, axis2 = compute_axis_svd(ca2)
#
# start1, end1 = get_endpoints(ca1, axis1)
# start2, end2 = get_endpoints(ca2, axis2)
#
# dist, point1, point2 = closest_distance_between_axes(start1, end1, start2, end2)
# angle = angle_between_vectors(axis1, axis2)
#
# # === Plotting ===
#
# fig = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot(111, projection='3d')
#
# # Helix backbones
# ax.plot(ca1[:, 0], ca1[:, 1], ca1[:, 2], color='skyblue', linewidth=2, label="Helix 1 Backbone")
# ax.plot(ca2[:, 0], ca2[:, 1], ca2[:, 2], color='lightcoral', linewidth=2, label="Helix 2 Backbone")
#
# # Principal axes
# ax.plot([start1[0], end1[0]], [start1[1], end1[1]], [start1[2], end1[2]],
#         color='blue', lw=3, label='Helix 1 Axis')
# ax.plot([start2[0], end2[0]], [start2[1], end2[1]], [start2[2], end2[2]],
#         color='red', lw=3, label='Helix 2 Axis')
#
# # Shortest distance line
# ax.plot([point1[0], point2[0]], [point1[1], point2[1]], [point1[2], point2[2]],
#         color='green', lw=3, label='Shortest Distance')
#
# # Distance label
# midpoint = (point1 + point2) / 2
# offset = np.array([2.0, 2.0, 2.0])
# ax.text(*(midpoint + offset), f"{dist:.2f} Å", color='green', fontsize=12, weight='bold')
#
# # Angle label
# angle_pos = (centroid1 + centroid2) / 2 + np.array([0, 0, 5.0])
# ax.text(*angle_pos, f"Angle: {angle:.1f}°", color='purple', fontsize=12, weight='bold')
#
# # Axes labels and legend
# ax.set_xlabel("X (Å)")
# ax.set_ylabel("Y (Å)")
# ax.set_zlabel("Z (Å)")
# ax.set_title("Helices with Principal Axes, Distance and Crossing Angle")
# ax.legend()
# plt.tight_layout()
# plt.show()

























###################################################################################################

###################################################################################################
###################### This is too perfect method to calculate the distance and angle ###########################
###################################################################################################

# import numpy as np
# import matplotlib.pyplot as plt
# import MDAnalysis as mda
#
# # File paths for helices
# pdb1 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/112ns-5a.pdb"
# pdb2 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/112ns-5b.pdb"
#
# # Load Cα atoms
# u1 = mda.Universe(pdb1)
# u2 = mda.Universe(pdb2)
# ca1 = u1.select_atoms("name CA").positions
# ca2 = u2.select_atoms("name CA").positions
#
# # === Consistent with trajectory code ===
#
# def unit_vector(v):
#     norm = np.linalg.norm(v)
#     return v / norm if norm > 1e-8 else np.zeros_like(v)
#
# def compute_axis_svd(coords):
#     centroid = coords.mean(axis=0)
#     centered = coords - centroid
#     _, _, vh = np.linalg.svd(centered, full_matrices=False)
#     axis = vh[0]  # same as used in trajectory analysis
#     return centroid, unit_vector(axis)
#
# def get_endpoints(coords, axis):
#     center = coords.mean(axis=0)
#     b = center + np.dot(coords[0] - center, axis) * axis
#     e = center + np.dot(coords[-1] - center, axis) * axis
#     return b, e
#
# def closest_distance_between_axes(b1, e1, b2, e2):
#     u = e1 - b1
#     v = e2 - b2
#     w0 = b1 - b2
#     a, b_, c = np.dot(u, u), np.dot(u, v), np.dot(v, v)
#     d, e_ = np.dot(u, w0), np.dot(v, w0)
#     denom = a * c - b_ * b_
#     s = (b_ * e_ - c * d) / denom if denom != 0 else 0
#     t = (a * e_ - b_ * d) / denom if denom != 0 else 0
#     s = np.clip(s, 0, 1)
#     t = np.clip(t, 0, 1)
#     p1 = b1 + s * u
#     p2 = b2 + t * v
#     dist = np.linalg.norm(p1 - p2)
#     return dist, p1, p2
#
# def angle_between_vectors(v1, v2):
#     v1_u = unit_vector(v1)
#     v2_u = unit_vector(v2)
#     dot = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
#     angle_deg = np.degrees(np.arccos(dot))
#     return 180 - angle_deg if angle_deg > 90 else angle_deg
#
# # === Compute axis and distance ===
#
# centroid1, axis1 = compute_axis_svd(ca1)
# centroid2, axis2 = compute_axis_svd(ca2)
#
# start1, end1 = get_endpoints(ca1, axis1)
# start2, end2 = get_endpoints(ca2, axis2)
#
# dist, point1, point2 = closest_distance_between_axes(start1, end1, start2, end2)
# angle = angle_between_vectors(axis1, axis2)
#
# # === Plot ===
#
# fig = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot(111, projection='3d')
#
# # Plot helix backbones
# ax.plot(ca1[:, 0], ca1[:, 1], ca1[:, 2], color='skyblue', linewidth=2, label="Helix 1 Backbone")
# ax.plot(ca2[:, 0], ca2[:, 1], ca2[:, 2], color='lightcoral', linewidth=2, label="Helix 2 Backbone")
#
# # Plot principal axes
# ax.plot([start1[0], end1[0]], [start1[1], end1[1]], [start1[2], end1[2]], color='blue', lw=3, label='Helix 1 Axis')
# ax.plot([start2[0], end2[0]], [start2[1], end2[1]], [start2[2], end2[2]], color='red', lw=3, label='Helix 2 Axis')
#
# # Plot shortest distance line
# ax.plot([point1[0], point2[0]], [point1[1], point2[1]], [point1[2], point2[2]], color='green', lw=3, label='Shortest Distance')
#
# # Annotate shortest distance
# midpoint = (point1 + point2) / 2
# offset = np.array([2.0, 2.0, 2.0])
# ax.text(*(midpoint + offset), f"{dist:.2f} Å", color='green', fontsize=12, weight='bold')
#
# # Annotate angle
# angle_pos = (centroid1 + centroid2) / 2 + np.array([0, 0, 5.0])
# ax.text(*angle_pos, f"Angle: {angle:.1f}°", color='purple', fontsize=12, weight='bold')
#
# # Axes
# ax.set_xlabel("X (Å)")
# ax.set_ylabel("Y (Å)")
# ax.set_zlabel("Z (Å)")
# ax.set_title("Helices with Principal Axes, Distance and Crossing Angle")
# ax.legend()
# plt.tight_layout()
# plt.show()

###################################################################################################
###################### This is too perfect method to calculate the distance and angle ###########################
###################################################################################################

############### use this code for the 6IYZ distance calculation between helices #####################

# import numpy as np
# import matplotlib.pyplot as plt
# import MDAnalysis as mda
#
# # File paths
# pdb1 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/60ns-5a.pdb"
# pdb2 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/60ns-5b.pdb"
#
# # Load and select CA atoms
# u1 = mda.Universe(pdb1)
# u2 = mda.Universe(pdb2)
# ca1 = u1.select_atoms("name CA").positions
# ca2 = u2.select_atoms("name CA").positions
#
# # === Utility functions ===
#
# def unit_vector(v):
#     norm = np.linalg.norm(v)
#     return v / norm if norm > 1e-8 else np.zeros_like(v)
#
# def compute_principal_axis(coords):
#     if len(coords) < 5:
#         raise ValueError("Helix too short for principal axis.")
#     p = coords[1:-2] + coords[3:] - 2 * coords[:-3]
#     p_mean = p.mean(axis=0)
#     centered = p - p_mean
#     inertia_tensor = np.dot(centered.T, centered)
#     eigvals, eigvecs = np.linalg.eigh(inertia_tensor)
#     axis = eigvecs[:, np.argmin(eigvals)]
#     return unit_vector(axis)
#
# def get_endpoints(coords, axis):
#     center = coords.mean(axis=0)
#     b = center + np.dot(coords[0] - center, axis) * axis
#     e = center + np.dot(coords[-1] - center, axis) * axis
#     return b, e
#
# def closest_distance_between_axes(b1, e1, b2, e2):
#     u = e1 - b1
#     v = e2 - b2
#     w0 = b1 - b2
#     a, b, c = np.dot(u, u), np.dot(u, v), np.dot(v, v)
#     d, e = np.dot(u, w0), np.dot(v, w0)
#     denom = a * c - b * b
#     s = (b * e - c * d) / denom if denom != 0 else 0
#     t = (a * e - b * d) / denom if denom != 0 else 0
#     s = np.clip(s, 0, 1)
#     t = np.clip(t, 0, 1)
#     p1 = b1 + s * u
#     p2 = b2 + t * v
#     dist = np.linalg.norm(p1 - p2)
#     return dist, p1, p2
#
# def angle_between_vectors(v1, v2):
#     v1_u = unit_vector(v1)
#     v2_u = unit_vector(v2)
#     dot = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
#     angle_deg = np.degrees(np.arccos(dot))
#     # Convert to acute angle
#     if angle_deg > 90:
#         angle_deg = 180 - angle_deg
#     return angle_deg
#
# # === Calculations ===
#
# axis1 = compute_principal_axis(ca1)
# axis2 = compute_principal_axis(ca2)
#
# start1, end1 = get_endpoints(ca1, axis1)
# start2, end2 = get_endpoints(ca2, axis2)
#
# dist, point1, point2 = closest_distance_between_axes(start1, end1, start2, end2)
# angle = angle_between_vectors(axis1, axis2)
#
# # === Plotting ===
# fig = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot(111, projection='3d')
#
# # Plot Cα backbone
# ax.plot(ca1[:, 0], ca1[:, 1], ca1[:, 2], color='skyblue', label='TM5A Backbone', linewidth=2)
# ax.plot(ca2[:, 0], ca2[:, 1], ca2[:, 2], color='lightcoral', label='TM5B Backbone', linewidth=2)
#
# # Plot principal axes
# ax.plot([start1[0], end1[0]], [start1[1], end1[1]], [start1[2], end1[2]], color='blue', lw=3, label='Helix 1 Axis')
# ax.plot([start2[0], end2[0]], [start2[1], end2[1]], [start2[2], end2[2]], color='red', lw=3, label='Helix 2 Axis')
#
# # Plot shortest distance line
# ax.plot([point1[0], point2[0]], [point1[1], point2[1]], [point1[2], point2[2]],
#         color='green', lw=3, label='Shortest Distance')
#
# # Annotate distance
# midpoint = (point1 + point2) / 2
# offset = np.array([2.0, 2.0, 2.0])
# ax.text(*(midpoint + offset), f"{dist:.2f} Å", color='green', fontsize=12, weight='bold')
#
# # Annotate angle
# # center1 = ca1.mean(axis=0)
# # center2 = ca2.mean(axis=0)
# # angle_pos = (center1 + center2) / 2 + np.array([0, 0, 5.0])
# # ax.text(*angle_pos, f"Angle: {angle:.1f}°", color='purple', fontsize=12, weight='bold')
#
# # Axes and legend
# ax.set_xlabel("X (Å)")
# ax.set_ylabel("Y (Å)")
# ax.set_zlabel("Z (Å)")
# ax.set_title("Helices with Principal Axes, Distance and Crossing Angle")
# ax.legend()
# plt.tight_layout()
# plt.show()























# import MDAnalysis as mda
# import numpy as np
# import matplotlib.pyplot as plt
# gro_file2 = "/media/supremeleader/Pantera/simulation/analysis_2025/6iyx/helix1.gro"
# xtc_file2 = "/media/supremeleader/Pantera/simulation/analysis_2025/6iyx/helix1.xtc"
#
# x1= mda.Universe(gro_file2, xtc_file2)
# helix_2={}
# for ts in x1.trajectory[::100]:
#     part1 = x1.select_atoms('name CA and resnum 23-39')
#     helix_2[ts.frame] = []
#     for atom in part1:
#         name = atom.name
#         resid = atom.resid
#         resname = atom.resname
#         x, y, z = atom.position
#         helix_2[ts.frame].append([resid,resname,x ,y, z])
# print(helix_2)

# (23-39) (51-74) (78-100) (113-135) (144-170) (184-202) (208-232)



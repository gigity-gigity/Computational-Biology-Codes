##############################################################################################################
##############################################################################################################
                                    # Code for average water molecule per bin #
##############################################################################################################
##############################################################################################################

import MDAnalysis as mda
import numpy as np

gro_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/pani.gro"
xtc_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/pani.xtc"
centroid_files = {
    "monomer1": "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/monomer1_centroid.xvg",
    "monomer2": "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/monomer2_centroid.xvg",
    "monomer3": "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/monomer3_centroid.xvg"
}

start_z_coordinate = 62.00
end_z_coordinate = 27.00
u = mda.Universe(gro_file, xtc_file)

# Dictionary to hold atom data for each frame
frame_data = {}
for ts in u.trajectory[::100]:
    atoms_in_range = u.select_atoms(f"prop z >= {end_z_coordinate} and prop z <= {start_z_coordinate}")
    frame_data[ts.frame * 10] = []
    for atom in atoms_in_range:
        name = atom.name
        resid = atom.resid
        x, y, z = atom.position
        frame_data[ts.frame * 10].append([resid, name, x, y, z])

# Load centroid data for each monomer
centroids = {}
for monomer, file_path in centroid_files.items():
    centroids[monomer] = {}
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith(("#", "@")):
                line = line.strip().split()
                frame = int(line[0])
                centroids[monomer][frame] = [float(line[1]) * 10, float(line[2]) * 10, float(line[3]) * 10]

# Function to check if atom is within cylinder
def central_cylinder(atom_coords, center_coords, radius, height):
    x_atom, y_atom, z_atom = atom_coords
    x_center, y_center, z_center = center_coords
    distance_xy = np.sqrt((x_atom - x_center) ** 2 + (y_atom - y_center) ** 2)
    within_radius = distance_xy <= radius
    within_height = (z_center - height / 2) <= z_atom <= (z_center + height / 2)
    return within_radius and within_height

# Analyze atoms inside the cylinder for each monomer
inside_cylinder = {monomer: {} for monomer in centroid_files}
radius = 5
height = 5

for monomer, centroid_data in centroids.items():
    for frame, center_coords in centroid_data.items():
        inside_cylinder[monomer][frame] = []
        if frame in frame_data:
            for atom_data in frame_data[frame]:
                resid, name, x, y, z = atom_data
                atom_coords = (x, y, z)
                if central_cylinder(atom_coords, center_coords, radius, height):
                    inside_cylinder[monomer][frame].append(atom_data)

# Save atom data to file
with open("/media/supremeleader/Pantera/simulation/analysis_2024/average_water_molecule/central_cylinder_atoms.txt", "w") as atom_file:
    for monomer, frames in inside_cylinder.items():
        atom_file.write(f"Monomer: {monomer}\n")
        for frame, atoms in frames.items():
            atom_file.write(f"Frame {frame}:\n")
            for atom_data in atoms:
                resid, name, x, y, z = atom_data
                atom_file.write(f"  ResID {resid} Atom {name} at ({x:.3f}, {y:.3f}, {z:.3f})\n")
        atom_file.write("\n")

# Calculate averages, standard deviations, and save results
average_results = {}
for monomer, frames in inside_cylinder.items():
    atom_counts = [len(atoms) for atoms in frames.values()]
    avg_atoms = np.mean(atom_counts)
    std_dev_atoms = np.std(atom_counts)
    average_results[monomer] = (avg_atoms, std_dev_atoms)
    # Print average and standard deviation for each monomer
    print(f"Average number of atoms per frame for {monomer}: {avg_atoms:.3f} (±{std_dev_atoms:.3f})")

# Calculate combined average and standard deviation across all monomers
all_atom_counts = [count for monomer, frames in inside_cylinder.items() for count in [len(atoms) for atoms in frames.values()]]
overall_avg = np.mean(all_atom_counts)
overall_std_dev = np.std(all_atom_counts)

# Print combined average and standard deviation for all monomers
print(f"Overall average number of atoms per frame for all monomers: {overall_avg:.3f} (±{overall_std_dev:.3f})")

# Save average and standard deviation data to file
with open("/media/supremeleader/Pantera/simulation/analysis_2024/average_water_molecule/central_cylinder_average.txt", "w") as avg_file:
    for monomer, (avg_atoms, std_dev_atoms) in average_results.items():
        avg_file.write(f"Average number of atoms per frame for {monomer}: {avg_atoms:.3f} (±{std_dev_atoms:.3f})\n")
    avg_file.write(f"Overall average number of atoms per frame for all monomers: {overall_avg:.3f} (±{overall_std_dev:.3f})\n")

#############################################################################################################
########################################### cylinder above and below ########################################
#############################################################################################################
# Parameters for above and below cylinders
radius = 5
height = 5
cylinder_count = 3  # Number of cylinders above and below

# Define function to check if atom is within a given cylinder
def check_cylinder(atom_coords, center_coords, radius, height):
    x_atom, y_atom, z_atom = atom_coords
    x_center, y_center, z_center = center_coords
    distance_xy = np.sqrt((x_atom - x_center) ** 2 + (y_atom - y_center) ** 2)
    within_radius = distance_xy <= radius
    within_height = (z_center - height / 2) <= z_atom <= (z_center + height / 2)
    return within_radius and within_height

# Initialize dictionaries to store results for above and below cylinders, with data separated by monomer
above_cylinders = {monomer: {f"above_cylinder{i+1}": {} for i in range(cylinder_count)} for monomer in centroid_files}
below_cylinders = {monomer: {f"below_cylinder{i+1}": {} for i in range(cylinder_count)} for monomer in centroid_files}

# Iterate through monomers and frames with centroids
for monomer, centroid_data in centroids.items():
    for frame, center_coords in centroid_data.items():
        z_center = center_coords[2]

        # Define cylinder positions above and below the central cylinder
        above_centers = [(center_coords[0], center_coords[1], z_center + 2.5 + i * height) for i in range(cylinder_count)]
        below_centers = [(center_coords[0], center_coords[1], z_center - 2.5 - i * height) for i in range(cylinder_count)]

        # Initialize lists to hold atoms inside each cylinder
        for i in range(cylinder_count):
            above_cylinders[monomer][f"above_cylinder{i+1}"].setdefault(frame, [])
            below_cylinders[monomer][f"below_cylinder{i+1}"].setdefault(frame, [])

        # Check each atom in frame_data to see if it is within any of the cylinders
        if frame in frame_data:
            for atom_data in frame_data[frame]:
                resid, name, x, y, z = atom_data
                atom_coords = (x, y, z)

                # Check above cylinders
                for i, center in enumerate(above_centers):
                    if check_cylinder(atom_coords, center, radius, height):
                        above_cylinders[monomer][f"above_cylinder{i+1}"][frame].append(atom_data)

                # Check below cylinders
                for i, center in enumerate(below_centers):
                    if check_cylinder(atom_coords, center, radius, height):
                        below_cylinders[monomer][f"below_cylinder{i+1}"][frame].append(atom_data)

# Function to save data and calculate statistics for cylinders, separated by monomer and cylinder
def save_cylinder_data(cylinder_data, filename_prefix):
    overall_avg_results = {}
    overall_avg_per_cylinder = {f"cylinder{i+1}": [] for i in range(cylinder_count)}

    for monomer, cylinders in cylinder_data.items():
        for cylinder, frames in cylinders.items():
            atom_counts = [len(atoms) for atoms in frames.values()]
            avg_atoms = np.mean(atom_counts)
            std_dev_atoms = np.std(atom_counts)
            overall_avg_results[(monomer, cylinder)] = (avg_atoms, std_dev_atoms)

            # Add to overall cylinder average across all monomers
            cylinder_index = int(cylinder.split("cylinder")[-1]) - 1
            overall_avg_per_cylinder[f"cylinder{cylinder_index+1}"].extend(atom_counts)

            # Save atoms to file
            with open(f"{filename_prefix}_{monomer}_{cylinder}.txt", "w") as file:
                file.write(f"Monomer: {monomer}\n")
                file.write(f"{cylinder}:\n")
                for frame, atoms in frames.items():
                    file.write(f"Frame {frame}:\n")
                    for atom_data in atoms:
                        resid, name, x, y, z = atom_data
                        file.write(f"  ResID {resid} Atom {name} at ({x:.3f}, {y:.3f}, {z:.3f})\n")
                file.write("\n")

    # Calculate and save per-cylinder and overall averages across all monomers
    with open(f"{filename_prefix}_averages.txt", "w") as file:
        total_atom_counts = []
        for (monomer, cylinder), (avg_atoms, std_dev_atoms) in overall_avg_results.items():
            file.write(f"Average number of atoms per frame for {monomer}, {cylinder}: {avg_atoms:.3f} (±{std_dev_atoms:.3f})\n")
            total_atom_counts.extend([avg_atoms] * len(cylinder_data[monomer][cylinder]))

        # Calculate overall average for each cylinder across all monomers
        for cylinder_name, counts in overall_avg_per_cylinder.items():
            avg_atoms_per_cylinder = np.mean(counts)
            std_dev_per_cylinder = np.std(counts)
            file.write(f"\nOverall average number of atoms per frame for all monomers, {cylinder_name}: {avg_atoms_per_cylinder:.3f} (±{std_dev_per_cylinder:.3f})\n")

        # Calculate overall average across all monomers and cylinders
        overall_avg_atoms = np.mean(total_atom_counts)
        overall_std_dev_atoms = np.std(total_atom_counts)
        file.write(f"\nOverall average number of atoms per frame for all monomers: {overall_avg_atoms:.3f} (±{overall_std_dev_atoms:.3f})\n")

# Save data for above and below cylinders, separated by monomer
save_cylinder_data(above_cylinders, "/media/supremeleader/Pantera/simulation/analysis_2024/average_water_molecule/above_cylinder_atoms")
save_cylinder_data(below_cylinders, "/media/supremeleader/Pantera/simulation/analysis_2024/average_water_molecule/below_cylinder_atoms")

############################# make figure ###############################################################
#########################################################################################################
import matplotlib.pyplot as plt
import re
import os

# Define file paths for the data files
above_file = "/media/supremeleader/Pantera/simulation/analysis_2024/average_water_molecule/above_cylinder_atoms_averages.txt"
below_file = "/media/supremeleader/Pantera/simulation/analysis_2024/average_water_molecule/below_cylinder_atoms_averages.txt"
central_file = "/media/supremeleader/Pantera/simulation/analysis_2024/average_water_molecule/central_cylinder_average.txt"

# Function to load data from a file into a dictionary
def load_data(file_path, is_central=False):
    data = {}
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if is_central:
                match = re.search(r"for (\w+): ([\d.]+) \(±([\d.]+)\)", line)
                if match:
                    monomer = match.group(1)
                    avg = float(match.group(2))
                    sd = float(match.group(3))
                    data[f"{monomer}_central"] = (avg, sd)
            else:
                if line.startswith("Average number of atoms"):
                    match = re.search(r"for (\w+),? (.*?): ([\d.]+) \(±([\d.]+)\)", line)
                    if match:
                        monomer = match.group(1)
                        cylinder = match.group(2).strip()
                        avg = float(match.group(3))
                        sd = float(match.group(4))
                        data[f"{monomer}_{cylinder}"] = (avg, sd)
                elif line.startswith("Overall average number of atoms"):
                    match = re.search(r"for all monomers, (.*?): ([\d.]+) \(±([\d.]+)\)", line)
                    if match:
                        cylinder = match.group(1).strip()
                        avg = float(match.group(2))
                        sd = float(match.group(3))
                        data[f"all_{cylinder}"] = (avg, sd)
    return data

# Load data from each file
above_data = load_data(above_file)
below_data = load_data(below_file)
central_data = load_data(central_file, is_central=True)

# Merge data dictionaries
combined_data = {}
combined_data.update(above_data)
combined_data.update(below_data)
combined_data.update(central_data)

# Function to plot histograms for a specified monomer or for all monomers
def plot_histograms(data, monomer):
    cylinders = [
        "above_cylinder3", "above_cylinder2", "above_cylinder1",
        "central",  # For central cylinder data
        "below_cylinder1", "below_cylinder2", "below_cylinder3"
    ]

    # Collect data for the specific monomer or for all monomers
    avg_values = []
    sd_values = []
    for cylinder in cylinders:
        key = f"{monomer}_{cylinder}" if monomer != "all" else f"all_{cylinder}"
        if key in data:
            avg, sd = data[key]
            avg_values.append(avg)
            sd_values.append(sd)
        else:
            avg_values.append(0)
            sd_values.append(0)

    # Define color scheme and label mapping
    color_scheme = [
        "lightblue", "deepskyblue", "dodgerblue", "red", "lightgreen", "mediumseagreen", "forestgreen"
    ]
    label_mapping = {
        "above_cylinder3": "Periplasm",
        "above_cylinder2": "Periplasmic side",
        "above_cylinder1": "Periplasmic side",
        "central": "Central cavity",
        "below_cylinder1": "Cytoplasmic side",
        "below_cylinder2": "Cytoplasmic side",
        "below_cylinder3": "Cytoplasm"
    }

    # Plot the histogram with error bars
    fig, ax = plt.subplots(figsize=(8, 5))
    bar_width = 0.7
    ax.bar(cylinders, avg_values, yerr=sd_values, capsize=5, color=color_scheme, width=bar_width)
    ax.set_xlabel("Channel coordinate", fontsize=12)
    ax.set_ylabel("Average number of water molecules", fontsize=12)
    ax.set_title(f"Hydration of the translocation pathway for 5WUC {monomer.capitalize()}", fontsize=14)

    # Set custom x-tick labels
    ax.set_xticks(range(len(cylinders)))
    ax.set_xticklabels([label_mapping[cylinder] for cylinder in cylinders], rotation=45)

    # Save the plot in both PDF and PNG formats
    output_dir = "/media/supremeleader/Pantera/simulation/analysis_2024/average_water_molecule"
    pdf_path = os.path.join(output_dir, f"{monomer}_hydration.pdf")
    png_path = os.path.join(output_dir, f"{monomer}_hydration.png")
    plt.tight_layout()
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight', pad_inches=0.1)
    plt.savefig(png_path, format='png', bbox_inches='tight', pad_inches=0.1)

    # Show the plot
    plt.show()

# Plot for each monomer individually
plot_histograms(combined_data, "monomer1")
plot_histograms(combined_data, "monomer2")
plot_histograms(combined_data, "monomer3")

###################################### trimer figure #####################################################
##########################################################################################################

# Function to load cylinder data into a dictionary
def load_cylinder_data(file_path, prefix):
    data = {}
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("Overall average number of atoms per frame"):
                match = re.search(r"for all monomers, (.*?): ([\d.]+) \(±([\d.]+)\)", line)
                if match:
                    cylinder = match.group(1).strip()
                    avg = float(match.group(2))
                    sd = float(match.group(3))
                    data[f"{prefix}_{cylinder}"] = (avg, sd)
    return data

# Define file paths for the data files
above_file = "/media/supremeleader/Pantera/simulation/analysis_2024/average_water_molecule/above_cylinder_atoms_averages.txt"
below_file = "/media/supremeleader/Pantera/simulation/analysis_2024/average_water_molecule/below_cylinder_atoms_averages.txt"
central_file = "/media/supremeleader/Pantera/simulation/analysis_2024/average_water_molecule/central_cylinder_average.txt"

# Load data from each file
above_data = load_cylinder_data(above_file, "above")
below_data = load_cylinder_data(below_file, "below")

# Load central cylinder data
central_data = {}
with open(central_file, 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith("Overall average number of atoms per frame"):
            match = re.search(r"for all monomers: ([\d.]+) \(±([\d.]+)\)", line)
            if match:
                avg = float(match.group(1))
                sd = float(match.group(2))
                central_data["central"] = (avg, sd)

# Prepare data for the histogram for all cylinders
cylinders = [
    "above_cylinder3", "above_cylinder2", "above_cylinder1",
    "central",
    "below_cylinder1", "below_cylinder2", "below_cylinder3"
]

avg_values = []
sd_values = []

# Collect averages and standard deviations for each specific cylinder
for cylinder in cylinders:
    if cylinder in above_data:
        avg, sd = above_data[cylinder]
        avg_values.append(avg)
        sd_values.append(sd)
    elif cylinder == "central":
        if "central" in central_data:
            avg, sd = central_data["central"]
            avg_values.append(avg)
            sd_values.append(sd)
    elif cylinder in below_data:
        avg, sd = below_data[cylinder]
        avg_values.append(avg)
        sd_values.append(sd)

# Define a unique color scheme for all cylinders histogram
all_color_scheme = [
    "lightblue",    # above_cylinder3
    "deepskyblue",  # above_cylinder2
    "dodgerblue",   # above_cylinder1
    "red",          # central
    "lightgreen",   # below_cylinder1
    "mediumseagreen", # below_cylinder2
    "forestgreen"   # below_cylinder3
]

label_mapping = {
        "above_cylinder3": "Periplasm",
        "above_cylinder2": "Periplasmic side",
        "above_cylinder1": "Periplasmic side",
        "central": "Central cavity",
        "below_cylinder1": "Cytoplasmic side",
        "below_cylinder2": "Cytoplasmic side",
        "below_cylinder3": "Cytoplasm"
    }

# Create a single histogram for all cylinders
fig, ax = plt.subplots(figsize=(8, 5))  # Increased figure size
bar_width = 0.7
ax.bar(cylinders, avg_values, yerr=sd_values, capsize=5, color=all_color_scheme, width=bar_width)
ax.set_xlabel("Channel coordinate", fontsize=12)
ax.set_ylabel("Average number of water molecules", fontsize=12)
ax.set_title("Hydration of the translocation pathway for 5WUC Trimer", fontsize=14)

# Set custom x-tick labels for all cylinders histogram
ax.set_xticks(cylinders)
ax.set_xticklabels([label_mapping[cylinder] for cylinder in cylinders], rotation=45, ha="right", fontsize=10)

# Adjust layout to ensure captions are within figure bounds
plt.tight_layout()

# Save the plot in both PDF and PNG formats
output_dir = "/media/supremeleader/Pantera/simulation/analysis_2024/average_water_molecule"
pdf_path = os.path.join(output_dir, "trimer_hydration.pdf")
png_path = os.path.join(output_dir, "trimer_hydration.png")
plt.savefig(pdf_path, format='pdf', bbox_inches="tight")
plt.savefig(png_path, format='png', bbox_inches="tight")

# Display the plot
plt.show()

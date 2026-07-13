# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
# import os
#
# # Load data
# df = pd.read_csv("/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wue/ring/5wue.csv")
#
# # Output path
# output_path = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wue/ring/top10_network_metrics_combined.png"
#
# # Choose label column
# label_col = "shared name" if "shared name" in df.columns else "name"
#
# # Metrics to plot
# metrics = ["Degree", "BetweennessCentrality", "ClosenessCentrality"]
#
# # Set plot style
# sns.set(style="whitegrid", context="talk")
#
# # Create a figure with subplots
# fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 6))
# palette = sns.color_palette("mako", 10)
#
# for ax, metric in zip(axes, metrics):
#     # Sort and select top 10
#     top10 = df.sort_values(by=metric, ascending=False).head(10)
#     sns.barplot(
#         data=top10,
#         x=metric,
#         y=label_col,
#         palette=palette,
#         ax=ax
#     )
#
#     ax.set_title(f"Top 10 by {metric}", fontsize=14)
#     ax.set_xlabel(metric)
#     ax.set_ylabel("")
#
# plt.tight_layout()
# plt.savefig(output_path, dpi=300)
# plt.show()
#
# print(f"Figure saved to: {output_path}")
#
#
#
# ############################################################################################################
# import pandas as pd
# import matplotlib.pyplot as plt
# from adjustText import adjust_text
# import os
#
# # Load CSV
# df = pd.read_csv("/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wue/ring/5wue.csv")
#
# # Output directory
# output_dir = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wue/ring/"
# label_col = "shared name" if "shared name" in df.columns else "name"
#
# # Thresholds for categorization
# degree_thresh = df["Degree"].median()
# bc_thresh = df["BetweennessCentrality"].median()
#
# # Classify into 4 categories
# def classify(row):
#     if row["Degree"] >= degree_thresh and row["BetweennessCentrality"] >= bc_thresh:
#         return "High Degree, High BC"
#     elif row["Degree"] >= degree_thresh and row["BetweennessCentrality"] < bc_thresh:
#         return "High Degree, Low BC"
#     elif row["Degree"] < degree_thresh and row["BetweennessCentrality"] >= bc_thresh:
#         return "Low Degree, High BC"
#     else:
#         return "Low Degree, Low BC"
#
# df["Category"] = df.apply(classify, axis=1)
#
# # Extract top 10 per category
# top10_per_category = {}
# highlight_mask = pd.Series(False, index=df.index)
#
# for cat in df["Category"].unique():
#     top10 = df[df["Category"] == cat].sort_values("BetweennessCentrality", ascending=False).head(10)
#     top10_per_category[cat] = top10
#     highlight_mask.loc[top10.index] = True
#
# # Print residue names with Degree and BC
# for cat, top10_df in top10_per_category.items():
#     print(f"\n🔹 Top 10 residues in category: {cat}")
#     for _, row in top10_df.iterrows():
#         print(f" - {row[label_col]} | Degree: {row['Degree']} | BC: {row['BetweennessCentrality']:.5f}")
#
# # Separate data
# highlighted_df = df[highlight_mask]
# background_df = df[~highlight_mask]
#
# # Colors
# highlight_palette = {
#     "High Degree, High BC": "red",
#     "High Degree, Low BC": "green",
#     "Low Degree, High BC": "blue",
#     "Low Degree, Low BC": "purple"
# }
# background_color = "#A0A0A0"
#
# # Plotting
# plt.figure(figsize=(10, 7))
# plt.style.use("seaborn-whitegrid")
#
# # Background
# plt.scatter(
#     background_df["Degree"],
#     background_df["BetweennessCentrality"],
#     color=background_color,
#     s=30,
#     label="_nolegend_"
# )
#
# # Highlighted points
# texts = []
# for cat, top10_df in top10_per_category.items():
#     plt.scatter(
#         top10_df["Degree"],
#         top10_df["BetweennessCentrality"],
#         label=cat,
#         s=80,
#         edgecolor="black",
#         color=highlight_palette[cat],
#         linewidth=0.6
#     )
#     for _, row in top10_df.iterrows():
#         texts.append(
#             plt.text(row["Degree"], row["BetweennessCentrality"], row[label_col], fontsize=7)
#         )
#
# # Improve label placement
# adjust_text(texts, arrowprops=dict(arrowstyle="->", color='black', lw=0.5))
#
# # Axis labels and title
# plt.title("Residue Network: Degree vs Betweenness Centrality", fontsize=14)
# plt.xlabel("Degree", fontsize=12)
# plt.ylabel("Betweenness Centrality", fontsize=12)
#
# # Legend with smaller font
# plt.legend(title="Category", fontsize=9, title_fontsize=10, loc="best")
# plt.tight_layout()
#
# # Save
# scatter_path = os.path.join(output_dir, "degree_vs_betweenness_top10_highlighted_clean.png")
# plt.savefig(scatter_path, dpi=300)
# plt.show()
#
# print(f"\n✅ Plot saved to: {scatter_path}")
#
#
# #############################################################################################################
# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
#
# # Load the data
# df = pd.read_csv("/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wue/ring/5wue.csv")
#
# # Check if residue column exists (adjust this to match your column name)
# residue_col = "Residue" if "Residue" in df.columns else "Node"
#
# # Set seaborn style
# sns.set(style="whitegrid")
#
# # Create figure and axes
# fig, axes = plt.subplots(1, 3, figsize=(16, 6))
#
# # Plot 1: Degree
# sns.boxplot(data=df, y="Degree", ax=axes[0], color="skyblue", fliersize=0)
# sns.stripplot(data=df, y="Degree", ax=axes[0], color="black", size=4, jitter=True, alpha=0.6)
# axes[0].set_title("Degree Distribution")
# axes[0].set_ylabel("Degree")
#
# # Plot 2: Betweenness Centrality
# sns.boxplot(data=df, y="BetweennessCentrality", ax=axes[1], color="salmon", fliersize=0)
# sns.stripplot(data=df, y="BetweennessCentrality", ax=axes[1], color="black", size=4, jitter=True, alpha=0.6)
# axes[1].set_title("Betweenness Centrality Distribution")
# axes[1].set_ylabel("Betweenness Centrality")
#
# # Plot 3: Closeness Centrality
# sns.boxplot(data=df, y="ClosenessCentrality", ax=axes[2], color="lightgreen", fliersize=0)
# sns.stripplot(data=df, y="ClosenessCentrality", ax=axes[2], color="black", size=4, jitter=True, alpha=0.6)
# axes[2].set_title("Closeness Centrality Distribution")
# axes[2].set_ylabel("Closeness Centrality")
#
# # Final layout adjustments
# plt.tight_layout()
#
# # Save the figure
# output_path = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wue/ring/centrality_boxplots.png"
# plt.savefig(output_path, dpi=300)
# plt.show()
#
# residue_col = "name"  # Column that holds residue names
#
# print("\n Top 10 Residues by Degree:")
# print(df[[residue_col, "Degree"]].nlargest(10, "Degree"))
#
# print("\n Top 10 Residues by Betweenness Centrality:")
# print(df[[residue_col, "BetweennessCentrality"]].nlargest(10, "BetweennessCentrality"))
#
# print("\n Top 10 Residues by Closeness Centrality:")
# print(df[[residue_col, "ClosenessCentrality"]].nlargest(10, "ClosenessCentrality"))
#
#
# print(f"\n✅ Plot saved to: {output_path}")

##################################################### Betweenness Centrality #############################################################################
# import pandas as pd
# import matplotlib.pyplot as plt
#
# # File paths
# file1 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wue/ring/5wue.csv"
# file2 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/ring/5wuc.csv"
# save_path = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/ring/betweenness_TM_comparison.png"
#
# # Read CSVs
# df1 = pd.read_csv(file1)
# df2 = pd.read_csv(file2)
#
# # Extract residue number from 'name' (e.g., GLY102 -> 102)
# df1['resnum'] = df1['name'].str.extract('(\d+)').astype(int)
# df2['resnum'] = df2['name'].str.extract('(\d+)').astype(int)
#
# # Define TM ranges for N and C triple helix bundles
# tm_regions_N = [(23, 39), (50, 73), (84, 98)]
# tm_regions_C = [(112, 138), (143, 170), (185, 202)]
#
# # Function to prepare data segments with gaps and labels
# def prepare_segments(df, tm_regions, tm_offset=0):
#     segments = []
#     ticks = []
#     x_labels = []
#     tm_labels = []
#     x_counter = 0
#
#     for idx, (start, end) in enumerate(tm_regions, start=1):
#         segment = df[(df['resnum'] >= start) & (df['resnum'] <= end)].copy()
#         segment = segment.sort_values(by='resnum')
#         segment['x'] = range(x_counter, x_counter + len(segment))
#         segments.append(segment)
#
#         # Store xtick labels
#         ticks.extend(segment['x'].tolist())
#         x_labels.extend(segment['name'].tolist())
#
#         # Store label position for helix (TM number adjusted by tm_offset)
#         tm_labels.append(((segment['x'].iloc[0] + segment['x'].iloc[-1]) / 2, f"TM{idx + tm_offset}"))
#
#         x_counter += len(segment) + 2  # Small gap of 2 between helices
#
#     return segments, ticks, x_labels, tm_labels
#
# # Prepare segments
# df1_segments_N, ticks_N, x_N, tm_labels_N = prepare_segments(df1, tm_regions_N, tm_offset=0)
# df2_segments_N, _, _, _ = prepare_segments(df2, tm_regions_N, tm_offset=0)
#
# df1_segments_C, ticks_C, x_C, tm_labels_C = prepare_segments(df1, tm_regions_C, tm_offset=3)
# df2_segments_C, _, _, _ = prepare_segments(df2, tm_regions_C, tm_offset=3)
#
# # Plotting
# fig, axs = plt.subplots(2, 1, figsize=(16, 10), sharey=True)
# colors = ['tab:blue', 'tab:orange']
#
# def plot_bundle(ax, segs1, segs2, ticks, labels, tm_labels, title):
#     # Plot once per dataset to avoid legend duplication
#     for i, (seg1, seg2) in enumerate(zip(segs1, segs2)):
#         ax.plot(seg1['x'], seg1['BetweennessCentrality'], marker='o', label='5WUE' if i == 0 else "", color=colors[0])
#         ax.plot(seg2['x'], seg2['BetweennessCentrality'], marker='s', label='5WUC' if i == 0 else "", color=colors[1])
#
#     ax.set_xticks(ticks)
#     ax.set_xticklabels(labels, fontsize=9, rotation=70)
#     ax.set_ylabel("Betweenness Centrality")
#     ax.set_title(title, fontsize=12)
#     for xpos, label in tm_labels:
#         ax.text(xpos, ax.get_ylim()[1]*0.95, label, ha='center', fontsize=9, fontweight='bold')
#     ax.legend()
#
# # N-bundle plot
# plot_bundle(axs[0], df1_segments_N, df2_segments_N, ticks_N, x_N, tm_labels_N, "N-triple Helix Bundle")
#
# # C-bundle plot
# plot_bundle(axs[1], df1_segments_C, df2_segments_C, ticks_C, x_C, tm_labels_C, "C-triple Helix Bundle")
#
# plt.tight_layout()
# plt.savefig(save_path, dpi=300)
# plt.close()
# print(f"Plot saved to: {save_path}")
# ################################################# Closeness Centrality ###################################################################
#
# import pandas as pd
# import matplotlib.pyplot as plt
#
# # File paths
# file1 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wue/ring/5wue.csv"
# file2 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/ring/5wuc.csv"
# save_path = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/ring/closeness_TM_comparison.png"
#
# # Read CSVs
# df1 = pd.read_csv(file1)
# df2 = pd.read_csv(file2)
#
# # Extract residue number from 'name' (e.g., GLY102 -> 102)
# df1['resnum'] = df1['name'].str.extract('(\d+)').astype(int)
# df2['resnum'] = df2['name'].str.extract('(\d+)').astype(int)
#
# # Define TM ranges for N and C triple helix bundles
# tm_regions_N = [(23, 39), (50, 73), (84, 98)]
# tm_regions_C = [(112, 138), (143, 170), (185, 202)]
#
# # Function to prepare data segments with gaps and labels
# def prepare_segments(df, tm_regions, label_offset=0):
#     segments = []
#     ticks = []
#     x_labels = []
#     tm_labels = []
#     x_counter = 0
#
#     for idx, (start, end) in enumerate(tm_regions, start=1):
#         segment = df[(df['resnum'] >= start) & (df['resnum'] <= end)].copy()
#         segment = segment.sort_values(by='resnum')
#         segment['x'] = range(x_counter, x_counter + len(segment))
#         segments.append(segment)
#
#         # Store xtick labels
#         ticks.extend(segment['x'].tolist())
#         x_labels.extend(segment['name'].tolist())
#
#         # Store label position for helix
#         tm_labels.append(((segment['x'].iloc[0] + segment['x'].iloc[-1]) / 2, f"TM{idx + label_offset}"))
#
#         x_counter += len(segment) + 2  # Small gap of 2 between helices
#
#     return segments, ticks, x_labels, tm_labels
#
# # Prepare segments
# df1_segments_N, ticks_N, x_N, tm_labels_N = prepare_segments(df1, tm_regions_N, label_offset=0)
# df2_segments_N, _, _, _ = prepare_segments(df2, tm_regions_N, label_offset=0)
#
# df1_segments_C, ticks_C, x_C, tm_labels_C = prepare_segments(df1, tm_regions_C, label_offset=3)
# df2_segments_C, _, _, _ = prepare_segments(df2, tm_regions_C, label_offset=3)
#
# # Plotting
# fig, axs = plt.subplots(2, 1, figsize=(16, 10), sharey=True)
# colors = ['tab:blue', 'tab:orange']
#
# def plot_bundle(ax, segs1, segs2, ticks, labels, tm_labels, title):
#     legend_shown = False
#     for seg1, seg2 in zip(segs1, segs2):
#         ax.plot(seg1['x'], seg1['ClosenessCentrality'], marker='o', label='5WUE' if not legend_shown else "", color=colors[0])
#         ax.plot(seg2['x'], seg2['ClosenessCentrality'], marker='s', label='5WUC' if not legend_shown else "", color=colors[1])
#         legend_shown = True
#
#     ax.set_xticks(ticks)
#     ax.set_xticklabels(labels, fontsize=9, rotation=70)
#     ax.set_ylabel("Closeness Centrality")
#     ax.set_title(title, fontsize=12)
#     for xpos, label in tm_labels:
#         ax.text(xpos, ax.get_ylim()[1]*0.95, label, ha='center', fontsize=9, fontweight='bold')
#     ax.legend()
#
# # N-bundle plot
# plot_bundle(axs[0], df1_segments_N, df2_segments_N, ticks_N, x_N, tm_labels_N, "N-triple Helix Bundle")
#
# # C-bundle plot
# plot_bundle(axs[1], df1_segments_C, df2_segments_C, ticks_C, x_C, tm_labels_C, "C-triple Helix Bundle")
#
# plt.tight_layout()
# plt.savefig(save_path, dpi=300)
# plt.close()
# print(f"Plot saved to: {save_path}")

################################################## DEGREE ########################################################
# import pandas as pd
# import matplotlib.pyplot as plt
#
# # File paths
# file1 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wue/ring/5wue.csv"
# file2 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/ring/5wuc.csv"
# save_path = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/ring/degree_TM_comparison.png"
#
# # Read CSVs
# df1 = pd.read_csv(file1)
# df2 = pd.read_csv(file2)
#
# # Extract residue number from 'name' (e.g., GLY102 -> 102)
# df1['resnum'] = df1['name'].str.extract('(\d+)').astype(int)
# df2['resnum'] = df2['name'].str.extract('(\d+)').astype(int)
#
# # Define TM ranges for N and C triple helix bundles
# tm_regions_N = [(23, 39), (50, 73), (84, 98)]
# tm_regions_C = [(112, 138), (143, 170), (185, 202)]
#
# # Function to prepare segments
# def prepare_segments(df, tm_regions, start_tm_number):
#     segments = []
#     ticks = []
#     x_labels = []
#     tm_labels = []
#     x_counter = 0
#
#     for idx, (start, end) in enumerate(tm_regions, start=start_tm_number):
#         segment = df[(df['resnum'] >= start) & (df['resnum'] <= end)].copy()
#         segment = segment.sort_values(by='resnum')
#         segment['x'] = range(x_counter, x_counter + len(segment))
#         segments.append(segment)
#
#         ticks.extend(segment['x'].tolist())
#         x_labels.extend(segment['name'].tolist())
#
#         tm_labels.append(((segment['x'].iloc[0] + segment['x'].iloc[-1]) / 2, f"TM{idx}"))
#
#         x_counter += len(segment) + 2  # gap
#
#     return segments, ticks, x_labels, tm_labels
#
# # Prepare segments
# df1_segments_N, ticks_N, x_N, tm_labels_N = prepare_segments(df1, tm_regions_N, start_tm_number=1)
# df2_segments_N, _, _, _ = prepare_segments(df2, tm_regions_N, start_tm_number=1)
#
# df1_segments_C, ticks_C, x_C, tm_labels_C = prepare_segments(df1, tm_regions_C, start_tm_number=4)
# df2_segments_C, _, _, _ = prepare_segments(df2, tm_regions_C, start_tm_number=4)
#
# # Plotting
# fig, axs = plt.subplots(2, 1, figsize=(16, 10), sharey=True)
# colors = ['tab:blue', 'tab:orange']
#
# def plot_bundle(ax, segs1, segs2, ticks, labels, tm_labels, title):
#     plotted_labels = set()
#
#     for seg1, seg2 in zip(segs1, segs2):
#         label1 = '5WUE' if '5WUE' not in plotted_labels else None
#         label2 = '5WUC' if '5WUC' not in plotted_labels else None
#
#         ax.plot(seg1['x'], seg1['Degree'], marker='o', color=colors[0], label=label1)
#         ax.plot(seg2['x'], seg2['Degree'], marker='s', color=colors[1], label=label2)
#
#         plotted_labels.update(filter(None, [label1, label2]))
#
#     ax.set_xticks(ticks)
#     ax.set_xticklabels(labels, fontsize=9, rotation=70)
#     ax.set_ylabel("Degree Centrality")
#     ax.set_title(title, fontsize=12)
#     for xpos, label in tm_labels:
#         ax.text(xpos, ax.get_ylim()[1]*0.95, label, ha='center', fontsize=9, fontweight='bold')
#     ax.legend()
#
# # N-bundle plot
# plot_bundle(axs[0], df1_segments_N, df2_segments_N, ticks_N, x_N, tm_labels_N, "N-triple Helix Bundle")
#
# # C-bundle plot
# plot_bundle(axs[1], df1_segments_C, df2_segments_C, ticks_C, x_C, tm_labels_C, "C-triple Helix Bundle")
#
# plt.tight_layout()
# plt.savefig(save_path, dpi=300)
# plt.close()
# print(f"Plot saved to: {save_path}")
#
# ###################################################### Clustering Coefficient ############################################################
#
# import pandas as pd
# import matplotlib.pyplot as plt
#
# # File paths
# file1 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wue/ring/5wue.csv"
# file2 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/ring/5wuc.csv"
# save_path = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/ring/clustering_TM_comparison.png"
#
# # Read CSVs
# df1 = pd.read_csv(file1)
# df2 = pd.read_csv(file2)
#
# # Extract residue number from 'name' (e.g., GLY102 -> 102)
# df1['resnum'] = df1['name'].str.extract('(\d+)').astype(int)
# df2['resnum'] = df2['name'].str.extract('(\d+)').astype(int)
#
# # Define TM ranges for N and C triple helix bundles
# tm_regions_N = [(23, 39), (50, 73), (84, 98)]
# tm_regions_C = [(112, 138), (143, 170), (185, 202)]
#
# # Function to prepare segments
# def prepare_segments(df, tm_regions, start_tm_number):
#     segments = []
#     ticks = []
#     x_labels = []
#     tm_labels = []
#     x_counter = 0
#
#     for idx, (start, end) in enumerate(tm_regions, start=start_tm_number):
#         segment = df[(df['resnum'] >= start) & (df['resnum'] <= end)].copy()
#         segment = segment.sort_values(by='resnum')
#         segment['x'] = range(x_counter, x_counter + len(segment))
#         segments.append(segment)
#
#         ticks.extend(segment['x'].tolist())
#         x_labels.extend(segment['name'].tolist())
#
#         tm_labels.append(((segment['x'].iloc[0] + segment['x'].iloc[-1]) / 2, f"TM{idx}"))
#
#         x_counter += len(segment) + 2  # gap between helices
#
#     return segments, ticks, x_labels, tm_labels
#
# # Prepare segments
# df1_segments_N, ticks_N, x_N, tm_labels_N = prepare_segments(df1, tm_regions_N, start_tm_number=1)
# df2_segments_N, _, _, _ = prepare_segments(df2, tm_regions_N, start_tm_number=1)
#
# df1_segments_C, ticks_C, x_C, tm_labels_C = prepare_segments(df1, tm_regions_C, start_tm_number=4)
# df2_segments_C, _, _, _ = prepare_segments(df2, tm_regions_C, start_tm_number=4)
#
# # Plotting
# fig, axs = plt.subplots(2, 1, figsize=(16, 10), sharey=True)
# colors = ['tab:blue', 'tab:orange']
#
# def plot_bundle(ax, segs1, segs2, ticks, labels, tm_labels, title):
#     plotted_labels = set()
#
#     for seg1, seg2 in zip(segs1, segs2):
#         label1 = '5WUE' if '5WUE' not in plotted_labels else None
#         label2 = '5WUC' if '5WUC' not in plotted_labels else None
#
#         ax.plot(seg1['x'], seg1['ClusteringCoefficient'], marker='o', color=colors[0], label=label1)
#         ax.plot(seg2['x'], seg2['ClusteringCoefficient'], marker='s', color=colors[1], label=label2)
#
#         plotted_labels.update(filter(None, [label1, label2]))
#
#     ax.set_xticks(ticks)
#     ax.set_xticklabels(labels, fontsize=9, rotation=70)
#     ax.set_ylabel("Clustering Coefficient")
#     ax.set_title(title, fontsize=12)
#     for xpos, label in tm_labels:
#         ax.text(xpos, ax.get_ylim()[1]*0.95, label, ha='center', fontsize=9, fontweight='bold')
#     ax.legend()
#
# # N-bundle plot
# plot_bundle(axs[0], df1_segments_N, df2_segments_N, ticks_N, x_N, tm_labels_N, "N-triple Helix Bundle")
#
# # C-bundle plot
# plot_bundle(axs[1], df1_segments_C, df2_segments_C, ticks_C, x_C, tm_labels_C, "C-triple Helix Bundle")
#
# plt.tight_layout()
# plt.savefig(save_path, dpi=300)
# plt.close()
# print(f"Plot saved to: {save_path}")
#
#
# ###################################################### Average Shortest Path Length ############################################################
#
# import pandas as pd
# import matplotlib.pyplot as plt
#
# # File paths
# file1 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wue/ring/5wue.csv"
# file2 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/ring/5wuc.csv"
# save_path = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/ring/aspl_TM_comparison.png"
#
# # Read CSVs
# df1 = pd.read_csv(file1)
# df2 = pd.read_csv(file2)
#
# # Extract residue number from 'name' (e.g., GLY102 -> 102)
# df1['resnum'] = df1['name'].str.extract('(\d+)').astype(int)
# df2['resnum'] = df2['name'].str.extract('(\d+)').astype(int)
#
# # Define TM regions
# tm_regions_N = [(23, 39), (50, 73), (84, 98)]
# tm_regions_C = [(112, 138), (143, 170), (185, 202)]
#
# # Function to prepare segments
# def prepare_segments(df, tm_regions, start_tm_number):
#     segments = []
#     ticks = []
#     x_labels = []
#     tm_labels = []
#     x_counter = 0
#
#     for idx, (start, end) in enumerate(tm_regions, start=start_tm_number):
#         segment = df[(df['resnum'] >= start) & (df['resnum'] <= end)].copy()
#         segment = segment.sort_values(by='resnum')
#         segment['x'] = range(x_counter, x_counter + len(segment))
#         segments.append(segment)
#
#         ticks.extend(segment['x'].tolist())
#         x_labels.extend(segment['name'].tolist())
#
#         tm_labels.append(((segment['x'].iloc[0] + segment['x'].iloc[-1]) / 2, f"TM{idx}"))
#
#         x_counter += len(segment) + 2  # gap between TM helices
#
#     return segments, ticks, x_labels, tm_labels
#
# # Prepare segments for both files
# df1_segments_N, ticks_N, x_N, tm_labels_N = prepare_segments(df1, tm_regions_N, start_tm_number=1)
# df2_segments_N, _, _, _ = prepare_segments(df2, tm_regions_N, start_tm_number=1)
#
# df1_segments_C, ticks_C, x_C, tm_labels_C = prepare_segments(df1, tm_regions_C, start_tm_number=4)
# df2_segments_C, _, _, _ = prepare_segments(df2, tm_regions_C, start_tm_number=4)
#
# # Plotting
# fig, axs = plt.subplots(2, 1, figsize=(16, 10), sharey=True)
# colors = ['tab:blue', 'tab:orange']
#
# def plot_bundle(ax, segs1, segs2, ticks, labels, tm_labels, title):
#     plotted_labels = set()
#
#     for seg1, seg2 in zip(segs1, segs2):
#         label1 = '5WUE' if '5WUE' not in plotted_labels else None
#         label2 = '5WUC' if '5WUC' not in plotted_labels else None
#
#         ax.plot(seg1['x'], seg1['AverageShortestPathLength'], marker='o', color=colors[0], label=label1)
#         ax.plot(seg2['x'], seg2['AverageShortestPathLength'], marker='s', color=colors[1], label=label2)
#
#         plotted_labels.update(filter(None, [label1, label2]))
#
#     ax.set_xticks(ticks)
#     ax.set_xticklabels(labels, fontsize=9, rotation=70)
#     ax.set_ylabel("Average Shortest Path Length")
#     ax.set_title(title, fontsize=12)
#     for xpos, label in tm_labels:
#         ax.text(xpos, ax.get_ylim()[1]*0.95, label, ha='center', fontsize=9, fontweight='bold')
#     ax.legend()
#
# # N-bundle plot
# plot_bundle(axs[0], df1_segments_N, df2_segments_N, ticks_N, x_N, tm_labels_N, "N-triple Helix Bundle")
#
# # C-bundle plot
# plot_bundle(axs[1], df1_segments_C, df2_segments_C, ticks_C, x_C, tm_labels_C, "C-triple Helix Bundle")
#
# plt.tight_layout()
# plt.savefig(save_path, dpi=300)
# plt.close()
# print(f"Plot saved to: {save_path}")

############################### euclidean distance #######################################3
# import os
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import MDAnalysis as mda
#
# # === File paths ===
# gro_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1_new.gro"
# xtc_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/chain1_new.xtc"
# output_dir = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/ring/"
# os.makedirs(output_dir, exist_ok=True)
#
# # === Load trajectory ===
# u = mda.Universe(gro_file, xtc_file)
#
# # === Atom selections ===
# k129_ca = u.select_atoms("resid 129 and name CA")
# k129_nz = u.select_atoms("resid 129 and name NZ")
#
# y29_ca = u.select_atoms("resid 29 and name CA")
# y29_oh = u.select_atoms("resid 29 and name OH")   # Tyr OH group
#
# s65_ca = u.select_atoms("resid 65 and name CA")
# s65_og = u.select_atoms("resid 65 and name OG")   # Ser OH group
#
# # Check that all selections returned exactly 1 atom
# sel_checks = [
#     ("K129 CA", k129_ca), ("K129 NZ", k129_nz),
#     ("Y29 CA", y29_ca), ("Y29 OH", y29_oh),
#     ("S65 CA", s65_ca), ("S65 OG", s65_og),
# ]
# for name, sel in sel_checks:
#     if len(sel) != 1:
#         raise ValueError(f"Selection for {name} returned {len(sel)} atoms — check resid/atom name.")
#
# # === Storage ===
# frames = []
# k129_y29_ca = []
# k129_y29_side = []
# k129_s65_ca = []
# k129_s65_side = []
# y29_s65_ca = []
# y29_s65_side = []
#
# # === Loop over trajectory every 100 frames ===
# for ts in u.trajectory[::100]:
#     frames.append(ts.frame)
#
#     # Positions
#     pos_kca = k129_ca.positions[0]
#     pos_knz = k129_nz.positions[0]
#     pos_yca = y29_ca.positions[0]
#     pos_yoh = y29_oh.positions[0]
#     pos_sca = s65_ca.positions[0]
#     pos_sog = s65_og.positions[0]
#
#     # Distances
#     d_kca_yca = np.linalg.norm(pos_kca - pos_yca)
#     d_knz_yoh = np.linalg.norm(pos_knz - pos_yoh)
#
#     d_kca_sca = np.linalg.norm(pos_kca - pos_sca)
#     d_knz_sog = np.linalg.norm(pos_knz - pos_sog)
#
#     d_yca_sca = np.linalg.norm(pos_yca - pos_sca)
#     d_yoh_sog = np.linalg.norm(pos_yoh - pos_sog)
#
#     # Append
#     k129_y29_ca.append(d_kca_yca)
#     k129_y29_side.append(d_knz_yoh)
#     k129_s65_ca.append(d_kca_sca)
#     k129_s65_side.append(d_knz_sog)
#     y29_s65_ca.append(d_yca_sca)
#     y29_s65_side.append(d_yoh_sog)
#
# # === Save to CSV ===
# df = pd.DataFrame({
#     "Frame": frames,
#     "K129_CA__Y29_CA_A": k129_y29_ca,
#     "K129_NZ__Y29_OH_A": k129_y29_side,
#     "K129_CA__S65_CA_A": k129_s65_ca,
#     "K129_NZ__S65_OG_A": k129_s65_side,
#     "Y29_CA__S65_CA_A": y29_s65_ca,
#     "Y29_OH__S65_OG_A": y29_s65_side
# })
# csv_path = os.path.join(output_dir, "K129_Y29_S65_all_distances_every100frames.csv")
# df.to_csv(csv_path, index=False)
# print(f"Saved distances CSV: {csv_path}")
#
# # === Plot ===
# plt.figure(figsize=(12, 10))
#
# # Subplot 1: CA-CA distances
# ax1 = plt.subplot(2, 1, 1)
# ax1.plot(frames, k129_y29_ca, label="K129 CA - Y29 CA", linewidth=1.5)
# ax1.plot(frames, k129_s65_ca, label="K129 CA - S65 CA", linewidth=1.5)
# ax1.plot(frames, y29_s65_ca, label="Y29 CA - S65 CA", linewidth=1.5)
# ax1.set_ylabel("Distance (Å)")
# ax1.set_title("Cα - Cα distances in 6IYX simulation")
# ax1.legend()
# ax1.grid(True, linestyle="--", alpha=0.5)
#
# # Subplot 2: Side-chain distances
# ax2 = plt.subplot(2, 1, 2)
# ax2.plot(frames, k129_y29_side, label="K129 NZ - Y29 OH", linewidth=1.5)
# ax2.plot(frames, k129_s65_side, label="K129 NZ - S65 OG", linewidth=1.5)
# ax2.plot(frames, y29_s65_side, label="Y29 OH - S65 OG", linewidth=1.5)
# ax2.set_xlabel("Frame")
# ax2.set_ylabel("Distance (Å)")
# ax2.set_title("Side-chain terminal atom distances in 6IYX simulation")
# ax2.legend()
# ax2.grid(True, linestyle="--", alpha=0.5)
#
# plt.tight_layout()
# plot_path = os.path.join(output_dir, "K129_Y29_S65_all_distances_every100frames.png")
# plt.savefig(plot_path, dpi=300)
# plt.show()
#
# print(f"Saved plot PNG: {plot_path}")



############################################### entire comparion ################################
# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
# import os
#
# # ===== User Inputs =====
# file1 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wue/ring/5wue.csv"
# file2 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/ring/5wuc.csv"
# file3 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/ring/6iyx.csv"
# file4 = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/ring/6iyz.csv"
# output_dir = "/media/supremeleader/Pantera/simulation/analysis_2024/ring"
#
# # ===== Load Data =====
# dfs = []
# for file, label in zip(
#     [file1, file2, file3, file4],
#     ["5wue", "5wuc", "6iyx", "6iyz"]
# ):
#     df = pd.read_csv(file)
#     df["Structure"] = label
#     dfs.append(df)
#
# df_all = pd.concat(dfs, ignore_index=True)
#
# # Determine label column safely
# possible_label_cols = ["shared name", "name", "Name"]
# label_col = next((col for col in possible_label_cols if col in df_all.columns), None)
# if label_col is None:
#     raise ValueError(f"None of {possible_label_cols} found in CSV columns.")
#
# # Metrics
# metrics = ["Degree", "BetweennessCentrality", "ClosenessCentrality"]
#
# # ===== Plot Setup =====
# sns.set(style="whitegrid", context="talk")
# palette = ["#a8e6cf", "#ffd3b6", "#aec6cf", "#ffb6c1"]
#
# fig, axes = plt.subplots(
#     nrows=2, ncols=len(metrics),
#     figsize=(18, 10)
# )
#
# # ===== Row 1: Boxplots =====
# for col_idx, metric in enumerate(metrics):
#     ax = axes[0, col_idx]
#     sns.boxplot(
#         data=df_all, x="Structure", y=metric,
#         palette=palette, ax=ax, showfliers=False
#     )
#     sns.stripplot(
#         data=df_all, x="Structure", y=metric,
#         color="black", size=3, jitter=True, alpha=0.3, ax=ax
#     )
#     ax.set_title(f"{metric} Distribution", fontsize=16, fontweight='bold')
#     ax.set_xlabel("Structure", fontsize=14, fontweight='bold')
#     ax.set_ylabel(metric, fontsize=14, fontweight='bold')
#     # safe legend removal
#     if ax.get_legend() is not None:
#         ax.get_legend().remove()
#
# # ===== Row 2: Grouped Top 5 residues =====
# for col_idx, metric in enumerate(metrics):
#     ax = axes[1, col_idx]
#
#     # Get top 5 residues per structure
#     top5_list = []
#     for structure in df_all["Structure"].unique():
#         sub_df = df_all[df_all["Structure"] == structure]
#         top5 = sub_df.nlargest(5, metric).copy()
#         top5_list.append(top5)
#
#     top5_data = pd.concat(top5_list, ignore_index=True)
#
#     # Residue labels grouped by structure
#     top5_data["Residue_Label"] = top5_data["Structure"] + "_" + top5_data[label_col]
#     order = top5_data["Residue_Label"]
#
#     sns.barplot(
#         data=top5_data,
#         x="Residue_Label", y=metric,
#         hue="Structure", palette=palette,
#         dodge=False, ax=ax,
#         order=order
#     )
#
#     ax.set_title(f"Top 5 Residues by {metric}", fontsize=14, fontweight='bold')
#     ax.set_xlabel("Residues (Grouped by Structure)", fontsize=12, fontweight='bold')
#     ax.set_ylabel(metric, fontsize=12, fontweight='bold')
#     ax.tick_params(axis='x', rotation=90, labelsize=8)
#     if ax.get_legend() is not None:
#         ax.get_legend().remove()
#
# plt.tight_layout()
# os.makedirs(output_dir, exist_ok=True)
# output_path = os.path.join(output_dir, "network_metrics_lightcolors_grouped_top5.png")
# plt.savefig(output_path, dpi=300)
# plt.show()
#
# print(f"Figure saved to: {output_path}")

################################ make freq files ###################################
     ######### distance between c-alpha atom ########

# import os
# import MDAnalysis as mda
# import numpy as np
# import pandas as pd
# from itertools import combinations
#
# # Input files
# gro_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wue/4_chain1.gro"
# xtc_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wue/4_chain1.xtc"
# output_folder = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wue/4_network"
#
# os.makedirs(output_folder, exist_ok=True)
#
# contacts_csv = os.path.join(output_folder, "contacts_heavyatom.csv")
# network_file = os.path.join(output_folder, "network_heavyatom.pdb_ringEdges")
#
# # Load universe
# u = mda.Universe(gro_file, xtc_file)
# residues = u.residues
#
# # Pre-select heavy atoms only (C, N, O, S)
# heavy_atoms = u.select_atoms("name C N O S")
#
# # Map residue → its heavy atoms
# res_to_heavy = {res.resid: res.atoms.select_atoms("name C N O S") for res in residues}
#
# all_contacts = []
#
# print("Processing trajectory (using heavy atom minimum distance)...")
#
# for ts in u.trajectory:
#     frame_num = ts.frame
#
#     # Loop over all unique residue pairs
#     for i, j in combinations(residues, 2):
#         heavy_i = res_to_heavy[i.resid]
#         heavy_j = res_to_heavy[j.resid]
#
#         # Skip if one residue has no heavy atoms (rare)
#         if len(heavy_i) == 0 or len(heavy_j) == 0:
#             continue
#
#         # Compute all inter-atomic distances between heavy atoms of residue i and j
#         diff = heavy_i.positions[:, None, :] - heavy_j.positions[None, :, :]
#         distances = np.linalg.norm(diff, axis=-1)
#         min_dist = np.min(distances)
#
#         # If any heavy atom pair within 4 Å, record
#         if min_dist <= 4.0:
#             res1_label = f"{i.resname}_{i.resid}"
#             res2_label = f"{j.resname}_{j.resid}"
#             all_contacts.append([frame_num, res1_label, res2_label, min_dist])
#
# print(f"Processed {len(u.trajectory)} frames. Writing to CSV...")
#
# # Save to CSV
# df = pd.DataFrame(all_contacts, columns=["frame", "residue1", "residue2", "min_distance"])
# df.to_csv(contacts_csv, index=False)
# print(f"Contacts (heavy-atom based) saved to: {contacts_csv}")
#
# # ---- Step 2: Check persistence (50% rule) ----
# print("Analyzing contact persistence...")
#
# n_frames = len(u.trajectory)
#
# # Normalize residue pair (order-independent)
# df["pair"] = df.apply(lambda x: tuple(sorted([x["residue1"], x["residue2"]])), axis=1)
#
# # Count unique frames per pair
# pair_frames = df.groupby("pair")["frame"].nunique().reset_index()
# pair_frames["fraction"] = pair_frames["frame"] / n_frames
#
# # Filter contacts present in >= 50% frames
# persistent_pairs = pair_frames[pair_frames["fraction"] >= 0.5]["pair"]
#
# # Compute average minimal distances for those pairs
# mean_distances = (
#     df[df["pair"].isin(persistent_pairs)]
#     .groupby("pair")["min_distance"]
#     .mean()
#     .reset_index()
# )
#
# # Write network edges file
# with open(network_file, "w") as f:
#     f.write("residue1,residue2,distance\n")
#     for _, row in mean_distances.iterrows():
#         res1, res2 = row["pair"]
#         f.write(f"{res1},{res2},{row['min_distance']:.3f}\n")
#
# print(f"Network edges written to: {network_file}")

##################################################################

# import pandas as pd
# import networkx as nx
# import matplotlib.pyplot as plt
#
# # === Step 1: Load your graph data ===
# file_path = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/4_graph/network_heavyatom.pdb_ringEdges"  # Path to your CSV file
# df = pd.read_csv(file_path)
#
# # === Step 2: Build the graph ===
# G = nx.Graph()
#
# # Optionally, filter by distance cutoff (e.g., ≤ 4 Å)
# distance_cutoff = 4.0
# for _, row in df.iterrows():
#     if row["distance"] <= distance_cutoff:
#         G.add_edge(row["residue1"], row["residue2"], weight=row["distance"])
#
# print(f"✅ Graph loaded with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
#
# # === Step 3: Detect communities ===
# communities = list(nx.algorithms.community.greedy_modularity_communities(G))
# print(f"✅ Number of communities found: {len(communities)}")
#
# # === Step 4: Assign a community ID to each node ===
# node_community = {}
# for i, comm in enumerate(communities):
#     for node in comm:
#         node_community[node] = i
#
# # === Step 5: Visualize the communities ===
# plt.figure(figsize=(12, 10))
#
# # Create a layout for better visualization
# pos = nx.spring_layout(G, seed=42, k=0.5)
#
# # Assign colors to communities
# colors = [node_community[n] for n in G.nodes()]
#
# # Draw nodes, edges, and labels
# nx.draw_networkx_nodes(G, pos, node_size=300, cmap=plt.cm.tab10, node_color=colors)
# nx.draw_networkx_edges(G, pos, alpha=0.3)
# nx.draw_networkx_labels(G, pos, font_size=8)
#
# plt.title("Protein Contact Network — Community Structure", fontsize=14)
# plt.axis("off")
# plt.tight_layout()
# plt.show()
#
# # === Step 6: Save community assignments ===
# community_data = []
# for i, comm in enumerate(communities, 1):
#     for node in comm:
#         community_data.append({"community_id": i, "residue": node})
#
# pd.DataFrame(community_data).to_csv("/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/4_graph/output_communities.csv", index=False)


########################################## analysis and figures everything ###########################
######################################## use this  ######################################################
##############################################################################################

# import pandas as pd
# import networkx as nx
# import os
# import matplotlib.pyplot as plt
# from networkx.algorithms import community
# import numpy as np
#
# # ============================================================
# # FILE PATHS
# # ============================================================
# input_csv = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/4_graph/network_heavyatom.pdb_ringEdges"
# output_dir = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/4_graph"
# os.makedirs(output_dir, exist_ok=True)
#
# # ============================================================
# # LOAD CSV
# # ============================================================
# df = pd.read_csv(input_csv)
# df.columns = [c.strip() for c in df.columns]
#
# # ============================================================
# # BUILD GRAPH
# # ============================================================
# G = nx.Graph()
# for _, row in df.iterrows():
#     res1, res2, dist = row["residue1"], row["residue2"], float(row["distance"])
#     G.add_edge(res1, res2, weight=dist)
#
# print(f"\nGraph loaded: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
#
# # ============================================================
# # CENTRALITY MEASURES
# # ============================================================
# degree_centrality = nx.degree_centrality(G)
# betweenness_centrality = nx.betweenness_centrality(G, weight="weight")
# closeness_centrality = nx.closeness_centrality(G, distance="weight")
# eigenvector_centrality = nx.eigenvector_centrality(G, weight="weight", max_iter=2000)
# clustering_coeff = nx.clustering(G, weight="weight")
#
# # ============================================================
# # COMMUNITY DETECTION
# # ============================================================
# communities = community.greedy_modularity_communities(G, weight="weight")
# community_of = {}
#
# for i, comm in enumerate(communities):
#     for node in comm:
#         community_of[node] = i + 1
#
# community_df = pd.DataFrame(
#     [(n, community_of[n]) for n in G.nodes()],
#     columns=["Residue", "Community"]
# )
# community_df.to_csv(os.path.join(output_dir, "network_communities.csv"), index=False)
#
# # ============================================================
# # CLIQUES
# # ============================================================
# cliques = list(nx.find_cliques(G))
# cliques_sorted = sorted(cliques, key=len, reverse=True)
#
# clique_df = pd.DataFrame({
#     "Clique_ID": range(1, len(cliques) + 1),
#     "Residues": [", ".join(c) for c in cliques],
#     "Size": [len(c) for c in cliques]
# })
# clique_df.to_csv(os.path.join(output_dir, "network_cliques.csv"), index=False)
#
# # ============================================================
# # N-TERMINAL / C-TERMINAL DETECTION
# # ============================================================
# aa_3to1 = {
#     "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
#     "GLU":"E","GLN":"Q","GLY":"G","HIS":"H","ILE":"I",
#     "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
#     "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
# }
#
# def residue_number(res):
#     return int(res.split("_")[-1])
#
# def residue_one_letter(res):
#     aa3 = res.split("_")[0]
#     num = res.split("_")[1]
#     return aa_3to1.get(aa3, "?") + num
#
# res_numbers = [residue_number(r) for r in G.nodes()]
# residues_sorted = sorted(G.nodes(), key=lambda x: residue_number(x))
#
# n_term = residues_sorted[0]
# c_term = residues_sorted[-1]
#
# n_term_label = residue_one_letter(n_term)
# c_term_label = residue_one_letter(c_term)
#
# print(f"N-terminal: {n_term_label}, C-terminal: {c_term_label}")
#
# # ============================================================
# # SHORTEST PATH
# # ============================================================
# try:
#     shortest_path_nodes = nx.shortest_path(G, source=n_term, target=c_term, weight="weight")
#     shortest_path_length = nx.shortest_path_length(G, source=n_term, target=c_term, weight="weight")
# except nx.NetworkXNoPath:
#     shortest_path_nodes = []
#     shortest_path_length = None
#
# # ============================================================
# # STORE ALL CENTRALITIES
# # ============================================================
# centrality_df = pd.DataFrame({
#     "Residue": list(G.nodes()),
#     "Label": [residue_one_letter(n) for n in G.nodes()],
#     "Degree_Centrality": [degree_centrality[n] for n in G.nodes()],
#     "Betweenness_Centrality": [betweenness_centrality[n] for n in G.nodes()],
#     "Closeness_Centrality": [closeness_centrality[n] for n in G.nodes()],
#     "Eigenvector_Centrality": [eigenvector_centrality[n] for n in G.nodes()],
#     "Clustering_Coefficient": [clustering_coeff[n] for n in G.nodes()],
#     "Community": [community_of[n] for n in G.nodes()]
# })
#
# centrality_df.to_csv(os.path.join(output_dir, "network_centrality_metrics.csv"), index=False)
#
# # ============================================================
# # PRINT TOP-5 PER METRIC
# # ============================================================
# def top5(metric_name):
#     print(f"\nTop 5 residues by {metric_name}:")
#     print(centrality_df.sort_values(metric_name, ascending=False)[["Residue", "Label", metric_name]].head(5))
#
# top5("Degree_Centrality")
# top5("Betweenness_Centrality")
# top5("Closeness_Centrality")
# top5("Eigenvector_Centrality")
# top5("Clustering_Coefficient")
#
# # ============================================================
# # SAVE SUMMARY
# # ============================================================
# with open(os.path.join(output_dir, "network_summary.txt"), "w") as f:
#     f.write(f"Nodes: {G.number_of_nodes()}\n")
#     f.write(f"Edges: {G.number_of_edges()}\n")
#     f.write(f"Communities: {len(communities)}\n")
#     f.write(f"N-terminal: {n_term_label}\n")
#     f.write(f"C-terminal: {c_term_label}\n")
#     f.write(f"Shortest path length: {shortest_path_length}\n")
#     f.write(f"Shortest path nodes: {shortest_path_nodes}\n")
#
# # ============================================================
# # LAYOUT FOR VISUALIZATION
# # ============================================================
# print("\nComputing spaced-out layout...")
# pos = nx.spring_layout(G, seed=42, k=1.2, iterations=300)
#
# # for labeling
# label_map = {n: residue_one_letter(n) for n in G.nodes()}
#
# def draw_shifted_labels(G, pos, labels, ax):
#     for n, (x, y) in pos.items():
#         ax.text(x, y + 0.04, labels[n], fontsize=8, ha="center", va="center")
#
# # ============================================================
# # PLOT 1 — COMMUNITIES
# # ============================================================
# plt.figure(figsize=(14, 12))
# ax = plt.gca()
#
# nx.draw_networkx_edges(G, pos, alpha=0.25, edge_color="gray", width=0.7)
#
# node_colors = [community_of[n] for n in G.nodes()]
# nx.draw_networkx_nodes(
#     G, pos,
#     node_color=node_colors,
#     cmap=plt.cm.tab20,
#     node_size=220,
#     edgecolors="black",
#     linewidths=0.6
# )
#
# draw_shifted_labels(G, pos, label_map, ax)
# plt.title("Amino Acid Interaction Network — Communities", fontsize=15)
# plt.axis("off")
# plt.tight_layout()
# plt.savefig(os.path.join(output_dir, "graph_communities_spaced.png"), dpi=300)
# plt.close()
#
# # ============================================================
# # PLOT 2 — CLIQUES
# # ============================================================
# plt.figure(figsize=(14, 12))
# ax = plt.gca()
#
# nx.draw_networkx_edges(G, pos, alpha=0.2, edge_color="gray")
# nx.draw_networkx_nodes(G, pos, node_color="lightgray", node_size=180, edgecolors="k")
#
# top_cliques = cliques_sorted[:3]
# colors = ["red", "blue", "green"]
#
# for idx, clique in enumerate(top_cliques):
#     nx.draw_networkx_nodes(
#         G, pos,
#         nodelist=clique,
#         node_color=colors[idx],
#         node_size=320,
#         edgecolors="black"
#     )
#     clique_edges = [(u, v) for u in clique for v in clique if G.has_edge(u, v)]
#     nx.draw_networkx_edges(G, pos, edgelist=clique_edges, edge_color=colors[idx], width=2.0)
#
# draw_shifted_labels(G, pos, label_map, ax)
# plt.title("Amino Acid Interaction Network — Largest Cliques", fontsize=15)
# plt.axis("off")
# plt.tight_layout()
# plt.savefig(os.path.join(output_dir, "graph_cliques_spaced.png"), dpi=300)
# plt.close()
#
# # ============================================================
# # PLOT 3 — SHORTEST PATH
# # ============================================================
# plt.figure(figsize=(14, 12))
# ax = plt.gca()
#
# nx.draw_networkx_edges(G, pos, alpha=0.15, edge_color="gray")
# nx.draw_networkx_nodes(G, pos, node_color="lightgray", node_size=160, edgecolors="k")
#
# if shortest_path_nodes:
#     nx.draw_networkx_nodes(
#         G, pos,
#         nodelist=shortest_path_nodes,
#         node_color="red",
#         node_size=350,
#         edgecolors="black"
#     )
#     edges_path = list(zip(shortest_path_nodes, shortest_path_nodes[1:]))
#     nx.draw_networkx_edges(G, pos, edgelist=edges_path, edge_color="red", width=3.0)
#
# draw_shifted_labels(G, pos, label_map, ax)
# plt.title(f"Shortest Path: {n_term_label} → {c_term_label}", fontsize=15)
# plt.axis("off")
# plt.tight_layout()
# plt.savefig(os.path.join(output_dir, "graph_shortest_path_spaced.png"), dpi=300)
# plt.close()
#
# print("\n✅ All analysis and visualizations completed successfully!")

##############################################################################################
                 # Find communities and there average degree of freedom
##############################################################################################

# import os
# import csv
# import networkx as nx
# import matplotlib.pyplot as plt
# from matplotlib import colormaps
# from community import community_louvain
#
# # -----------------------------------------
# # FILE PATHS
# # -----------------------------------------
# input_csv = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/4_graph/network_heavyatom.pdb_ringEdges"
# output_dir = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/4_graph/"
#
# os.makedirs(output_dir, exist_ok=True)
#
# # -----------------------------------------
# # LOAD GRAPH
# # -----------------------------------------
# G = nx.Graph()
# with open(input_csv, "r") as f:
#     for row in csv.reader(f):
#         if len(row) >= 2:
#             G.add_edge(row[0].strip(), row[1].strip())
#
# print(f"Loaded graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
#
# # -----------------------------------------
# # COMMUNITY DETECTION
# # -----------------------------------------
# partition = community_louvain.best_partition(G)
#
# communities = {}
# for node, cid in partition.items():
#     communities.setdefault(cid, []).append(node)
#
# print(f"Detected {len(communities)} communities")
#
# # -----------------------------------------
# # AVERAGE DEGREE PER COMMUNITY
# # -----------------------------------------
# avg_degree = {}
# for cid, nodes in communities.items():
#     degs = [G.degree(n) for n in nodes]
#     avg_degree[cid] = sum(degs) / len(degs)
#
# with open(os.path.join(output_dir, "community_avg_degree.txt"), "w") as f:
#     for cid in sorted(avg_degree.keys()):
#         f.write(f"Community {cid} | Nodes: {len(communities[cid])} | Avg Degree: {avg_degree[cid]:.3f}\n")
#
# print("Average degree per community saved.")
#
# # -----------------------------------------
# # NODE LAYOUT — MUCH BETTER SEPARATION
# # -----------------------------------------
# # Increase k to increase spacing; more iterations = smoother layout
# pos = nx.spring_layout(
#     G,
#     k=1.8,              # ⬅ increase spacing between nodes
#     iterations=500,     # ⬅ more stable layout
#     seed=42
# )
# # -----------------------------------------
# # DEGREE DETAILS PER COMMUNITY
# # -----------------------------------------
# degree_details_path = os.path.join(output_dir, "community_degree_details.txt")
#
# with open(degree_details_path, "w") as f:
#     for cid in sorted(communities.keys()):
#         nodes = communities[cid]
#         degrees = {n: G.degree(n) for n in nodes}
#
#         total_degree = sum(degrees.values())
#         avg_degree_val = total_degree / len(nodes)
#
#         f.write(f"===== Community {cid} =====\n")
#         f.write(f"Total Nodes: {len(nodes)}\n")
#         f.write(f"Total Degree: {total_degree}\n")
#         f.write(f"Average Degree: {avg_degree_val:.3f}\n")
#         f.write("Node Degrees:\n")
#
#         # Write each node degree
#         for node, deg in sorted(degrees.items(), key=lambda x: x[0]):
#             f.write(f"  {node}: {deg}\n")
#
#         f.write("\n")
#
# print(f"Degree details per community saved to: {degree_details_path}")
# # -----------------------------------------
# # PLOT WITH COMMUNITY COLORS (now includes avg degree in legend)
# # -----------------------------------------
# plt.figure(figsize=(20, 16))
#
# color_map = colormaps.get_cmap("tab20")   # modern API
# sorted_cids = sorted(communities.keys())  # sorted community IDs
#
# # Draw nodes by community
# for idx, cid in enumerate(sorted_cids):
#     nodes = communities[cid]
#
#     # ✅ NEW: Legend label contains the average degree
#     label = f"Community {cid} (Avg Deg: {avg_degree[cid]:.2f})"
#
#     nx.draw_networkx_nodes(
#         G,
#         pos,
#         nodelist=nodes,
#         node_size=450,
#         node_color=[color_map(idx % 20)],
#         label=label,
#         alpha=0.9
#     )
#
# # Draw edges
# nx.draw_networkx_edges(G, pos, width=0.8, alpha=0.4)
# nx.draw_networkx_labels(G, pos, font_size=7)
#
# # -----------------------------------------
# # FIXED LEGEND — SORTED, SPACED, WITH AVG DEGREE
# # -----------------------------------------
# plt.legend(
#     title="Communities (with Avg Degree)",
#     fontsize=10,
#     title_fontsize=12,
#     loc="center left",
#     bbox_to_anchor=(1.02, 0.5),
#     borderpad=1.2
# )
#
# plt.title("Amino Acid Interaction Network for 6IYZ", fontsize=20)
# plt.axis("off")
#
# output_plot = os.path.join(output_dir, "network_communities.png")
# plt.savefig(output_plot, dpi=300, bbox_inches="tight")
# plt.close()
#
# print(f"Community graph saved at: {output_plot}")
#
#
# ##############################################################################################
# ###############################################################################################
#
# ##############################################################################################
#                  # Find communities and there average degree of freedom
# ##############################################################################################
# # ---------------------------------------------------------
# # FIND MAXIMAL CLIQUES OF SIZE BETWEEN 6 AND 8
# # ---------------------------------------------------------
# all_cliques = list(nx.find_cliques(G))     # All maximal cliques
# cliques = [c for c in all_cliques if 5 <= len(c) <= 6]
#
# print(f"Detected {len(cliques)} maximal cliques of size 6–8")
#
#
#
#
# # ---------------------------------------------------------
# # PLOT MAXIMAL CLIQUES (SIZE 6–8 ONLY)
# # ---------------------------------------------------------
# plt.figure(figsize=(20, 16))
#
# clique_color_map = colormaps.get_cmap("tab20")
#
# legend_handles = []
#
# for idx, clique in enumerate(cliques):
#     color = clique_color_map(idx % 20)
#
#     nx.draw_networkx_nodes(
#         G,
#         pos,
#         nodelist=clique,
#         node_size=500,
#         node_color=[color],
#         alpha=0.85
#     )
#
#     legend_handles.append(
#         plt.Line2D(
#             [0], [0],
#             marker='o',
#             color='white',
#             label=f"Clique {idx} (Size {len(clique)})",
#             markerfacecolor=color,
#             markersize=12
#         )
#     )
#
# nx.draw_networkx_edges(G, pos, alpha=0.3)
# nx.draw_networkx_labels(G, pos, font_size=7)
#
#
# plt.title("Amino Acid Interaction Network for 6IYZ", fontsize=20)
# plt.axis("off")
#
# clique_plot = os.path.join(output_dir, "network_cliques.png")
# plt.savefig(clique_plot, dpi=300, bbox_inches="tight")
# plt.close()
#
# print(f"Clique graph saved at: {clique_plot}")

#################### highest degree for each community of the interaction network ##################
############### this is this small code which will help to calculate the residue with highest degree of freedom for each community######
########## this part of the code is very helpful to identify these residues and how they ####################

# import os
# from Bio.PDB import PDBParser, NeighborSearch
# import networkx as nx
# import matplotlib.pyplot as plt
#
# # ----------------------------
# # USER INPUT
# # ----------------------------
# pdb_path = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/4_graph/6IYX.pdb"
# cutoff = 4.0  # Angstrom cutoff
# output_png = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/4_graph/amino_acid_interaction_network.png"
#
# # ----------------------------
# # One-letter amino acid mapping
# # ----------------------------
# aa3to1 = {
#     'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
#     'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
#     'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
#     'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
#     'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
# }
#
# # ----------------------------
# # Load PDB Structure
# # ----------------------------
# parser = PDBParser(QUIET=True)
# structure = parser.get_structure("protein", pdb_path)
#
# # Extract heavy atoms
# atoms = [atom for atom in structure.get_atoms() if atom.element != 'H']
#
# # Setup neighbor search
# ns = NeighborSearch(atoms)
#
# # Create graph
# G = nx.Graph()
#
# # Store residue → node name mapping
# res_to_node = {}
#
# # ----------------------------
# # Build Nodes
# # ----------------------------
# for model in structure:
#     for chain in model:
#         for residue in chain:
#             if residue.get_resname() in aa3to1:
#                 resname = aa3to1[residue.get_resname()]
#                 resid = residue.get_id()[1]
#                 node_name = f"{resname}:{resid}"
#
#                 G.add_node(node_name, chain=chain.id, res_name=resname, resid=resid)
#                 res_to_node[(chain.id, residue.get_id())] = node_name
#
# # ----------------------------
# # Build Edges (heavy-atom < 4 Å)
# # ----------------------------
# residue_list = list(res_to_node.keys())
#
# for i in range(len(residue_list)):
#     for j in range(i + 1, len(residue_list)):
#         (chainA, idA) = residue_list[i]
#         (chainB, idB) = residue_list[j]
#
#         resA = structure[0][chainA][idA]
#         resB = structure[0][chainB][idB]
#
#         # Get heavy atoms from each residue
#         atomsA = [atom for atom in resA if atom.element != 'H']
#         atomsB = [atom for atom in resB if atom.element != 'H']
#
#         found_contact = False
#         for atomA in atomsA:
#             if found_contact:
#                 break
#             neighbors = ns.search(atomA.coord, cutoff)
#             for atomB in neighbors:
#                 if atomB.get_parent() == resB:
#                     found_contact = True
#                     break
#
#         if found_contact:
#             resA_name = res_to_node[(chainA, idA)]
#             resB_name = res_to_node[(chainB, idB)]
#             G.add_edge(resA_name, resB_name)
#
# # ----------------------------
# # Plot Graph (single color)
# # ----------------------------
# plt.figure(figsize=(12, 10))
# pos = nx.spring_layout(G, seed=42)
#
# nx.draw_networkx_nodes(G, pos, node_size=500, node_color="skyblue")
# nx.draw_networkx_edges(G, pos, alpha=0.4)
# nx.draw_networkx_labels(G, pos, font_size=7)
#
# plt.title("Amino Acid Interaction Network (<4Å heavy-atom contacts)")
# plt.axis("off")
# plt.tight_layout()
#
# plt.savefig(output_png, dpi=600)
# plt.show()
#
# print(f"\nNetwork created successfully!")
# print(f"Nodes: {len(G.nodes)}")
# print(f"Edges: {len(G.edges)}")
# print(f"Graph saved to: {output_png}")

##################################################################################


# import os
# import csv
# import networkx as nx
# import matplotlib.pyplot as plt
# from community import community_louvain
#
# # ---------------------------------------------------
# # FILE PATHS
# # ---------------------------------------------------
# input_csv = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/4_graph/network_heavyatom.pdb_ringEdges"
# output_dir = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/4_graph"
#
# os.makedirs(output_dir, exist_ok=True)
#
# # ---------------------------------------------------
# # LOAD GRAPH
# # ---------------------------------------------------
# G = nx.Graph()
# with open(input_csv, "r") as f:
#     for row in csv.reader(f):
#         if len(row) >= 2:
#             G.add_edge(row[0].strip(), row[1].strip())
#
# print(f"Loaded graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
#
# # ---------------------------------------------------
# # COMMUNITY DETECTION
# # ---------------------------------------------------
# partition = community_louvain.best_partition(G)
#
# communities = {}
# for node, cid in partition.items():
#     communities.setdefault(cid, []).append(node)
#
# print(f"Detected {len(communities)} communities")
#
# # ---------------------------------------------------
# # GLOBAL LAYOUT (same layout for all community plots)
# # ---------------------------------------------------
# pos = nx.spring_layout(G, k=1.5, iterations=400, seed=42)
#
# # ---------------------------------------------------
# # GENERATE SEPARATE FIGURES FOR EACH COMMUNITY
# # ---------------------------------------------------
# for cid, nodes in sorted(communities.items()):
#     subG = G.subgraph(nodes).copy()
#
#     plt.figure(figsize=(12, 10))
#
#     # Draw nodes - dark color
#     nx.draw_networkx_nodes(
#         subG,
#         pos,
#         nodelist=nodes,
#         node_color="#0d1a26",   # Dark navy color
#         node_size=650,
#         alpha=0.95
#     )
#
#     # Draw edges
#     nx.draw_networkx_edges(
#         subG,
#         pos,
#         width=1.3,
#         alpha=0.75,
#         edge_color="gray"
#     )
#
#     # Draw labels with font size 14, bold
#     nx.draw_networkx_labels(
#         subG,
#         pos,
#         font_size=14,
#         font_weight="bold",
#         font_color="black"
#     )
#
#     plt.title(
#         f"Community {cid} — Node Count: {len(nodes)}",
#         fontsize=18,
#         fontweight="bold"
#     )
#     plt.axis("off")
#
#     # Save file
#     out_path = os.path.join(output_dir, f"community_{cid}.png")
#     plt.savefig(out_path, dpi=300, bbox_inches="tight")
#     plt.close()
#
#     print(f"Saved: {out_path}")
#
# print("\nAll community figures generated successfully.")


################################################################################################
################################################################################################

# import os
# from Bio.PDB import PDBParser, NeighborSearch
# import networkx as nx
# import matplotlib.pyplot as plt
# import numpy as np
#
# # ----------------------------
# # USER INPUT
# # ----------------------------
# pdb_path = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/4_graph/6IYX.pdb"
# cutoff = 4.0  # Angstrom cutoff
# output_png = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/4_graph/amino_acid_interaction_network.png"
# atom_output = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/4_graph/interaction_atom_level.txt"
#
# # ----------------------------
# # One-letter amino acid mapping
# # ----------------------------
# aa3to1 = {
#     'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
#     'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
#     'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
#     'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
#     'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
# }
#
# # ----------------------------
# # Load PDB Structure
# # ----------------------------
# parser = PDBParser(QUIET=True)
# structure = parser.get_structure("protein", pdb_path)
#
# # Extract heavy atoms
# atoms = [atom for atom in structure.get_atoms() if atom.element != 'H']
#
# # Setup neighbor search
# ns = NeighborSearch(atoms)
#
# # Create graph
# G = nx.Graph()
#
# # Store residue → node name mapping
# res_to_node = {}
#
# # ----------------------------
# # Build Nodes
# # ----------------------------
# for model in structure:
#     for chain in model:
#         for residue in chain:
#             if residue.get_resname() in aa3to1:
#                 resname = aa3to1[residue.get_resname()]
#                 resid = residue.get_id()[1]
#                 node_name = f"{resname}:{resid}"
#
#                 G.add_node(node_name, chain=chain.id, res_name=resname, resid=resid)
#                 res_to_node[(chain.id, residue.get_id())] = node_name
#
# # ----------------------------
# # Build Edges AND record atom-level contacts
# # ----------------------------
# atom_contacts = []  # store detailed atom-level interaction info
#
# residue_list = list(res_to_node.keys())
#
# for i in range(len(residue_list)):
#     for j in range(i + 1, len(residue_list)):
#         (chainA, idA) = residue_list[i]
#         (chainB, idB) = residue_list[j]
#
#         resA = structure[0][chainA][idA]
#         resB = structure[0][chainB][idB]
#
#         # Heavy atoms only
#         atomsA = [atom for atom in resA if atom.element != 'H']
#         atomsB = [atom for atom in resB if atom.element != 'H']
#
#         has_contact = False
#
#         # Compare all heavy atoms and record closest atomic contacts
#         for atomA in atomsA:
#             for atomB in atomsB:
#                 d = np.linalg.norm(atomA.coord - atomB.coord)
#                 if d < cutoff:
#                     has_contact = True
#                     atom_contacts.append((
#                         res_to_node[(chainA, idA)],
#                         atomA.get_name(),
#                         res_to_node[(chainB, idB)],
#                         atomB.get_name(),
#                         d
#                     ))
#
#         # Add residue-level edge if ANY atom pair < cutoff
#         if has_contact:
#             resA_name = res_to_node[(chainA, idA)]
#             resB_name = res_to_node[(chainB, idB)]
#             G.add_edge(resA_name, resB_name)
#
# # ----------------------------
# # SAVE ATOM-LEVEL CONTACT FILE
# # ----------------------------
# with open(atom_output, "w") as f:
#     for resA, atomA, resB, atomB, dist in atom_contacts:
#         f.write(
#             f"{resA}({atomA})  —  {resB}({atomB})   distance = {dist:.2f} Å\n"
#         )
#
# print(f"\nAtom-level contact info saved to: {atom_output}")
#
# # ----------------------------
# # Plot Graph (single color)
# # ----------------------------
# plt.figure(figsize=(12, 10))
# pos = nx.spring_layout(G, seed=42)
#
# nx.draw_networkx_nodes(G, pos, node_size=500, node_color="skyblue")
# nx.draw_networkx_edges(G, pos, alpha=0.4)
# nx.draw_networkx_labels(G, pos, font_size=7)
#
# plt.title("Amino Acid Interaction Network (<4Å heavy-atom contacts)")
# plt.axis("off")
# plt.tight_layout()
#
# plt.savefig(output_png, dpi=600)
# plt.show()
#
# print(f"\nNetwork created successfully!")
# print(f"Nodes: {len(G.nodes)}")
# print(f"Edges: {len(G.edges)}")
# print(f"Graph saved to: {output_png}")

################### remove i and i+1 and backbone interactions ######################

# import re
# from pathlib import Path
#
# input_path = Path("/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/4_graph/interaction_atom_level.txt")
# output_path = input_path.with_name("interaction_atom_level_filtered.txt")
#
# # Regex: identify atom names inside parentheses and residue numbers
# paren_atom_re = re.compile(r'\(([^\)]+)\)')
# resnum_re = re.compile(r':(\d+)\(')
#
# def norm(atom):
#     """Normalize atom name: uppercase, remove non-alphanumerics."""
#     a = atom.upper()
#     a = re.sub(r'[^A-Z0-9]', '', a)
#     return a
#
# removed = 0
# total = 0
#
# # Backbone atoms
# BACKBONE = {"N", "CA", "C", "O"}
#
# with input_path.open("r") as fin, output_path.open("w") as fout:
#     for line in fin:
#         total += 1
#
#         atoms = paren_atom_re.findall(line)
#         resnums = resnum_re.findall(line)
#
#         # Ensure valid parsing
#         if len(atoms) == 2 and len(resnums) == 2:
#             a1 = norm(atoms[0])
#             a2 = norm(atoms[1])
#
#             r1 = int(resnums[0])
#             r2 = int(resnums[1])
#
#             # ---------------------------------------------------
#             # RULE 1: Remove interactions between residues i & i+1
#             # ---------------------------------------------------
#             if abs(r1 - r2) == 1:
#                 removed += 1
#                 continue
#
#             # ---------------------------------------------------
#             # RULE 2: Remove specific unwanted atom pairs
#             # ---------------------------------------------------
#             remove_pairs = {
#                 ("CA", "CA"),
#                 ("C", "N"), ("N", "C"),
#                 ("CA", "N"), ("N", "CA"),
#                 ("C", "CA"), ("CA", "C"),
#                 ("O", "N"), ("N", "O"),
#                 ("O", "CA"), ("CA", "O"),
#             }
#
#             if (a1, a2) in remove_pairs:
#                 removed += 1
#                 continue
#
#             # ---------------------------------------------------
#             # RULE 3: Remove backbone–backbone interactions
#             # ---------------------------------------------------
#             if a1 in BACKBONE and a2 in BACKBONE:
#                 removed += 1
#                 continue
#
#         # If none of the removal conditions matched → keep line
#         fout.write(line)
#
# print(f"Processed {total} lines.")
# print(f"Removed {removed} lines (i/i+1, unwanted atom pairs, backbone-backbone).")
# print("Filtered output saved as:", output_path)

######################################## remove i and i+1 interactions #################

######################## remove redundant interaction types and change it into a graph ##############################


# import re
# from pathlib import Path
#
# # Input/output files
# input_file = Path("/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/4_graph/interaction_atom_level_filtered.txt")
# output_file = input_file.with_name("residue_level_graph.txt")
#
# # Regex to capture residues (e.g., Q:8)
# residue_pattern = re.compile(r'([A-Z]:\d+)')
#
# pairs = set()  # use a set to avoid duplicates
#
# with input_file.open("r") as f:
#     for line in f:
#         # extract both residues in the line
#         res = residue_pattern.findall(line)
#
#         if len(res) == 2:
#             r1, r2 = res
#
#             # sort so "E:11, Q:8" and "Q:8, E:11" are treated the same
#             pair = tuple(sorted([r1, r2]))
#
#             pairs.add(pair)
#
# # Write output
# with output_file.open("w") as f:
#     for r1, r2 in sorted(pairs):
#         f.write(f"{r1} , {r2}\n")
#
# print(f"✓ Converted to residue-level graph")
# print(f"✓ Unique edges: {len(pairs)}")
# print(f"Saved to: {output_file}")

################## remove node redundancy ######################################
# import re
#
# file_path = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/4_graph/residue_level_graph.txt"
# output_path = file_path.replace(".txt", "_CLEANED.txt")
#
# # Regex for residue like Q:8 or L:211
# res_pattern = re.compile(r"([A-Z]):\s*(\d+)")
#
# clean_edges = set()
# unique_nodes = set()
#
# with open(file_path, "r") as f:
#     for line in f:
#         # Extract both residues
#         matches = res_pattern.findall(line)
#
#         if len(matches) != 2:
#             continue
#
#         (aa1, num1), (aa2, num2) = matches
#
#         # Normalize formatting
#         r1 = f"{aa1}:{int(num1)}"
#         r2 = f"{aa2}:{int(num2)}"
#
#         # build unique edge: sorted ensures no duplication (A,B) == (B,A)
#         edge = tuple(sorted([r1, r2]))
#         clean_edges.add(edge)
#
#         unique_nodes.update([r1, r2])
#
# # Write clean edges
# with open(output_path, "w") as f:
#     for r1, r2 in sorted(clean_edges):
#         f.write(f"{r1},{r2}\n")
#
# print("\n========== CLEANING REPORT ==========")
# print(f"Original file: {file_path}")
# print(f"Cleaned output: {output_path}")
# print(f"Unique residues (nodes): {len(unique_nodes)}")
# print(f"Unique edges: {len(clean_edges)}")
# print("=====================================\n")






############################################# new method for distance calculation #########################################
###########################################################################################################################


# import MDAnalysis as mda
# import numpy as np
# import os
#
# # ==============================================
# # USER INPUT
# # ==============================================
# gro = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/new_5wuc.gro"
# xtc = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/new_5wuc.xtc"
# output_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/new_graph/all_frames_contacts.txt"
#
# # Ensure output folder exists
# os.makedirs(os.path.dirname(output_file), exist_ok=True)
#
# # ==============================================
# # ONE-LETTER AMINO ACID MAP
# # HSD will be converted to H
# # ==============================================
# three_to_one = {
#     'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C',
#     'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'HSD':'H',
#     'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F',
#     'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'
# }
#
# # ==============================================
# # LOAD TRAJECTORY
# # ==============================================
# print("Loading trajectory...")
# u = mda.Universe(gro, xtc)
# residues = u.residues
#
# # Pre-cache atoms per residue for speed
# res_to_atoms = {res.resid: res.atoms for res in residues}
# res_to_resname = {res.resid: res.resname for res in residues}
#
# # ==============================================
# # PARAMETERS
# # ==============================================
# cutoff = 4.0               # Å
# interval_ns = 1.0
# interval_ps = interval_ns * 1000.0
# last_frame_time = -1.0
#
# # ==============================================
# # MAIN PROCESSING
# # ==============================================
# print("Starting calculations...")
#
# with open(output_file, "w") as out:
#     out.write("ALL FRAME CONTACTS (≤4 Å)\n")
#     out.write("=============================================\n")
#
#     for ts in u.trajectory:
#         current_ps = ts.time
#
#         # Only process frames every 1 ns
#         if current_ps - last_frame_time < interval_ps:
#             continue
#
#         last_frame_time = current_ps
#         frame_ns = int(round(current_ps / 1000.0))
#
#         out.write(f"\nFrame : {frame_ns} ns\n")
#         out.write("Res1\tRes2\tAtom1\tAtom2\tDistance\n")
#
#         for i, res1 in enumerate(residues):
#             r1 = res1.resid
#             resname1 = res_to_resname[r1]
#             aa1 = three_to_one.get(resname1, 'X')  # fallback to X if unknown
#             atoms1 = res_to_atoms[r1]
#
#             for j in range(i + 1, len(residues)):
#                 res2 = residues[j]
#                 r2 = res2.resid
#                 resname2 = res_to_resname[r2]
#                 aa2 = three_to_one.get(resname2, 'X')
#                 atoms2 = res_to_atoms[r2]
#
#                 # Exclude only i+1 and i-1
#                 if r2 == r1 + 1 or r2 == r1 - 1:
#                     continue
#
#                 # Vectorized positions
#                 pos1 = atoms1.positions
#                 pos2 = atoms2.positions
#                 diff = pos1[:, None, :] - pos2[None, :, :]
#                 dist2 = np.sum(diff * diff, axis=2)
#                 dist = np.sqrt(dist2)
#
#                 idx1, idx2 = np.where(dist <= cutoff)
#
#                 for a1_idx, a2_idx in zip(idx1, idx2):
#                     atom1 = atoms1[a1_idx]
#                     atom2 = atoms2[a2_idx]
#                     d = dist[a1_idx, a2_idx]
#
#                     out.write(f"{r1}-{aa1}\t{r2}-{aa2}\t{atom1.name}\t{atom2.name}\t{d:.3f}\n")
#
#         print(f"Processed frame {frame_ns} ns")
#
# print(f"\nDone. Output saved to:\n{output_file}")

#########################################################################################

############################## present 50% of time #########################################
############################################################################################


# import os
#
# # ---------------- USER INPUT ----------------
# input_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/new_graph/all_frames_contacts.txt"
# output_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/new_graph/persistent_atom_atom_interactions.txt"
#
# THRESHOLD_FRACTION = 0.5   # 50%
# # -------------------------------------------------
#
#
# # Checks
# if not os.path.exists(input_file):
#     raise SystemExit(f"Input file does not exist:\n{input_file}")
#
#
# print("\nReading contact file...\n")
#
# # Key = (Res1, Res2, Atom1, Atom2)
# # Value = set(frames where present)
# atom_pair_frames = {}
#
# frames_found = 0
# current_frame = None
#
# with open(input_file, "r") as f:
#     for line in f:
#         line = line.strip()
#
#         if not line:
#             continue
#
#         # Detect frame line: "Frame : X ns"
#         if line.startswith("Frame :"):
#             frames_found += 1
#             current_frame = frames_found
#             continue
#
#         # Skip headers and decorations
#         if line.startswith("Res1") or line.startswith("ALL FRAME") or line.startswith("="):
#             continue
#
#         # Parse interaction line:
#         # Res1  Res2  Atom1  Atom2  Distance
#         parts = line.split()
#         if len(parts) < 5:
#             continue
#
#         Res1 = parts[0]      # e.g. "7-L"
#         Res2 = parts[1]      # e.g. "9-L"
#         Atom1 = parts[2]     # e.g. "O"
#         Atom2 = parts[3]     # e.g. "N"
#
#         key = (Res1, Res2, Atom1, Atom2)
#
#         if key not in atom_pair_frames:
#             atom_pair_frames[key] = set()
#
#         atom_pair_frames[key].add(current_frame)
#
#
# # -------------------------------------------------
# # Compute persistent interactions
# # -------------------------------------------------
# print(f"Total frames detected = {frames_found}")
#
# threshold_count = int(frames_found * THRESHOLD_FRACTION)
#
# persistent = []
# for key, frame_set in atom_pair_frames.items():
#     count = len(frame_set)
#     if count >= threshold_count:
#         persistent.append((key, count))
#
# # Sort by number of frames matched
# persistent.sort(key=lambda x: x[1], reverse=True)
#
# # -------------------------------------------------
# # Write results
# # -------------------------------------------------
# os.makedirs(os.path.dirname(output_file), exist_ok=True)
#
# with open(output_file, "w") as out:
#     out.write("PERSISTENT ATOM–ATOM INTERACTIONS (≥50% FRAMES)\n")
#     out.write("=================================================\n")
#     out.write(f"Total frames: {frames_found}\n")
#     out.write(f"Threshold : {THRESHOLD_FRACTION*100:.0f}% → Minimum frames = {threshold_count}\n\n")
#     out.write("Res1\tRes2\tAtom1\tAtom2\tFrames_present\tFraction\n")
#
#     for (Res1, Res2, Atom1, Atom2), count in persistent:
#         out.write(f"{Res1}\t{Res2}\t{Atom1}\t{Atom2}\t{count}\t{count/frames_found:.3f}\n")
#
# print("\nDone! Output written to:")
# print(output_file)

############################################################################################################

################################### remove redundancy ############################################
##################################################################################################
#
# import os
#
# # ------------ USER INPUT ------------------
# input_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/new_graph/persistent_atom_atom_interactions.txt"
# output_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/new_graph/unique_residue_pairs.txt"
# # ------------------------------------------
#
#
# if not os.path.exists(input_file):
#     raise SystemExit(f"Input file does not exist:\n{input_file}")
#
# unique_pairs = set()
#
# with open(input_file, "r") as f:
#     for line in f:
#         line = line.strip()
#         if not line or line.startswith("Res1") or line.startswith("PERSISTENT") or line.startswith("=") or line.startswith("Total"):
#             continue
#
#         parts = line.split()
#         if len(parts) < 2:
#             continue
#
#         Res1 = parts[0]
#         Res2 = parts[1]
#
#         # Normalize pair so A-B == B-A (remove duplicates)
#         if Res1 <= Res2:
#             pair = (Res1, Res2)
#         else:
#             pair = (Res2, Res1)
#
#         unique_pairs.add(pair)
#
# # Sort alphabetically or numerically as needed
# unique_pairs = sorted(unique_pairs, key=lambda x: (x[0], x[1]))
#
# # Write output
# with open(output_file, "w") as out:
#     out.write("UNIQUE RESIDUE–RESIDUE PAIRS\n")
#     out.write("=============================\n\n")
#     out.write("Res1\tRes2\n")
#
#     for Res1, Res2 in unique_pairs:
#         out.write(f"{Res1}\t{Res2}\n")
#
# print("\nDone! Unique residue pairs written to:")
# print(output_file)


#########################################################################################################

# import networkx as nx
# import matplotlib.pyplot as plt
# import os
#
# # ================= USER INPUT =================
# input_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/4_graph/new_result/unique_residue_pairs.txt"
# output_dir = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/4_graph/new_result/"
#
# os.makedirs(output_dir, exist_ok=True)
#
# # =============================================
# # LOAD NETWORK
# # =============================================
# G = nx.Graph()
#
# with open(input_file) as f:
#     for line in f:
#         line = line.strip()
#         if not line or line.startswith("Res") or line.startswith("="):
#             continue
#         r1, r2 = line.split()
#         G.add_edge(r1, r2)
#
# print("Nodes:", G.number_of_nodes())
# print("Edges:", G.number_of_edges())
#
# # =============================================
# # COMMUNITY DETECTION
# # =============================================
# communities = list(nx.algorithms.community.greedy_modularity_communities(G))
#
# # Write clean community list
# community_file = os.path.join(output_dir, "communities_clean.txt")
# with open(community_file, "w") as f:
#     for i, comm in enumerate(communities, 1):
#         f.write(f"Community {i} (size = {len(comm)}):\n")
#         f.write(", ".join(sorted(comm)) + "\n\n")
#
# print("Community list written to:", community_file)
#
# # =============================================
# # IMPROVED LAYOUT (NO OVERLAP)
# # =============================================
#
# # Large k = more repulsion = better spacing
# pos = nx.spring_layout(
#     G,
#     k=1.5,           # increase spacing between nodes
#     iterations=200,  # more relaxation
#     seed=42
# )
#
# # Color mapping per community
# node_color_map = {}
# for i, comm in enumerate(communities):
#     for node in comm:
#         node_color_map[node] = i
#
# node_colors = [node_color_map[n] for n in G.nodes()]
#
# # =============================================
# # PLOT
# # =============================================
# plt.figure(figsize=(14, 12))
#
# nx.draw_networkx_nodes(
#     G, pos,
#     node_size=700,
#     node_color=node_colors,
#     cmap=plt.cm.tab10,
#     alpha=0.9
# )
#
# nx.draw_networkx_edges(
#     G, pos,
#     alpha=0.4,
#     width=1.5
# )
#
# nx.draw_networkx_labels(
#     G, pos,
#     font_size=9,
#     font_weight="bold"
# )
#
# plt.title("Amino Acid Interaction Network – Community Structure", fontsize=14)
# plt.axis("off")
#
# fig_path = os.path.join(output_dir, "interaction_network_communities_no_overlap.png")
# plt.savefig(fig_path, dpi=300, bbox_inches="tight")
# plt.show()
#
# print("Figure saved to:", fig_path)

########################################################################################################################
########################################################################################################################

# import networkx as nx
# import matplotlib.pyplot as plt
# import os
# import numpy as np
#
# # ================= USER INPUT =================
# input_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/new_graph/unique_residue_pairs.txt"
# output_dir = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/new_graph/"
#
# os.makedirs(output_dir, exist_ok=True)
#
# # =============================================
# # LOAD NETWORK
# # =============================================
# G = nx.Graph()
#
# with open(input_file) as f:
#     for line in f:
#         line = line.strip()
#         if not line or line.startswith("Res") or line.startswith("="):
#             continue
#         r1, r2 = line.split()
#         G.add_edge(r1, r2)
#
# print("Nodes:", G.number_of_nodes())
# print("Edges:", G.number_of_edges())
#
# # =============================================
# # DEGREE OF EACH NODE
# # =============================================
# degree_dict = dict(G.degree())
#
# # =============================================
# # COMMUNITY DETECTION
# # =============================================
# communities = list(nx.algorithms.community.greedy_modularity_communities(G))
#
# # =============================================
# # WRITE COMMUNITY LIST + AVERAGE DEGREE
# # =============================================
# community_file = os.path.join(output_dir, "communities_clean.txt")
#
# with open(community_file, "w") as f:
#     f.write("COMMUNITY ANALYSIS\n")
#     f.write("=============================\n\n")
#
#     for i, comm in enumerate(communities, 1):
#         degrees = [degree_dict[node] for node in comm]
#         avg_degree = np.mean(degrees)
#
#         f.write(f"Community {i}\n")
#         f.write(f"Size = {len(comm)}\n")
#         f.write(f"Average degree of freedom = {avg_degree:.2f}\n")
#         f.write("Residues:\n")
#         f.write(", ".join(sorted(comm)) + "\n\n")
#
# print("Community list with average degree written to:", community_file)
#
# # =============================================
# # IMPROVED LAYOUT (NO OVERLAP)
# # =============================================
# pos = nx.spring_layout(
#     G,
#     k=1.5,
#     iterations=200,
#     seed=42
# )
#
# # Color mapping per community
# node_color_map = {}
# for i, comm in enumerate(communities):
#     for node in comm:
#         node_color_map[node] = i
#
# node_colors = [node_color_map[n] for n in G.nodes()]
#
# # =============================================
# # PLOT
# # =============================================
# plt.figure(figsize=(14, 12))
#
# nx.draw_networkx_nodes(
#     G, pos,
#     node_size=700,
#     node_color=node_colors,
#     cmap=plt.cm.tab10,
#     alpha=0.9
# )
#
# nx.draw_networkx_edges(
#     G, pos,
#     alpha=0.4,
#     width=1.5
# )
#
# nx.draw_networkx_labels(
#     G, pos,
#     font_size=9,
#     font_weight="bold"
# )
#
# plt.title("Amino Acid Interaction Network – Community Structure", fontsize=14)
# plt.axis("off")
#
# fig_path = os.path.join(output_dir, "interaction_network_communities_no_overlap.png")
# plt.savefig(fig_path, dpi=300, bbox_inches="tight")
# plt.show()
#
# print("Figure saved to:", fig_path)

######################################## top 5 residues from community with highest degree ##########################
from collections import defaultdict

# ---------------- USER INPUT ----------------
network = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/new_graph/unique_residue_pairs.txt"
community_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/new_graph/communities_clean.txt"
output_file = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wuc/production_run/new_graph/top5_residues_per_community.txt"

TOP_N = 5
# -------------------------------------------


# ---------- STEP 1: Compute degrees ----------
degree = defaultdict(int)

with open(network) as f:
    for line in f:
        if not line.strip():
            continue
        r1, r2 = line.split()
        degree[r1] += 1
        degree[r2] += 1


# ---------- STEP 2: Parse communities ----------
communities = {}
current_comm = None
reading_residues = False

with open(community_file) as f:
    for line in f:
        line = line.strip()

        if line.startswith("Community"):
            current_comm = line
            communities[current_comm] = []
            reading_residues = False

        elif line.startswith("Residues:"):
            reading_residues = True

        elif reading_residues:
            # stop if a new section starts
            if line.startswith("Community") or line == "":
                reading_residues = False
            else:
                residues = [r.strip() for r in line.split(",") if r.strip()]
                communities[current_comm].extend(residues)


# ---------- STEP 3: Top residues per community ----------
with open(output_file, "w") as out:
    out.write("TOP 5 RESIDUES BY DEGREE PER COMMUNITY\n")
    out.write("=====================================\n\n")

    for comm, residues in communities.items():
        res_deg = [(r, degree.get(r, 0)) for r in residues]
        res_deg.sort(key=lambda x: x[1], reverse=True)

        out.write(f"{comm}\n")
        for r, d in res_deg[:TOP_N]:
            out.write(f"  {r}\tdegree = {d}\n")
        out.write("\n")

print(f"✔ Correct output written to:\n{output_file}")


































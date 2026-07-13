#################### analysis of dat file and figure drawing from VMD water analysis ###############
#
# import re
# import os
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# from collections import defaultdict
#
# # -------------------------------------------------------
# # FILE PATHS — change these to match your system
# # -------------------------------------------------------
# input_file = "/media/supremeleader/Pantera/simulation/lipid_simulation/6IZF/charmm-gui-7977903380/gromacs/6izf/water_analysis_per_frame.dat"
# output_dir = os.path.dirname(os.path.abspath(input_file))
#
# residue_file   = os.path.join(output_dir, "residue_water_contacts.dat")
# conserved_file = os.path.join(output_dir, "conserved_interactions.dat")
# summary_file   = os.path.join(output_dir, "analysis_summary.txt")
# plot_png       = os.path.join(output_dir, "water_analysis_publication.png")
# plot_pdf       = os.path.join(output_dir, "water_analysis_publication.pdf")
#
# # -------------------------------------------------------
# # AA three to single letter
# # -------------------------------------------------------
# aa3to1 = {
#     'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
#     'GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
#     'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
#     'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
#     'HSD':'H','HSE':'H','HSP':'H'
# }
#
# # -------------------------------------------------------
# # STEP 1 — Parse water_analysis_per_frame.dat
# # -------------------------------------------------------
# print("Parsing input file...")
# frame_data      = {}
# current_frame   = None
# current_monomer = None
#
# with open(input_file) as f:
#     for line in f:
#         line = line.rstrip()
#
#         frame_match = re.match(r'^Frame (\d+):', line)
#         if frame_match:
#             current_frame = int(frame_match.group(1))
#             frame_data[current_frame] = {}
#             continue
#
#         monomer_match = re.match(r'\s+Monomer (\d+): (\d+) water molecules', line)
#         if monomer_match:
#             current_monomer = int(monomer_match.group(1))
#             frame_data[current_frame][current_monomer] = {
#                 'water_count' : int(monomer_match.group(2)),
#                 'water_list'  : [],
#                 'interactions': []
#             }
#             continue
#
#         wat_match = re.match(r'\s+List of water molecules: (.+)', line)
#         if wat_match and current_frame and current_monomer:
#             wat_str = wat_match.group(1).strip()
#             if wat_str != 'none':
#                 frame_data[current_frame][current_monomer]['water_list'] = [
#                     int(w) for w in wat_str.split()
#                 ]
#             continue
#
#         int_match = re.match(
#             r'\s+ResID_(\d+) \((\w+)\) atom (\w+) -- HOH_(\d+) distance=([\d.]+)A type=(\w+)', line)
#         if int_match and current_frame and current_monomer:
#             frame_data[current_frame][current_monomer]['interactions'].append({
#                 'resid'   : int(int_match.group(1)),
#                 'resname' : int_match.group(2),
#                 'atom'    : int_match.group(3),
#                 'hoh'     : int(int_match.group(4)),
#                 'distance': float(int_match.group(5)),
#                 'type'    : int_match.group(6)
#             })
#
# total_frames   = len(frame_data)
# monomers_found = sorted(set(m for f in frame_data.values() for m in f.keys()))
# print(f"Total frames parsed : {total_frames}")
# print(f"Monomers found      : {monomers_found}")
#
# # -------------------------------------------------------
# # STEP 2 — Water counts per monomer and trimer
# # -------------------------------------------------------
# monomer_counts = {}
# monomer_avg    = {}
# monomer_std    = {}
#
# for monomer in [1, 2, 3]:
#     counts = [frame_data[f][monomer]['water_count']
#               for f in sorted(frame_data) if monomer in frame_data[f]]
#     monomer_counts[monomer] = counts
#     monomer_avg[monomer]    = np.mean(counts)
#     monomer_std[monomer]    = np.std(counts)
#
# trimer_counts = [
#     sum(frame_data[f][m]['water_count'] for m in frame_data[f])
#     for f in sorted(frame_data)
# ]
# avg_trimer = np.mean(trimer_counts)
# std_trimer = np.std(trimer_counts)
#
# # -------------------------------------------------------
# # STEP 3 — Contact count per residue-atom per monomer
# # -------------------------------------------------------
# print("Calculating residue-atom contact percentages...")
# contact_count = defaultdict(lambda: defaultdict(int))
#
# for frame, monomers in frame_data.items():
#     for monomer, data in monomers.items():
#         seen_this_frame = set()
#         for interaction in data['interactions']:
#             key = (interaction['resid'],
#                    interaction['resname'],
#                    interaction['atom'],
#                    interaction['type'])
#             if key not in seen_this_frame:
#                 contact_count[monomer][key] += 1
#                 seen_this_frame.add(key)
#
# # -------------------------------------------------------
# # STEP 4 — Save residue contact file
# # -------------------------------------------------------
# monomer_residue_data = {}
# log = open(summary_file, "w")
#
# def write(msg=""):
#     log.write(msg + "\n")
#
# write("=" * 70)
# write("WATER ANALYSIS SUMMARY")
# write("=" * 70)
# write(f"Input file   : {input_file}")
# write(f"Total frames : {total_frames}")
# write(f"Monomers     : {monomers_found}")
# write()
# write("=" * 70)
# write("ANALYSIS 1 — Average water count for entire trimer")
# write("=" * 70)
# write(f"Average : {avg_trimer:.2f} +/- {std_trimer:.2f}")
# write(f"Min     : {min(trimer_counts)}")
# write(f"Max     : {max(trimer_counts)}")
# write()
# write("=" * 70)
# write("ANALYSIS 2 — Average water count per monomer")
# write("=" * 70)
# for monomer in [1, 2, 3]:
#     write(f"Monomer {monomer}: {monomer_avg[monomer]:.2f} +/- {monomer_std[monomer]:.2f}  "
#           f"(min={min(monomer_counts[monomer])}  max={max(monomer_counts[monomer])})")
# write()
# write("=" * 70)
# write("ANALYSIS 3 — Residue atom contact percentage per monomer")
# write("=" * 70)
#
# outfile_residue = open(residue_file, "w")
# outfile_residue.write("# Residue-atom water contact analysis\n")
# outfile_residue.write(f"# Total frames: {total_frames}\n")
# outfile_residue.write(f"# {'Monomer':<10} {'ResID':<8} {'ResName':<10} {'Atom':<8} "
#                       f"{'Type':<22} {'Count':<8} {'Percentage'}\n")
# outfile_residue.write("=" * 80 + "\n")
#
# for monomer in [1, 2, 3]:
#     frames_with_monomer  = sum(1 for f in frame_data if monomer in frame_data[f])
#     monomer_residue_data[monomer] = {}
#     sorted_contacts = sorted(contact_count[monomer].items(),
#                              key=lambda x: x[1], reverse=True)
#     write(f"\nMonomer {monomer} (frames: {frames_with_monomer}):")
#     write(f"  {'ResID':<8} {'ResName':<10} {'Atom':<8} {'Type':<22} {'Count':<8} %")
#     write(f"  {'-'*65}")
#     outfile_residue.write(f"\nMonomer {monomer}:\n")
#
#     for (resid, resname, atom, itype), count in sorted_contacts:
#         pct = 100.0 * count / frames_with_monomer
#         write(f"  {resid:<8} {resname:<10} {atom:<8} {itype:<22} {count:<8} {pct:.1f}%")
#         outfile_residue.write(
#             f"  Monomer_{monomer}  ResID_{resid:<6} {resname:<10} {atom:<8} "
#             f"{itype:<22} {count:<8} {pct:.1f}%\n"
#         )
#         monomer_residue_data[monomer][(resid, resname, atom)] = pct
#
# outfile_residue.close()
#
# # -------------------------------------------------------
# # STEP 5 — Conserved interactions across all 3 monomers
# # -------------------------------------------------------
# write()
# write("=" * 70)
# write("ANALYSIS 4 — Residue atoms in ALL 3 monomers")
# write("=" * 70)
#
# common_keys = (set(monomer_residue_data[1].keys()) &
#                set(monomer_residue_data[2].keys()) &
#                set(monomer_residue_data[3].keys()))
#
# conserved_list = []
# for (resid, resname, atom) in sorted(common_keys, key=lambda x: x[0]):
#     p1  = monomer_residue_data[1].get((resid, resname, atom), 0)
#     p2  = monomer_residue_data[2].get((resid, resname, atom), 0)
#     p3  = monomer_residue_data[3].get((resid, resname, atom), 0)
#     avg = np.mean([p1, p2, p3])
#     conserved_list.append((resid, resname, atom, p1, p2, p3, avg))
#
# conserved_list.sort(key=lambda x: x[6], reverse=True)
#
# write(f"Total conserved residue-atom pairs: {len(conserved_list)}")
# write()
# write(f"{'ResID':<8} {'ResName':<10} {'Atom':<8} {'Mono1%':<10} {'Mono2%':<10} {'Mono3%':<10} Average%")
# write("-" * 68)
#
# outfile_conserved = open(conserved_file, "w")
# outfile_conserved.write("# Residue atoms interacting with water in ALL 3 monomers\n")
# outfile_conserved.write(f"# Total frames: {total_frames}\n")
# outfile_conserved.write(f"{'ResID':<8} {'ResName':<10} {'Atom':<8} "
#                         f"{'Mono1%':<10} {'Mono2%':<10} {'Mono3%':<10} {'Average%'}\n")
# outfile_conserved.write("=" * 68 + "\n")
#
# for (resid, resname, atom, p1, p2, p3, avg) in conserved_list:
#     line = f"{resid:<8} {resname:<10} {atom:<8} {p1:<10.1f} {p2:<10.1f} {p3:<10.1f} {avg:.1f}%"
#     write(line)
#     outfile_conserved.write(line + "\n")
#
# outfile_conserved.close()
# write()
# write(f"residue_water_contacts.dat  -> {residue_file}")
# write(f"conserved_interactions.dat  -> {conserved_file}")
# log.close()
# print(f"Analysis complete! Summary saved to {summary_file}")
#
# # -------------------------------------------------------
# # STEP 6 — Publication figure
# # -------------------------------------------------------
# print("Generating publication figure...")
#
# top20 = conserved_list[:20]
#
# plt.rcParams.update({
#     'font.family'       : 'DejaVu Sans',
#     'font.size'         : 13,
#     'axes.titlesize'    : 15,
#     'axes.titleweight'  : 'bold',
#     'axes.labelsize'    : 13,
#     'axes.labelweight'  : 'bold',
#     'xtick.labelsize'   : 12,
#     'ytick.labelsize'   : 12,
#     'legend.fontsize'   : 12,
#     'legend.frameon'    : False,
#     'axes.linewidth'    : 1.3,
#     'xtick.major.width' : 1.3,
#     'ytick.major.width' : 1.3,
#     'xtick.major.size'  : 5,
#     'ytick.major.size'  : 5,
#     'axes.spines.top'   : False,
#     'axes.spines.right' : False,
#     'figure.facecolor'  : 'white',
#     'axes.facecolor'    : 'white',
# })
#
# C1   = '#0072B2'
# C2   = '#D55E00'
# C3   = '#009E73'
# CAVG = '#CC0000'
#
# fig = plt.figure(figsize=(18, 13))
# gs  = gridspec.GridSpec(
#     2, 2, figure=fig,
#     hspace=0.38, wspace=0.28,
#     height_ratios=[1, 1.05]
# )
# ax1 = fig.add_subplot(gs[0, 0])
# ax2 = fig.add_subplot(gs[0, 1])
# ax3 = fig.add_subplot(gs[1, :])
#
# frames_x = list(sorted(frame_data.keys()))
# window   = max(1, total_frames // 200)
#
# # --- Plot A ---
# smoothed = np.convolve(trimer_counts, np.ones(window)/window, mode='same')
# ax1.plot(frames_x, trimer_counts, color=C1, linewidth=0.6, alpha=0.35, zorder=2)
# ax1.plot(frames_x, smoothed, color=C1, linewidth=2.0, alpha=0.95,
#          zorder=3, label='Rolling mean')
# ax1.axhline(avg_trimer, color=CAVG, linestyle='--', linewidth=2.0, zorder=4,
#             label=f'Mean = {avg_trimer:.1f} \u00b1 {std_trimer:.1f}')
# ymin_a = min(trimer_counts) * 0.97
# ymax_a = max(trimer_counts) * 1.03
# ax1.set_ylim(ymin_a, ymax_a)
# ax1.set_xlim(frames_x[0], frames_x[-1])
# ax1.fill_between(frames_x, trimer_counts, ymin_a, alpha=0.08, color=C1, zorder=1)
# ax1.set_xlabel('Frame')
# ax1.set_ylabel('Number of water molecules')
# ax1.set_title('Total pore water — trimer')
# ax1.legend(loc='upper right', handlelength=1.8)
# ax1.grid(True, alpha=0.2, linestyle=':', zorder=0)
#
# # --- Plot B ---
# colors_mono = {1: C1, 2: C2, 3: C3}
# for monomer in [1, 2, 3]:
#     cnt     = monomer_counts[monomer]
#     smooth2 = np.convolve(cnt, np.ones(window)/window, mode='same')
#     avg     = monomer_avg[monomer]
#     color   = colors_mono[monomer]
#     ax2.plot(list(range(len(cnt))), cnt,
#              color=color, linewidth=0.5, alpha=0.25, zorder=2)
#     ax2.plot(list(range(len(smooth2))), smooth2,
#              color=color, linewidth=2.0, alpha=0.95, zorder=3,
#              label=f'Monomer {monomer}  (avg {avg:.1f})')
#
# ymin_b = min(min(v) for v in monomer_counts.values()) * 0.96
# ymax_b = max(max(v) for v in monomer_counts.values()) * 1.04
# ax2.set_ylim(ymin_b, ymax_b)
# ax2.set_xlim(0, max(len(v) for v in monomer_counts.values()) - 1)
# ax2.set_xlabel('Frame')
# ax2.set_ylabel('Number of water molecules')
# ax2.set_title('Pore water per monomer per frame')
# ax2.legend(loc='upper right', handlelength=1.8)
# ax2.grid(True, alpha=0.2, linestyle=':', zorder=0)
#
# # --- Plot C ---
# if top20:
#     labels = []
#     for rid, rn, at, p1, p2, p3, avg in top20:
#         s = aa3to1.get(rn.upper(), rn[0])
#         labels.append(f'{s}{rid}\n{at}')
#
#     avgs       = [avg for _, _, _, p1, p2, p3, avg in top20]
#     x          = np.arange(len(labels))
#     norm       = plt.Normalize(vmin=min(avgs)*0.95, vmax=100)
#     cmap       = plt.cm.Blues
#     bar_colors = [cmap(norm(v) * 0.75 + 0.25) for v in avgs]
#
#     bars = ax3.bar(x, avgs, color=bar_colors, edgecolor='#2C4A7C',
#                    linewidth=0.7, width=0.65, zorder=2)
#
#     for bar, val in zip(bars, avgs):
#         ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.6,
#                  f'{val:.1f}%', ha='center', va='bottom',
#                  fontsize=10.5, fontweight='bold', color='#2C4A7C')
#
#     ax3.axhline(50, color='#E06C00', linestyle='--', linewidth=1.6,
#                 alpha=0.8, zorder=3, label='50% threshold')
#     ax3.axhline(80, color='#CC0000', linestyle=':', linewidth=1.6,
#                 alpha=0.8, zorder=3, label='80% threshold')
#     ax3.set_xticks(x)
#     ax3.set_xticklabels(labels, fontsize=12, ha='center', linespacing=1.4)
#     ax3.set_ylabel('Average contact percentage (%)')
#     ax3.set_title('Top 20 conserved residue\u2013atom contacts with water across all 3 monomers')
#     ax3.set_xlim(-0.5, len(labels) - 0.5)
#     ax3.set_ylim(0, 115)
#     ax3.grid(True, alpha=0.2, linestyle=':', axis='y', zorder=0)
#     ax3.legend(loc='lower left', handlelength=2.0)
#
#     sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
#     sm.set_array([])
#     cbar = fig.colorbar(sm, ax=ax3, orientation='vertical',
#                         fraction=0.015, pad=0.01, shrink=0.85)
#     cbar.set_label('Contact %', fontsize=11)
#     cbar.ax.tick_params(labelsize=10)
#
# # --- Panel labels ---
# for ax, label in zip([ax1, ax2, ax3], ['A', 'B', 'C']):
#     ax.text(-0.07, 1.04, label, transform=ax.transAxes,
#             fontsize=20, fontweight='bold', va='top', color='black')
#
# plt.savefig(plot_png, dpi=600, bbox_inches='tight',
#             facecolor='white', format='png')
# plt.savefig(plot_pdf, bbox_inches='tight',
#             facecolor='white', format='pdf')
# plt.close()
#
# print(f"\nAll done! Output files:")
# print(f"  1. {summary_file}")
# print(f"  2. {residue_file}")
# print(f"  3. {conserved_file}")
# print(f"  4. {plot_png}")
# print(f"  5. {plot_pdf}")


###################################### interaction analysis ##########################

# import os
# import re
# from collections import defaultdict
#
# # ─────────────────────────────────────────────
# # File paths
# # ─────────────────────────────────────────────
# CRYSTAL_FILES = {
#     "5wuc": "/media/supremeleader/Pantera/simulation/static_structure_HOH/crystal_water_5wuc_interactions.dat",
#     "5wue": "/media/supremeleader/Pantera/simulation/static_structure_HOH/crystal_water_5wue_interactions.dat",
#     "6iyx": "/media/supremeleader/Pantera/simulation/static_structure_HOH/crystal_water_6iyx_interactions.dat",
#     "6iyz": "/media/supremeleader/Pantera/simulation/static_structure_HOH/crystal_water_6iyz_interactions.dat",
# }
#
# SIM_FILES = {
#     "5wuc": "/media/supremeleader/Pantera/simulation/static_structure_HOH/simulated_conserved_5wuc_interactions.dat",
#     "5wue": "/media/supremeleader/Pantera/simulation/static_structure_HOH/simulated_conserved_5wue_interactions.dat",
#     "6iyx": "/media/supremeleader/Pantera/simulation/static_structure_HOH/simulated_conserved_6iyx_interactions.dat",
#     "6iyz": "/media/supremeleader/Pantera/simulation/static_structure_HOH/simulated_conserved_6iyz_interactions.dat",
# }
#
# THRESHOLD = 50.0  # % threshold for simulated interactions
#
# # ─────────────────────────────────────────────
# # Parsers
# # ─────────────────────────────────────────────
#
# def parse_crystal(filepath):
#     interactions = {}
#     with open(filepath) as f:
#         for line in f:
#             line = line.strip()
#             if not line or line.startswith("#"):
#                 continue
#             parts = line.split()
#             if len(parts) < 6:
#                 continue
#             resid_raw = parts[0]          # e.g. ResID_45
#             resname   = parts[1].upper()  # e.g. THR
#             atom      = parts[2].upper()  # e.g. O
#             itype     = parts[5]          # e.g. Backbone_HBond
#
#             # extract numeric resid
#             match = re.search(r"(\d+)", resid_raw)
#             if not match:
#                 continue
#             resid = int(match.group(1))
#
#             key = (resid, resname, atom)
#             interactions[key] = itype
#     return interactions
#
#
# def parse_simulated(filepath, threshold=THRESHOLD):
#
#     all_sim = {}
#     conserved_sim = {}
#
#     with open(filepath) as f:
#         for line in f:
#             line = line.strip()
#             if not line or line.startswith("#"):
#                 continue
#             # skip header line
#             if line.startswith("ResID") or line.startswith("==="):
#                 continue
#             parts = line.split()
#             if len(parts) < 7:
#                 continue
#             try:
#                 resid   = int(parts[0])
#                 resname = parts[1].upper()
#                 atom    = parts[2].upper()
#                 mono1   = float(parts[3])
#                 mono2   = float(parts[4])
#                 mono3   = float(parts[5])
#                 avg_str = parts[6].replace("%", "")
#                 avg     = float(avg_str)
#             except ValueError:
#                 continue
#
#             key = (resid, resname, atom)
#             all_sim[key] = (mono1, mono2, mono3, avg)
#
#             # conserved if average > threshold OR any single monomer > threshold
#             if avg > threshold or mono1 > threshold or mono2 > threshold or mono3 > threshold:
#                 conserved_sim[key] = (mono1, mono2, mono3, avg)
#
#     return all_sim, conserved_sim
#
#
# # ─────────────────────────────────────────────
# # Comparison
# # ─────────────────────────────────────────────
#
# def compare(crystal_interactions, conserved_sim, all_sim):
#
#     crystal_keys = set(crystal_interactions.keys())
#     sim_keys     = set(conserved_sim.keys())
#
#     common       = crystal_keys & sim_keys
#     only_crystal = crystal_keys - sim_keys
#     only_sim     = sim_keys     - crystal_keys
#
#     return common, only_crystal, only_sim
#
#
# # ─────────────────────────────────────────────
# # Reporting
# # ─────────────────────────────────────────────
#
# def format_sim_pct(vals):
#     m1, m2, m3, avg = vals
#     return f"Mono1={m1:.1f}%  Mono2={m2:.1f}%  Mono3={m3:.1f}%  Avg={avg:.1f}%"
#
#
# def write_report(structure_id, crystal_interactions, conserved_sim, all_sim, outfile):
#     common, only_crystal, only_sim = compare(crystal_interactions, conserved_sim, all_sim)
#
#     lines = []
#     sep  = "=" * 80
#     sep2 = "-" * 80
#
#     lines.append(sep)
#     lines.append(f"  STRUCTURE: {structure_id.upper()}")
#     lines.append(f"  Crystal interactions   : {len(crystal_interactions)}")
#     lines.append(f"  Simulated conserved (>{THRESHOLD}% threshold): {len(conserved_sim)}")
#     lines.append(sep)
#
#     # ── 1) COMMON ────────────────────────────────────────────────
#     lines.append(f"\n{'─'*80}")
#     lines.append(f"  [1] COMMON INTERACTIONS  (in crystal AND in simulation >{THRESHOLD}%)")
#     lines.append(f"      Count: {len(common)}")
#     lines.append(f"{'─'*80}")
#     lines.append(f"  {'ResID':<8} {'ResName':<9} {'Atom':<8} {'Crystal_Type':<22} {'Simulation %'}")
#     lines.append(f"  {'-'*7} {'-'*8} {'-'*7} {'-'*21} {'-'*46}")
#     for key in sorted(common):
#         resid, resname, atom = key
#         itype = crystal_interactions[key]
#         sim_v = all_sim.get(key, (0,0,0,0))
#         lines.append(f"  {resid:<8} {resname:<9} {atom:<8} {itype:<22} {format_sim_pct(sim_v)}")
#
#     # ── 2) ONLY IN CRYSTAL (not conserved in simulation) ─────────
#     lines.append(f"\n{'─'*80}")
#     lines.append(f"  [2] ONLY IN CRYSTAL  (not conserved in simulation >{THRESHOLD}%)")
#     lines.append(f"      Count: {len(only_crystal)}")
#     lines.append(f"{'─'*80}")
#     lines.append(f"  {'ResID':<8} {'ResName':<9} {'Atom':<8} {'Crystal_Type':<22} {'Simulation % (if any)'}")
#     lines.append(f"  {'-'*7} {'-'*8} {'-'*7} {'-'*21} {'-'*46}")
#     for key in sorted(only_crystal):
#         resid, resname, atom = key
#         itype = crystal_interactions[key]
#         if key in all_sim:
#             sim_v = all_sim[key]
#             sim_str = format_sim_pct(sim_v) + "  [below threshold]"
#         else:
#             sim_str = "NOT FOUND in simulation"
#         lines.append(f"  {resid:<8} {resname:<9} {atom:<8} {itype:<22} {sim_str}")
#
#     # ── 3) ONLY IN SIMULATION (not in crystal) ───────────────────
#     lines.append(f"\n{'─'*80}")
#     lines.append(f"  [3] ONLY IN SIMULATION  (conserved >{THRESHOLD}% but absent from crystal)")
#     lines.append(f"      Count: {len(only_sim)}")
#     lines.append(f"{'─'*80}")
#     lines.append(f"  {'ResID':<8} {'ResName':<9} {'Atom':<8} {'Simulation %'}")
#     lines.append(f"  {'-'*7} {'-'*8} {'-'*7} {'-'*46}")
#     for key in sorted(only_sim):
#         resid, resname, atom = key
#         sim_v = all_sim[key]
#         lines.append(f"  {resid:<8} {resname:<9} {atom:<8} {format_sim_pct(sim_v)}")
#
#     lines.append("")
#
#     text = "\n".join(lines)
#     outfile.write(text + "\n")
#     return text
#
#
# # ─────────────────────────────────────────────
# # Summary table: reads output file and appends
# # ─────────────────────────────────────────────
#
# def parse_and_append_summary(output_path, threshold):
#     """
#     Reads the already-written output file, extracts per-structure counts
#     from the lines written by write_report(), then appends a summary table.
#
#     Lines it looks for (written by write_report):
#         "  STRUCTURE: 5WUC"
#         "  Crystal interactions   : 40"
#         "      Count: 23"   ← appears 3 times per structure:
#                                first  = common
#                                second = only_crystal
#                                third  = only_sim
#     """
#
#     with open(output_path, "r") as f:
#         content = f.read()
#
#     # Split into per-structure blocks on the STRUCTURE: marker
#     # Each block starts at "  STRUCTURE: XXXX"
#     blocks = re.split(r"(?=\s{2}STRUCTURE:\s)", content)
#
#     rows = []  # list of (sid, crystal, common, only_crystal, only_sim)
#
#     for block in blocks:
#         sid_match = re.search(r"STRUCTURE:\s+(\S+)", block)
#         if not sid_match:
#             continue
#         sid = sid_match.group(1).upper()
#
#         # Crystal interactions count
#         cry_match = re.search(r"Crystal interactions\s*:\s*(\d+)", block)
#         if not cry_match:
#             continue
#         crystal_total = int(cry_match.group(1))
#
#         # The three "Count: N" lines appear in order: common, only_crystal, only_sim
#         counts = re.findall(r"Count:\s*(\d+)", block)
#         if len(counts) < 3:
#             continue
#         common_n    = int(counts[0])
#         only_cry_n  = int(counts[1])
#         only_sim_n  = int(counts[2])
#
#         rows.append((sid, crystal_total, common_n, only_cry_n, only_sim_n))
#
#     if not rows:
#         print("[WARN] Summary parser found no structure blocks — table not written.")
#         return
#
#     # Build the table
#     sep  = "=" * 65
#     sep2 = "-" * 65
#     lines = [
#         "",
#         sep,
#         "  QUICK SUMMARY OF RESULTS",
#         sep,
#         f"  {'Structure':<15} {'Crystal':>10} {'Common':>10} {'Only Crystal':>14} {'Only Simulation':>17}",
#         f"  {'-'*13} {'-'*10} {'-'*10} {'-'*14} {'-'*17}",
#     ]
#     for sid, cry_total, common_n, only_cry_n, only_sim_n in rows:
#         lines.append(
#             f"  {sid:<15} {cry_total:>10} {common_n:>10} {only_cry_n:>14} {only_sim_n:>17}"
#         )
#     lines += [
#         sep,
#         f"  Threshold: >{threshold}% (Average OR any single monomer) to be considered conserved in simulation",
#         sep,
#         "",
#     ]
#
#     table_text = "\n".join(lines)
#
#     # Append to the file
#     with open(output_path, "a") as f:
#         f.write(table_text + "\n")
#
#     print(table_text)
#     print(f"✅  Summary table appended to: {output_path}")
#
#
# # ─────────────────────────────────────────────
# # Main
# # ─────────────────────────────────────────────
#
# def main():
#     output_path = "/media/supremeleader/Pantera/simulation/static_structure_HOH/water_interaction_comparison.txt"
#
#     with open(output_path, "w") as outf:
#         header = (
#             "WATER INTERACTION COMPARISON: Crystal vs Simulated\n"
#             f"Threshold for simulated conservation: >{THRESHOLD}% in Average OR any monomer\n"
#             "Crystal = monomer | Simulated = trimer (3 monomers)\n"
#             "Matching key: (ResID, ResName, Atom)  [water identity ignored]\n\n"
#         )
#         outf.write(header)
#         print(header, end="")
#
#         for sid in ["5wuc", "5wue", "6iyx", "6iyz"]:
#             cry_file = CRYSTAL_FILES[sid]
#             sim_file = SIM_FILES[sid]
#
#             if not os.path.exists(cry_file):
#                 print(f"[WARN] Missing crystal file: {cry_file}")
#                 continue
#             if not os.path.exists(sim_file):
#                 print(f"[WARN] Missing simulated file: {sim_file}")
#                 continue
#
#             crystal_interactions = parse_crystal(cry_file)
#             all_sim, conserved_sim = parse_simulated(sim_file, THRESHOLD)
#
#             report = write_report(sid, crystal_interactions, conserved_sim, all_sim, outf)
#             print(report)
#
#     # ── After file is fully written, read it back and append summary ──
#     parse_and_append_summary(output_path, THRESHOLD)
#
#
# if __name__ == "__main__":
#     main()

############################ comparision between structures ######################

########################## WATER MOLECULES AT TRIMER INTERFACE #####################
# import numpy as np
#
# water_counts = []
# frames = []
#
# with open("/media/supremeleader/Pantera/simulation/2024_simulation_analysis/5wue/trimer_water.dat", "r") as f:
#     for line in f:
#         if line.startswith("Frame") and "Trimer interface" in line:
#             frame = int(line.split("Frame")[1].split(":")[0].strip())
#             count = int(line.split("Waters=")[1].split("|")[0].strip())
#             water_counts.append(count)
#             frames.append(frame)
#
# water_counts = np.array(water_counts)
# frames = np.array(frames)
#
# max_idx = np.argmax(water_counts)
#
# print(f"Total frames analyzed  : {len(water_counts)}")
# print(f"Average water molecules: {np.mean(water_counts):.2f}")
# print(f"Standard deviation     : {np.std(water_counts):.2f}")
# print(f"Min                    : {np.min(water_counts)}")
# print(f"Max                    : {np.max(water_counts)}")
# print(f"Max water frame        : Frame {frames[max_idx]} with {water_counts[max_idx]} waters")


# ###################### water coordinating mimicing K+ or Na+ coordination ####################
# #################################################################################
# import numpy as np
# import MDAnalysis as mda
# from MDAnalysis.analysis import distances
# from collections import defaultdict
# import matplotlib.pyplot as plt
# import os
#
# # ══════════════════════════════════════════════════════════════
# #  TRIC Channel — Chemically Rigorous Ion Binding Site Analysis
# #  pH 7 corrected — only chemically valid coordinating atoms
# #  K+ range: 2.7-3.3 Å | Na+ range: 2.3-2.7 Å
# # ══════════════════════════════════════════════════════════════
#
# # ── Output directory ──────────────────────────────────────────
# outdir = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/coordination"
# os.makedirs(outdir, exist_ok=True)
# print(f"Output directory: {outdir}")
#
# # ── Load trajectory ───────────────────────────────────────────
# u = mda.Universe(
#     "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/aligned.gro",
#     "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/aligned.xtc"
# )
#
# # ── Monomer boundaries ────────────────────────────────────────
# prot     = u.select_atoms("protein")
# per_mono = len(prot) // 3
# frags = {
#     0: u.select_atoms(f"protein and index {prot.indices[0]}:{prot.indices[per_mono-1]}"),
#     1: u.select_atoms(f"protein and index {prot.indices[per_mono]}:{prot.indices[2*per_mono-1]}"),
#     2: u.select_atoms(f"protein and index {prot.indices[2*per_mono]}:{prot.indices[3*per_mono-1]}")
# }
# waters = u.select_atoms("resname TIP3 and name OH2")
#
# print(f"Protein atoms    : {len(prot)}")
# print(f"Atoms per monomer: {per_mono}")
# print(f"TIP3 waters      : {len(waters)}")
#
# # ══════════════════════════════════════════════════════════════
# #  CHEMICALLY VALID COORDINATING ATOMS AT pH 7
# #  Based on literature:
# #  - Backbone N (amide) → EXCLUDED — protonated at pH 7, no lone pair
# #  - LYS NZ             → EXCLUDED — positively charged at pH 7
# #  - ARG NE/NH1/NH2     → EXCLUDED — positively charged at pH 7
# #  - ASN ND2, GLN NE2   → EXCLUDED — too weak for ion coordination
# #  - HSD NE2            → INCLUDED — unprotonated N at pH 7 (CHARMM HSD)
# #  - HSE ND1            → INCLUDED — unprotonated N at pH 7 (CHARMM HSE)
# #  - All backbone O     → INCLUDED — primary K+ coordination (KcsA)
# #  - All sidechain O    → INCLUDED — important for Na+ coordination
# # ══════════════════════════════════════════════════════════════
#
# # Backbone carbonyl oxygen — PRIMARY coordination for K+
# MC_O = {"O"}
#
# # Sidechain oxygens — valid at pH 7
# SC_O = {
#     "OG",              # SER — hydroxyl
#     "OG1",             # THR — hydroxyl
#     "OD1", "OD2",      # ASP — carboxylate (negatively charged at pH 7) ✅ strong
#     "OE1", "OE2",      # GLU — carboxylate (negatively charged at pH 7) ✅ strong
#     "OD1",             # ASN — amide carbonyl O
#     "OE1",             # GLN — amide carbonyl O
#     "OH",              # TYR — phenolic hydroxyl
# }
#
# # Histidine N — ONLY unprotonated N is valid
# # In CHARMM force field:
# #   HSD = ND1 protonated, NE2 FREE → NE2 coordinates
# #   HSE = NE2 protonated, ND1 FREE → ND1 coordinates
# #   HSP = both protonated → NEITHER coordinates
# HIS_N_MAP = {
#     "HSD": "NE2",   # unprotonated N in HSD
#     "HSE": "ND1",   # unprotonated N in HSE
# }
#
# # Combined valid atoms (O only — His N handled separately)
# ALL_VALID_O = MC_O | SC_O
#
# # ── Parameters ────────────────────────────────────────────────
# # K+ coordination range (KcsA crystal structure: 2.7-3.1 Å avg 2.85 Å)
# K_MIN        = 2.7
# K_MAX        = 3.3
#
# # Na+ coordination range (NavAb crystal structure: ~2.3-2.5 Å)
# NA_MIN       = 2.3
# NA_MAX       = 2.7
#
# # Water-water H-bond cutoff (unchanged)
# WAT_WAT_CUTOFF = 3.5
#
# # Cylinder parameters
# CYL_RADIUS   = 8.0
# CYL_RADIUS2  = CYL_RADIUS ** 2
# Z_HALF       = 11.0
#
# # Minimum protein atoms required for ion-like coordination
# # (at least 2 protein O atoms must be within range)
# MIN_PROT_COORD = 2
#
# # K+ coordination number range
# K_COORD_MIN  = 6
# K_COORD_MAX  = 8
#
# # Na+ coordination number range
# NA_COORD_MIN = 4
# NA_COORD_MAX = 6
#
# STRIDE       = 50
#
# print(f"\nCylinder       : radius={CYL_RADIUS}Å  Z=±{Z_HALF}Å")
# print(f"K+  range      : {K_MIN}-{K_MAX}Å  coord={K_COORD_MIN}-{K_COORD_MAX}")
# print(f"Na+ range      : {NA_MIN}-{NA_MAX}Å  coord={NA_COORD_MIN}-{NA_COORD_MAX}")
# print(f"Water-water    : <{WAT_WAT_CUTOFF}Å")
# print(f"Min prot atoms : {MIN_PROT_COORD}")
# print(f"Stride         : every {STRIDE} frames")
# print(f"Total frames   : {len(u.trajectory)}")
# print(f"Frames analyzed: {len(u.trajectory)//STRIDE}")
#
# # ── Output files ──────────────────────────────────────────────
# k_out   = open(f"{outdir}/K_coordination.dat",   "w")
# na_out  = open(f"{outdir}/Na_coordination.dat",  "w")
# all_out = open(f"{outdir}/all_coordination.dat", "w")
# sum_out = open(f"{outdir}/coordination_summary.dat", "w")
#
# header = (
#     "# Frame  Frag  WaterResID  "
#     "ProtO_in_range  WatWat  TotalCoord  "
#     "IonType  Coordinating_atoms\n"
# )
# k_out.write(  f"# K+  coordination events (prot O: {K_MIN}-{K_MAX}Å, total coord: {K_COORD_MIN}-{K_COORD_MAX})\n" + header)
# na_out.write( f"# Na+ coordination events (prot O: {NA_MIN}-{NA_MAX}Å, total coord: {NA_COORD_MIN}-{NA_COORD_MAX})\n" + header)
# all_out.write(f"# All coordination events\n" + header)
#
# # ── Storage ───────────────────────────────────────────────────
# k_events  = []
# na_events = []
# all_events = []
#
# k_dists   = []   # all K+ range distances
# na_dists  = []   # all Na+ range distances
# ww_dists  = []   # all water-water distances
#
# residue_k_count  = defaultdict(int)
# residue_na_count = defaultdict(int)
#
# total_frames = 0
#
# print("\n══════════════════════════════════════════════")
# print("  Starting chemically rigorous analysis...")
# print("══════════════════════════════════════════════")
#
# for i, ts in enumerate(u.trajectory[::STRIDE]):
#     frame = ts.frame
#
#     if i % 100 == 0:
#         print(f"  Frame {frame}/{len(u.trajectory)}  ({i*STRIDE*100//len(u.trajectory)}% done)")
#
#     total_frames += 1
#
#     # ── Build list of valid coordinating protein atoms ────────
#     # Do this once per frame for efficiency
#     coord_pos      = []
#     coord_names    = []
#     coord_resids   = []
#     coord_resnames = []
#     coord_types    = []   # "MC-O", "SC-O", "HIS-N"
#
#     for fid, fg in frags.items():
#         for atom in fg.atoms:
#             include = False
#             atype   = None
#
#             # Check oxygen atoms
#             if atom.name in ALL_VALID_O:
#                 include = True
#                 atype   = "MC-O" if atom.name == "O" else "SC-O"
#
#             # Check Histidine N — only unprotonated N
#             elif atom.resname in HIS_N_MAP:
#                 valid_n = HIS_N_MAP[atom.resname]
#                 if atom.name == valid_n:
#                     include = True
#                     atype   = "HIS-N"
#
#             if include:
#                 coord_pos.append(atom.position)
#                 coord_names.append(atom.name)
#                 coord_resids.append(atom.resid)
#                 coord_resnames.append(atom.resname)
#                 coord_types.append(atype)
#
#     coord_pos = np.array(coord_pos)
#
#     # ── Per monomer cylinder ──────────────────────────────────
#     for frag_id, frag in frags.items():
#
#         # Cylinder center from this monomer
#         com        = frag.center_of_mass()
#         cx, cy, cz = com
#
#         # Waters inside cylinder
#         wp      = waters.positions
#         in_cyl  = ((wp[:,0]-cx)**2 + (wp[:,1]-cy)**2) <= CYL_RADIUS2
#         in_z    = (wp[:,2] >= cz - Z_HALF) & (wp[:,2] <= cz + Z_HALF)
#         cyl_wat = waters[in_cyl & in_z]
#
#         if len(cyl_wat) == 0:
#             continue
#
#         # Distance matrix: cylinder waters vs valid protein atoms
#         d_prot = distances.distance_array(cyl_wat.positions, coord_pos)
#
#         # Distance matrix: cylinder waters vs all waters
#         d_wat  = distances.distance_array(cyl_wat.positions, waters.positions)
#
#         # ── Per water analysis ────────────────────────────────
#         for wi in range(len(cyl_wat)):
#             wat_resid = cyl_wat[wi].resid
#
#             # ── Count protein atoms in K+ range ──────────────
#             k_prot_count  = 0
#             k_prot_res    = []
#
#             # ── Count protein atoms in Na+ range ─────────────
#             na_prot_count = 0
#             na_prot_res   = []
#
#             for pi, (d, name, resid, resname, atype) in enumerate(zip(
#                     d_prot[wi], coord_names, coord_resids,
#                     coord_resnames, coord_types)):
#
#                 # K+ range: 2.7 - 3.3 Å
#                 if K_MIN <= d <= K_MAX:
#                     k_prot_count += 1
#                     k_dists.append(d)
#                     k_prot_res.append(f"{resname}{resid}-{name}({atype})")
#                     residue_k_count[resid] += 1
#
#                 # Na+ range: 2.3 - 2.7 Å
#                 if NA_MIN <= d <= NA_MAX:
#                     na_prot_count += 1
#                     na_dists.append(d)
#                     na_prot_res.append(f"{resname}{resid}-{name}({atype})")
#                     residue_na_count[resid] += 1
#
#             # ── Water-water coordination (< 3.5 Å) ───────────
#             ww_count  = 0
#             ww_res    = []
#
#             for wj in range(len(waters)):
#                 if waters[wj].resid == wat_resid:
#                     continue
#                 if d_wat[wi][wj] < WAT_WAT_CUTOFF:
#                     ww_count += 1
#                     ww_dists.append(d_wat[wi][wj])
#                     ww_res.append(f"TIP3{waters[wj].resid}(WAT)")
#
#             # ── Total coordination numbers ────────────────────
#             k_total  = k_prot_count  + ww_count
#             na_total = na_prot_count + ww_count
#
#             # ── Ion type classification ───────────────────────
#             # STRICT criteria:
#             # K+:  >= MIN_PROT_COORD protein O in K+ range
#             #      AND total coord in K+ range (6-8)
#             # Na+: >= MIN_PROT_COORD protein O in Na+ range
#             #      AND total coord in Na+ range (4-6)
#             # Note: Na+ range is subset of K+ range so check Na+ first
#
#             is_na = (na_prot_count >= MIN_PROT_COORD and
#                      NA_COORD_MIN <= na_total <= NA_COORD_MAX)
#
#             is_k  = (k_prot_count  >= MIN_PROT_COORD and
#                      K_COORD_MIN  <= k_total  <= K_COORD_MAX)
#
#             if not is_na and not is_k:
#                 continue
#
#             # Determine ion type
#             if is_na and is_k:
#                 ion_type = "K+/Na+"   # ambiguous — overlapping range
#             elif is_k:
#                 ion_type = "K+"
#             else:
#                 ion_type = "Na+"
#
#             # Build coordinating atom string
#             if ion_type in ("K+", "K+/Na+"):
#                 coord_str = "|".join((k_prot_res + ww_res)[:12])
#                 prot_count = k_prot_count
#                 total_coord = k_total
#             else:
#                 coord_str = "|".join((na_prot_res + ww_res)[:12])
#                 prot_count = na_prot_count
#                 total_coord = na_total
#
#             line = (
#                 f"{frame:<8} {frag_id:<6} {wat_resid:<12} "
#                 f"{prot_count:<16} {ww_count:<8} {total_coord:<12} "
#                 f"{ion_type:<12} {coord_str}\n"
#             )
#
#             all_out.write(line)
#
#             event = {
#                 "frame":      frame,
#                 "frag":       frag_id,
#                 "wat":        wat_resid,
#                 "k_prot":     k_prot_count,
#                 "na_prot":    na_prot_count,
#                 "ww":         ww_count,
#                 "k_total":    k_total,
#                 "na_total":   na_total,
#                 "ion_type":   ion_type,
#                 "k_res":      k_prot_res,
#                 "na_res":     na_prot_res,
#                 "ww_res":     ww_res
#             }
#
#             if ion_type == "K+":
#                 k_out.write(line)
#                 k_events.append(event)
#             elif ion_type == "Na+":
#                 na_out.write(line)
#                 na_events.append(event)
#             else:  # K+/Na+ ambiguous
#                 k_out.write(line)
#                 na_out.write(line)
#                 k_events.append(event)
#                 na_events.append(event)
#
#             all_events.append(event)
#
# # ── Close files ───────────────────────────────────────────────
# k_out.close()
# na_out.close()
# all_out.close()
#
# print(f"\n✅ Analysis complete!")
# print(f"   Frames analyzed : {total_frames}")
# print(f"   K+  events      : {len(k_events)}")
# print(f"   Na+ events      : {len(na_events)}")
#
# # ── Summary ───────────────────────────────────────────────────
# resname_map = {atom.resid: atom.resname for atom in prot}
#
# summary = []
# summary.append("═"*70)
# summary.append("  TRIC CHANNEL — ION BINDING SITE ANALYSIS (pH 7 corrected)")
# summary.append("═"*70)
# summary.append(f"  Frames analyzed        : {total_frames}")
# summary.append(f"  Stride                 : every {STRIDE} frames")
# summary.append(f"  K+  events (coord {K_COORD_MIN}-{K_COORD_MAX}) : {len(k_events)}")
# summary.append(f"  Na+ events (coord {NA_COORD_MIN}-{NA_COORD_MAX}) : {len(na_events)}")
# summary.append(f"\n  Chemistry corrections applied:")
# summary.append(f"    ❌ Backbone N removed  (protonated at pH 7)")
# summary.append(f"    ❌ LYS NZ removed      (positively charged at pH 7)")
# summary.append(f"    ❌ ARG NE/NH1/NH2 removed (positively charged at pH 7)")
# summary.append(f"    ✅ HSD NE2 kept        (unprotonated at pH 7)")
# summary.append(f"    ✅ HSE ND1 kept        (unprotonated at pH 7)")
# summary.append(f"    ✅ All backbone O kept  (primary K+ coordination)")
# summary.append(f"    ✅ All sidechain O kept (Na+ coordination)")
#
# # Per monomer
# summary.append(f"\n{'═'*70}")
# summary.append("  PER MONOMER SUMMARY")
# summary.append("═"*70)
#
# for fid in range(3):
#     k_frag  = [e for e in k_events  if e["frag"] == fid]
#     na_frag = [e for e in na_events if e["frag"] == fid]
#
#     summary.append(f"\n  Monomer {fid+1} (Fragment {fid}):")
#     summary.append(f"    K+  events : {len(k_frag)}")
#     summary.append(f"    Na+ events : {len(na_frag)}")
#
#     if k_frag:
#         rc = defaultdict(int)
#         for e in k_frag:
#             for r in e["k_res"]:
#                 rc[r.split("(")[0]] += 1
#         top = sorted(rc.items(), key=lambda x: x[1], reverse=True)[:5]
#         summary.append(f"    Top K+  residues: {[r[0] for r in top]}")
#
#     if na_frag:
#         rc = defaultdict(int)
#         for e in na_frag:
#             for r in e["na_res"]:
#                 rc[r.split("(")[0]] += 1
#         top = sorted(rc.items(), key=lambda x: x[1], reverse=True)[:5]
#         summary.append(f"    Top Na+ residues: {[r[0] for r in top]}")
#
# # Top K+ residues overall
# summary.append(f"\n{'═'*70}")
# summary.append("  TOP 15 RESIDUES IN K+ COORDINATION SITES")
# summary.append("═"*70)
# summary.append(f"{'ResID':<8} {'ResName':<8} {'Count':<10} {'Freq%'}")
# summary.append("─"*40)
# sorted_k = sorted(residue_k_count.items(), key=lambda x: x[1], reverse=True)[:15]
# for resid, count in sorted_k:
#     pct = 100.0 * count / (total_frames * 3)
#     summary.append(f"{resid:<8} {resname_map.get(resid,'???'):<8} {count:<10} {pct:.2f}%")
#
# # Top Na+ residues overall
# summary.append(f"\n{'═'*70}")
# summary.append("  TOP 15 RESIDUES IN Na+ COORDINATION SITES")
# summary.append("═"*70)
# summary.append(f"{'ResID':<8} {'ResName':<8} {'Count':<10} {'Freq%'}")
# summary.append("─"*40)
# sorted_na = sorted(residue_na_count.items(), key=lambda x: x[1], reverse=True)[:15]
# for resid, count in sorted_na:
#     pct = 100.0 * count / (total_frames * 3)
#     summary.append(f"{resid:<8} {resname_map.get(resid,'???'):<8} {count:<10} {pct:.2f}%")
#
# # Top 20 highest coordination events
# summary.append(f"\n{'═'*70}")
# summary.append("  TOP 20 K+ COORDINATION EVENTS")
# summary.append("═"*70)
# summary.append(f"{'Frame':<8} {'Frag':<6} {'WatID':<8} {'ProtO':<7} {'W-W':<6} {'Total':<7} {'IonType'}")
# summary.append("─"*60)
# top_k = sorted(k_events, key=lambda x: x["k_total"], reverse=True)[:20]
# for ev in top_k:
#     summary.append(
#         f"{ev['frame']:<8} {ev['frag']:<6} {ev['wat']:<8} "
#         f"{ev['k_prot']:<7} {ev['ww']:<6} {ev['k_total']:<7} {ev['ion_type']}"
#     )
#     summary.append(f"         Residues: {', '.join(ev['k_res'][:5])}")
#
# summary.append(f"\n{'═'*70}")
# summary.append("  TOP 20 Na+ COORDINATION EVENTS")
# summary.append("═"*70)
# summary.append(f"{'Frame':<8} {'Frag':<6} {'WatID':<8} {'ProtO':<7} {'W-W':<6} {'Total':<7} {'IonType'}")
# summary.append("─"*60)
# top_na = sorted(na_events, key=lambda x: x["na_total"], reverse=True)[:20]
# for ev in top_na:
#     summary.append(
#         f"{ev['frame']:<8} {ev['frag']:<6} {ev['wat']:<8} "
#         f"{ev['na_prot']:<7} {ev['ww']:<6} {ev['na_total']:<7} {ev['ion_type']}"
#     )
#     summary.append(f"         Residues: {', '.join(ev['na_res'][:5])}")
#
# # Print and save summary
# for line in summary:
#     print(line)
#
# sum_out.write("\n".join(summary))
# sum_out.close()
#
# # ── Plots ─────────────────────────────────────────────────────
# fig, axes = plt.subplots(2, 3, figsize=(18, 10))
#
# # Plot 1: K+ distance distribution
# if len(k_dists) > 0:
#     bins = np.arange(2.5, 3.5, 0.05)
#     hist, edges = np.histogram(k_dists, bins=bins, density=True)
#     centers = (edges[:-1] + edges[1:]) / 2
#     axes[0][0].bar(centers, hist, width=0.04, color="green", alpha=0.7, edgecolor="black")
#     axes[0][0].axvline(x=2.85, color="red", linestyle="--", linewidth=2, label="KcsA avg (2.85Å)")
#     axes[0][0].axvspan(2.7, 3.3, alpha=0.1, color="green", label="K+ range")
#     axes[0][0].set_title("K+ coordination distances\n(protein O atoms)")
#     axes[0][0].set_xlabel("Distance (Å)")
#     axes[0][0].set_ylabel("Density")
#     axes[0][0].legend(fontsize=8)
#
# # Plot 2: Na+ distance distribution
# if len(na_dists) > 0:
#     bins = np.arange(2.0, 3.0, 0.05)
#     hist, edges = np.histogram(na_dists, bins=bins, density=True)
#     centers = (edges[:-1] + edges[1:]) / 2
#     axes[0][1].bar(centers, hist, width=0.04, color="orange", alpha=0.7, edgecolor="black")
#     axes[0][1].axvline(x=2.4, color="red", linestyle="--", linewidth=2, label="NavAb avg (2.4Å)")
#     axes[0][1].axvspan(2.3, 2.7, alpha=0.1, color="orange", label="Na+ range")
#     axes[0][1].set_title("Na+ coordination distances\n(protein O atoms)")
#     axes[0][1].set_xlabel("Distance (Å)")
#     axes[0][1].set_ylabel("Density")
#     axes[0][1].legend(fontsize=8)
#
# # Plot 3: Water-water distances
# if len(ww_dists) > 0:
#     bins = np.arange(2.0, 3.8, 0.05)
#     hist, edges = np.histogram(ww_dists, bins=bins, density=True)
#     centers = (edges[:-1] + edges[1:]) / 2
#     axes[0][2].bar(centers, hist, width=0.04, color="blue", alpha=0.7, edgecolor="black")
#     axes[0][2].axvline(x=2.8, color="red", linestyle="--", linewidth=2, label="H-bond peak (~2.8Å)")
#     axes[0][2].set_title("Water-water distances\n(H-bond coordination)")
#     axes[0][2].set_xlabel("Distance (Å)")
#     axes[0][2].set_ylabel("Density")
#     axes[0][2].legend(fontsize=8)
#
# # Plot 4: K+ events per monomer
# k_per_frag  = [len([e for e in k_events  if e["frag"]==f]) for f in range(3)]
# na_per_frag = [len([e for e in na_events if e["frag"]==f]) for f in range(3)]
# x = np.arange(3); w = 0.35
# axes[1][0].bar(x-w/2, k_per_frag,  w, color="green",  alpha=0.8, label="K+")
# axes[1][0].bar(x+w/2, na_per_frag, w, color="orange", alpha=0.8, label="Na+")
# axes[1][0].set_xticks(x)
# axes[1][0].set_xticklabels(["Monomer 1", "Monomer 2", "Monomer 3"])
# axes[1][0].set_title("Ion-like events per monomer")
# axes[1][0].set_ylabel("Event count")
# axes[1][0].legend()
#
# # Plot 5: Top K+ residues
# if sorted_k:
#     res_labels = [f"{resname_map.get(r,'???')}{r}" for r, _ in sorted_k[:10]]
#     res_counts = [c for _, c in sorted_k[:10]]
#     axes[1][1].barh(res_labels, res_counts, color="green", alpha=0.8)
#     axes[1][1].set_title("Top residues in K+ sites")
#     axes[1][1].set_xlabel("Contact count")
#     axes[1][1].invert_yaxis()
#
# # Plot 6: Top Na+ residues
# if sorted_na:
#     res_labels = [f"{resname_map.get(r,'???')}{r}" for r, _ in sorted_na[:10]]
#     res_counts = [c for _, c in sorted_na[:10]]
#     axes[1][2].barh(res_labels, res_counts, color="orange", alpha=0.8)
#     axes[1][2].set_title("Top residues in Na+ sites")
#     axes[1][2].set_xlabel("Contact count")
#     axes[1][2].invert_yaxis()
#
# plt.suptitle(
#     "TRIC Channel — Chemically Rigorous Ion Binding Analysis (pH 7)\n"
#     f"K+: {K_MIN}-{K_MAX}Å coord {K_COORD_MIN}-{K_COORD_MAX} | "
#     f"Na+: {NA_MIN}-{NA_MAX}Å coord {NA_COORD_MIN}-{NA_COORD_MAX} | "
#     f"Stride={STRIDE}",
#     fontsize=11, fontweight="bold"
# )
# plt.tight_layout()
# plt.savefig(f"{outdir}/ion_coordination_pH7_corrected.png", dpi=150)
# plt.show()
#
# print(f"\n{'═'*70}")
# print(f"  OUTPUT FILES in {outdir}/")
# print(f"{'═'*70}")
# print(f"  K_coordination.dat              — K+  events only")
# print(f"  Na_coordination.dat             — Na+ events only")
# print(f"  all_coordination.dat            — all events")
# print(f"  coordination_summary.dat        — full statistics")
# print(f"  ion_coordination_pH7_corrected.png — plots")
# print(f"{'═'*70}")






####################################
##################################
################# publication quality fig ########################

# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.ticker as ticker
# from matplotlib import rcParams
# from collections import defaultdict
# import MDAnalysis as mda
# from MDAnalysis.analysis import distances
# import os
#
# # ── Paths ─────────────────────────────────────────────────────
# gro_file  = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/aligned.gro"
# xtc_file  = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/aligned.xtc"
# infile    = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/coordination/all_coordination.dat"
# outdir    = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/coordination"
# outdat    = f"{outdir}/all_coordination_with_distances.dat"
# outfig    = f"{outdir}/6IYZ_cation_binding_sites_REAL.png"
# os.makedirs(outdir, exist_ok=True)
#
# # ── Parameters ────────────────────────────────────────────────
# K_MIN        = 2.7;  K_MAX  = 3.3
# NA_MIN       = 2.3;  NA_MAX = 2.7
# WAT_CUTOFF   = 3.5
#
# # ── Load trajectory ───────────────────────────────────────────
# print("Loading trajectory...")
# u = mda.Universe(gro_file, xtc_file)
# print(f"Total frames: {len(u.trajectory)}")
#
# # ── Monomer boundaries ────────────────────────────────────────
# prot     = u.select_atoms("protein")
# per_mono = len(prot) // 3
# frag_idx = {
#     0: (prot.indices[0],            prot.indices[per_mono-1]),
#     1: (prot.indices[per_mono],     prot.indices[2*per_mono-1]),
#     2: (prot.indices[2*per_mono],   prot.indices[3*per_mono-1])
# }
# print(f"Monomer boundaries: {frag_idx}")
#
# # ── Parse dat file ────────────────────────────────────────────
# print("\nParsing all_coordination.dat...")
#
# records = []
# with open(infile, "r") as f:
#     for line in f:
#         line = line.strip()
#         if line.startswith("#") or line == "":
#             continue
#         parts = line.split()
#         if len(parts) < 7:
#             continue
#         try:
#             frame      = int(parts[0])
#             frag       = int(parts[1])
#             wat_resid  = int(parts[2])
#             prot_count = int(parts[3])
#             ww_count   = int(parts[4])
#             total_coord= int(parts[5])
#             ion_type   = parts[6]
#             coord_str  = parts[7] if len(parts) > 7 else ""
#         except:
#             continue
#
#         coord_atoms = [c for c in coord_str.split("|") if c != ""]
#         records.append({
#             "frame":      frame,
#             "frag":       frag,
#             "wat_resid":  wat_resid,
#             "prot_count": prot_count,
#             "ww_count":   ww_count,
#             "total":      total_coord,
#             "ion_type":   ion_type,
#             "coord_atoms":coord_atoms
#         })
#
# print(f"Total records to process: {len(records)}")
#
# # ── Get unique frames to visit ────────────────────────────────
# unique_frames = sorted(set(r["frame"] for r in records))
# print(f"Unique frames to visit: {len(unique_frames)}")
#
# # ── Build frame → record index mapping ───────────────────────
# from collections import defaultdict
# frame_to_records = defaultdict(list)
# for i, rec in enumerate(records):
#     frame_to_records[rec["frame"]].append(i)
#
# # ── Storage for distances and results ────────────────────────
# k_prot_dists  = []   # real K+  protein O distances
# na_prot_dists = []   # real Na+ protein O distances
# ww_dists      = []   # real water-water distances
#
# k_res_count   = defaultdict(int)
# na_res_count  = defaultdict(int)
# k_per_frag    = [0, 0, 0]
# na_per_frag   = [0, 0, 0]
#
# # store enriched records
# enriched = []
#
# # ── Output file header ────────────────────────────────────────
# fout = open(outdat, "w")
# fout.write(
#     "# Frame  Frag  WaterResID  ProtO_in_range  WatWat  TotalCoord  "
#     "IonType  Atom:Distance(A)|...\n"
# )
#
# # ── Main loop — visit each unique frame ───────────────────────
# print("\nCalculating real distances...")
# processed = 0
#
# for frame_num in unique_frames:
#     # go to this frame
#     u.trajectory[frame_num]
#
#     # get all water positions once
#     all_waters = u.select_atoms("resname TIP3 and name OH2")
#
#     recs = frame_to_records[frame_num]
#
#     for ridx in recs:
#         rec      = records[ridx]
#         frag_id  = rec["frag"]
#         wat_resid= rec["wat_resid"]
#         ion_type = rec["ion_type"]
#
#         # get central water position
#         wat_sel = u.select_atoms(f"resid {wat_resid} and resname TIP3 and name OH2")
#         if len(wat_sel) == 0:
#             continue
#         wat_pos = wat_sel.positions[0]
#
#         # monomer selection
#         f_start, f_end = frag_idx[frag_id]
#
#         atom_dist_list = []
#         real_prot_count = 0
#         real_ww_count   = 0
#
#         for atom_str in rec["coord_atoms"]:
#             atom_str = atom_str.strip()
#             if atom_str == "":
#                 continue
#
#             if "WAT" in atom_str or "TIP3" in atom_str:
#                 # water-water interaction
#                 # extract resid from TIP3XXXXX format
#                 try:
#                     wresid = int(atom_str.replace("TIP3","").replace("(WAT)","").split("-")[0])
#                     other_wat = u.select_atoms(f"resid {wresid} and resname TIP3 and name OH2")
#                     if len(other_wat) > 0:
#                         d = float(np.linalg.norm(wat_pos - other_wat.positions[0]))
#                         atom_dist_list.append(f"{atom_str}:{d:.3f}")
#                         ww_dists.append(d)
#                         real_ww_count += 1
#                 except:
#                     atom_dist_list.append(f"{atom_str}:?")
#
#             else:
#                 # protein atom interaction
#                 # parse: RESNAMEresid-ATOMNAME(TYPE)
#                 # e.g.: SER65-OG(SC-O)  or  HSD61-NE2(HIS-N)
#                 try:
#                     # split off type
#                     base = atom_str.split("(")[0]   # e.g. SER65-OG
#                     parts2 = base.split("-")
#                     # resname+resid = parts2[0], atom = parts2[1]
#                     resname_resid = parts2[0]        # e.g. SER65
#                     atom_name     = parts2[1]        # e.g. OG
#
#                     # extract resname and resid
#                     # resname = leading letters, resid = trailing digits
#                     import re
#                     m = re.match(r"([A-Za-z]+)(\d+)", resname_resid)
#                     if not m:
#                         atom_dist_list.append(f"{atom_str}:?")
#                         continue
#                     resname = m.group(1)
#                     resid   = int(m.group(2))
#
#                     # select atom from correct monomer
#                     prot_atom = u.select_atoms(
#                         f"resid {resid} and name {atom_name} "
#                         f"and index {f_start}:{f_end}"
#                     )
#                     if len(prot_atom) == 0:
#                         atom_dist_list.append(f"{atom_str}:?")
#                         continue
#
#                     d = float(np.linalg.norm(wat_pos - prot_atom.positions[0]))
#                     atom_dist_list.append(f"{atom_str}:{d:.3f}")
#
#                     # classify into K+ or Na+ distance bins
#                     if ion_type in ("K+", "K+/Na+"):
#                         if K_MIN <= d <= K_MAX:
#                             k_prot_dists.append(d)
#                             k_res_count[base] += 1
#                             real_prot_count += 1
#                     if ion_type in ("Na+", "K+/Na+"):
#                         if NA_MIN <= d <= NA_MAX:
#                             na_prot_dists.append(d)
#                             na_res_count[base] += 1
#                             real_prot_count += 1
#
#                 except Exception as e:
#                     atom_dist_list.append(f"{atom_str}:?")
#
#         # per frag counts
#         if ion_type in ("K+", "K+/Na+"):
#             k_per_frag[frag_id] += 1
#         if ion_type in ("Na+", "K+/Na+"):
#             na_per_frag[frag_id] += 1
#
#         # write to output file
#         dist_str = "|".join(atom_dist_list)
#         fout.write(
#             f"{rec['frame']:<8} {frag_id:<6} {wat_resid:<12} "
#             f"{real_prot_count:<16} {real_ww_count:<8} {rec['total']:<12} "
#             f"{ion_type:<12} {dist_str}\n"
#         )
#
#         enriched.append({
#             "frame":    rec["frame"],
#             "frag":     frag_id,
#             "wat":      wat_resid,
#             "ion_type": ion_type,
#             "dists":    atom_dist_list
#         })
#
#     processed += 1
#     if processed % 100 == 0:
#         print(f"  Processed {processed}/{len(unique_frames)} frames  "
#               f"K+ dists: {len(k_prot_dists)}  "
#               f"Na+ dists: {len(na_prot_dists)}")
#
# fout.close()
#
# # save distance arrays for future use
# np.save(f"{outdir}/k_distances_REAL.npy",  np.array(k_prot_dists))
# np.save(f"{outdir}/na_distances_REAL.npy", np.array(na_prot_dists))
# np.save(f"{outdir}/ww_distances_REAL.npy", np.array(ww_dists))
#
# print(f"\n✅ Analysis complete!")
# print(f"   K+  real distances : {len(k_prot_dists)}")
# print(f"   Na+ real distances : {len(na_prot_dists)}")
# print(f"   W-W real distances : {len(ww_dists)}")
# print(f"   Output dat file    : {outdat}")
#
# # ── Top residues ──────────────────────────────────────────────
# sorted_k  = sorted(k_res_count.items(),
#                    key=lambda x: x[1], reverse=True)[:10]
# sorted_na = sorted(na_res_count.items(),
#                    key=lambda x: x[1], reverse=True)[:10]
#
# k_dists  = np.array(k_prot_dists)
# na_dists = np.array(na_prot_dists)
# ww_dists = np.array(ww_dists)
#
# # ══════════════════════════════════════════════════════════════
# #  PUBLICATION QUALITY FIGURE
# # ══════════════════════════════════════════════════════════════
# rcParams['font.family']       = 'DejaVu Sans'
# rcParams['font.weight']       = 'bold'
# rcParams['axes.labelweight']  = 'bold'
# rcParams['axes.titleweight']  = 'bold'
# rcParams['axes.linewidth']    = 1.5
# rcParams['xtick.major.width'] = 1.5
# rcParams['ytick.major.width'] = 1.5
# rcParams['xtick.labelsize']   = 11
# rcParams['ytick.labelsize']   = 11
# rcParams['axes.labelsize']    = 12
# rcParams['axes.titlesize']    = 12
# rcParams['legend.fontsize']   = 10
# rcParams['legend.framealpha'] = 0.9
# rcParams['legend.edgecolor']  = 'black'
# rcParams['figure.dpi']        = 150
# rcParams['savefig.dpi']       = 300
#
# C_K    = "#2ecc71"
# C_NA   = "#e67e22"
# C_WW   = "#3498db"
# C_LINE = "#e74c3c"
#
# fig = plt.figure(figsize=(22, 14))
# gs  = fig.add_gridspec(
#     2, 3,
#     left=0.07, right=0.97,
#     top=0.88,  bottom=0.08,
#     hspace=0.45, wspace=0.38
# )
# axes = [[fig.add_subplot(gs[r, c]) for c in range(3)] for r in range(2)]
#
# def style_ax(ax):
#     for spine in ax.spines.values():
#         spine.set_linewidth(1.5)
#     ax.tick_params(width=1.5, length=5)
#     ax.yaxis.set_tick_params(labelsize=11)
#     ax.xaxis.set_tick_params(labelsize=11)
#
# # ── Panel 1: K+ distances — WIDER RANGE ──────────────────────
# ax = axes[0][0]
# if len(k_dists) > 0:
#     bins = np.arange(2.5, 3.6, 0.05)
#     ax.hist(k_dists, bins=bins, color=C_K, alpha=0.85,
#             edgecolor="white", linewidth=0.8, zorder=3)
#     ax.axvline(x=2.85, color=C_LINE, linestyle="--",
#                linewidth=2.0, zorder=4, label="KcsA avg (2.85 Å)")
#     ax.axvspan(K_MIN, K_MAX, alpha=0.1, color=C_K,
#                label=f"K$^+$ range ({K_MIN}–{K_MAX} Å)", zorder=1)
# ax.set_xlim(2.45, 3.6)
# ax.set_xlabel("Distance (Å)", fontweight="bold")
# ax.set_ylabel("Frequency (count)", fontweight="bold")
# ax.set_title(r"$\mathbf{K^+}$ coordination distances"
#              + "\n(protein O atoms)", fontweight="bold")
# ax.legend(loc="upper right", prop={"weight": "bold", "size": 10})
# ax.yaxis.set_major_formatter(
#     ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
# ax.grid(axis="y", linestyle="--", alpha=0.4, zorder=0)
# style_ax(ax)
#
# # ── Panel 2: Na+ distances — WIDER RANGE ─────────────────────
# ax = axes[0][1]
# if len(na_dists) > 0:
#     bins = np.arange(2.0, 3.0, 0.05)
#     ax.hist(na_dists, bins=bins, color=C_NA, alpha=0.85,
#             edgecolor="white", linewidth=0.8, zorder=3)
#     ax.axvline(x=2.4, color=C_LINE, linestyle="--",
#                linewidth=2.0, zorder=4, label="NavAb avg (2.4 Å)")
#     ax.axvspan(NA_MIN, NA_MAX, alpha=0.1, color=C_NA,
#                label=f"Na$^+$ range ({NA_MIN}–{NA_MAX} Å)", zorder=1)
# ax.set_xlim(1.95, 3.0)
# ax.set_xlabel("Distance (Å)", fontweight="bold")
# ax.set_ylabel("Frequency (count)", fontweight="bold")
# ax.set_title(r"$\mathbf{Na^+}$ coordination distances"
#              + "\n(protein O atoms)", fontweight="bold")
# ax.legend(loc="upper right", prop={"weight": "bold", "size": 10})
# ax.yaxis.set_major_formatter(
#     ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
# ax.grid(axis="y", linestyle="--", alpha=0.4, zorder=0)
# style_ax(ax)
#
# # ── Panel 3: Water-water — WIDER RANGE ───────────────────────
# ax = axes[0][2]
# if len(ww_dists) > 0:
#     bins = np.arange(2.3, WAT_CUTOFF + 0.05, 0.05)
#     ax.hist(ww_dists, bins=bins, color=C_WW, alpha=0.85,
#             edgecolor="white", linewidth=0.8, zorder=3)
#     ax.axvline(x=2.8, color=C_LINE, linestyle="--",
#                linewidth=2.0, zorder=4, label="H-bond peak (~2.8 Å)")
# ax.set_xlim(2.25, WAT_CUTOFF + 0.1)
# ax.set_xlabel("Distance (Å)", fontweight="bold")
# ax.set_ylabel("Frequency (count)", fontweight="bold")
# ax.set_title("Water–Water distances"
#              + "\n(H-bond coordination)", fontweight="bold")
# ax.legend(loc="upper right", prop={"weight": "bold", "size": 10})
# ax.yaxis.set_major_formatter(
#     ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
# ax.grid(axis="y", linestyle="--", alpha=0.4, zorder=0)
# style_ax(ax)
#
# # ── Panel 4: Events per monomer ──────────────────────────────
# ax = axes[1][0]
# x = np.arange(3); w = 0.32
# bk  = ax.bar(x-w/2, k_per_frag,  w, color=C_K,  alpha=0.85,
#              edgecolor="black", linewidth=1.0,
#              label=r"$\mathbf{K^+}$",  zorder=3)
# bna = ax.bar(x+w/2, na_per_frag, w, color=C_NA, alpha=0.85,
#              edgecolor="black", linewidth=1.0,
#              label=r"$\mathbf{Na^+}$", zorder=3)
# for bar in bk:
#     h = bar.get_height()
#     ax.text(bar.get_x()+bar.get_width()/2, h*0.92,
#             f"{int(h):,}", ha="center", va="top",
#             fontsize=10, fontweight="bold", color="white")
# for bar in bna:
#     h = bar.get_height()
#     ax.text(bar.get_x()+bar.get_width()/2, h*0.92,
#             f"{int(h):,}", ha="center", va="top",
#             fontsize=10, fontweight="bold", color="white")
# ax.set_xticks(x)
# ax.set_xticklabels(["Monomer 1", "Monomer 2", "Monomer 3"],
#                    fontweight="bold", fontsize=11)
# ax.set_title("Ion-like events per monomer", fontweight="bold")
# ax.set_ylabel("Event count", fontweight="bold")
# ax.legend(loc="upper right", prop={"weight": "bold", "size": 11})
# ax.yaxis.set_major_formatter(
#     ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
# ax.set_ylim(0, max(k_per_frag) * 1.15)
# ax.grid(axis="y", linestyle="--", alpha=0.4, zorder=0)
# style_ax(ax)
#
# # ── Panel 5: Top K+ residues ──────────────────────────────────
# ax = axes[1][1]
# if sorted_k:
#     labels = [r[0] for r in sorted_k]
#     counts = [r[1] for r in sorted_k]
#     y_pos  = np.arange(len(labels))
#     bars   = ax.barh(y_pos, counts, color=C_K, alpha=0.85,
#                      edgecolor="black", linewidth=0.8, zorder=3)
#     ax.set_yticks(y_pos)
#     ax.set_yticklabels(labels, fontweight="bold", fontsize=10)
#     ax.invert_yaxis()
#     for bar, cnt in zip(bars, counts):
#         ax.text(bar.get_width()*0.96,
#                 bar.get_y()+bar.get_height()/2,
#                 f"{cnt:,}", ha="right", va="center",
#                 fontsize=9, fontweight="bold", color="white")
#     ax.set_xlim(0, max(counts)*1.05)
# ax.set_title(r"Top residues in $\bf{K^+}$ sites", fontweight="bold")
# ax.set_xlabel("Contact count", fontweight="bold")
# ax.xaxis.set_major_formatter(
#     ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
# ax.grid(axis="x", linestyle="--", alpha=0.4, zorder=0)
# style_ax(ax)
#
# # ── Panel 6: Top Na+ residues ─────────────────────────────────
# ax = axes[1][2]
# if sorted_na:
#     labels = [r[0] for r in sorted_na]
#     counts = [r[1] for r in sorted_na]
#     y_pos  = np.arange(len(labels))
#     bars   = ax.barh(y_pos, counts, color=C_NA, alpha=0.85,
#                      edgecolor="black", linewidth=0.8, zorder=3)
#     ax.set_yticks(y_pos)
#     ax.set_yticklabels(labels, fontweight="bold", fontsize=10)
#     ax.invert_yaxis()
#     for bar, cnt in zip(bars, counts):
#         ax.text(bar.get_width()*0.96,
#                 bar.get_y()+bar.get_height()/2,
#                 f"{cnt:,}", ha="right", va="center",
#                 fontsize=9, fontweight="bold", color="white")
#     ax.set_xlim(0, max(counts)*1.05)
# ax.set_title(r"Top residues in $\bf{Na^+}$ sites", fontweight="bold")
# ax.set_xlabel("Contact count", fontweight="bold")
# ax.xaxis.set_major_formatter(
#     ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
# ax.grid(axis="x", linestyle="--", alpha=0.4, zorder=0)
# style_ax(ax)
#
# # ── Main title ────────────────────────────────────────────────
# fig.suptitle(
#     r"Monovalent Cation Binding Site Mimics in 6IYZ Structure"
#     "\n"
#     r"$\mathbf{K^+}$ coord 6–8  |  $\mathbf{Na^+}$ coord 4–6 ",
#     fontsize=14, fontweight="bold", y=0.97
# )
#
# plt.savefig(outfig, dpi=300, bbox_inches="tight",
#             facecolor="white", edgecolor="none")
# plt.show()
#
# print(f"\n✅ Publication figure saved : {outfig}")
# print(f"✅ Enriched dat file saved  : {outdat}")
# print(f"✅ K+  distances saved      : {outdir}/k_distances_REAL.npy")
# print(f"✅ Na+ distances saved      : {outdir}/na_distances_REAL.npy")
# print(f"✅ W-W distances saved      : {outdir}/ww_distances_REAL.npy")
#

####################################
##################################

################# verification ########################
import re
from collections import defaultdict
#
# # ── Paths ─────────────────────────────────────────────────────
# infile = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/coordination/all_coordination_with_distances.dat"
# outfile = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/coordination/verification_report.dat"
#
# # ── Parameters ────────────────────────────────────────────────
# K_MIN  = 2.7;  K_MAX  = 3.3
# NA_MIN = 2.3;  NA_MAX = 2.7
# WAT_CUTOFF = 3.5
#
# # ── Counters ──────────────────────────────────────────────────
# total_records     = 0
# total_prot_pairs  = 0
# total_ww_pairs    = 0
#
# # K+ verification
# k_records         = 0
# k_prot_verified   = 0   # distance within K+ range ✅
# k_prot_mismatch   = 0   # distance outside K+ range ❌
# k_ww_verified     = 0
# k_ww_mismatch     = 0
#
# # Na+ verification
# na_records        = 0
# na_prot_verified  = 0
# na_prot_mismatch  = 0
# na_ww_verified    = 0
# na_ww_mismatch    = 0
#
# # store mismatches for report
# mismatches = []
#
# # distance collections for stats
# k_dists_ok    = []
# k_dists_bad   = []
# na_dists_ok   = []
# na_dists_bad  = []
# ww_dists_ok   = []
# ww_dists_bad  = []
#
# print("Reading and verifying all_coordination_with_distances.dat...")
# print("═"*60)
#
# with open(infile, "r") as f:
#     for line in f:
#         line = line.strip()
#         if line.startswith("#") or line == "":
#             continue
#
#         parts = line.split()
#         if len(parts) < 7:
#             continue
#
#         try:
#             frame      = int(parts[0])
#             frag       = int(parts[1])
#             wat_resid  = int(parts[2])
#             ion_type   = parts[6]
#             coord_str  = parts[7] if len(parts) > 7 else ""
#         except:
#             continue
#
#         total_records += 1
#         if ion_type in ("K+",    "K+/Na+"): k_records  += 1
#         if ion_type in ("Na+",   "K+/Na+"): na_records += 1
#
#         coord_atoms = [c for c in coord_str.split("|") if c != ""]
#
#         for atom_str in coord_atoms:
#             atom_str = atom_str.strip()
#             if atom_str == "" or atom_str == "?":
#                 continue
#
#             # parse distance — format: ATOM:distance  e.g. SER65-OG(SC-O):2.834
#             if ":" not in atom_str:
#                 continue
#
#             atom_part, dist_part = atom_str.rsplit(":", 1)
#
#             if dist_part == "?":
#                 continue
#
#             try:
#                 d = float(dist_part)
#             except:
#                 continue
#
#             is_water = "WAT" in atom_part or "TIP3" in atom_part
#
#             # ── Verify water-water ────────────────────────────
#             if is_water:
#                 total_ww_pairs += 1
#                 if d <= WAT_CUTOFF:
#                     k_ww_verified  += 1 if ion_type in ("K+",  "K+/Na+") else 0
#                     na_ww_verified += 1 if ion_type in ("Na+", "K+/Na+") else 0
#                     ww_dists_ok.append(d)
#                 else:
#                     k_ww_mismatch  += 1 if ion_type in ("K+",  "K+/Na+") else 0
#                     na_ww_mismatch += 1 if ion_type in ("Na+", "K+/Na+") else 0
#                     ww_dists_bad.append(d)
#                     mismatches.append({
#                         "frame": frame, "frag": frag,
#                         "wat": wat_resid, "ion_type": ion_type,
#                         "atom": atom_part, "dist": d,
#                         "expected": f"<={WAT_CUTOFF}",
#                         "type": "WAT-WAT"
#                     })
#
#             # ── Verify protein atom ───────────────────────────
#             else:
#                 total_prot_pairs += 1
#
#                 # K+ verification
#                 if ion_type in ("K+", "K+/Na+"):
#                     if K_MIN <= d <= K_MAX:
#                         k_prot_verified += 1
#                         k_dists_ok.append(d)
#                     else:
#                         k_prot_mismatch += 1
#                         k_dists_bad.append(d)
#                         mismatches.append({
#                             "frame": frame, "frag": frag,
#                             "wat": wat_resid, "ion_type": ion_type,
#                             "atom": atom_part, "dist": d,
#                             "expected": f"{K_MIN}-{K_MAX}",
#                             "type": "K+-PROT"
#                         })
#
#                 # Na+ verification
#                 if ion_type in ("Na+", "K+/Na+"):
#                     if NA_MIN <= d <= NA_MAX:
#                         na_prot_verified += 1
#                         na_dists_ok.append(d)
#                     else:
#                         na_prot_mismatch += 1
#                         na_dists_bad.append(d)
#                         mismatches.append({
#                             "frame": frame, "frag": frag,
#                             "wat": wat_resid, "ion_type": ion_type,
#                             "atom": atom_part, "dist": d,
#                             "expected": f"{NA_MIN}-{NA_MAX}",
#                             "type": "Na+-PROT"
#                         })
#
# # ── Write verification report ─────────────────────────────────
# import numpy as np
#
# with open(outfile, "w") as f:
#
#     f.write("═"*65 + "\n")
#     f.write("  VERIFICATION REPORT — all_coordination_with_distances.dat\n")
#     f.write("═"*65 + "\n\n")
#
#     f.write(f"  K+  distance range  : {K_MIN} - {K_MAX} Å\n")
#     f.write(f"  Na+ distance range  : {NA_MIN} - {NA_MAX} Å\n")
#     f.write(f"  Water-water cutoff  : <= {WAT_CUTOFF} Å\n\n")
#
#     f.write(f"  Total records       : {total_records}\n")
#     f.write(f"  K+  records         : {k_records}\n")
#     f.write(f"  Na+ records         : {na_records}\n")
#     f.write(f"  Total protein pairs : {total_prot_pairs}\n")
#     f.write(f"  Total W-W pairs     : {total_ww_pairs}\n\n")
#
#     f.write("─"*65 + "\n")
#     f.write("  K+ PROTEIN DISTANCE VERIFICATION\n")
#     f.write("─"*65 + "\n")
#     k_total_prot = k_prot_verified + k_prot_mismatch
#     k_pass_pct   = 100*k_prot_verified/k_total_prot if k_total_prot > 0 else 0
#     f.write(f"  Pairs checked       : {k_total_prot}\n")
#     f.write(f"  ✅ Within range     : {k_prot_verified} ({k_pass_pct:.1f}%)\n")
#     f.write(f"  ❌ Outside range    : {k_prot_mismatch} ({100-k_pass_pct:.1f}%)\n")
#     if k_dists_ok:
#         f.write(f"  Mean dist (valid)  : {np.mean(k_dists_ok):.3f} ± {np.std(k_dists_ok):.3f} Å\n")
#         f.write(f"  Min/Max (valid)    : {np.min(k_dists_ok):.3f} / {np.max(k_dists_ok):.3f} Å\n")
#     if k_dists_bad:
#         f.write(f"  Mean dist (bad)    : {np.mean(k_dists_bad):.3f} ± {np.std(k_dists_bad):.3f} Å\n")
#
#     f.write("\n" + "─"*65 + "\n")
#     f.write("  Na+ PROTEIN DISTANCE VERIFICATION\n")
#     f.write("─"*65 + "\n")
#     na_total_prot = na_prot_verified + na_prot_mismatch
#     na_pass_pct   = 100*na_prot_verified/na_total_prot if na_total_prot > 0 else 0
#     f.write(f"  Pairs checked       : {na_total_prot}\n")
#     f.write(f"  ✅ Within range     : {na_prot_verified} ({na_pass_pct:.1f}%)\n")
#     f.write(f"  ❌ Outside range    : {na_prot_mismatch} ({100-na_pass_pct:.1f}%)\n")
#     if na_dists_ok:
#         f.write(f"  Mean dist (valid)  : {np.mean(na_dists_ok):.3f} ± {np.std(na_dists_ok):.3f} Å\n")
#         f.write(f"  Min/Max (valid)    : {np.min(na_dists_ok):.3f} / {np.max(na_dists_ok):.3f} Å\n")
#     if na_dists_bad:
#         f.write(f"  Mean dist (bad)    : {np.mean(na_dists_bad):.3f} ± {np.std(na_dists_bad):.3f} Å\n")
#
#     f.write("\n" + "─"*65 + "\n")
#     f.write("  WATER-WATER DISTANCE VERIFICATION\n")
#     f.write("─"*65 + "\n")
#     ww_total    = k_ww_verified + k_ww_mismatch
#     ww_pass_pct = 100*k_ww_verified/ww_total if ww_total > 0 else 0
#     f.write(f"  Pairs checked       : {ww_total}\n")
#     f.write(f"  ✅ Within cutoff    : {k_ww_verified} ({ww_pass_pct:.1f}%)\n")
#     f.write(f"  ❌ Outside cutoff   : {k_ww_mismatch} ({100-ww_pass_pct:.1f}%)\n")
#     if ww_dists_ok:
#         f.write(f"  Mean dist (valid)  : {np.mean(ww_dists_ok):.3f} ± {np.std(ww_dists_ok):.3f} Å\n")
#         f.write(f"  Min/Max (valid)    : {np.min(ww_dists_ok):.3f} / {np.max(ww_dists_ok):.3f} Å\n")
#
#     # ── Overall verdict ───────────────────────────────────────
#     f.write("\n" + "═"*65 + "\n")
#     f.write("  OVERALL VERDICT\n")
#     f.write("═"*65 + "\n")
#
#     overall_verified = k_prot_verified + na_prot_verified + k_ww_verified
#     overall_total    = (k_prot_verified + k_prot_mismatch +
#                         na_prot_verified + na_prot_mismatch +
#                         k_ww_verified + k_ww_mismatch)
#     overall_pct = 100*overall_verified/overall_total if overall_total > 0 else 0
#
#     f.write(f"  Total pairs verified : {overall_total}\n")
#     f.write(f"  Overall pass rate    : {overall_pct:.2f}%\n\n")
#
#     if overall_pct >= 95:
#         f.write("  ✅ ANALYSIS VERIFIED — previous analysis is CORRECT\n")
#         f.write("     >95% of distances match expected ranges\n")
#         f.write("     Results are suitable for publication\n")
#     elif overall_pct >= 85:
#         f.write("  ⚠️  MOSTLY VERIFIED — minor discrepancies present\n")
#         f.write("     85-95% of distances match — acceptable but review mismatches\n")
#     else:
#         f.write("  ❌ VERIFICATION FAILED — significant mismatches found\n")
#         f.write("     <85% of distances match — review analysis parameters\n")
#
#     # ── Top 20 mismatches ─────────────────────────────────────
#     if mismatches:
#         f.write(f"\n{'─'*65}\n")
#         f.write(f"  FIRST 20 MISMATCHES (out of {len(mismatches)} total)\n")
#         f.write(f"{'─'*65}\n")
#         f.write(f"{'Frame':<8} {'Frag':<6} {'WatID':<8} {'Type':<12} "
#                 f"{'Atom':<25} {'Dist(Å)':<10} {'Expected'}\n")
#         f.write("─"*65 + "\n")
#         for m in mismatches[:20]:
#             f.write(
#                 f"{m['frame']:<8} {m['frag']:<6} {m['wat']:<8} "
#                 f"{m['type']:<12} {m['atom']:<25} "
#                 f"{m['dist']:<10.3f} {m['expected']}\n"
#             )
#
# # ── Print summary to console ──────────────────────────────────
# print(f"  Total records       : {total_records}")
# print(f"  K+  records         : {k_records}")
# print(f"  Na+ records         : {na_records}")
# print(f"\n  K+  PROTEIN:")
# print(f"    ✅ Verified        : {k_prot_verified} ({k_pass_pct:.1f}%)")
# print(f"    ❌ Mismatches      : {k_prot_mismatch}")
# if k_dists_ok:
#     print(f"    Mean distance     : {np.mean(k_dists_ok):.3f} ± {np.std(k_dists_ok):.3f} Å")
#
# print(f"\n  Na+ PROTEIN:")
# print(f"    ✅ Verified        : {na_prot_verified} ({na_pass_pct:.1f}%)")
# print(f"    ❌ Mismatches      : {na_prot_mismatch}")
# if na_dists_ok:
#     print(f"    Mean distance     : {np.mean(na_dists_ok):.3f} ± {np.std(na_dists_ok):.3f} Å")
#
# print(f"\n  WATER-WATER:")
# print(f"    ✅ Verified        : {k_ww_verified} ({ww_pass_pct:.1f}%)")
# print(f"    ❌ Mismatches      : {k_ww_mismatch}")
# if ww_dists_ok:
#     print(f"    Mean distance     : {np.mean(ww_dists_ok):.3f} ± {np.std(ww_dists_ok):.3f} Å")
#
# print(f"\n  Total mismatches    : {len(mismatches)}")
# print(f"  Overall pass rate   : {overall_pct:.2f}%")
#
# if overall_pct >= 95:
#     print(f"\n  ✅ ANALYSIS VERIFIED — results suitable for publication!")
# elif overall_pct >= 85:
#     print(f"\n  ⚠️  MOSTLY VERIFIED — review mismatches before publication")
# else:
#     print(f"\n  ❌ VERIFICATION FAILED — significant issues found")
#
# print(f"\n✅ Full report saved: {outfile}")
#

##############################################################################################################################
############################### Residence time ###############################################################################
##############################################################################################################################
# """
# PART ONE OF THE ANALYSIS - FIND WATER INSIDE THE TRIMER PORE OF THE CHANNEL AT ANY GIVEN TIME
# Pore Water Count Analysis — TRIC Channel
# =========================================
# For every frame:
#   1. Calculate protein COM
#   2. Define cylinder: radius=5A, Z = COM ± 10A
#   3. Count TIP3 oxygen (OH2) atoms inside cylinder ONLY
#   4. Save time(ps) and water count to output file
# """
#
# import MDAnalysis as mda
# import numpy as np
# import os
#
# # ─────────────────────────────────────────────────────────────────────────────
# # PATHS
# # ─────────────────────────────────────────────────────────────────────────────
# GRO     = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/center.gro"
# XTC     = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/center.xtc"
# OUT_DIR = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/residence_time/"
#
# os.makedirs(OUT_DIR, exist_ok=True)
#
# OUT_FILE = os.path.join(OUT_DIR, "pore_water_count.dat")
#
# # ─────────────────────────────────────────────────────────────────────────────
# # PARAMETERS
# # ─────────────────────────────────────────────────────────────────────────────
# RADIUS_A     = 5.0    # Cylinder radius (Angstroms)
# Z_HALF       = 10.0   # Z half-window (Angstroms) — ±10 A from protein COM
# TIME_STEP_PS = 10.0   # ps per frame (dt=0.002 x nstxout-compressed=5000)
#
# # ─────────────────────────────────────────────────────────────────────────────
# # LOAD
# # ─────────────────────────────────────────────────────────────────────────────
# print("Loading topology and trajectory...")
# u        = mda.Universe(GRO, XTC)
# n_frames = len(u.trajectory)
# print(f"  Total frames        : {n_frames}")
# print(f"  Time per frame      : {TIME_STEP_PS} ps")
# print(f"  Total sim time      : {n_frames * TIME_STEP_PS:.1f} ps")
# print(f"  Cylinder radius     : {RADIUS_A} A")
# print(f"  Z window            : COM +/- {Z_HALF} A")
#
# # Pre-select atom groups (done once — faster than re-selecting every frame)
# protein  = u.select_atoms("protein")
# tip3_oxy = u.select_atoms("resname TIP3 and name OH2")
#
# print(f"  Protein atoms       : {len(protein)}")
# print(f"  TIP3 OH2 atoms      : {len(tip3_oxy)}")
#
# # ─────────────────────────────────────────────────────────────────────────────
# # MAIN LOOP
# # ─────────────────────────────────────────────────────────────────────────────
# print(f"\nRunning analysis...")
#
# with open(OUT_FILE, 'w') as fout:
#
#     # Header
#     fout.write("# Pore Water Count — TRIC Channel (6IYX)\n")
#     fout.write(f"# Cylinder: radius={RADIUS_A}A  Z=COM+/-{Z_HALF}A\n")
#     fout.write(f"# TIP3 oxygen (OH2) only\n")
#     fout.write(f"# Columns: Time(ps)   N_waters   Water_ResIDs\n")
#     fout.write("#\n")
#
#     for frame_idx, ts in enumerate(u.trajectory):
#
#         # Time in ps
#         time_ps = frame_idx * TIME_STEP_PS
#
#         # Protein COM at this frame
#         com        = protein.center_of_mass()
#         cx, cy, cz = com[0], com[1], com[2]
#         zmin       = cz - Z_HALF
#         zmax       = cz + Z_HALF
#
#         # TIP3 oxygen positions
#         pos  = tip3_oxy.positions
#         dx   = pos[:, 0] - cx
#         dy   = pos[:, 1] - cy
#         dz   = pos[:, 2]
#
#         # Cylinder selection: r <= 5A AND z within range
#         in_cyl   = (dx**2 + dy**2 <= RADIUS_A**2) & (dz >= zmin) & (dz <= zmax)
#         n_waters = int(np.sum(in_cyl))
#
#         # Get residue IDs of waters inside cylinder
#         pore_oxy    = tip3_oxy[in_cyl]
#         water_ids   = sorted(set(pore_oxy.resids.tolist()))
#         ids_str     = ",".join(str(r) for r in water_ids) if water_ids else "none"
#
#         # Write to file
#         fout.write(f"{time_ps:10.1f}   {n_waters:5d}   {ids_str}\n")
#
#         # Progress every 1000 frames
#         if frame_idx % 1000 == 0 or frame_idx == n_frames - 1:
#             print(f"  Frame {frame_idx+1:6d}/{n_frames} | "
#                   f"{time_ps:10.1f} ps | "
#                   f"Waters in pore: {n_waters}")
#
# print(f"\nDone!")
# print(f"Output saved to: {OUT_FILE}")
#
#
#
# """
# PART TWO CALCULATION OF THE RESIDENCE TIME OF WATER INSIDE THE TRIMER FRAME
# Water Residence Time Analysis — TRIC Channel Pore
# ===================================================
# Reads pore_water_count.dat and computes:
#   1. Total unique water molecules that entered the pore
#   2. Average number of waters per frame
#   3. Consecutive residence times (broken if water leaves and re-enters)
#   4. Average residence time
#   5. Top 5 longest residence times
# """
#
# import os
# from collections import defaultdict
#
# # ─────────────────────────────────────────────────────────────────────────────
# # PATHS
# # ─────────────────────────────────────────────────────────────────────────────
# IN_FILE  = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/residence_time/pore_water_count.dat"
# OUT_FILE = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/residence_time/residence_time_summary.dat"
#
# TIME_STEP_PS = 10.0   # ps per frame
#
# # ─────────────────────────────────────────────────────────────────────────────
# # READ INPUT FILE
# # ─────────────────────────────────────────────────────────────────────────────
# print("Reading input file...")
#
# times       = []   # list of time values (ps)
# water_sets  = []   # list of sets of water resids per frame
#
# with open(IN_FILE, 'r') as f:
#     for line in f:
#         line = line.strip()
#         if not line or line.startswith('#'):
#             continue
#         parts = line.split()
#         time_ps = float(parts[0])
#         n_wat   = int(parts[1])
#         if n_wat == 0 or len(parts) < 3:
#             waters = set()
#         else:
#             waters = set(int(x) for x in parts[2].split(','))
#         times.append(time_ps)
#         water_sets.append(waters)
#
# n_frames = len(times)
# print(f"  Frames read         : {n_frames}")
# print(f"  Time range          : {times[0]:.1f} – {times[-1]:.1f} ps")
#
# # ─────────────────────────────────────────────────────────────────────────────
# # 1. UNIQUE WATER MOLECULES & AVERAGE COUNT
# # ─────────────────────────────────────────────────────────────────────────────
# all_unique_waters = set()
# total_water_count = 0
#
# for ws in water_sets:
#     all_unique_waters.update(ws)
#     total_water_count += len(ws)
#
# n_unique     = len(all_unique_waters)
# avg_per_frame = total_water_count / n_frames
#
# print(f"\n  Unique water molecules that entered pore : {n_unique}")
# print(f"  Average waters per frame                 : {avg_per_frame:.2f}")
#
# # ─────────────────────────────────────────────────────────────────────────────
# # 2. CONSECUTIVE RESIDENCE TIMES
# # For each water molecule, find all consecutive runs (unbroken stretches)
# # A run breaks when the water is absent for even one frame
# # Each run = one residence event
# # ─────────────────────────────────────────────────────────────────────────────
# print("\nCalculating consecutive residence times...")
#
# # Build presence list per water: list of frame indices where it is present
# water_presence = defaultdict(list)
# for fi, ws in enumerate(water_sets):
#     for wid in ws:
#         water_presence[wid].append(fi)
#
# # For each water find consecutive runs
# # A consecutive run = frame indices with no gap (difference = 1)
# residence_events = []   # list of (water_id, start_time_ps, end_time_ps, duration_ps)
#
# for wid, frame_list in water_presence.items():
#     frame_list = sorted(frame_list)
#
#     # Find consecutive groups
#     run_start = frame_list[0]
#     prev      = frame_list[0]
#
#     for fi in frame_list[1:]:
#         if fi == prev + 1:
#             # Still consecutive
#             prev = fi
#         else:
#             # Gap found — save completed run
#             duration = (prev - run_start + 1) * TIME_STEP_PS
#             start_t  = times[run_start]
#             end_t    = times[prev]
#             residence_events.append((wid, start_t, end_t, duration))
#             # Start new run
#             run_start = fi
#             prev      = fi
#
#     # Save last run
#     duration = (prev - run_start + 1) * TIME_STEP_PS
#     start_t  = times[run_start]
#     end_t    = times[prev]
#     residence_events.append((wid, start_t, end_t, duration))
#
# n_events = len(residence_events)
# print(f"  Total residence events (consecutive runs): {n_events}")
#
# # ─────────────────────────────────────────────────────────────────────────────
# # 3. STATISTICS
# # ─────────────────────────────────────────────────────────────────────────────
# durations   = [e[3] for e in residence_events]
# avg_res     = sum(durations) / len(durations)
# max_res     = max(durations)
# min_res     = min(durations)
#
# # Median
# sorted_dur  = sorted(durations)
# mid         = len(sorted_dur) // 2
# median_res  = sorted_dur[mid] if len(sorted_dur) % 2 != 0 \
#               else (sorted_dur[mid-1] + sorted_dur[mid]) / 2
#
# # Top 5 longest
# top5 = sorted(residence_events, key=lambda x: x[3], reverse=True)[:5]
#
# # Waters that re-entered (had more than 1 residence event)
# reentry_count = defaultdict(int)
# for wid, _, _, _ in residence_events:
#     reentry_count[wid] += 1
# re_entrants = {wid: cnt for wid, cnt in reentry_count.items() if cnt > 1}
#
# print(f"  Average residence time : {avg_res:.1f} ps")
# print(f"  Median residence time  : {median_res:.1f} ps")
# print(f"  Longest residence time : {max_res:.1f} ps")
# print(f"  Shortest residence time: {min_res:.1f} ps")
# print(f"  Waters that re-entered : {len(re_entrants)}")
#
# # ─────────────────────────────────────────────────────────────────────────────
# # 4. WRITE SUMMARY FILE
# # ─────────────────────────────────────────────────────────────────────────────
# print(f"\nWriting summary to {OUT_FILE}...")
#
# with open(OUT_FILE, 'w') as f:
#
#     f.write("=" * 65 + "\n")
#     f.write("  WATER RESIDENCE TIME SUMMARY — TRIC Pore (6IYX)\n")
#     f.write("=" * 65 + "\n\n")
#
#     # ── Section 1 ────────────────────────────────────────────────
#     f.write("─" * 65 + "\n")
#     f.write("  1. GENERAL STATISTICS\n")
#     f.write("─" * 65 + "\n")
#     f.write(f"  Total frames analyzed          : {n_frames}\n")
#     f.write(f"  Time range                     : {times[0]:.1f} – {times[-1]:.1f} ps\n")
#     f.write(f"  Time step                      : {TIME_STEP_PS:.1f} ps\n")
#     f.write(f"  Cylinder definition            : r=5.0A, Z=COM+/-10.0A\n\n")
#     f.write(f"  Unique water molecules (total) : {n_unique}\n")
#     f.write(f"  Average waters per frame       : {avg_per_frame:.2f}\n")
#     f.write(f"  Total residence events         : {n_events}\n")
#     f.write(f"  Waters that re-entered pore    : {len(re_entrants)}\n\n")
#
#     # ── Section 2 ────────────────────────────────────────────────
#     f.write("─" * 65 + "\n")
#     f.write("  2. RESIDENCE TIME STATISTICS\n")
#     f.write("─" * 65 + "\n")
#     f.write(f"  Average residence time         : {avg_res:.1f} ps\n")
#     f.write(f"  Median residence time          : {median_res:.1f} ps\n")
#     f.write(f"  Longest residence time         : {max_res:.1f} ps\n")
#     f.write(f"  Shortest residence time        : {min_res:.1f} ps\n\n")
#
#     # Distribution
#     bins = [10, 50, 100, 500, 1000, 5000, 10000, float('inf')]
#     bin_labels = ["10", "50", "100", "500", "1000", "5000", "10000", ">10000"]
#     f.write(f"  Residence time distribution:\n")
#     f.write(f"  {'Range (ps)':<20} {'Count':<10} {'Percent'}\n")
#     f.write(f"  {'-'*18} {'-'*8} {'-'*8}\n")
#     prev_b = 0
#     for b, label in zip(bins, bin_labels):
#         count = sum(1 for d in durations if prev_b < d <= b)
#         pct   = count / n_events * 100
#         range_str = f"{prev_b}–{label}"
#         f.write(f"  {range_str:<20} {count:<10} {pct:.1f}%\n")
#         prev_b = b
#
#     # ── Section 3 ────────────────────────────────────────────────
#     f.write(f"\n{'─'*65}\n")
#     f.write("  3. TOP 5 LONGEST RESIDENCE EVENTS\n")
#     f.write("─" * 65 + "\n")
#     f.write(f"  {'Rank':<6} {'WaterID':<10} {'Start(ps)':<12} "
#             f"{'End(ps)':<12} {'Duration(ps)':<14} {'Duration(ns)'}\n")
#     f.write(f"  {'-'*5} {'-'*9} {'-'*11} {'-'*11} {'-'*13} {'-'*12}\n")
#     for rank, (wid, st, et, dur) in enumerate(top5, 1):
#         f.write(f"  {rank:<6} {wid:<10} {st:<12.1f} {et:<12.1f} "
#                 f"{dur:<14.1f} {dur/1000:.3f}\n")
#
#     # ── Section 4 ────────────────────────────────────────────────
#     f.write(f"\n{'─'*65}\n")
#     f.write("  4. WATERS THAT RE-ENTERED THE PORE\n")
#     f.write("─" * 65 + "\n")
#     if re_entrants:
#         f.write(f"  {'WaterID':<10} {'N_entries':<12} {'Total_time(ps)'}\n")
#         f.write(f"  {'-'*9} {'-'*11} {'-'*14}\n")
#         for wid, cnt in sorted(re_entrants.items(),
#                                 key=lambda x: x[1], reverse=True):
#             total_t = sum(e[3] for e in residence_events if e[0] == wid)
#             f.write(f"  {wid:<10} {cnt:<12} {total_t:.1f}\n")
#     else:
#         f.write("  No water molecules re-entered the pore.\n")
#
#     # ── Section 5 — Full event list ───────────────────────────────
#     f.write(f"\n{'─'*65}\n")
#     f.write("  5. ALL RESIDENCE EVENTS (sorted by duration, longest first)\n")
#     f.write("─" * 65 + "\n")
#     f.write(f"  {'WaterID':<10} {'Start(ps)':<12} {'End(ps)':<12} "
#             f"{'Duration(ps)':<14} {'Duration(ns)'}\n")
#     f.write(f"  {'-'*9} {'-'*11} {'-'*11} {'-'*13} {'-'*12}\n")
#     for wid, st, et, dur in sorted(residence_events,
#                                     key=lambda x: x[3], reverse=True):
#         f.write(f"  {wid:<10} {st:<12.1f} {et:<12.1f} "
#                 f"{dur:<14.1f} {dur/1000:.3f}\n")
#
#     f.write(f"\n{'='*65}\n")
#     f.write("  END OF SUMMARY\n")
#     f.write("=" * 65 + "\n")
#
# print(f"\nDone! Summary saved to:\n  {OUT_FILE}")
# ##############################################################################################################################
# ############################### Residence time ###############################################################################
# ##############################################################################################################################

############ fig residence time #######
# """
# Publication-grade comparison plot: 6IYX vs 6IYZ
# Residence Time Inside the Trimer Pore — REMADE
# """
#
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
# import matplotlib.gridspec as gridspec
# import matplotlib.ticker as ticker
# import numpy as np
# import re
# import os
#
# # ─────────────────────────────────────────────────────────────────────────────
# # PATHS
# # ─────────────────────────────────────────────────────────────────────────────
# SUMMARY_6IYX = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/residence_time/residence_time_summary.dat"
# SUMMARY_6IYZ = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/residence_time/residence_time_summary.dat"
# COUNT_6IYX   = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyx/production_run/residence_time/pore_water_count.dat"
# COUNT_6IYZ   = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/6iyz/production_run/residence_time/pore_water_count.dat"
# OUT_DIR      = "/media/supremeleader/Pantera/simulation/2024_simulation_analysis/"
# OUT_FILE     = os.path.join(OUT_DIR, "pore_water_comparison_v2.png")
#
# # ─────────────────────────────────────────────────────────────────────────────
# # PARSE FUNCTIONS
# # ─────────────────────────────────────────────────────────────────────────────
# def parse_summary(filepath):
#     stats = {}
#     with open(filepath, 'r') as f:
#         content = f.read()
#     patterns = {
#         'avg_per_frame'   : r'Average waters per frame\s*:\s*([\d.]+)',
#         'avg_residence'   : r'Average residence time\s*:\s*([\d.]+)',
#         'median_residence': r'Median residence time\s*:\s*([\d.]+)',
#         'longest'         : r'Longest residence time\s*:\s*([\d.]+)',
#         'shortest'        : r'Shortest residence time\s*:\s*([\d.]+)',
#         'n_unique'        : r'Unique water molecules \(total\)\s*:\s*(\d+)',
#         'n_events'        : r'Total residence events\s*:\s*(\d+)',
#     }
#     for key, pat in patterns.items():
#         m = re.search(pat, content)
#         stats[key] = float(m.group(1)) if m else 0.0
#     return stats
#
# def parse_count(filepath):
#     times, counts = [], []
#     with open(filepath, 'r') as f:
#         for line in f:
#             line = line.strip()
#             if not line or line.startswith('#'):
#                 continue
#             parts = line.split()
#             times.append(float(parts[0]))
#             counts.append(int(parts[1]))
#     return np.array(times), np.array(counts)
#
# def parse_all_events(filepath):
#     durations = []
#     in_section5 = False
#     with open(filepath, 'r') as f:
#         for line in f:
#             if '5. ALL RESIDENCE EVENTS' in line:
#                 in_section5 = True
#                 continue
#             if in_section5:
#                 if line.strip().startswith('=') or 'END' in line:
#                     break
#                 if line.strip().startswith('-') or line.strip().startswith('Water'):
#                     continue
#                 parts = line.strip().split()
#                 if len(parts) >= 4:
#                     try:
#                         durations.append(float(parts[3]))
#                     except:
#                         pass
#     return np.array(durations)
#
# print("Parsing data...")
# sx  = parse_summary(SUMMARY_6IYX)
# sz  = parse_summary(SUMMARY_6IYZ)
# _, counts_x = parse_count(COUNT_6IYX)
# _, counts_z = parse_count(COUNT_6IYZ)
# dur_x = parse_all_events(SUMMARY_6IYX)
# dur_z = parse_all_events(SUMMARY_6IYZ)
#
# print(f"  6IYX: avg_waters={sx['avg_per_frame']:.2f}  avg_res={sx['avg_residence']:.1f} ps  longest={sx['longest']:.0f} ps")
# print(f"  6IYZ: avg_waters={sz['avg_per_frame']:.2f}  avg_res={sz['avg_residence']:.1f} ps  longest={sz['longest']:.0f} ps")
#
# # ─────────────────────────────────────────────────────────────────────────────
# # STYLE
# # ─────────────────────────────────────────────────────────────────────────────
# plt.rcParams.update({
#     'font.family'       : 'DejaVu Sans',
#     'font.size'         : 11,
#     'axes.linewidth'    : 1.2,
#     'axes.spines.top'   : False,
#     'axes.spines.right' : False,
#     'xtick.direction'   : 'out',
#     'ytick.direction'   : 'out',
#     'xtick.major.width' : 1.2,
#     'ytick.major.width' : 1.2,
#     'xtick.major.size'  : 4,
#     'ytick.major.size'  : 4,
#     'pdf.fonttype'      : 42,
# })
#
# COL_X  = '#1565C0'   # blue  — 6IYX closed
# COL_Z  = '#C62828'   # red   — 6IYZ open
# BARCOL_X = '#4A90D9'
# BARCOL_Z = '#E57373'
#
# # ─────────────────────────────────────────────────────────────────────────────
# # FIGURE — 1 row × 3 panels
# #   A: Average waters per frame (bar)
# #   B: Average residence time (bar)
# #   C: Residence time distribution (box + strip — SELF EXPLANATORY)
# # ─────────────────────────────────────────────────────────────────────────────
# fig = plt.figure(figsize=(14, 5.5))
# fig.patch.set_facecolor('white')
#
# gs = gridspec.GridSpec(1, 3, figure=fig,
#                        left=0.07, right=0.97,
#                        top=0.87, bottom=0.14,
#                        wspace=0.38)
#
# fig.suptitle("Residence Time Inside the Trimer Pore",
#              fontsize=15, fontweight='bold', y=0.98,
#              color='#1B2D4F')
#
# # Shared x-tick labels and positions
# structs  = ['6IYX\n(Closed,\nCa²⁺-bound)', '6IYZ\n(Open,\nCa²⁺-free)']
# x_pos    = [0, 1]
# bar_w    = 0.5
# colors   = [BARCOL_X, BARCOL_Z]
# edgecols = [COL_X, COL_Z]
#
# def style_bar_ax(ax, title, ylabel, values, colors, edgecols, structs):
#     """Draw a clean bar chart on ax."""
#     bars = ax.bar(x_pos, values, bar_w,
#                   color=colors, edgecolor=edgecols,
#                   linewidth=1.4, zorder=3, alpha=0.88)
#     for bar, val in zip(bars, values):
#         ax.text(bar.get_x() + bar.get_width()/2,
#                 bar.get_height() + max(values)*0.02,
#                 f'{val:.2f}' if val < 100 else f'{val:.1f}',
#                 ha='center', va='bottom',
#                 fontsize=11, fontweight='bold',
#                 color='#1B2D4F')
#     ax.set_xticks(x_pos)
#     ax.set_xticklabels(structs, fontsize=9.5)
#     ax.set_ylabel(ylabel, fontsize=11, fontweight='bold', color='#1B2D4F')
#     ax.set_ylim(0, max(values) * 1.55)
#     ax.tick_params(axis='x', bottom=False)
#     ax.grid(axis='y', linestyle='--', alpha=0.3, zorder=0)
#     ax.set_title(title, loc='left', fontsize=13,
#                  fontweight='bold', pad=6, color='#1B2D4F')
#
# # ─────────────────────────────────────────────────────────────────────────────
# # PANEL A — Average waters per frame
# # ─────────────────────────────────────────────────────────────────────────────
# ax_a = fig.add_subplot(gs[0])
# style_bar_ax(ax_a,
#              title='A',
#              ylabel='Average waters per frame',
#              values=[sx['avg_per_frame'], sz['avg_per_frame']],
#              colors=colors, edgecols=edgecols,
#              structs=structs)
#
# # ─────────────────────────────────────────────────────────────────────────────
# # PANEL B — Average residence time (ps)
# # ─────────────────────────────────────────────────────────────────────────────
# ax_b = fig.add_subplot(gs[1])
# avg_res_vals = [sx['avg_residence'], sz['avg_residence']]
# style_bar_ax(ax_b,
#              title='B',
#              ylabel='Average residence time (ps)',
#              values=avg_res_vals,
#              colors=colors, edgecols=edgecols,
#              structs=structs)
#
# # Median dashed line inside each bar — white outline for visibility
# for xi, (med, col) in enumerate(
#         zip([sx['median_residence'], sz['median_residence']], edgecols)):
#     ax_b.plot([xi - bar_w/2 + 0.04, xi + bar_w/2 - 0.04],
#               [med, med],
#               color='white', linewidth=3.5, zorder=5, solid_capstyle='round')
#     ax_b.plot([xi - bar_w/2 + 0.04, xi + bar_w/2 - 0.04],
#               [med, med],
#               color=col, linewidth=2.0, zorder=6,
#               linestyle='--', solid_capstyle='round')
#     # Label the median value next to the line
#     ax_b.text(xi + bar_w/2 - 0.05, med + max(avg_res_vals)*0.015,
#               f'Median: {med:.1f} ps',
#               ha='right', va='bottom',
#               fontsize=7.5, color=col, style='italic')
#
# # Longest residence time — clean text annotation at TOP of figure
# for xi, (lval, col) in enumerate(
#         zip([sx['longest'], sz['longest']], edgecols)):
#     ax_b.text(xi, max(avg_res_vals) * 1.28,
#               f'Longest event: {lval/1000:.2f} ns',
#               ha='center', va='center',
#               fontsize=8, color=col, fontweight='bold',
#               bbox=dict(boxstyle='round,pad=0.4',
#                         facecolor='white', edgecolor=col,
#                         alpha=0.95, linewidth=1.2))
#
# # ─────────────────────────────────────────────────────────────────────────────
# # PANEL C — Residence time distribution
# #   Use log-binned bar chart (not histogram) with CLEAR labels
# #   Bins: 10ps, 20-50ps, 50-100ps, 100-500ps, 500-1000ps, >1000ps
# #   This is self-explanatory with named x-axis bins
# # ─────────────────────────────────────────────────────────────────────────────
# ax_c = fig.add_subplot(gs[2])
#
# # Define named bins — human readable
# bin_edges  = [0, 10, 50, 100, 500, 1000, 5000, np.inf]
# bin_labels = ['≤10', '10–50', '50–100', '100–500',
#               '500–1000', '1–5 ns', '>5 ns']
# n_bins     = len(bin_labels)
# bx         = np.arange(n_bins)
# bar_bw     = 0.35
#
# def bin_durations(durations, bin_edges):
#     counts = []
#     for i in range(len(bin_edges)-1):
#         lo, hi = bin_edges[i], bin_edges[i+1]
#         if i == 0:
#             c = np.sum(durations <= hi)
#         else:
#             c = np.sum((durations > lo) & (durations <= hi))
#         counts.append(c)
#     return np.array(counts)
#
# cnt_x = bin_durations(dur_x, bin_edges) if len(dur_x) > 0 else np.zeros(n_bins)
# cnt_z = bin_durations(dur_z, bin_edges) if len(dur_z) > 0 else np.zeros(n_bins)
#
# # Normalise to percentage of total events for fair comparison
# pct_x = cnt_x / cnt_x.sum() * 100 if cnt_x.sum() > 0 else cnt_x
# pct_z = cnt_z / cnt_z.sum() * 100 if cnt_z.sum() > 0 else cnt_z
#
# bars_x = ax_c.bar(bx - bar_bw/2, cnt_x, bar_bw,
#                   color=BARCOL_X, edgecolor=COL_X,
#                   linewidth=1.2, alpha=0.88,
#                   label='6IYX (Closed, Ca²⁺-bound)', zorder=3)
# bars_z = ax_c.bar(bx + bar_bw/2, cnt_z, bar_bw,
#                   color=BARCOL_Z, edgecolor=COL_Z,
#                   linewidth=1.2, alpha=0.88,
#                   label='6IYZ (Open, Ca²⁺-free)', zorder=3)
#
# ax_c.set_xticks(bx)
# ax_c.set_xticklabels(bin_labels, fontsize=8.5, rotation=30, ha='right')
# ax_c.set_xlabel('Residence time (ps)', fontsize=11,
#                 fontweight='bold', labelpad=4)
# ax_c.set_ylabel('Number of residence events', fontsize=11,
#                 fontweight='bold', color='#1B2D4F')
# ax_c.set_title('C', loc='left', fontsize=13,
#                fontweight='bold', pad=6, color='#1B2D4F')
# ax_c.grid(axis='y', linestyle='--', alpha=0.3, zorder=0)
# ax_c.tick_params(axis='x', bottom=False)
#
# # Legend inside panel — top right, no overlap with bars
# patch_x = mpatches.Patch(facecolor=BARCOL_X, edgecolor=COL_X,
#                           linewidth=1.2, alpha=0.88,
#                           label='6IYX (Closed, Ca²⁺-bound)')
# patch_z = mpatches.Patch(facecolor=BARCOL_Z, edgecolor=COL_Z,
#                           linewidth=1.2, alpha=0.88,
#                           label='6IYZ (Open, Ca²⁺-free)')
# ax_c.legend(handles=[patch_x, patch_z],
#             fontsize=8.5, loc='upper right',
#             framealpha=0.92, edgecolor='#CCCCCC',
#             handlelength=1.2, handletextpad=0.5)
#
# # Subtitle explaining the panel
# ax_c.text(0.5, -0.22,
#           'Each bar shows number of water molecules with that residence time.',
#           transform=ax_c.transAxes,
#           ha='center', va='top', fontsize=8,
#           color='#555555', style='italic')
#
# # ─────────────────────────────────────────────────────────────────────────────
# # SHARED FOOTER NOTE
# # ─────────────────────────────────────────────────────────────────────────────
# # fig.text(0.5, 0.01,
# #          'Cylinder definition: radius = 5 Å, Z = protein COM ± 10 Å  |  '
# #          'TIP3P water oxygen (OH2) only  |  Time step = 10 ps',
# #          ha='center', va='bottom', fontsize=8,
# #          color='#888888', style='italic')
#
# # ─────────────────────────────────────────────────────────────────────────────
# # SAVE
# # ─────────────────────────────────────────────────────────────────────────────
# plt.savefig(OUT_FILE, dpi=300, bbox_inches='tight',
#             facecolor='white', edgecolor='none')
# print(f"\nFigure saved to:\n  {OUT_FILE}")
# plt.close()
# ############ fig residence time #######
#
##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hydration profile along the TRIC translocation pathway (trimer average), 5WUC,
with AUTOCORRELATION-AWARE error bars (block averaging / Flyvbjerg-Petersen).

Monomer center = per-frame geometric center of the monomer's protein atoms
  == VMD: measure center [atomselect "protein and fragment N"].

Bins (radius 8 A, contiguous): central [+/-2.5]; above 1-3 and below 1-3, height 5 A.
Sampling: every 250 ps over the last 500 ns. Water = TIP3 OH2.

Error bars (ERROR_MODE):
  'block_sem'  (default, recommended): per monomer/bin SEM from block averaging
               (accounts for time autocorrelation); trimer SEM = propagated over monomers.
  'monomer_std': std across the 3 monomer means (old behaviour).
  'frame_std'  : raw frame-to-frame std (a SPREAD, not an error on the mean).
  'sem_naive'  : std/sqrt(N_frames)  -> UNDERESTIMATES (ignores autocorrelation).

Requires: MDAnalysis, numpy, matplotlib
"""

import os
import numpy as np
import MDAnalysis as mda
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIG ============================
GRO = "/media/supremeleader/Pantera/simulation/lipid_simulation/6IZF/charmm-gui-7977903380/gromacs/6izf/center.gro"
XTC = "/media/supremeleader/Pantera/simulation/lipid_simulation/6IZF/charmm-gui-7977903380/gromacs/6izf/center.xtc"
OUTDIR = "/media/supremeleader/Pantera/simulation/apo_halo/water_analysis"

WATER_SEL = "resname TIP3 and name OH2"
N_CHAINS  = 3
RADIUS    = 8.0
BIN_H     = 5.0
N_ABOVE   = 3
N_BELOW   = 3
SAMPLE_PS = 250.0
LAST_NS   = 500.0
PER_MONOMER_FIGS = True
ERROR_MODE = "block_sem"          # 'block_sem' | 'monomer_std' | 'frame_std' | 'sem_naive'
BLOCKING_DIAGNOSTIC = True         # save SEM-vs-blocksize curve for monomer1/central

# x-axis bin labels, one per bar, in plotting order LEFT -> RIGHT.
# Bars run +z (top/SR-lumen side) on the LEFT  to  -z (cytoplasm) on the RIGHT.
# These are just cosmetic annotations; edit any string. (SR lumen = eukaryotic SR TRIC, e.g. 6IYX/6IZF.)
LABELS = [
    "SR lumen",            # above 3  (farthest +z)
    "Lumenal\nside",       # above 2
    "Lumenal\nside",       # above 1
    "Central\ncavity",     # central
    "Cytoplasmic\nside",   # below 1
    "Cytoplasmic\nside",   # below 2
    "Cytoplasm",           # below 3  (farthest -z)
]

AA = ("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL "
      "HSD HSE HSP HID HIE HIP HISD HISE HISP CYX CYM ASH GLH LYN ARN").split()
# ================================================================

R2 = RADIUS * RADIUS
OFFSETS = np.concatenate([np.arange(N_ABOVE, 0, -1), [0], -np.arange(1, N_BELOW + 1)]) * BIN_H
NBIN = len(OFFSETS)
if len(LABELS) != NBIN:
    raise ValueError(f"LABELS has {len(LABELS)} entries but there are {NBIN} bins.")


def make_colors():
    cols, bl, gr = [], plt.cm.Blues, plt.cm.Greens
    for off in OFFSETS:
        k = int(round(off / BIN_H))
        if k == 0:    cols.append("#D62728")
        elif k > 0:   cols.append(bl(0.30 + 0.55 * (k / N_ABOVE)))
        else:         cols.append(gr(0.30 + 0.55 * (abs(k) / N_BELOW)))
    return cols


COLORS = make_colors()


def block_sem(x, min_blocks=8):
    """Flyvbjerg-Petersen block-averaging SEM.
    Returns (plateau_sem, blocksizes, sems, sem_errs).
    Repeatedly pairwise-averages the series; the SEM estimate rises then plateaus.
    The plateau (taken as the max SEM among levels with >= min_blocks blocks) is the
    autocorrelation-corrected standard error of the mean."""
    x = np.asarray(x, float)
    n0 = len(x)
    if n0 < 4:
        return (x.std(ddof=1) / np.sqrt(n0)) if n0 > 1 else 0.0, [], [], []
    sems, errs, bsizes, nblocks = [], [], [], []
    data, bs = x.copy(), 1
    while len(data) >= 4:
        m = len(data)
        var = data.var(ddof=1)
        sem = np.sqrt(var / m)
        sems.append(sem); errs.append(sem / np.sqrt(2 * (m - 1)))
        bsizes.append(bs); nblocks.append(m)
        if m % 2:
            data = data[:-1]
        data = 0.5 * (data[0::2] + data[1::2])
        bs *= 2
    cand = [s for s, nb in zip(sems, nblocks) if nb >= min_blocks]
    plateau = max(cand) if cand else max(sems)
    return plateau, bsizes, sems, errs


def monomer_atomgroups(u):
    res = u.select_atoms("name CA and (resname " + " ".join(AA) + ")").residues
    n = len(res)
    if n == 0 or n % N_CHAINS != 0:
        raise ValueError(f"{n} protein residues not divisible by {N_CHAINS}.")
    per = n // N_CHAINS
    return [res[k * per:(k + 1) * per].atoms for k in range(N_CHAINS)]


def frame_window(u):
    n = len(u.trajectory)
    t0 = u.trajectory[0].time
    dt = (u.trajectory[1].time - t0) if n > 1 else SAMPLE_PS
    last = u.trajectory[-1].time
    step = max(1, int(round(SAMPLE_PS / dt))) if dt > 0 else 1
    start_time = last - LAST_NS * 1000.0
    if start_time < t0:
        print(f"  ! trajectory < {LAST_NS} ns (total {(last-t0)/1000:.1f} ns); using all.")
        start_time = t0
    return max(0, int(np.ceil((start_time - t0) / dt))), step, dt


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    u = mda.Universe(GRO, XTC)
    monomers = monomer_atomgroups(u)
    water = u.select_atoms(WATER_SEL)
    if len(water) == 0:
        raise ValueError(f"No water from '{WATER_SEL}'.")
    start_frame, step, dt = frame_window(u)

    counts, times = [], []
    for ts in u.trajectory[start_frame::step]:
        wpos = water.positions
        fc = np.zeros((N_CHAINS, NBIN), dtype=int)
        for m, ag in enumerate(monomers):
            c = ag.center_of_geometry()
            dx = wpos[:, 0] - c[0]; dy = wpos[:, 1] - c[1]
            in_r = (dx * dx + dy * dy) < R2
            dz = wpos[in_r, 2] - c[2]
            for b, off in enumerate(OFFSETS):
                fc[m, b] = np.count_nonzero((dz >= off - BIN_H / 2) & (dz < off + BIN_H / 2))
        counts.append(fc); times.append(ts.time)

    counts = np.array(counts)                 # (n_frames, n_mon, n_bin)
    nfr = counts.shape[0]
    print(f"{len(water)} waters, {N_CHAINS} monomers. {nfr} frames "
          f"{times[0]/1000:.1f}-{times[-1]/1000:.1f} ns (every {step*dt:.0f} ps). "
          f"error mode = {ERROR_MODE}")

    mon_mean = counts.mean(axis=0)            # (n_mon, n_bin)
    trimer_mean = mon_mean.mean(axis=0)

    # per-monomer & trimer errors per chosen mode
    mon_block_sem = np.array([[block_sem(counts[:, m, b])[0] for b in range(NBIN)]
                              for m in range(N_CHAINS)])
    if ERROR_MODE == "block_sem":
        mon_err = mon_block_sem
        trimer_err = np.sqrt((mon_block_sem ** 2).sum(axis=0)) / N_CHAINS   # propagated
    elif ERROR_MODE == "frame_std":
        mon_err = counts.std(axis=0, ddof=1)
        trimer_err = mon_mean.std(axis=0, ddof=1)
    elif ERROR_MODE == "sem_naive":
        mon_err = counts.std(axis=0, ddof=1) / np.sqrt(nfr)
        trimer_err = np.sqrt((mon_err ** 2).sum(axis=0)) / N_CHAINS
    else:  # monomer_std
        mon_err = counts.std(axis=0, ddof=1)
        trimer_err = mon_mean.std(axis=0, ddof=1)

    # ---- save data (always include block-SEM for reference) ----
    with open(os.path.join(OUTDIR, "hydration_averages.dat"), "w") as f:
        f.write(f"# radius={RADIUS} bin_h={BIN_H} window=last {LAST_NS}ns every {SAMPLE_PS}ps "
                f"nframes={nfr} error_mode={ERROR_MODE}\n")
        f.write("# bin label mean err  block_sem  frame_std\n")
        for m in range(N_CHAINS):
            f.write(f"# monomer {m+1}\n")
            for b in range(NBIN):
                f.write(f"{b} {LABELS[b].replace(chr(10),' ')} {mon_mean[m,b]:.3f} "
                        f"{mon_err[m,b]:.3f} {mon_block_sem[m,b]:.3f} "
                        f"{counts[:,m,b].std(ddof=1):.3f}\n")
        f.write("# TRIMER\n")
        for b in range(NBIN):
            f.write(f"{b} {LABELS[b].replace(chr(10),' ')} {trimer_mean[b]:.3f} {trimer_err[b]:.3f}\n")

    # ---- blocking diagnostic (verify the plateau) ----
    if BLOCKING_DIAGNOSTIC:
        cb = int(np.where(OFFSETS == 0)[0][0])
        _, bsizes, sems, errs = block_sem(counts[:, 0, cb])
        if bsizes:
            fig, ax = plt.subplots(figsize=(6, 4))
            tb = np.array(bsizes) * step * dt / 1000.0     # block length in ns
            ax.errorbar(tb, sems, yerr=errs, fmt="o-", capsize=3, color="#444")
            ax.set_xscale("log")
            ax.set_xlabel("Block length (ns)"); ax.set_ylabel("SEM estimate (waters)")
            ax.set_title("Block-averaging convergence (monomer 1, central bin)")
            fig.tight_layout()
            for ext in ("png", "pdf"):
                fig.savefig(os.path.join(OUTDIR, f"blocking_curve.{ext}"), dpi=300, bbox_inches="tight")
            plt.close(fig)
            print("  saved blocking_curve (check that SEM plateaus)")

    # ---- histograms ----
    def histogram(mean_vec, err_vec, title, outname):
        plt.rcParams.update({"font.family": "sans-serif", "font.size": 12})
        fig, ax = plt.subplots(figsize=(8, 5))
        x = np.arange(NBIN)
        ax.bar(x, mean_vec, yerr=err_vec, capsize=5, color=COLORS,
               width=0.72, edgecolor="0.25", linewidth=0.6)
        ax.set_xlabel("Channel coordinate")
        ax.set_ylabel("Average number of water molecules")
        ax.set_title(title)
        ax.set_xticks(x); ax.set_xticklabels(LABELS, rotation=30, ha="right")
        ax.margins(x=0.01)
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(OUTDIR, f"{outname}.{ext}"), dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  saved {outname}")

    histogram(trimer_mean, trimer_err,
              f"Hydration of the translocation pathway \u2013 6IZF trimer (last {LAST_NS:.0f} ns)",
              "hydration_trimer")
    if PER_MONOMER_FIGS:
        for m in range(N_CHAINS):
            histogram(mon_mean[m], mon_err[m],
                      f"Hydration of the translocation pathway \u2013 6IZF monomer {m+1}",
                      f"hydration_monomer{m+1}")
    print(f"\nAll outputs in: {OUTDIR}")


if __name__ == "__main__":
    main()
#

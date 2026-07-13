# import os
#
# # ─── Paths ────────────────────────────────────────────────────────────────────
# INPUT_FILE  = "/home/supremeleader/Desktop/article/seq_analysis/bacteria_seq/8_gram-ve_percentage_conservation.txt"
# OUTPUT_FILE = "/home/supremeleader/Desktop/article/seq_analysis/bacteria_seq/9_gram-ve_group_conservation.txt"
#
# # ─── Amino acid groups (Proline excluded) ────────────────────────────────────
# GROUPS = {
#     "GSACT" : ["G", "S", "A", "C", "T"],
#     "YFW"   : ["Y", "F", "W"],
#     "NQ"    : ["N", "Q"],
#     "DE"    : ["D", "E"],
#     "KRH"   : ["K", "R", "H"],
#     "LIVM"  : ["L", "I", "V", "M"],
#     "P"     : ["P"],
# }
#
# # Fixed column order from the input file
# AA_ORDER = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
#             "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-"]
#
# # Pre-build index lookup for speed
# AA_INDEX = {aa: i for i, aa in enumerate(AA_ORDER)}
#
#
# def sum_group(values, aa_list):
#     """Sum percentage values for a list of amino acids."""
#     return sum(values[AA_INDEX[aa]] for aa in aa_list)
#
#
# def format_group_line(gen_num, res_label, values):
#     """Build the output line for one residue position."""
#     parts = [gen_num, res_label]
#     for grp_name, aa_list in GROUPS.items():
#         total = sum_group(values, aa_list)
#         parts.append(f" {grp_name}= {round(total, 2)}")
#     return ",".join(parts)
#
#
# def main():
#     print(f"Reading : {INPUT_FILE}")
#
#     with open(INPUT_FILE) as fh:
#         raw_lines = [l.rstrip("\n") for l in fh]
#
#     out_lines  = []
#     data_rows  = 0
#     header_written = False   # tracks whether we've written first header yet
#
#     # Output header template (no AA columns — just positional labels)
#     HEADER     = "generic_numbering,residues_numbering,,,,,,,,"
#     BLANK      = ",,,,,,,,"
#
#     for line in raw_lines:
#
#         # ── Skip the input file's own header line ────────────────────────────
#         if line.startswith("generic_n"):
#             # Don't echo the input header — we write our own
#             continue
#
#         # ── Blank / separator rows ────────────────────────────────────────────
#         if set(line.replace(",", "").strip()) == set():   # all commas or empty
#             continue   # we'll emit blanks ourselves around TM headers
#
#         fields = line.split(",")
#
#         # ── TM helix header row  e.g. "TM1,,,,,..." ──────────────────────────
#         if fields[0].startswith("TM"):
#             tm_name = fields[0]
#             # Reprint full header block before every TM section
#             if header_written:
#                 out_lines.append("")          # blank separator
#             out_lines.append(HEADER)
#             out_lines.append(BLANK)
#             out_lines.append(tm_name + "," * (len(BLANK.split(",")) - 1))
#             header_written = True
#             continue
#
#         # ── Data row  e.g. "1.36,D10,6.13,0.00,..." ─────────────────────────
#         # Must have at least 2 label fields + 21 AA value fields
#         if len(fields) < 23:
#             continue
#
#         gen_num   = fields[0].strip()
#         res_label = fields[1].strip()
#
#         # Skip rows where first field doesn't look like a generic number
#         if not gen_num or not any(c.isdigit() for c in gen_num):
#             continue
#
#         try:
#             values = [float(fields[i + 2]) for i in range(len(AA_ORDER))]
#         except (ValueError, IndexError):
#             continue
#
#         out_lines.append(format_group_line(gen_num, res_label, values))
#         data_rows += 1
#
#     # ── Write output ──────────────────────────────────────────────────────────
#     os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
#     with open(OUTPUT_FILE, "w") as fh:
#         fh.write("\n".join(out_lines) + "\n")
#
#     print(f"  Data rows written : {data_rows}")
#     print(f"  Output            : {OUTPUT_FILE}")
#     print("Done.")
#
#
# if __name__ == "__main__":
#     main()



########################################################################################################
########################## overall code for bacterial cluster analysis #####################
#########################################################################################################
# """
# batch_conservation.py
# ──────────────────────
# For every cluster alignment file in the input directory, calculates:
#   1. % conservation of each amino acid at every TM position
#   2. Group conservation (GSACT, YFW, NQ, DE, KRH, LIVM, P)
#
# Output per cluster (saved in same directory):
#   clusterX_percentage_conservation.txt
#   clusterX_group_conservation.txt
#
# Directory : /home/supremeleader/Desktop/article/phylogenetic_analysis/gram_negative/cluster_sequence
# """
#
# import os
# import glob
# from collections import Counter
#
# # ─── Directory ────────────────────────────────────────────────────────────────
# CLUSTER_DIR = "/home/supremeleader/Desktop/article/phylogenetic_analysis/gram_positive/cluster_sequence"
#
# # ─── Column order ─────────────────────────────────────────────────────────────
# AA_ORDER = ["A","C","D","E","F","G","H","I","K","L",
#             "M","N","P","Q","R","S","T","V","W","Y","-"]
#
# # ─── Amino acid groups ────────────────────────────────────────────────────────
# GROUPS = {
#     "GSACT" : ["G","S","A","C","T"],
#     "YFW"   : ["Y","F","W"],
#     "NQ"    : ["N","Q"],
#     "DE"    : ["D","E"],
#     "KRH"   : ["K","R","H"],
#     "LIVM"  : ["L","I","V","M"],
#     "P"     : ["P"],
# }
#
# AA_INDEX = {aa: i for i, aa in enumerate(AA_ORDER)}
#
# # ─── Positions: (generic_numbering, residue_label) ───────────────────────────
# # Listed in exact MSA column order — section headers use (TM_name, None)
# POSITIONS = [
#     # TM1
#     ("TM1",  None),
#     ("1.36", "D10"),  ("1.37", "I11"),  ("1.38", "F12"),
#     ("1.39", "N13"),  ("1.4",  "Y14"),  ("1.41", "I15"),
#     ("1.42", "G16"),  ("1.43", "I17"),  ("1.44", "V18"),
#     ("1.45", "A19"),  ("1.46", "F20"),  ("1.47", "A21"),
#     ("1.48", "I22"),  ("1.49", "S23"),  ("1.5",  "G24"),
#     ("1.51", "A25"),  ("1.52", "I26"),  ("1.53", "K27"),
#     ("1.54", "A28"),  ("1.55", "V29"),  ("1.56", "K30"),
#     # TM2
#     ("TM2",  None),
#     ("2.36", "L35"),  ("2.37", "L36"),  ("2.38", "G37"),
#     ("2.39", "V38"),  ("2.4",  "L39"),  ("2.41", "V40"),
#     ("2.42", "L41"),  ("2.43", "G42"),  ("2.44", "F43"),
#     ("2.45", "S44"),  ("2.46", "T45"),  ("2.47", "A46"),
#     ("2.48", "L47"),  ("2.49", "G48"),  ("2.5",  "G49"),
#     ("2.51", "G50"),  ("2.52", "I51"),  ("2.53", "I52"),
#     ("2.54", "S53"),  ("2.55", "N54"),  ("2.56", "L55"),
#     ("2.57", "L56"),  ("2.58", "L57"),  ("2.59", "G58"),
#     ("2.6",  "K59"),
#     # TM3
#     ("TM3",  None),
#     ("3.45", "L65"),  ("3.46", "I66"),  ("3.47", "Y67"),
#     ("3.48", "Y68"),  ("3.49", "P69"),  ("3.5",  "Y70"),
#     ("3.51", "P71"),  ("3.52", "I72"),  ("3.53", "T73"),
#     ("3.54", "A74"),  ("3.55", "F75"),  ("3.56", "L76"),
#     ("3.57", "A77"),  ("3.58", "S78"),  ("3.59", "L79"),
#     ("3.6",  "A80"),  ("3.61", "T81"),  ("3.62", "F82"),
#     ("3.63", "V83"),  ("3.64", "F84"),
#     # TM4
#     ("TM4",  None),
#     ("4.43", "G92"),  ("4.44", "K93"),  ("4.45", "P94"),
#     ("4.46", "L95"),  ("4.47", "L96"),  ("4.48", "Y97"),
#     ("4.49", "A98"),  ("4.5",  "D99"),  ("4.51", "A100"),
#     ("4.52", "I101"), ("4.53", "G102"), ("4.54", "L103"),
#     ("4.55", "G104"), ("4.56", "A105"), ("4.57", "F106"),
#     ("4.58", "A107"), ("4.59", "S108"), ("4.6",  "S109"),
#     ("4.61", "G110"), ("4.62", "A111"), ("4.63", "S112"),
#     ("4.64", "L113"), ("4.65", "A114"), ("4.66", "Y115"),
#     ("4.67", "S116"),
#     # TM5
#     ("TM5",  None),
#     ("5.36", "V121"), ("5.37", "I122"), ("5.38", "L123"),
#     ("5.39", "V124"), ("5.4",  "V125"), ("5.41", "I126"),
#     ("5.42", "V127"), ("5.43", "G128"), ("5.44", "A129"),
#     ("5.45", "I130"), ("5.46", "T131"), ("5.47", "A132"),
#     ("5.48", "V133"), ("5.49", "G134"), ("5.5",  "G135"),
#     ("5.51", "G136"), ("5.52", "V137"), ("5.53", "I138"),
#     ("5.54", "R139"), ("5.55", "D140"), ("5.56", "I141"),
#     ("5.57", "L142"), ("5.58", "S143"), ("5.59", "N144"),
#     ("5.6",  "E145"),
#     # TM6
#     ("TM6",  None),
#     ("6.5",  "Y155"), ("6.51", "A156"), ("6.52", "T157"),
#     ("6.53", "T158"), ("6.54", "A159"), ("6.55", "V160"),
#     ("6.56", "I161"), ("6.57", "G162"), ("6.58", "S163"),
#     ("6.59", "F164"), ("6.6",  "V165"), ("6.61", "Y166"),
#     ("6.62", "F167"), ("6.63", "I168"), ("6.64", "A169"),
#     ("6.65", "S170"), ("6.66", "D171"), ("6.67", "L172"),
#     ("6.68", "S173"), ("6.69", "V174"),
#     # TM7
#     ("TM7",  None),
#     ("7.38", "A179"), ("7.39", "L180"), ("7.4",  "I181"),
#     ("7.41", "V182"), ("7.42", "S183"), ("7.43", "F184"),
#     ("7.44", "L185"), ("7.45", "I186"), ("7.46", "T187"),
#     ("7.47", "L188"), ("7.48", "I189"), ("7.49", "L190"),
#     ("7.5",  "R191"), ("7.51", "I192"), ("7.52", "L193"),
#     ("7.53", "A194"), ("7.54", "M195"), ("7.55", "E196"),
#     ("7.56", "L197"), ("7.57", "K198"), ("7.58", "W199"),
# ]
#
# BLANK  = "," * (len(AA_ORDER) + 1)
# HEADER = "generic_numbering,residues_numbering," + ",".join(AA_ORDER)
#
# GRP_BLANK  = ",,,,,,,,"
# GRP_HEADER = "generic_numbering,residues_numbering,,,,,,,,"
#
#
# # ─── Step 1: calculate % conservation ────────────────────────────────────────
#
# def calc_pct_conservation(input_file):
#     """Read FASTA alignment, return list of output lines (% conservation)."""
#     seqs = []
#     current = []
#     with open(input_file) as fh:
#         for line in fh:
#             line = line.rstrip()
#             if line.startswith(">"):
#                 if current:
#                     seqs.append("".join(current))
#                 current = []
#             else:
#                 current.append(line)
#     if current:
#         seqs.append("".join(current))
#
#     n_seqs = len(seqs)
#     out    = [HEADER, BLANK]
#     col    = 0
#     first_tm = True
#
#     for gen_num, res_label in POSITIONS:
#         if res_label is None:
#             if not first_tm:
#                 out.append(BLANK)
#                 out.append(HEADER)
#                 out.append(BLANK)
#             first_tm = False
#             out.append(gen_num + "," * (len(AA_ORDER) + 1))
#             continue
#
#         counts = Counter(seq[col].upper() for seq in seqs)
#         pct    = [f"{counts.get(aa, 0) / n_seqs * 100:.2f}" for aa in AA_ORDER]
#         out.append(f"{gen_num},{res_label}," + ",".join(pct))
#         col += 1
#
#     out.append(BLANK)
#     return n_seqs, col, out
#
#
# # ─── Step 2: calculate group conservation from % conservation lines ───────────
#
# def calc_group_conservation(pct_lines):
#     """Take % conservation lines, return group conservation output lines."""
#     out           = []
#     header_written = False
#
#     for line in pct_lines:
#         if line.startswith("generic_n"):
#             continue
#         if set(line.replace(",", "").strip()) == set():
#             continue
#
#         fields = line.split(",")
#
#         if fields[0].startswith("TM"):
#             if header_written:
#                 out.append("")
#             out.append(GRP_HEADER)
#             out.append(GRP_BLANK)
#             out.append(fields[0] + "," * (len(GRP_BLANK.split(",")) - 1))
#             header_written = True
#             continue
#
#         if len(fields) < 23:
#             continue
#
#         gen_num   = fields[0].strip()
#         res_label = fields[1].strip()
#
#         if not gen_num or not any(c.isdigit() for c in gen_num):
#             continue
#
#         try:
#             values = [float(fields[i + 2]) for i in range(len(AA_ORDER))]
#         except (ValueError, IndexError):
#             continue
#
#         parts = [gen_num, res_label]
#         for grp_name, aa_list in GROUPS.items():
#             total = sum(values[AA_INDEX[aa]] for aa in aa_list)
#             parts.append(f" {grp_name}= {round(total, 2)}")
#         out.append(",".join(parts))
#
#     return out
#
#
# # ─── Main: process all cluster files ─────────────────────────────────────────
#
# def main():
#     pattern = os.path.join(CLUSTER_DIR, "cluster*_alignment.txt")
#     files   = sorted(glob.glob(pattern))
#
#     if not files:
#         print(f"No cluster alignment files found in:\n  {CLUSTER_DIR}")
#         return
#
#     print(f"Found {len(files)} alignment file(s):\n")
#
#     for aln_file in files:
#         basename = os.path.basename(aln_file)                        # cluster1_alignment.txt
#         prefix   = basename.replace("_alignment.txt", "")            # cluster1
#         pct_out  = os.path.join(CLUSTER_DIR, f"{prefix}_percentage_conservation.txt")
#         grp_out  = os.path.join(CLUSTER_DIR, f"{prefix}_group_conservation.txt")
#
#         print(f"  Processing : {basename}")
#
#         # Step 1 — % conservation
#         n_seqs, n_rows, pct_lines = calc_pct_conservation(aln_file)
#         with open(pct_out, "w") as fh:
#             fh.write("\n".join(pct_lines) + "\n")
#         print(f"    Sequences : {n_seqs}  |  Positions : {n_rows}")
#         print(f"    → {os.path.basename(pct_out)}")
#
#         # Step 2 — group conservation
#         grp_lines = calc_group_conservation(pct_lines)
#         with open(grp_out, "w") as fh:
#             fh.write("\n".join(grp_lines) + "\n")
#         print(f"    → {os.path.basename(grp_out)}")
#
#     print(f"\nDone. All outputs saved to:\n  {CLUSTER_DIR}")
#
#
# if __name__ == "__main__":
#     main()
#
#
# ########################################################################################################
# ########################## overall code for bacterial cluster analysis #####################
# #########################################################################################################

# ########################################################################################################
# ########################## overall code for eukaryote cluster analysis #####################
# #########################################################################################################

"""
batch_euk_conservation.py
──────────────────────────
For every cluster alignment file in the input directory, calculates:
  1. % conservation of each amino acid at every TM position
  2. Group conservation (GSACT, YFW, NQ, DE, KRH, LIVM, P)

Output per cluster (saved in same directory):
  clusterX_percentage_conservation.txt
  clusterX_group_conservation.txt

Directory : /home/supremeleader/Desktop/article/phylogenetic_analysis/eukaryote/cluster_sequence
"""

import os
import glob
from collections import Counter

# ─── Directory ────────────────────────────────────────────────────────────────
CLUSTER_DIR = "/home/supremeleader/Desktop/article/phylogenetic_analysis/eukaryote/cluster_sequence"

# ─── Column order ─────────────────────────────────────────────────────────────
AA_ORDER = ["A","C","D","E","F","G","H","I","K","L",
            "M","N","P","Q","R","S","T","V","W","Y","-"]

# ─── Amino acid groups ────────────────────────────────────────────────────────
GROUPS = {
    "GSACT" : ["G","S","A","C","T"],
    "YFW"   : ["Y","F","W"],
    "NQ"    : ["N","Q"],
    "DE"    : ["D","E"],
    "KRH"   : ["K","R","H"],
    "LIVM"  : ["L","I","V","M"],
    "P"     : ["P"],
}

AA_INDEX = {aa: i for i, aa in enumerate(AA_ORDER)}

# ─── Positions: (generic_numbering, residue_label) ───────────────────────────
# Listed in exact MSA column order — section headers use (TM_name, None)
POSITIONS = [
    # TM1
    ("TM1",   None),
    ("1.48",  "P23"),  ("1.49",  "L24"),  ("1.50",  "F25"),
    ("1.51",  "D26"),  ("1.52",  "A27"),  ("1.53",  "A28"),
    ("1.54",  "Y29"),  ("1.55",  "F30"),  ("1.56",  "I31"),
    ("1.57",  "V32"),  ("1.58",  "S33"),  ("1.59",  "V34"),
    ("1.60",  "L35"),  ("1.61",  "Y36"),  ("1.62",  "L37"),
    ("1.63",  "Y38"),  ("1.64",  "L39"),  ("1.65",  "E40"),
    ("1.66",  "P41"),
    # TM2
    ("TM2",   None),
    ("2.38",  "P51"),  ("2.39",  "F52"),  ("2.40",  "A53"),
    ("2.41",  "S54"),  ("2.42",  "W55"),  ("2.43",  "L56"),
    ("2.44",  "C57"),  ("2.45",  "A58"),  ("2.46",  "M59"),
    ("2.47",  "L60"),  ("2.48",  "H61"),  ("2.49",  "C62"),
    ("2.50",  "F63"),  ("2.51",  "G64"),  ("2.52",  "S65"),
    ("2.53",  "Y66"),  ("2.54",  "I67"),  ("2.55",  "L68"),
    ("2.56",  "A69"),  ("2.57",  "D70"),  ("2.58",  "L71"),
    ("2.59",  "L72"),  ("2.60",  "L73"),  ("2.61",  "G74"),
    # TM3
    ("TM3",   None),
    ("3.37",  "S82"),  ("3.38",  "N83"),  ("3.39",  "N84"),
    ("3.40",  "S85"),  ("3.41",  "S86"),  ("3.42",  "V87"),
    ("3.43",  "I88"),  ("3.44",  "L89"),  ("3.45",  "A90"),
    ("3.46",  "T91"),  ("3.47",  "A92"),  ("3.48",  "V93"),
    ("3.49",  "W94"),  ("3.50",  "Y95"),  ("3.51",  "L96"),
    ("3.52",  "I97"),  ("3.53",  "F98"),  ("3.54",  "F99"),
    ("3.55",  "C100"),
    # TM4
    ("TM4",   None),
    ("4.37",  "P113"), ("4.38",  "V114"), ("4.39",  "K115"),
    ("4.40",  "L116"), ("4.41",  "I117"), ("4.42",  "F118"),
    ("4.43",  "V119"), ("4.44",  "A120"), ("4.45",  "M121"),
    ("4.46",  "K122"), ("4.47",  "E123"), ("4.48",  "V124"),
    ("4.49",  "V125"), ("4.50",  "R126"), ("4.51",  "V127"),
    ("4.52",  "R128"), ("4.53",  "K129"), ("4.54",  "I130"),
    ("4.55",  "A131"), ("4.56",  "A132"), ("4.57",  "G133"),
    ("4.58",  "V134"), ("4.59",  "H135"), ("4.60",  "H136"),
    ("4.61",  "A137"), ("4.62",  "H138"), ("4.63",  "H139"),
    ("4.64",  "Q140"),
    # TM5
    ("TM5",   None),
    ("5.44",  "F146"), ("5.45",  "I147"), ("5.46",  "M148"),
    ("5.47",  "M149"), ("5.48",  "A150"), ("5.49",  "T151"),
    ("5.50",  "G152"), ("5.51",  "W153"), ("5.52",  "V154"),
    ("5.53",  "K155"), ("5.54",  "G156"), ("5.55",  "S157"),
    ("5.56",  "G158"), ("5.57",  "V159"), ("5.58",  "A160"),
    ("5.59",  "L161"), ("5.60",  "M162"), ("5.61",  "S163"),
    ("5.62",  "N164"), ("5.63",  "F165"), ("5.64",  "E166"),
    ("5.65",  "Q167"), ("5.66",  "L168"), ("5.67",  "L169"),
    ("5.68",  "R170"),
    # TM6
    ("TM6",   None),
    ("6.46",  "S184"), ("6.47",  "F185"), ("6.48",  "P186"),
    ("6.49",  "T187"), ("6.50",  "K188"), ("6.51",  "A189"),
    ("6.52",  "S190"), ("6.53",  "L191"), ("6.54",  "Y192"),
    ("6.55",  "G193"), ("6.56",  "T194"), ("6.57",  "V195"),
    ("6.58",  "L196"), ("6.59",  "F197"), ("6.60",  "T198"),
    ("6.61",  "L199"), ("6.62",  "Q200"), ("6.63",  "Q201"),
    ("6.64",  "T202"),
    # TM7
    ("TM7",   None),
    ("7.46",  "S208"), ("7.47",  "E209"), ("7.48",  "A210"),
    ("7.49",  "N211"), ("7.50",  "L212"), ("7.51",  "V213"),
    ("7.52",  "F214"), ("7.53",  "F215"), ("7.54",  "F216"),
    ("7.55",  "T217"), ("7.56",  "M218"), ("7.57",  "F219"),
    ("7.58",  "M220"), ("7.59",  "I221"), ("7.60",  "V222"),
    ("7.61",  "C223"), ("7.62",  "K224"), ("7.63",  "V225"),
    ("7.64",  "F226"), ("7.65",  "M227"), ("7.66",  "T228"),
    ("7.67",  "A229"), ("7.68",  "T230"), ("7.69",  "H231"),
]

BLANK  = "," * (len(AA_ORDER) + 1)
HEADER = "generic_numbering,residues_numbering," + ",".join(AA_ORDER)

GRP_BLANK  = ",,,,,,,,"
GRP_HEADER = "generic_numbering,residues_numbering,,,,,,,,"


# ─── Step 1: calculate % conservation ────────────────────────────────────────

def calc_pct_conservation(input_file):
    """Read FASTA alignment, return (n_seqs, n_cols, list of output lines)."""
    seqs = []
    current = []
    with open(input_file) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current:
                    seqs.append("".join(current))
                current = []
            else:
                current.append(line)
    if current:
        seqs.append("".join(current))

    n_seqs   = len(seqs)
    out      = [HEADER, BLANK]
    col      = 0
    first_tm = True

    for gen_num, res_label in POSITIONS:
        if res_label is None:
            if not first_tm:
                out.append(BLANK)
                out.append(HEADER)
                out.append(BLANK)
            first_tm = False
            out.append(gen_num + "," * (len(AA_ORDER) + 1))
            continue

        counts = Counter(seq[col].upper() for seq in seqs)
        pct    = [f"{counts.get(aa, 0) / n_seqs * 100:.2f}" for aa in AA_ORDER]
        out.append(f"{gen_num},{res_label}," + ",".join(pct))
        col += 1

    out.append(BLANK)
    return n_seqs, col, out


# ─── Step 2: calculate group conservation from % conservation lines ───────────

def calc_group_conservation(pct_lines):
    """Take % conservation lines, return group conservation output lines."""
    out            = []
    header_written = False

    for line in pct_lines:
        if line.startswith("generic_n"):
            continue
        if set(line.replace(",", "").strip()) == set():
            continue

        fields = line.split(",")

        if fields[0].startswith("TM"):
            if header_written:
                out.append("")
            out.append(GRP_HEADER)
            out.append(GRP_BLANK)
            out.append(fields[0] + "," * (len(GRP_BLANK.split(",")) - 1))
            header_written = True
            continue

        if len(fields) < 23:
            continue

        gen_num   = fields[0].strip()
        res_label = fields[1].strip()

        if not gen_num or not any(c.isdigit() for c in gen_num):
            continue

        try:
            values = [float(fields[i + 2]) for i in range(len(AA_ORDER))]
        except (ValueError, IndexError):
            continue

        parts = [gen_num, res_label]
        for grp_name, aa_list in GROUPS.items():
            total = sum(values[AA_INDEX[aa]] for aa in aa_list)
            parts.append(f" {grp_name}= {round(total, 2)}")
        out.append(",".join(parts))

    return out


# ─── Main: process all cluster files ─────────────────────────────────────────

def main():
    pattern = os.path.join(CLUSTER_DIR, "cluster*_alignment.txt")
    files   = sorted(glob.glob(pattern))

    if not files:
        print(f"No cluster alignment files found in:\n  {CLUSTER_DIR}")
        return

    print(f"Found {len(files)} alignment file(s):\n")

    for aln_file in files:
        basename = os.path.basename(aln_file)                         # cluster1_alignment.txt
        prefix   = basename.replace("_alignment.txt", "")             # cluster1
        pct_out  = os.path.join(CLUSTER_DIR, f"{prefix}_percentage_conservation.txt")
        grp_out  = os.path.join(CLUSTER_DIR, f"{prefix}_group_conservation.txt")

        print(f"  Processing : {basename}")

        # Step 1 — % conservation
        n_seqs, n_rows, pct_lines = calc_pct_conservation(aln_file)
        with open(pct_out, "w") as fh:
            fh.write("\n".join(pct_lines) + "\n")
        print(f"    Sequences : {n_seqs}  |  Positions : {n_rows}")
        print(f"    → {os.path.basename(pct_out)}")

        # Step 2 — group conservation
        grp_lines = calc_group_conservation(pct_lines)
        with open(grp_out, "w") as fh:
            fh.write("\n".join(grp_lines) + "\n")
        print(f"    → {os.path.basename(grp_out)}")

    print(f"\nDone. All outputs saved to:\n  {CLUSTER_DIR}")


if __name__ == "__main__":
    main()
    
# ########################################################################################################
# ########################## overall code for eukaryote cluster analysis #####################
# #########################################################################################################

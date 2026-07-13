"""
prokaryote_motif_analysis.py
─────────────────────────────
Reads a prokaryotic TM-only FASTA alignment and extracts
functional residue combinations for each sequence, then
summarises their frequencies.

Reference structure : 5WUC
Input  : /home/supremeleader/Desktop/article/seq_analysis/prokaryote_analysis/1_prokar_overall_analysis/prokaryote_alignment.txt
Output : /home/supremeleader/Desktop/article/database/
            prokaryote_motif_per_sequence.txt   (file 1)
            prokaryote_motif_summary.txt         (file 2)

Helix boundaries (PDB residue numbering, loops removed):
    TM1  10-30   offset 10   cols   0-20
    TM2  35-59   offset 14   cols  21-45
    TM3  65-84   offset 19   cols  46-65
    TM4  92-116  offset 26   cols  66-90
    TM5 121-145  offset 30   cols  91-115
    TM6 155-174  offset 39   cols 116-135
    TM7 179-199  offset 43   cols 136-156

Key residue → MSA column (0-based):
    N13  → col   3      (N-THB gate,        TM1)
    F20  → col  10      (N-THB constriction, TM1)
    G48  → col  34      (N-THB kink,         TM2)
    G49  → col  35      (N-THB kink,         TM2)
    G50  → col  36      (N-THB kink,         TM2)
    S53  → col  39      (N-THB gate,         TM2)
    N54  → col  40      (N-THB gate,         TM2)
    Y68  → col  49      (N-THB gate,         TM3)
    D99  → col  73      (C-THB gate,         TM4)
    F106 → col  80      (C-THB constriction, TM4)
    G134 → col 104      (C-THB kink,         TM5)
    G135 → col 105      (C-THB kink,         TM5)
    G136 → col 106      (C-THB kink,         TM5)
    R139 → col 109      (C-THB gate,         TM5)
    D140 → col 110      (C-THB gate,         TM5)
    Y155 → col 116      (C-THB gate,         TM6)
"""

# import os
# from collections import Counter
#
# # ─── Paths ────────────────────────────────────────────────────────────────────
# INPUT_FILE  = "/home/supremeleader/Desktop/article/phylogenetic_analysis/gram_pos_tree/cluster_sequence/cluster4_alignment.txt"
# OUTPUT_DIR  = "/home/supremeleader/Desktop/article/database"
# OUT_FILE1   = os.path.join(OUTPUT_DIR, "G-P_cluster4_per_seq.txt")
# OUT_FILE2   = os.path.join(OUTPUT_DIR, "G-P_cluster4_summary.txt")
#
# # ─── MSA column indices (0-based) for each key residue ───────────────────────
# COL = {
#     # N-THB
#     "N13":  3,
#     "F20":  10,
#     "G48":  34,
#     "G49":  35,
#     "G50":  36,
#     "S53":  39,
#     "N54":  40,
#     "Y68":  49,
#     # C-THB
#     "D99":  73,
#     "F106": 80,
#     "G134": 104,
#     "G135": 105,
#     "G136": 106,
#     "R139": 109,
#     "D140": 110,
#     "Y155": 116,
# }
#
#
# # ─── Read FASTA ───────────────────────────────────────────────────────────────
# def read_fasta(path):
#     """Return list of (header, sequence) tuples."""
#     entries, current_hdr, current_seq = [], None, []
#     with open(path) as fh:
#         for line in fh:
#             line = line.rstrip()
#             if line.startswith(">"):
#                 if current_hdr is not None:
#                     entries.append((current_hdr, "".join(current_seq)))
#                 current_hdr, current_seq = line[1:].strip(), []
#             else:
#                 current_seq.append(line.upper())
#     if current_hdr is not None:
#         entries.append((current_hdr, "".join(current_seq)))
#     return entries
#
#
# # ─── Extract residue at a given column ───────────────────────────────────────
# def aa(seq, col):
#     """Return amino acid at MSA column (col), or '-' if out of bounds."""
#     return seq[col] if col < len(seq) else "-"
#
#
# # ─── Build motifs for one sequence ───────────────────────────────────────────
# def extract_motifs(seq):
#     """Return dict with all six motif strings for this sequence."""
#     return {
#         "nthb_kink":         aa(seq, COL["G48"])  + aa(seq, COL["G49"])  + aa(seq, COL["G50"]),
#         "nthb_gate":         aa(seq, COL["N13"])  + aa(seq, COL["S53"])  + aa(seq, COL["N54"])  + aa(seq, COL["Y68"]),
#         "nthb_constriction": aa(seq, COL["F20"]),
#         "cthb_kink":         aa(seq, COL["G134"]) + aa(seq, COL["G135"]) + aa(seq, COL["G136"]),
#         "cthb_gate":         aa(seq, COL["D99"])  + aa(seq, COL["R139"]) + aa(seq, COL["D140"]) + aa(seq, COL["Y155"]),
#         "cthb_constriction": aa(seq, COL["F106"]),
#     }
#
#
# # ─── Main ─────────────────────────────────────────────────────────────────────
# def main():
#     print(f"Reading : {INPUT_FILE}")
#     entries = read_fasta(INPUT_FILE)
#     n       = len(entries)
#     print(f"  {n} sequences found.")
#
#     os.makedirs(OUTPUT_DIR, exist_ok=True)
#
#     # ── Counters for summary ──────────────────────────────────────────────────
#     cnt_nthb_kink         = Counter()
#     cnt_nthb_gate         = Counter()
#     cnt_nthb_constr       = Counter()
#     cnt_cthb_kink         = Counter()
#     cnt_cthb_gate         = Counter()
#     cnt_cthb_constr       = Counter()
#     cnt_kink_pair         = Counter()   # nthb_kink + cthb_kink combined
#     cnt_gate_pair         = Counter()   # nthb_gate + cthb_gate combined
#     cnt_constr_pair       = Counter()   # nthb_constr + cthb_constr combined
#
#     # ── File 1: per-sequence table ────────────────────────────────────────────
#     header = (
#         "sequence_id\t"
#         "nthb_kink(G48,G49,G50)\t"
#         "nthb_gate(N13,S53,N54,Y68)\t"
#         "nthb_constriction(F20)\t"
#         "cthb_kink(G134,G135,G136)\t"
#         "cthb_gate(D99,R139,D140,Y155)\t"
#         "cthb_constriction(F106)"
#     )
#
#     lines1 = [header]
#     for hdr, seq in entries:
#         m = extract_motifs(seq)
#
#         # Accumulate counts
#         cnt_nthb_kink[m["nthb_kink"]]                         += 1
#         cnt_nthb_gate[m["nthb_gate"]]                         += 1
#         cnt_nthb_constr[m["nthb_constriction"]]                += 1
#         cnt_cthb_kink[m["cthb_kink"]]                         += 1
#         cnt_cthb_gate[m["cthb_gate"]]                         += 1
#         cnt_cthb_constr[m["cthb_constriction"]]                += 1
#         cnt_kink_pair[f"{m['nthb_kink']} | {m['cthb_kink']}"] += 1
#         cnt_gate_pair[f"{m['nthb_gate']} | {m['cthb_gate']}"] += 1
#         cnt_constr_pair[f"{m['nthb_constriction']}{m['cthb_constriction']}"] += 1
#
#         lines1.append(
#             f"{hdr}\t"
#             f"{m['nthb_kink']}\t"
#             f"{m['nthb_gate']}\t"
#             f"{m['nthb_constriction']}\t"
#             f"{m['cthb_kink']}\t"
#             f"{m['cthb_gate']}\t"
#             f"{m['cthb_constriction']}"
#         )
#
#     with open(OUT_FILE1, "w") as fh:
#         fh.write("\n".join(lines1) + "\n")
#     print(f"  File 1 written : {OUT_FILE1}")
#
#     # ── File 2: summary ───────────────────────────────────────────────────────
#     def fmt_counter(counter, total, title, pos_label):
#         """Format a counter block as readable text."""
#         lines = [f"\n{title}", f"  Positions : {pos_label}"]
#         for combo, count in counter.most_common():
#             pct = count / total * 100
#             lines.append(f"    {combo:<20}  {count:6d}  ({pct:.2f}%)")
#         return "\n".join(lines)
#
#     lines2 = []
#     lines2.append("=" * 70)
#     lines2.append("PROKARYOTIC TRIC — FUNCTIONAL MOTIF SUMMARY")
#     lines2.append(f"Total sequences analysed : {n}")
#     lines2.append("=" * 70)
#
#     # ── Part 1: individual motif frequencies ─────────────────────────────────
#     lines2.append("\n" + "─" * 70)
#     lines2.append("PART 1 — Individual Motif Frequencies")
#     lines2.append("─" * 70)
#
#     lines2.append(fmt_counter(cnt_nthb_constr, n,
#         "1) N-THB Constriction  [F20 / TM1]",
#         "F20"))
#     lines2.append(fmt_counter(cnt_cthb_constr, n,
#         "2) C-THB Constriction  [F106 / TM4]",
#         "F106"))
#     lines2.append(fmt_counter(cnt_nthb_kink, n,
#         "3) N-THB Kink  [G48-G49-G50 / TM2]",
#         "G48, G49, G50"))
#     lines2.append(fmt_counter(cnt_cthb_kink, n,
#         "4) C-THB Kink  [G134-G135-G136 / TM5]",
#         "G134, G135, G136"))
#     lines2.append(fmt_counter(cnt_nthb_gate, n,
#         "5) N-THB Gate  [N13-S53-N54-Y68 / TM1,TM2,TM3]",
#         "N13, S53, N54, Y68"))
#     lines2.append(fmt_counter(cnt_cthb_gate, n,
#         "6) C-THB Gate  [D99-R139-D140-Y155 / TM4,TM5,TM6]",
#         "D99, R139, D140, Y155"))
#
#     # ── Part 2: paired (N-THB + C-THB) combination frequencies ───────────────
#     lines2.append("\n" + "─" * 70)
#     lines2.append("PART 2 — N-THB + C-THB Paired Combination Frequencies")
#     lines2.append("─" * 70)
#
#     lines2.append(fmt_counter(cnt_constr_pair, n,
#         "1) Constriction pair  [F20 | F106]  (format: N-THB C-THB)",
#         "F20, F106"))
#     lines2.append(fmt_counter(cnt_kink_pair, n,
#         "2) Kink pair  [G48-50 | G134-136]  (format: N-THB kink | C-THB kink)",
#         "G48-G50, G134-G136"))
#     lines2.append(fmt_counter(cnt_gate_pair, n,
#         "3) Gate pair  [N13-S53-N54-Y68 | D99-R139-D140-Y155]  (format: N-THB gate | C-THB gate)",
#         "N13-S53-N54-Y68, D99-R139-D140-Y155"))
#
#     lines2.append("\n" + "=" * 70 + "\n")
#
#     with open(OUT_FILE2, "w") as fh:
#         fh.write("\n".join(lines2) + "\n")
#     print(f"  File 2 written : {OUT_FILE2}")
#     print("Done.")
#
#
# if __name__ == "__main__":
#     main()


######################### eukaryotes combination #####################

# """
# eukaryote_motif_analysis.py
# ────────────────────────────
# Reads a eukaryotic TM-only FASTA alignment and extracts
# functional residue combinations for each sequence, then
# summarises their frequencies.
#
# Reference structure : 6IYX
# Input  : /home/supremeleader/Desktop/article/seq_analysis/eukaryote_analysis/eukaryotic_overall_analysis/1_euk_alignment.txt
# Output : /home/supremeleader/Desktop/article/database/
#             eukaryote_motif_per_sequence.txt   (file 1)
#             eukaryote_motif_summary.txt         (file 2)
#
# Helix boundaries (PDB residue numbering, loops removed):
#     TM1   23-41   offset 23   cols   0-18
#     TM2   51-74   offset 32   cols  19-42
#     TM3   82-100  offset 39   cols  43-61
#     TM4  113-140  offset 51   cols  62-89
#     TM5  146-170  offset 56   cols  90-114
#     TM6  184-202  offset 69   cols 115-133
#     TM7  208-231  offset 74   cols 134-157
#
# Key residue → MSA column (0-based):
#     Y29  → col   6   (Gating,          TM1)
#     G64  → col  32   (TM2 kink N-THB,  TM2)
#     S65  → col  33   (TM2 kink N-THB + Gating, TM2)
#     Y66  → col  34   (TM2 kink N-THB,  TM2)
#     K129 → col  78   (Gating + Constriction, TM4)
#     H136 → col  85   (Constriction,    TM4)
#     G156 → col 100   (TM5 kink C-THB,  TM5)
#     S157 → col 101   (TM5 kink C-THB,  TM5)
#     G158 → col 102   (TM5 kink C-THB,  TM5)
#
# Functional combinations:
#     Gating            : Y29, S65, K129         → e.g. YSK
#     Constriction      : K129, H136             → e.g. KH
#     TM2 kink (N-THB)  : G64, S65, Y66          → e.g. GSY
#     TM5 kink (C-THB)  : G156, S157, G158       → e.g. GSG
#
# Paired combinations (N-THB + C-THB):
#     Kink pair         : TM2 kink | TM5 kink
#     Constriction pair : K129 | H136   (single residues from TM4)
# """
#
# import os
# from collections import Counter
#
# # ─── Paths ────────────────────────────────────────────────────────────────────
# INPUT_FILE = "/home/supremeleader/Desktop/article/phylogenetic_analysis/eukaryote_tree/cluster_sequence/6_cluster6_alignment.txt"
# OUTPUT_DIR = "/home/supremeleader/Desktop/article/database"
# OUT_FILE1  = os.path.join(OUTPUT_DIR, "euk_cluster6_per_sequence.txt")
# OUT_FILE2  = os.path.join(OUTPUT_DIR, "euk_cluster6_summary.txt")
#
# # ─── MSA column indices (0-based) ─────────────────────────────────────────────
# COL = {
#     "Y29":  6,
#     "G64":  32,
#     "S65":  33,
#     "Y66":  34,
#     "K129": 78,
#     "H136": 85,
#     "G156": 100,
#     "S157": 101,
#     "G158": 102,
# }
#
#
# # ─── Read FASTA ───────────────────────────────────────────────────────────────
# def read_fasta(path):
#     """Return list of (header, sequence) tuples."""
#     entries, current_hdr, current_seq = [], None, []
#     with open(path) as fh:
#         for line in fh:
#             line = line.rstrip()
#             if line.startswith(">"):
#                 if current_hdr is not None:
#                     entries.append((current_hdr, "".join(current_seq)))
#                 current_hdr, current_seq = line[1:].strip(), []
#             else:
#                 current_seq.append(line.upper())
#     if current_hdr is not None:
#         entries.append((current_hdr, "".join(current_seq)))
#     return entries
#
#
# # ─── Extract residue at a given column ───────────────────────────────────────
# def aa(seq, col):
#     return seq[col] if col < len(seq) else "-"
#
#
# # ─── Build motifs for one sequence ───────────────────────────────────────────
# def extract_motifs(seq):
#     return {
#         "gating":           aa(seq, COL["Y29"])  + aa(seq, COL["S65"])  + aa(seq, COL["K129"]),
#         "constriction":     aa(seq, COL["K129"]) + aa(seq, COL["H136"]),
#         "tm2_kink":         aa(seq, COL["G64"])  + aa(seq, COL["S65"])  + aa(seq, COL["Y66"]),
#         "tm5_kink":         aa(seq, COL["G156"]) + aa(seq, COL["S157"]) + aa(seq, COL["G158"]),
#     }
#
#
# # ─── Main ─────────────────────────────────────────────────────────────────────
# def main():
#     print(f"Reading : {INPUT_FILE}")
#     entries = read_fasta(INPUT_FILE)
#     n       = len(entries)
#     print(f"  {n} sequences found.")
#
#     os.makedirs(OUTPUT_DIR, exist_ok=True)
#
#     # ── Counters ──────────────────────────────────────────────────────────────
#     cnt_gating      = Counter()
#     cnt_constr      = Counter()
#     cnt_tm2_kink    = Counter()
#     cnt_tm5_kink    = Counter()
#     cnt_kink_pair   = Counter()   # tm2_kink | tm5_kink
#     cnt_constr_pair = Counter()   # K129 | H136 (already in constriction, but shown as pair)
#
#     # ── File 1: per-sequence ──────────────────────────────────────────────────
#     header = (
#         "sequence_id\t"
#         "gating(Y29,S65,K129)\t"
#         "constriction(K129,H136)\t"
#         "tm2_kink_nthb(G64,S65,Y66)\t"
#         "tm5_kink_cthb(G156,S157,G158)"
#     )
#
#     lines1 = [header]
#     for hdr, seq in entries:
#         m = extract_motifs(seq)
#
#         cnt_gating[m["gating"]]                                    += 1
#         cnt_constr[m["constriction"]]                              += 1
#         cnt_tm2_kink[m["tm2_kink"]]                                += 1
#         cnt_tm5_kink[m["tm5_kink"]]                                += 1
#         cnt_kink_pair[f"{m['tm2_kink']} | {m['tm5_kink']}"]       += 1
#         cnt_constr_pair[f"{m['constriction'][0]} | {m['constriction'][1]}"] += 1
#
#         lines1.append(
#             f"{hdr}\t"
#             f"{m['gating']}\t"
#             f"{m['constriction']}\t"
#             f"{m['tm2_kink']}\t"
#             f"{m['tm5_kink']}"
#         )
#
#     with open(OUT_FILE1, "w") as fh:
#         fh.write("\n".join(lines1) + "\n")
#     print(f"  File 1 written : {OUT_FILE1}")
#
#     # ── File 2: summary ───────────────────────────────────────────────────────
#     def fmt_counter(counter, total, title, pos_label):
#         lines = [f"\n{title}", f"  Positions : {pos_label}"]
#         for combo, count in counter.most_common():
#             pct = count / total * 100
#             lines.append(f"    {combo:<25}  {count:6d}  ({pct:.2f}%)")
#         return "\n".join(lines)
#
#     lines2 = []
#     lines2.append("=" * 70)
#     lines2.append("EUKARYOTIC TRIC — FUNCTIONAL MOTIF SUMMARY")
#     lines2.append(f"Total sequences analysed : {n}")
#     lines2.append("=" * 70)
#
#     # Part 1 — individual motif frequencies
#     lines2.append("\n" + "─" * 70)
#     lines2.append("PART 1 — Individual Motif Frequencies")
#     lines2.append("─" * 70)
#
#     lines2.append(fmt_counter(cnt_gating, n,
#         "1) Gating  [Y29-S65-K129 / TM1, TM2, TM4]",
#         "Y29, S65, K129"))
#     lines2.append(fmt_counter(cnt_constr, n,
#         "2) Constriction  [K129-H136 / TM4]",
#         "K129, H136"))
#     lines2.append(fmt_counter(cnt_tm2_kink, n,
#         "3) TM2 Kink N-THB  [G64-S65-Y66 / TM2]",
#         "G64, S65, Y66"))
#     lines2.append(fmt_counter(cnt_tm5_kink, n,
#         "4) TM5 Kink C-THB  [G156-S157-G158 / TM5]",
#         "G156, S157, G158"))
#
#     # Part 2 — paired combination frequencies
#     lines2.append("\n" + "─" * 70)
#     lines2.append("PART 2 — Paired Combination Frequencies")
#     lines2.append("─" * 70)
#
#     lines2.append(fmt_counter(cnt_kink_pair, n,
#         "1) Kink pair  [TM2 kink | TM5 kink]  (format: N-THB | C-THB)",
#         "G64-S65-Y66, G156-S157-G158"))
#     lines2.append(fmt_counter(cnt_constr_pair, n,
#         "2) Constriction pair  [K129 | H136]  (format: K129 residue | H136 residue)",
#         "K129, H136"))
#
#     lines2.append("\n" + "=" * 70 + "\n")
#
#     with open(OUT_FILE2, "w") as fh:
#         fh.write("\n".join(lines2) + "\n")
#     print(f"  File 2 written : {OUT_FILE2}")
#     print("Done.")
#
#
# if __name__ == "__main__":
#     main()
#

################### seperate sequences #####################3
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
# FULL_SEQ  = "/home/supremeleader/Desktop/article/full_sequence/eukaryote_full_length.fa"
# TRICA_IDS = "/home/supremeleader/Desktop/Paper1_data/paper1/eukaryotes/TRIC-A_id.txt"
# TRICB_IDS = "/home/supremeleader/Desktop/Paper1_data/paper1/eukaryotes/TRIC-B_id.txt"
# OUT_DIR   = "/home/supremeleader/Desktop/article/full_sequence/"
#
# # ── Read ID files ─────────────────────────────────────────────────────────────
# def read_ids(path):
#     ids = set()
#     with open(path) as fh:
#         for line in fh:
#             line = line.strip()
#             if line:
#                 ids.add(line)
#     return ids
#
# trica_ids = read_ids(TRICA_IDS)
# tricb_ids = read_ids(TRICB_IDS)
# print(f"TRIC-A IDs loaded : {len(trica_ids)}")
# print(f"TRIC-B IDs loaded : {len(tricb_ids)}")
#
# # ── Read full FASTA ───────────────────────────────────────────────────────────
# def read_fasta(path):
#     records = {}
#     current_id, current_seq = None, []
#     with open(path) as fh:
#         for line in fh:
#             line = line.rstrip()
#             if line.startswith(">"):
#                 if current_id:
#                     records[current_id] = "".join(current_seq)
#                 # take just the first word after > as ID
#                 current_id  = line[1:].split()[0]
#                 current_seq = []
#             else:
#                 current_seq.append(line)
#     if current_id:
#         records[current_id] = "".join(current_seq)
#     return records
#
# sequences = read_fasta(FULL_SEQ)
# print(f"Total sequences    : {len(sequences)}")
#
# # ── Match and write ───────────────────────────────────────────────────────────
# def write_fasta(seq_dict, filepath):
#     with open(filepath, "w") as fh:
#         for sid, seq in seq_dict.items():
#             fh.write(f">{sid}\n{seq}\n")
#
# # match: check if any ID from the id-file is a substring of the fasta header id
# # (handles cases like sp|Q9H6F2|... vs plain Q9H6F2)
# def match_ids(sequences, id_set):
#     matched = {}
#     for seq_id, seq in sequences.items():
#         for target in id_set:
#             if target in seq_id or seq_id in target:
#                 matched[seq_id] = seq
#                 break
#     return matched
#
# trica_seqs = match_ids(sequences, trica_ids)
# tricb_seqs = match_ids(sequences, tricb_ids)
#
# out_a = OUT_DIR + "TRIC-A_sequences.fa"
# out_b = OUT_DIR + "TRIC-B_sequences.fa"
#
# write_fasta(trica_seqs, out_a)
# write_fasta(tricb_seqs, out_b)
#
# print(f"\nTRIC-A sequences written : {len(trica_seqs)}  →  {out_a}")
# print(f"TRIC-B sequences written : {len(tricb_seqs)}  →  {out_b}")
#
# # ── Report unmatched ──────────────────────────────────────────────────────────
# unmatched_a = trica_ids - set(trica_seqs.keys())
# unmatched_b = tricb_ids - set(tricb_seqs.keys())
# if unmatched_a:
#     print(f"\nWARNING: {len(unmatched_a)} TRIC-A IDs not found in FASTA:")
#     for u in sorted(unmatched_a): print(f"  {u}")
# if unmatched_b:
#     print(f"\nWARNING: {len(unmatched_b)} TRIC-B IDs not found in FASTA:")
#     for u in sorted(unmatched_b): print(f"  {u}")



############################################################################################
########################## Pairwise sequence identity ######################################
############################################################################################
# """
# tric_exact_pairwise.py
# ──────────────────────────────────────────────────────────────────────────────
# Calculates EXACT average pairwise % identity and % similarity for all
# within-group and between-group comparisons.
#
# 8 x 8 MATRIX:
#   Pro | Euk | TRIC-A | TRIC-B | Gram-N | Gram-P | Archaea | Bacteria
#
# Strategy (fast + 100% exact — no sampling):
#   For each comparison (GroupA vs GroupB):
#     ─ Split the SMALLER group into individual single-sequence .fa files
#     ─ Run: needle  seq_i.fa  full_GroupB.fa  →  one output file per query
#     ─ needle aligns seq_i against EVERY sequence in GroupB in one call
#     ─ Parse ALL alignment blocks from the output (reads "# Identity:" lines)
#     ─ Skip 100.0% self-alignments (same as your original code)
#     ─ Average all collected identity% and similarity% values
#
#   Example: Archaea (44 seq) vs Prokaryote (3238 seq)
#     → split Archaea into 44 individual files
#     → run 44 needle calls, each producing 3238 alignment blocks
#     → parse 44 × 3238 = 142,472 pairs  (completes in minutes)
#
#   For within-group (e.g. Bacteria vs Bacteria):
#     → split ALL sequences into individual files
#     → run each seq_i against the full group file
#     → self-alignments (100.0%) are automatically skipped
#
# Output structure
# ──────────────────
#   pair_wise_alignment/
#       Archaea_vs_Pro/
#           individual_seqs/     ← one .fa file per split sequence
#           results/             ← one .needle file per query sequence
#           all_pairs.txt        ← seq_A  seq_B  identity%  similarity%
#           summary.txt          ← avg identity, avg similarity, n pairs, time
#       Pro_vs_Pro/
#           ...
#       summary_table.txt        ← full 8×8 matrix (human-readable)
#       summary_table.csv        ← same, importable to Excel
#
# Usage
# ──────
#   python3 tric_exact_pairwise.py
#
# Requirements
# ─────────────
#   EMBOSS needle in PATH (already confirmed on your system)
# """
#
# import os
# import sys
# import csv
# import subprocess
# import time
# from datetime import datetime
#
# # ─────────────────────────────────────────────────────────────────────────────
# # CONFIGURATION
# # ─────────────────────────────────────────────────────────────────────────────
# SEQ_DIR = "/home/supremeleader/Desktop/article/sequences"
# OUT_DIR = "/media/supremeleader/Pantera/pairwise_alignment"
#
# FILES = {
#     "Pro"     : os.path.join(SEQ_DIR, "prokaryote_full_length.fa"),
#     "Euk"     : os.path.join(SEQ_DIR, "eukaryote_full_length.fa"),
#     "TRICA"   : os.path.join(SEQ_DIR, "TRIC-A_sequences.fa"),
#     "TRICB"   : os.path.join(SEQ_DIR, "TRIC-B_sequences.fa"),
#     "GramN"   : os.path.join(SEQ_DIR, "gram-ve_full_length.fa"),
#     "GramP"   : os.path.join(SEQ_DIR, "gram+ve_full_length.fa"),
#     "Archaea" : os.path.join(SEQ_DIR, "archaea_full_length.fa"),
#     "Bacteria": os.path.join(SEQ_DIR, "bacteria_full_length.fa"),
# }
#
# DISPLAY = {
#     "Pro"     : "Prokaryote",
#     "Euk"     : "Eukaryote",
#     "TRICA"   : "TRIC-A",
#     "TRICB"   : "TRIC-B",
#     "GramN"   : "Gram-Neg",
#     "GramP"   : "Gram-Pos",
#     "Archaea" : "Archaea",
#     "Bacteria": "Bacteria",
# }
#
# GROUPS = ["Pro", "Euk", "TRICA", "TRICB", "GramN", "GramP", "Archaea", "Bacteria"]
#
#
#
# # ─────────────────────────────────────────────────────────────────────────────
# # COUNT SEQUENCES IN A FASTA FILE (fast, no full load)
# # ─────────────────────────────────────────────────────────────────────────────
# def count_fasta(path):
#     n = 0
#     with open(path) as f:
#         for line in f:
#             if line.startswith(">"):
#                 n += 1
#     return n
#
# # ─────────────────────────────────────────────────────────────────────────────
# # READ FASTA
# # ─────────────────────────────────────────────────────────────────────────────
# def read_fasta(path):
#     """Returns list of (seq_id, sequence) tuples."""
#     records = []
#     cur_id, cur_seq = None, []
#     with open(path) as fh:
#         for line in fh:
#             line = line.rstrip()
#             if line.startswith(">"):
#                 if cur_id is not None:
#                     records.append((cur_id, "".join(cur_seq)))
#                 cur_id  = line[1:].split()[0]
#                 cur_seq = []
#             else:
#                 cur_seq.append(line.upper())
#     if cur_id is not None:
#         records.append((cur_id, "".join(cur_seq)))
#     return records
#
# # ─────────────────────────────────────────────────────────────────────────────
# # SPLIT FASTA INTO INDIVIDUAL FILES
# # ─────────────────────────────────────────────────────────────────────────────
# def split_fasta(records, out_dir):
#     """
#     Write each sequence as its own .fa file into out_dir.
#     Returns list of (seq_id, filepath) tuples.
#     """
#     os.makedirs(out_dir, exist_ok=True)
#     paths = []
#     for seq_id, seq in records:
#         safe_id = seq_id.replace("|", "_").replace("/", "_").replace(" ", "_")
#         fname   = os.path.join(out_dir, f"{safe_id}.fa")
#         with open(fname, "w") as f:
#             f.write(f">{seq_id}\n{seq}\n")
#         paths.append((seq_id, fname))
#     return paths
#
# # ─────────────────────────────────────────────────────────────────────────────
# # PARSE NEEDLE OUTPUT — same logic as your calculate_average_sd
# # Handles multiple alignment blocks in one file (one-vs-all output)
# # ─────────────────────────────────────────────────────────────────────────────
# def parse_needle_file(needle_file):
#     """
#     Parses a needle output file containing one or many alignment blocks.
#     Returns list of (seq_b_id, identity_pct, similarity_pct).
#     Skips 100.0% self-alignments exactly like your original code.
#     """
#     results  = []
#     cur_seqb = None
#     cur_id   = None
#     cur_sim  = None
#
#     with open(needle_file) as fh:
#         for line in fh:
#
#             # Seq B name at start of block
#             if line.startswith("# 2:"):
#                 cur_seqb = line.split("# 2:")[-1].strip()
#                 cur_id   = None
#                 cur_sim  = None
#
#             # Identity line — skip self-alignments (100.0%)
#             elif line.startswith("# Identity:") and "100.0%" not in line:
#                 cur_id = float(line.split("(")[-1].split("%")[0])
#
#             # Similarity line — skip self-alignments
#             elif line.startswith("# Similarity:") and "100.0%" not in line:
#                 cur_sim = float(line.split("(")[-1].split("%")[0])
#
#             # End of block separator — collect result if both values found
#             elif line.startswith("#=====") and cur_seqb is not None:
#                 if cur_id is not None and cur_sim is not None:
#                     results.append((cur_seqb, cur_id, cur_sim))
#                 cur_id  = None
#                 cur_sim = None
#
#     return results
#
# # ─────────────────────────────────────────────────────────────────────────────
# # RUN NEEDLE: single_seq.fa  vs  group_file.fa
# # ─────────────────────────────────────────────────────────────────────────────
# def run_needle(seq_a_path, group_b_path, out_path):
#     cmd = (
#         f"needle "
#         f"-asequence {seq_a_path} "
#         f"-bsequence {group_b_path} "
#         f"-outfile {out_path} "
#         f"-gapopen 10 "
#         f"-gapextend 0.5 "
#         f"-auto"
#     )
#     ret = subprocess.run(cmd, shell=True, capture_output=True)
#     return ret.returncode == 0
#
# # ─────────────────────────────────────────────────────────────────────────────
# # CORE COMPARISON FUNCTION
# # ─────────────────────────────────────────────────────────────────────────────
# def compare_groups(label_split, label_search, file_split, file_search, out_folder):
#     """
#     Splits file_split into individual sequences.
#     Runs needle for each one against entire file_search.
#     Parses results and writes output.
#     Returns (avg_identity, avg_similarity, n_pairs).
#     """
#     ind_dir     = os.path.join(out_folder, "individual_seqs")
#     results_dir = os.path.join(out_folder, "results")
#     os.makedirs(out_folder,  exist_ok=True)
#     os.makedirs(ind_dir,     exist_ok=True)
#     os.makedirs(results_dir, exist_ok=True)
#
#     all_pairs_file = os.path.join(out_folder, "all_pairs.txt")
#     summary_file   = os.path.join(out_folder, "summary.txt")
#
#     print(f"  Loading {label_split} sequences to split...")
#     records_split = read_fasta(file_split)
#     print(f"  {label_split}: {len(records_split)} sequences → splitting into individual files")
#     ind_files = split_fasta(records_split, ind_dir)
#
#     n_search   = count_fasta(file_search)
#     total_est  = len(ind_files) * n_search
#     print(f"  {label_search}: {n_search} sequences (search target)")
#     print(f"  Estimated pairs: ~{total_est:,}")
#
#     sum_id  = 0.0
#     sum_sim = 0.0
#     n_pairs = 0
#     n_fail  = 0
#     t_start = time.time()
#
#     with open(all_pairs_file, "w") as pf:
#         pf.write(f"# {label_split} vs {label_search}\n")
#         pf.write(f"# Date: {datetime.now()}\n")
#         pf.write(f"{'Seq_A':<35} {'Seq_B':<35} {'Identity%':>10} {'Similarity%':>12}\n")
#         pf.write("-" * 96 + "\n")
#
#         for i, (seq_id_a, ind_fa) in enumerate(ind_files):
#             safe_id    = os.path.basename(ind_fa).replace(".fa", "")
#             needle_out = os.path.join(results_dir, f"{safe_id}.needle")
#
#             # Run needle: seq_i vs entire group file
#             ok = run_needle(ind_fa, file_search, needle_out)
#             if not ok:
#                 print(f"\n  WARNING: needle failed for {seq_id_a}")
#                 n_fail += 1
#                 continue
#
#             # Parse all alignment blocks from this output file
#             block_results = parse_needle_file(needle_out)
#
#             for seq_id_b, pct_id, pct_sim in block_results:
#                 sum_id  += pct_id
#                 sum_sim += pct_sim
#                 n_pairs += 1
#                 pf.write(f"{seq_id_a:<35} {seq_id_b:<35} {pct_id:>10.2f} {pct_sim:>12.2f}\n")
#
#             # Progress line
#             elapsed  = time.time() - t_start
#             rate     = (i + 1) / elapsed if elapsed > 0 else 0
#             eta      = (len(ind_files) - i - 1) / rate / 60 if rate > 0 else 0
#             cur_avg  = sum_id / n_pairs if n_pairs > 0 else 0
#             cur_sim  = sum_sim / n_pairs if n_pairs > 0 else 0
#             print(f"  [{i+1:>5}/{len(ind_files)}]  "
#                   f"pairs={n_pairs:>8,}  "
#                   f"avg_id={cur_avg:.2f}%  avg_sim={cur_sim:.2f}%  "
#                   f"ETA={eta:.1f} min     ",
#                   end="\r")
#
#     print()  # newline after \r progress
#
#     avg_id  = sum_id  / n_pairs if n_pairs > 0 else 0.0
#     avg_sim = sum_sim / n_pairs if n_pairs > 0 else 0.0
#     elapsed = time.time() - t_start
#
#     with open(summary_file, "w") as sf:
#         sf.write(f"Comparison        : {label_split} vs {label_search}\n")
#         sf.write(f"Date              : {datetime.now()}\n")
#         sf.write(f"Sequences split   : {len(records_split)} ({label_split})\n")
#         sf.write(f"Sequences searched: {n_search} ({label_search})\n")
#         sf.write(f"Total pairs parsed: {n_pairs:,}\n")
#         sf.write(f"Failed queries    : {n_fail:,}\n")
#         sf.write(f"Avg Identity      : {avg_id:.2f} %\n")
#         sf.write(f"Avg Similarity    : {avg_sim:.2f} %\n")
#         sf.write(f"Time taken        : {elapsed/60:.1f} min\n")
#
#     print(f"  ✓  {label_split} vs {label_search}:  "
#           f"identity={avg_id:.2f}%  similarity={avg_sim:.2f}%  "
#           f"({n_pairs:,} pairs,  {elapsed/60:.1f} min)")
#
#     return avg_id, avg_sim, n_pairs
#
# # ─────────────────────────────────────────────────────────────────────────────
# # WRITE SUMMARY MATRIX
# # ─────────────────────────────────────────────────────────────────────────────
# def write_matrix(results, out_dir):
#     col_w = 24
#
#     txt_path = os.path.join(out_dir, "summary_table.txt")
#     with open(txt_path, "w") as f:
#         f.write("TRIC PAIRWISE SEQUENCE IDENTITY [SIMILARITY] MATRIX\n")
#         f.write("Format: avg_identity% [avg_similarity%]\n")
#         f.write(f"Date  : {datetime.now()}\n\n")
#         f.write(f"{'':22}")
#         for g in GROUPS:
#             f.write(f"{DISPLAY[g]:>{col_w}}")
#         f.write("\n")
#         f.write("─" * (22 + col_w * len(GROUPS)) + "\n")
#         for row in GROUPS:
#             f.write(f"{DISPLAY[row]:<22}")
#             for col in GROUPS:
#                 val = results.get((row, col)) or results.get((col, row))
#                 cell = f"{val[0]:.2f} [{val[1]:.2f}]" if val else "—"
#                 f.write(f"{cell:>{col_w}}")
#             f.write("\n")
#
#     csv_path = os.path.join(out_dir, "summary_table.csv")
#     with open(csv_path, "w", newline="") as cf:
#         writer = csv.writer(cf)
#         writer.writerow([""] + [DISPLAY[g] for g in GROUPS])
#         for row in GROUPS:
#             cells = [DISPLAY[row]]
#             for col in GROUPS:
#                 val = results.get((row, col)) or results.get((col, row))
#                 cells.append(f"{val[0]:.2f} [{val[1]:.2f}]" if val else "")
#             writer.writerow(cells)
#
#     print("\n" + "=" * 80)
#     print("CURRENT 8×8 MATRIX — Identity% [Similarity%]")
#     print("=" * 80)
#     with open(txt_path) as f:
#         print(f.read())
#
# # ─────────────────────────────────────────────────────────────────────────────
# # MAIN
# # ─────────────────────────────────────────────────────────────────────────────
# def main():
#     print("=" * 70)
#     print("TRIC Exact Pairwise Identity / Similarity — 8×8 Matrix")
#     print(f"Started: {datetime.now()}")
#     print("=" * 70)
#
#     ret = subprocess.run(["which", "needle"], capture_output=True, text=True)
#     if ret.returncode != 0:
#         print("ERROR: needle not found. Install: sudo apt-get install emboss")
#         sys.exit(1)
#     print(f"✓ needle: {ret.stdout.strip()}\n")
#
#     os.makedirs(OUT_DIR, exist_ok=True)
#
#     print("Input files:")
#     seq_counts = {}
#     for label, fpath in FILES.items():
#         if os.path.exists(fpath):
#             n = count_fasta(fpath)
#             seq_counts[label] = n
#             print(f"  {label:10}: {n:6,} sequences  {fpath}")
#         else:
#             seq_counts[label] = 0
#             print(f"  {label:10}: MISSING  {fpath}")
#     print()
#
#     # Upper triangle + diagonal = 36 comparisons for 8 groups
#     comparisons = [
#         (GROUPS[i], GROUPS[j])
#         for i in range(len(GROUPS))
#         for j in range(i, len(GROUPS))
#     ]
#     print(f"Total comparisons: {len(comparisons)}\n")
#
#     results     = {}
#     grand_start = time.time()
#
#     for idx, (ga, gb) in enumerate(comparisons, 1):
#         print(f"\n{'─'*60}")
#         print(f"[{idx}/{len(comparisons)}]  {ga}  vs  {gb}")
#         print(f"{'─'*60}")
#
#         if seq_counts.get(ga, 0) == 0 or seq_counts.get(gb, 0) == 0:
#             print("  SKIPPED — missing file or no sequences")
#             continue
#
#         folder_name = f"{ga}_vs_{gb}"
#         out_folder  = os.path.join(OUT_DIR, folder_name)
#
#         # For within-group: split ga, search against ga file
#         # For between-group: split SMALLER group, search against LARGER
#         if ga == gb:
#             file_split  = FILES[ga]
#             file_search = FILES[gb]
#             lbl_split   = ga
#             lbl_search  = gb
#         else:
#             na = seq_counts[ga]
#             nb = seq_counts[gb]
#             if na <= nb:
#                 file_split, file_search = FILES[ga], FILES[gb]
#                 lbl_split,  lbl_search  = ga, gb
#             else:
#                 file_split, file_search = FILES[gb], FILES[ga]
#                 lbl_split,  lbl_search  = gb, ga
#             print(f"  Splitting {lbl_split} ({min(na,nb):,}) "
#                   f"against {lbl_search} ({max(na,nb):,})")
#
#         avg_id, avg_sim, n_pairs = compare_groups(
#             lbl_split, lbl_search,
#             file_split, file_search,
#             out_folder
#         )
#         results[(ga, gb)] = (avg_id, avg_sim)
#
#         # Save checkpoint after every comparison
#         write_matrix(results, OUT_DIR)
#
#     write_matrix(results, OUT_DIR)
#     total_time = (time.time() - grand_start) / 60
#     print(f"\n{'='*70}")
#     print(f"ALL DONE in {total_time:.1f} minutes")
#     print(f"Results: {OUT_DIR}")
#     print(f"{'='*70}")
#
#
# if __name__ == "__main__":
#     main()

#################################################################################################################################
########################################### download structure from alphafold database ########################
#################################################################################################################################
# import os
# import requests
#
# INPUT = "/home/supremeleader/Desktop/article/uniprot_ids.txt"
# OUTPUT_DIR = "/home/supremeleader/alphafold_structures"
# NOT_FOUND = os.path.join(OUTPUT_DIR, "not_found.txt")
#
# os.makedirs(OUTPUT_DIR, exist_ok=True)
#
# with open(INPUT, 'r') as f:
#     ids = [line.strip() for line in f if line.strip()]
#
# not_found = []
#
# for i, uid in enumerate(ids, 1):
#     print(f"[{i}/{len(ids)}] Processing {uid}...")
#
#     # Step 1 — query the API to get the actual PDB download URL
#     api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uid}"
#     try:
#         api_resp = requests.get(api_url, timeout=30)
#         if api_resp.status_code != 200:
#             print(f"  ✗ Not found in AlphaFold: {uid} (HTTP {api_resp.status_code})")
#             not_found.append(uid)
#             continue
#
#         data = api_resp.json()
#         pdb_url = data[0]['pdbUrl']
#
#         # Step 2 — download the PDB file
#         pdb_resp = requests.get(pdb_url, timeout=30)
#         if pdb_resp.status_code == 200:
#             out_file = os.path.join(OUTPUT_DIR, f"{uid}.pdb")
#             with open(out_file, 'wb') as f:
#                 f.write(pdb_resp.content)
#             print(f"  ✓ Saved {uid}.pdb")
#         else:
#             print(f"  ✗ Failed to download PDB for {uid} (HTTP {pdb_resp.status_code})")
#             not_found.append(uid)
#
#     except Exception as e:
#         print(f"  ✗ Error for {uid}: {e}")
#         not_found.append(uid)
#
# with open(NOT_FOUND, 'w') as f:
#     f.write('\n'.join(not_found))
#
# print(f"\nDone.")
# print(f"Downloaded : {len(ids) - len(not_found)} structures")
# print(f"Not found  : {len(not_found)} — saved to {NOT_FOUND}")
#################################################################################################################################
########################################### download structure from alphafold database ########################
#################################################################################################################################

#################################################################################################################################
########################################### make the list of sequences which were not downloaded and save it in a file###########
#################################################################################################################################
#!/usr/bin/env python3
"""
find_missing_sequences.py
--------------------------
Finds sequences that could NOT be downloaded from the AlphaFold database
by comparing:
  1. Uniprot IDs already downloaded (PDB files in alphafold_structures/)
  2. Full list of Uniprot IDs (uniprot_ids.txt)
  3. Full sequence file (entire_sequence.txt) in FASTA format

Saves the missing sequences to: /home/supremeleader/modelling_sequence.txt

Usage:
    python3 find_missing_sequences.py
"""

import os
import re

# ─────────────────────────────────────────────
#  PATHS — edit if needed
# ─────────────────────────────────────────────

ALPHAFOLD_STRUCTURES_DIR = "/home/supremeleader/alphafold_structures"
UNIPROT_IDS_FILE         = "/home/supremeleader/Desktop/article/uniprot_ids.txt"
SEQUENCES_FILE           = "/home/supremeleader/Desktop/article/sequences/entire_sequence.txt"
OUTPUT_FILE              = "/home/supremeleader/modelling_sequence.txt"

# ─────────────────────────────────────────────
#  STEP 1 — Get IDs already downloaded (from PDB filenames)
# ─────────────────────────────────────────────

print("Step 1: Scanning downloaded AlphaFold structures...")

if not os.path.isdir(ALPHAFOLD_STRUCTURES_DIR):
    print(f"  ERROR: Directory not found: {ALPHAFOLD_STRUCTURES_DIR}")
    exit(1)

downloaded_ids = set()
for filename in os.listdir(ALPHAFOLD_STRUCTURES_DIR):
    if filename.endswith(".pdb"):
        # Filenames are simply "UNIPROTID.pdb" e.g. "A0A011UMQ6.pdb"
        uniprot_id = filename.replace(".pdb", "").strip().upper()
        downloaded_ids.add(uniprot_id)

print(f"  Found {len(downloaded_ids)} downloaded structures.")

# ─────────────────────────────────────────────
#  STEP 2 — Get full list of Uniprot IDs
# ─────────────────────────────────────────────

print("Step 2: Reading full Uniprot ID list...")

if not os.path.isfile(UNIPROT_IDS_FILE):
    print(f"  ERROR: File not found: {UNIPROT_IDS_FILE}")
    exit(1)

all_ids = set()
with open(UNIPROT_IDS_FILE, "r") as f:
    for line in f:
        uid = line.strip().upper()
        if uid:
            all_ids.add(uid)

print(f"  Found {len(all_ids)} total Uniprot IDs.")

# ─────────────────────────────────────────────
#  STEP 3 — Find missing IDs
# ─────────────────────────────────────────────

print("Step 3: Finding missing IDs (not downloaded from AlphaFold)...")

missing_ids = all_ids - downloaded_ids

print(f"  Total IDs         : {len(all_ids)}")
print(f"  Downloaded        : {len(downloaded_ids)}")
print(f"  Missing (to model): {len(missing_ids)}")

if not missing_ids:
    print("\n  All sequences were already downloaded! Nothing to model.")
    exit(0)

# ─────────────────────────────────────────────
#  STEP 4 — Parse the full sequence FASTA file
# ─────────────────────────────────────────────

print("Step 4: Parsing sequence file...")

if not os.path.isfile(SEQUENCES_FILE):
    print(f"  ERROR: File not found: {SEQUENCES_FILE}")
    exit(1)

# Parse FASTA — handles both standard FASTA and plain sequence lists
sequences = {}  # { uniprot_id : (header, sequence) }

with open(SEQUENCES_FILE, "r") as f:
    content = f.read().strip()

if content.startswith(">"):
    # Standard FASTA format
    current_header = None
    current_seq_lines = []

    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            # Save previous entry
            if current_header is not None:
                seq = "".join(current_seq_lines).upper()
                # Try to extract Uniprot ID from header
                # Handles: >P12345, >sp|P12345|NAME, >tr|P12345|NAME, >P12345 description
                uid_match = re.search(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', current_header)
                if uid_match:
                    uid = uid_match.group(0).upper()
                else:
                    # fallback: use first word after >
                    uid = current_header.split()[0].lstrip(">").split("|")[1] \
                        if "|" in current_header else current_header.split()[0].lstrip(">")
                    uid = uid.upper()
                sequences[uid] = (current_header, seq)

            current_header = line
            current_seq_lines = []
        else:
            current_seq_lines.append(line)

    # Don't forget the last entry
    if current_header is not None and current_seq_lines:
        seq = "".join(current_seq_lines).upper()
        uid_match = re.search(r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', current_header)
        if uid_match:
            uid = uid_match.group(0).upper()
        else:
            uid = current_header.split()[0].lstrip(">").split("|")[1] \
                if "|" in current_header else current_header.split()[0].lstrip(">")
            uid = uid.upper()
        sequences[uid] = (current_header, seq)

else:
    # Plain list — match line number to ID order in uniprot_ids.txt
    print("  Detected plain sequence list. Matching by order to uniprot_ids.txt...")
    all_ids_ordered = []
    with open(UNIPROT_IDS_FILE, "r") as f:
        for line in f:
            uid = line.strip().upper()
            if uid:
                all_ids_ordered.append(uid)

    for i, line in enumerate(content.splitlines()):
        seq = line.strip().upper()
        if seq and i < len(all_ids_ordered):
            uid = all_ids_ordered[i]
            sequences[uid] = (f">{uid}", seq)

print(f"  Parsed {len(sequences)} sequences from file.")

# ─────────────────────────────────────────────
#  STEP 5 — Write missing sequences to output
# ─────────────────────────────────────────────

print("Step 5: Writing missing sequences to output file...")

found_count = 0
not_found_ids = []

with open(OUTPUT_FILE, "w") as out:
    for uid in sorted(missing_ids):
        if uid in sequences:
            header, seq = sequences[uid]
            # Ensure header starts with >
            if not header.startswith(">"):
                header = f">{header}"
            out.write(f"{header}\n{seq}\n")
            found_count += 1
        else:
            not_found_ids.append(uid)

# ─────────────────────────────────────────────
#  SUMMARY
# ─────────────────────────────────────────────

print("\n" + "=" * 50)
print("  SUMMARY")
print("=" * 50)
print(f"  Total IDs in list        : {len(all_ids)}")
print(f"  Already downloaded       : {len(downloaded_ids)}")
print(f"  Missing (needed)         : {len(missing_ids)}")
print(f"  Sequences found & saved  : {found_count}")
print(f"  IDs with no sequence     : {len(not_found_ids)}")
print(f"\n  Output saved to: {OUTPUT_FILE}")

if not_found_ids:
    print(f"\n  WARNING: These {len(not_found_ids)} IDs had no matching sequence in the file:")
    for uid in not_found_ids:
        print(f"    - {uid}")

print("\nDone! You can now run:")
print(f"  colabfold_batch {OUTPUT_FILE} /home/supremeleader/alpha_models/models/")

















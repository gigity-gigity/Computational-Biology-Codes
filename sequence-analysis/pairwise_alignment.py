import sys
import os
import subprocess
import pandas as pd
import re
from collections import Counter
file_directory = "/home/supremeleader/first_paper_data"
bacteria_csv_file = "/home/supremeleader/first_paper_data/updated_bacteria.csv"
eukaryotes_csv_file = "/home/supremeleader/first_paper_data/database_eukaryote_data.csv"
# gram_positive=["terrabacteria"]
# gram_negative=["Proteobacteria", "FCB group","unclassified Bacteria","Spirochaetes","Fusobacteria","PVC group","Deferribacteres","Acidobacteria", "Aquificae", "Calditrichaeota"]
# archaea=["TACK group","Euryarchaeota"]
# archaea_id=[]
# gram_positive_id=[]
# gram_negative_id=[]
# eukaryota_id=[]

# with open(bacteria_csv_file) as handle2:
#     lines=handle2.readlines()
#     for line in lines:
#         if any(re.search(gram_positive,line,re.IGNORECASE)for gram_positive in gram_positive):
#             ids=line.split(",")
#             gram_positive_id.append(ids[1])
# print(len(gram_positive_id))
# with open(bacteria_csv_file) as handle:
#     lines=handle.readlines()
#     for line in lines:
#         if any(re.search(gram_negative,line,re.IGNORECASE)for gram_negative in gram_negative):
#             ids=line.split(",")
#             gram_negative_id.append(ids[1])
# print(len(gram_negative_id))
# with open(bacteria_csv_file) as handle:
#     lines=handle.readlines()
#     for line in lines:
#         if any(re.search(archaea,line,re.IGNORECASE)for archaea in archaea):
#             ids=line.split(",")
#             archaea_id.append(ids[1])
# print(len(archaea_id))
# with open(eukaryotes_csv_file)as handle:
#     for i in handle:
#         x=i.split(",")
#         eukaryota_id.append(x[1])
# print(len(eukaryota_id))


# df=pd.read_csv(eukaryotes_csv_file)
# sequence={}
# for ids in eukaryota_id:
#     fasta_sequence = df.loc[df['Entry'].str.contains(ids, case=False), "sequence"].dropna()
#     if not fasta_sequence.empty:
#         fasta_sequence=fasta_sequence.iloc[0]
#         sequence[ids] = fasta_sequence
#     output_filename = f"eukaryote_full_length.fa"
#     with open(output_filename, 'w') as f:
#         for key, value in sequence.items():
#             output_line = f">{key}\n{value}\n"
#             f.write(output_line)
# print(len(sequence))
# print(sequence)

# def split_fasta(input_file):
#     with open(input_file, 'r') as infile:
#         current_sequence = ''
#         for line in infile:
#             if line.startswith('>'):
#                 if current_sequence:
#                     write_sequence_to_file(sequence_name, current_sequence)
#                 sequence_name = line.strip().lstrip('>')
#                 current_sequence = ''
#             else:
#                 current_sequence += line.strip()
#         if current_sequence:
#             write_sequence_to_file(sequence_name, current_sequence)
#
# def write_sequence_to_file(sequence_name, sequence):
#     output_file = f'{sequence_name}.fa'
#     with open(output_file, 'w') as outfile:
#         outfile.write(f'>{sequence_name}\n{sequence}')
#
# def make_individual_fasta_file_and_delete_original_file(fasta_sequence_location):
#     files=os.listdir(fasta_sequence_location)
#     current_directory = os.getcwd()
#     for file in files:
#         cluster=file.split(".")[0]
#         source_directory = os.path.join(current_directory,cluster)
#         sequenceA = os.path.join(source_directory, 'sequenceA')
#         fasta_file_to_copy = cluster + ".fa"
#         fasta_file_location = os.path.join(sequenceA, fasta_file_to_copy)
#         os.chdir(sequenceA)
#         split_fasta(fasta_file_location)
#         os.remove(fasta_file_location)
#         os.chdir(source_directory)
# make_individual_fasta_file_and_delete_original_file('/home/supremeleader/first_paper_data/archaea_seq_full_length')


import matplotlib.pyplot as plt
import numpy as np

# Define the range for x
x = np.linspace(-3, 3, 400)

# Define the functions
y1 = 2 * x**2
y2 = x**3

# Plot the functions
plt.figure(figsize=(10, 6))
plt.plot(x, y1, label=r'$y = 2x^2$')
plt.plot(x, y2, label=r'$y = x^3$')

# Adding labels and title
plt.xlabel('x')
plt.ylabel('y')
plt.title('Plot of $y = 2x^2$ and $y = x^3$')
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.grid(color='gray', linestyle='--', linewidth=0.5)
plt.legend()

# Show the plot
plt.show()



















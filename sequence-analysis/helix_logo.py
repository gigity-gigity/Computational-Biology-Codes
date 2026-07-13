import os
import pandas as pd
import shutil

##################################### Make cluster folder #############################################################
current_working_directory=os.getcwd()
def make_directory_for_sequence_logo(path_to_id_folder):
    files=os.listdir(path_to_id_folder)
    os.makedirs("eukaryotes_helix_files")
    for file in files:
        cluster_name=file.split('_')[0]
        directory_name=os.path.join(current_working_directory,"eukaryotes_helix_files")
        cluster_folder=os.path.join(directory_name,cluster_name)
        os.makedirs(cluster_folder)
make_directory_for_sequence_logo("/home/supremeleader/PycharmProjects/Programming/pairwise_alignment/folder")
##################################### Give path to this function ######################################################

##################################### Get helices for helix logo ######################################################
# fasta_sequence = df.loc[df['fasta_id'].str.contains(itol_ids, case=False), "sequence"].dropna()
def get_helix_sequence_from_cluster_id(path_to_cluster_id):
    df = pd.read_csv("/home/supremeleader/PycharmProjects/Programming/pairwise_alignment/database_eukaryote_data.csv")
    helices=['TM1','TM2','TM3','TM4','TM5','TM6','TM7'] #have to name this list for desired helices
    sequence={}
    with open(path_to_cluster_id)as handle:
        for helix in helices:
            handle.seek(0)
            for i in handle:
                entries=i.strip().split("_")[-1]
                helix_data = df.loc[df['fasta_id'].str.contains(entries, case=False), helix].dropna()
                sequence[entries]=helix_data.iloc[0] if not helix_data.empty else None
            with open (f'{helix}_helix_sequences.fa','w')as lines:
                for key,value in sequence.items():
                    write_lines=f'>{key}\n{value}\n'
                    lines.write(write_lines)
            sequence={}
# get_helix_sequence_from_cluster_id("/home/supremeleader/PycharmProjects/Programming/pairwise_alignment/cluster3_id.txt")
############################# No need to give path to this function it will be given below ################################

################################# Get individual helix sequence for each cluster ##########################################

def get_helix_file_for_each_cluster(path_to_id_containing_folder):
    path_to_cluster_files=os.listdir(path_to_id_containing_folder)
    eukaryotes_helix_files=os.path.join(current_working_directory,'eukaryotes_helix_files')
    for cluster in path_to_cluster_files:
        path_of_cluster_id_files=os.path.join(path_to_id_containing_folder,cluster)
        cluster_directory=cluster.split("_")[0]
        path_to_each_cluster=os.path.join(eukaryotes_helix_files,cluster_directory)
        os.chdir(path_to_each_cluster)
        get_helix_sequence_from_cluster_id(path_of_cluster_id_files)
        os.chdir(current_working_directory)
get_helix_file_for_each_cluster("/home/supremeleader/PycharmProjects/Programming/pairwise_alignment/folder")
############################# Give path to this function  #################################################################


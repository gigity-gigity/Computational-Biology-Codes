import pandas as pd
import shutil
import os
import matplotlib.pyplot as plt

################################ Make cluster specific analysis folder (No need to run it) ######################
working_directory=os.getcwd()
def make_cluster_specific_analysis_folder(itol_id_folder_path):
    os.chdir(itol_id_folder_path)
    files=os.listdir(itol_id_folder_path)
    for cluster_ids in files:
        cluster=cluster_ids.split("_")
        os.mkdir(cluster[0])

# make_cluster_specific_analysis_folder('/home/supremeleader/PycharmProjects/Programming/pairwise_alignment/folder')
####################################### No need to use this function ##################################################

################################# function to get sequence file using cluster_ids################################
###########################################(No need to use this function)###################
def get_fasta_seq_of_phylogenetic_cluster(phylogeny_cluster):
    sequence={}
    cluster=phylogeny_cluster.split("_")
    cluster_name=cluster[0]
    current_directory=os.getcwd()
    folder_directory=os.path.join(current_directory,'folder')
    phylogeny_cluster_name=os.path.join(folder_directory,phylogeny_cluster)
    df=pd.read_csv("database_eukaryote_data.csv")
    with open(phylogeny_cluster_name)as handle:
        for ids in handle:
            itol_ids = ids.strip()
            fasta_sequence = df.loc[df['fasta_id'].str.contains(itol_ids, case=False), "sequence"].dropna()
            if not fasta_sequence.empty:
                fasta_sequence = fasta_sequence.iloc[0]
                sequence[itol_ids] = fasta_sequence
    output_filename = f"{cluster_name}.fa"
    if os.path.isdir('fasta_sequences')==False:
        os.mkdir("fasta_sequences")
    os.chdir("fasta_sequences")
    with open(output_filename, 'w') as f:
        for key, value in sequence.items():
            output_line = f">{key}\n{value}\n"
            f.write(output_line)
    os.chdir("..")

# get_fasta_seq_of_phylogenetic_cluster('cluster3_id.txt')
# ############################## No need to use this  ###########################################################
##################################################################################################################

############################# Make fasta sequence file for every cluster ##########################################
###########################################(have to use this function)###################
def fasta_sequence_file_for_every_cluster(itol_id_folder_path):
    os.chdir(itol_id_folder_path)
    os.system('pwd')
    files = os.listdir(".")
    os.chdir("..")
    for file in files:
        get_fasta_seq_of_phylogenetic_cluster(file)

fasta_sequence_file_for_every_cluster("/home/supremeleader/PycharmProjects/Programming/pairwise_alignment/folder")
#################################### have to use this as well###############################
##################################################################################################################

################# Distribute the sequences to the specific folder ################################################
###########################################(have to use this function)###################
def copy_fasta_sequences(id_folder_name):
    files = os.listdir(id_folder_name)
    source_file="/home/supremeleader/PycharmProjects/Programming/pairwise_alignment/fasta_sequences"
    cluster_folder=os.getcwd()
    for file in files:
        folder_name=file.split("_")[0]
        # print(folder_name[0])
        output_folder=os.path.join(cluster_folder,folder_name)
        if not os.path.exists(output_folder):
            shutil.copytree(source_file,output_folder)

copy_fasta_sequences('/home/supremeleader/PycharmProjects/Programming/pairwise_alignment/folder')
#################################### have to use this as well###############################
#################################################################################################################

##################### Make folder for sequenceB files and copy them there also make result folder ###############
###########################################(have to use this function)###################
def make_seqB_folder_and_result_folder(id_folder_name):
    files = os.listdir(id_folder_name)
    path_to_current_folder= os.getcwd()
    for file in files:
        folder_name=file.split("_")[0]
        src_directory=os.path.join(path_to_current_folder,folder_name)
        sequenceb_path=os.path.join(src_directory,'sequenceB')
        try:
            shutil.copytree(src_directory, sequenceb_path)
        except Exception as e:
            print(f"Error during copytree: {e}")
        result_folder = os.path.join(src_directory, "result")
        os.makedirs(result_folder)

make_seqB_folder_and_result_folder('/home/supremeleader/PycharmProjects/Programming/pairwise_alignment/folder')
#################################### have to use this as well###############################
#################################################################################################################

########## make sequenceA folder and copy cluster specific fasta file to the folder #############################
###########################################(have to use this function)###################
def make_seqB_folder_and_copy_cluster_specific_file(fasta_sequence_location):
    files=os.listdir(fasta_sequence_location)
    current_directory=os.getcwd()
    for file in files:
        ids=file.split(".")[0]
        source_directory=os.path.join(current_directory,ids)
        sequenceA=os.path.join(source_directory,'sequenceA')
        os.makedirs(sequenceA)
        fasta_file_to_copy=ids+".fa"
        fasta_file_location=os.path.join(fasta_sequence_location,fasta_file_to_copy)
        shutil.copy(fasta_file_location,sequenceA)

make_seqB_folder_and_copy_cluster_specific_file('/home/supremeleader/PycharmProjects/Programming/pairwise_alignment/fasta_sequences')
#################################### have to use this as well###############################
#################################################################################################################

########################### make individual fasta file from a single large fasta file #################################
 ########################### ( do not give path to this function ######################################################
def split_fasta(input_file):
    with open(input_file, 'r') as infile:
        current_sequence = ''
        for line in infile:
            if line.startswith('>'):
                if current_sequence:
                    write_sequence_to_file(sequence_name, current_sequence)
                sequence_name = line.strip().lstrip('>')
                current_sequence = ''
            else:
                current_sequence += line.strip()
        if current_sequence:
            write_sequence_to_file(sequence_name, current_sequence)

def write_sequence_to_file(sequence_name, sequence):
    output_file = f'{sequence_name}.fa'
    with open(output_file, 'w') as outfile:
        outfile.write(f'>{sequence_name}\n{sequence}')
########################################## do not use this function ######################################################


########################### In sequenceA folder fasta file to make individual fasta sequence and delete the original one #############
###########################################(have to use this function)###################
def make_individual_fasta_file_and_delete_original_file(fasta_sequence_location):
    files=os.listdir(fasta_sequence_location)

    current_directory = os.getcwd()
    for file in files:
        cluster=file.split(".")[0]
        source_directory = os.path.join(current_directory,cluster)
        sequenceA = os.path.join(source_directory, 'sequenceA')
        fasta_file_to_copy = cluster + ".fa"
        fasta_file_location = os.path.join(sequenceA, fasta_file_to_copy)
        os.chdir(sequenceA)
        split_fasta(fasta_file_location)
        os.remove(fasta_file_location)
        os.chdir(source_directory)
make_individual_fasta_file_and_delete_original_file('/home/supremeleader/PycharmProjects/Programming/pairwise_alignment/fasta_sequences')
###################################### Have to give path to this function ##########################################################

################################################################################################################################################
                                  # Making folder for cluster specific pairwise sequence analysis #
################################################################################################################################################

                         #Running pairwise sequence alignment and analysis#
                    ############################################################


################################ Run pairwise sequence alignment ######################################################################
###########################################(have to use this function)###################
def run_pairwise_Sequence_alignment(fasta_sequence_location):
    files=os.listdir(fasta_sequence_location)
    os.chdir(working_directory)
    current_directory = os.getcwd()
    for file in files:
        folder_name = file.split(".")[0]
        directory_for_alignment=os.path.join(current_directory,folder_name)
        sequencea_directory=os.path.join(directory_for_alignment,"sequenceA")
        sequencea_list=os.listdir(sequencea_directory)
        sequenceb_directory=os.path.join(directory_for_alignment,"sequenceB")
        sequenceb_list = os.listdir(sequenceb_directory)
        result_directory=os.path.join(directory_for_alignment,"result")
        for sequencea in sequencea_list:
            sequenceA_name=sequencea.split(".")[0]
            for sequenceb in sequenceb_list:
                sequenceb_name=sequenceb.split(".")[0]
                cmd="needle -asequence " +str(sequencea_directory)+"/"+str(sequencea)+ " -bsequence " +str(sequenceb_directory)+"/"+ str(sequenceb) + " -outfile " + str(result_directory) +"/"+str (sequenceA_name) +"-"+ str(sequenceb_name)+".txt" +" -auto"
                os.system(cmd)
run_pairwise_Sequence_alignment('/home/supremeleader/PycharmProjects/Programming/pairwise_alignment/fasta_sequences')
############################# Have to give fucntion the path of the fasta_sequence_folder#####################################################


                        # Combine alignment result files form each result folder#
                    ################################################################

############################# Combine cluster specific alignment result files ###################################################
###########################################(have to use this function)###########################################################
def combine_alignment_result_file(path_to_fasta_files):
    files=os.listdir(path_to_fasta_files)
    current_directory=os.getcwd()
    file_path=[]
    for file in files:
        cluster_name=file.split(".")[0]
        cluster_path=os.path.join(current_directory,cluster_name)
        result_directory=os.path.join(cluster_path,"result")
        result_files=os.listdir(result_directory)
        for x in files:
            cluster_specific=x.split('.')[0]
            for result_name in result_files:
                if cluster_specific in result_name:
                    file_path.append(os.path.join(result_directory,result_name))
            with open(f'{result_directory}/{cluster_specific}_all_result.txt','w') as output_file:
                for file_path in file_path:
                    with open(file_path,'r')as input_file:
                        for line in input_file:
                            output_file.write(line)
            file_path=[]
combine_alignment_result_file('/home/supremeleader/PycharmProjects/Programming/pairwise_alignment/fasta_sequences')
#################################### Have to give path to this function ###############################################

                       # Calculate Average/SD sequence identity and similarity#
                     ################################################################

####################### Calculate Average/SD of sequence identity and specificity  #################################
###########################################(have to use this function)##############################################
def calculate_average_sd(file_path):
    identity=[]
    similarity=[]
    with open(file_path)as handle:
        for line in handle:
            if line.startswith("# Identity:"):
                if "100.0%" not in line:
                    uncurated=line.split("(")[-1]
                    curated=uncurated.split("%")[0]
                    identity.append(curated)
            if line.startswith("# Similarity:"):
                if "100.0%" not in line:
                    uncurated=line.split("(")[-1]
                    curated=uncurated.split("%")[0]
                    similarity.append(curated)
        identity_score=[eval(i) for i in identity]
        similarity_score=[eval(i)for i in similarity]
        sequence_identity=round(sum(identity_score)/len(identity_score),2)
        similarity_score=round(sum(similarity_score)/len(similarity_score),2)
        scores= str(sequence_identity)+"("+str(similarity_score)+")"
        return scores
# calculate_average_sd("cluster3_all_result.txt") # No need to give path here, function will be used further down

def calculate_Average_and_SD_sequence_identity_similarity_of_clusters(path_to_fasta_files):
    files=os.listdir(path_to_fasta_files)
    current_directory=os.getcwd()
    cluster_alignment_result=[]
    pairwise_alignment_dict={}
    for file in files:
        cluster_folder=file.split(".")[0]
        cluster = os.path.join(current_directory,cluster_folder)
        result_folder = os.path.join(cluster,'result')
        result_files=os.listdir(result_folder)
        for result_file in result_files:
            if result_file.endswith("_all_result.txt"):
                cluster_alignment_result.append(os.path.join(result_folder,result_file))
        for path in cluster_alignment_result:
            name=path.split("/")[-1]
            naming=name.split("_")[0]
            alignment_score_name=cluster_folder+'-'+naming
            alignment_score=calculate_average_sd(path)
            pairwise_alignment_dict[alignment_score_name]=alignment_score
        cluster_alignment_result=[]
    ############################## Make table #####################################
    rows_column=[]
    for clusters in files:
        cluster=clusters.split(".")[0]
        rows_column.append(cluster)
    df=pd.DataFrame(index=rows_column,columns=rows_column)
    for key,value in pairwise_alignment_dict.items():
        cluster1,cluster2=key.split('-')
        crossponding_value=str(value)
        df.at[cluster1,cluster2] = crossponding_value
    df.to_csv('unordered.csv')
    df = pd.read_csv('unordered.csv',index_col=0)
    cluster_numbers = [int(cluster.replace('cluster', '')) for cluster in df.columns]
    desired_order = ['cluster' + str(cluster) for cluster in sorted(cluster_numbers)]
    df_reordered = df.reindex(index=desired_order, columns=desired_order)
    df_reordered.to_csv('eukaryotic_cluster_specific_pairwise_seq_analysis.csv')
    print(df_reordered)
    df = pd.DataFrame(df_reordered)
    fig,ax = plt.subplots(figsize=(20,10))
    ax.axis('off')
    table = ax.table(cellText=df.values, colLabels=df.columns, rowLabels=df.index, loc='center', cellLoc='center', colWidths=[0.15]*len(df.columns))
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    plt.show()
calculate_Average_and_SD_sequence_identity_similarity_of_clusters('/home/supremeleader/PycharmProjects/Programming/pairwise_alignment/fasta_sequences')
#################################### Have to give path to this function ###############################################

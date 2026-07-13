import os
import shutil
path='/home/supremeleader/modelling_code'
folder_name=[]
with open("ids_bacteria.txt")as handle:
	handle1=handle.readlines()
for ids in handle1:
	iwd=ids.strip()
	folder_name.append(iwd)
for folder in folder_name:
	curr_dir= os.getcwd()
	os.mkdir(folder)
	os.chdir(folder)
	###########################file to add in the folder#########################################
	f=open("align2d.py",'w') #create align2d.py file here u want to add in the folder
	model1=open("model-single.py",'w') #create model-single.py file here u want to add in the folder

	#############################################################################################
	shutil.copy(curr_dir+"/5wue.pdb",".")
	shutil.copy(curr_dir+"/open_state"+"/"+folder+".ali",".")

	########################### making align2d.py file ##########################################
	targetfilename= folder
	align_code="""from modeller import *
env = Environ()
aln = Alignment(env)
mdl = Model(env, file='5wue', model_segment=('FIRST:A','LAST:A'))
aln.append_model(mdl, align_codes='5wueA', atom_files='5wue.pdb')
aln.append(file='"""+targetfilename+""".ali', align_codes='"""+targetfilename+"""')
aln.align2d(max_gap_length=50)
aln.write(file='"""+targetfilename+"""-5wueA.ali', alignment_format='PIR')
aln.write(file='"""+targetfilename+"""-5wue.pap', alignment_format='PAP')
"""
    ############################ making single-model.py file #####################################
	
	model="""from modeller import *
from modeller.automodel import *
env = Environ()
a = AutoModel(env, alnfile='"""+targetfilename +"""-5wueA.ali',knowns='5wueA', sequence='"""+targetfilename +"""',assess_methods=(assess.DOPE,assess.GA341))
a.starting_model = 1
a.ending_model = 5
a.make()
"""
    #################### writing and closing the files #############################################
	f.write("%s\n"%align_code)
	f.close()

	model1.write("%s\n"%model)
	model1.close()

	os.chdir(curr_dir)
#print(os.listdir('/home/supremeleader/modelling_code')) 
#############################################################################################################
# using read() to work with file
# with open ("A0A4D6HMN6.fas") as handle:
# 	handle1=handle.read()
# name=handle1.split("\n")[1][0:]
# sequence=''.join (handle1.split('\n')[1:])
# print("id={0}".format(name))
# print("seq={0}".format(sequence))

#########################################################

#using readline() to loop through the lines in a file format
# sequence =''
# with open("A0A4D6HMN6.fas") as fh:
# 	name=fh.readline()[1:-1]
# 	for line in fh:
# 		sequence += line.replace("\n","")
# print("id = {0}".format(name))
# print("seq = {0}".format(sequence))

#########################################################
#calculating charge of a protein
# sequence = ''
# charge=-0.002
# aa_charge={'C':-0.45,'D': -0.999, 'E':-0.998, 'H': -0.91, 'K':1, 'R':1, 'Y': -0.01 }
# with open("A0A4D6HMN6.fas")as fh:
# 	fh.readline()
# 	for line in fh:
# 		sequence += line[:-1].upper()
# print(sequence)
# for aa in sequence:
# 	charge += aa_charge.get(aa,0)
# print(charge)
############################################################
#writing into the file
# sequence = ''
# charge=-0.002
# aa_charge={'C':-0.45,'D': -0.999, 'E':-0.998, 'H': -0.91, 'K':1, 'R':1, 'Y': -0.01 }
# with open("A0A4D6HMN6.fas")as fh:
# 	fh.readline()
# 	for line in fh:
# 		sequence += line[:-1].upper()
# for aa in sequence:
# 	charge += aa_charge.get(aa,0)
# with open("charge.txt","w")as handle:
# 	handle.write(str(charge))
###############################################################
#working with csv file
# total_len=0
# with open("sample.csv")as fh:
# 	next(fh)
# 	for n,line in enumerate(fh):
# 		data=line.split(",")
# 		total_len += int(data[1])
# print(n)
# print(total_len/n)
##############################################################
#using csv module
# import csv
# total_len=0
# lines=csv.reader(open("sample.csv"))
# next(lines)
# for n,line in enumerate(lines):
# 	total_len += int(line[1])
# print(total_len/n)
##############################################################
#import csv
# data=list(csv.reader(open("sample.csv")))
# print(data[0][2])
# print(data[1][1])
##############################################################
# import os
# cwd= os.getcwd()
# print(cwd)
# os.chdir('/home/supremeleader/Desktop')
# cwd= os.getcwd()
# print(cwd)
# f=os.listdir('/home/supremeleader/Desktop')
# print(f)
# print(os.path.isfile("5wucout.pdb"))
# os.chdir('/home/supremeleader/')
# #os.rename('/home/supremeleader/5wuc.pdb','/home/supremeleader/Desktop/5wucout.pdb')
# os.chdir('/home/supremeleader/modelling_code')
# print(os.listdir('/home/supremeleader/modelling_code'))
# #os.mkdir('/home/supremeleader/modelling_code/modle_folder1')
# print(os.path.join('/home/supremeleader/modelling_code/model_folder','/home/supremeleader/modelling_code/model_folder1'))
# print(os.path.exists(os.path.join(os.getcwd(),'sample.csv')))
#############################################################################################################################
#creating multiple folder with name from list
import os
import shutil
path='/home/supremeleader/modelling_code'
folder_name=['prashant','reddy', 'pranav']
ids=["prashant.ali","reddy.ali","pranav.ali"]
for folder in folder_name:
	curr_dir= os.getcwd()
	os.mkdir(folder)
	os.chdir(folder)
	f=open("align.py",'w')
	shutil.copy(curr_dir+"/5wuc.pdb",".")
	shutil.copy(curr_dir+"/alignfiles"+"/"+folder+".ali",".")
	###########################align.py making #####################################
	inputpdb="'5wuc'"
	inputpdbchains= "'5wuc'"
	inputpdbfilename= "'5wuc'"
	targetalifile= "'5wuc'"
	targetid= "'5wuc'"
	outputpir= "'5wuc'"
	outputpap= "'5wuc'"
	align_code="""from modeller import *
env = Environ()
aln = Alignment(env)
mdl = Model(env, file="""+inputpdb+""", model_segment=('FIRST:A','LAST:A'))
aln.append_model(mdl, align_codes="""+inputpdbchains+""", atom_files="""+inputpdbfilename+""")
aln.append(file="""+targetalifile+""", align_codes="""+targetid+""")
aln.align2d(max_gap_length=50)
aln.write(file="""+outputpir+""", alignment_format='PIR')
aln.write(file="""+outputpap+""", alignment_format='PAP')
"""
	f.write("%s\n"%align_code)
	f.close()
	#################################################################
	#modelling=""" """
	os.chdir(curr_dir)
print(os.listdir('/home/supremeleader/modelling_code')) 


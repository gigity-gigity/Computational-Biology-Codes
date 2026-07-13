import os
import sys
p=[]
name_dir=os.listdir("/home/supremeleader/modelling_code/bacteria_model")
os.chdir("/home/supremeleader/modelling_code/bacteria_model")
for foldername in name_dir:
	os.chdir(foldername)
	files=os.listdir('.')
	for k in files:
		if k.endswith(".pap"):
			with open(k)as handle:
				handle1=handle.readlines()
			for s in handle1:
				p.append(s)
			p.append("\n###############################################################################################\n\n")
	os.chdir("/home/supremeleader/modelling_code/bacteria_model")
os.chdir("/home/supremeleader/modelling_code")
with open("alignment.txt","w")as f:
	for items in p:
		f.writelines('%s'%items)
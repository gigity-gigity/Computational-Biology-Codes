import os
import sys
p=[]
name_dir=os.listdir("/home/supremeleader/modelling_code/open_state_bacteria")
os.chdir("/home/supremeleader/modelling_code/open_state_bacteria")
for foldername in name_dir:
	os.chdir(foldername)
	files=os.listdir('.')
	for k in files:
		if k.endswith("single.log"):
			with open(k)as handle:
				handle1=handle.readlines()
			for s in handle1:
				if ".pdb" in s:
					p.append(s)

# 				p.append(s)
			p.append("\n###############################################################################################\n\n")
	os.chdir("/home/supremeleader/modelling_code/open_state_bacteria")
os.chdir("/home/supremeleader/modelling_code")
with open("open_model_soap.txt","w")as f:
	for items in p:
		f.writelines('%s'%items)
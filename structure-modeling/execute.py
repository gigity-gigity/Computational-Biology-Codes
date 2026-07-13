import os
import sys
name_dir=os.listdir("/home/supremeleader/modelling_code/open_state_bacteria")
os.chdir("/home/supremeleader/modelling_code/open_state_bacteria")
for foldername in name_dir:
	os.chdir(foldername)
	#os.system("mod10.0 align2d.py" )
	os.system("mod10.0 model-single.py")
	os.chdir("/home/supremeleader/modelling_code/open_state_bacteria")
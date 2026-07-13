starting_number_pot = int(52393)
ending_number_pot = int(52453)
starting_number_cla = int(52454)
ending_number_cla = int(52531)
potassium_command_index_file="echo \"a {}\ndel 0-16\nq\" | /home/supremeleader/softwares/executable/gromacs22/bin/gmx make_ndx -f /home/supremeleader/error/paramater_files/short_run.gro -o /media/supremeleader/Pantera/simulation/2025/5egi/potassium_permeation/pot_{}.ndx"
chloride_command_index_file="echo \"a {}\ndel 0-16\nq\" | /home/supremeleader/softwares/executable/gromacs22/bin/gmx make_ndx -f /home/supremeleader/error/paramater_files/short_run.gro -o /media/supremeleader/Pantera/simulation/2025/5egi/chloride_permeation/chl_{}.ndx"
potassium_com_xvg_file = "/home/supremeleader/softwares/executable/gromacs22/bin/gmx traj -f /home/supremeleader/error/paramater_files/center.xtc -s /home/supremeleader/error/paramater_files/short_run.gro -n /media/supremeleader/Pantera/simulation/2025/5egi/potassium_permeation/pot_{}.ndx -ox /media/supremeleader/Pantera/simulation/2025/5egi/potassium_permeation/pot_{}.xvg -dt 10 -x no -y no -b 0 -e 500000"
chloride_com_xvg_file = "/home/supremeleader/softwares/executable/gromacs22/bin/gmx traj -f /home/supremeleader/error/paramater_files/center.xtc  -s /home/supremeleader/error/paramater_files/short_run.gro -n /media/supremeleader/Pantera/simulation/2025/5egi/chloride_permeation/chl_{}.ndx -ox /media/supremeleader/Pantera/simulation/2025/5egi/chloride_permeation/chl_{}.xvg -dt 10 -x no -y no -b 0 -e 500000"

import sys
import os
import subprocess
potassium_index_file = []
chloride_index_file = []
potassium_xvg_file = []
chloride_xvg_file = []
for i in range (starting_number_pot, ending_number_pot+1):
    command = potassium_command_index_file.format(i,i)
    potassium_index_file.append(command)
for i in range (starting_number_pot, ending_number_pot+1):
    command = potassium_com_xvg_file.format(i,i)
    potassium_xvg_file.append(command)

for i in range (starting_number_cla, ending_number_cla+1):
    command = chloride_command_index_file.format(i,i)
    chloride_index_file.append(command)
for i in range (starting_number_cla, ending_number_cla+1):
    command = chloride_com_xvg_file.format(i,i)
    chloride_xvg_file.append(command)
#
# print(potassium_index_file)
# print(len(potassium_index_file))
# print(chloride_index_file)
# print(len(chloride_index_file))
#
# print(potassium_xvg_file)
# print(len(potassium_xvg_file))
# print(chloride_xvg_file)
# print(len(chloride_xvg_file))

for command in chloride_xvg_file:
    subprocess.run(command,shell=True)
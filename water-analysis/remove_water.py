import MDAnalysis as mda
import numpy as np


def remove_water_around_glycerol(input_gro, output_gro, glycerol_resname='POT', water_resname='TIP3',
                                 cutoff_distance=5):
    u = mda.Universe(input_gro)

    glycerol_atoms = u.select_atoms(f'resname {glycerol_resname}')
    water_atoms = u.select_atoms(f'resname {water_resname}')

    water_molecules_to_remove = set()

    for glycerol_atom in glycerol_atoms:
        for water_atom in water_atoms:
            distance = np.linalg.norm(glycerol_atom.position - water_atom.position)
            if distance <= cutoff_distance:
                water_molecules_to_remove.add(water_atom.resindex)

    u.atoms = u.atoms[np.invert(np.isin(u.atoms.resindices, list(water_molecules_to_remove)))]

    u.atoms.write(output_gro)


remove_water_around_glycerol(
    '/home/supremeleader/6izf_umb/charmm-gui-8305938907/gromacs/step5_input.gro',
    '/home/supremeleader/6izf_umb/charmm-gui-8305938907/gromacs/final_umb.gro')


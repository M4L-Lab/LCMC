"""Console script for geneticsro."""
import sys
import click
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.io import jsonio
from ase.build import sort
from mcdft.mcdft import MCDFT
from mcdft.calculators import vasp_calculator, emt_calculator
import os
import json

def run_mcdft(atoms,N_steps,traj,Temps):
    calc = emt_calculator(atoms)
    mcdft = MCDFT(atoms, calc, N_steps, traj)
    mcdft.build_traj(Temps)
    
 

@click.command()
def main(arg=None):
    with open("mcdft_parameters.json", "r") as jsonfile:
        data = json.load(jsonfile)

    dir = os.getcwd()
    atoms=read(data['structure_file'])
    N_steps=data['N_steps']
    Temps=data['Temps']
    traj= Trajectory(os.path.join(dir,data['traj_file']), "w")
    run_mcdft(atoms,N_steps,traj,Temps)
    return 0

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover

"""Console script for geneticsro."""
import sys
import click
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.io import jsonio
from ase.build import sort
from mcdft.mcdft import MCDFT
from mcdft.calculators import vasp_calculator, emt_calculator
from mcdft.fast_string import FastString
from ase.ga.utilities import get_nndist
import os
import json


def get_fast_string_maker(atoms):
    natoms=atoms.copy()
    natoms=natoms*(3,3,3)
    dm = natoms.get_all_distances(mic=True)
    nndist = get_nndist(natoms, dm) + 0.2
    fast_string_maker=FastString(dm,nndist)
    return fast_string_maker


def run_mcdft(atoms,N_steps,traj,Temps,compare=False, chain_length=10, max_try_for_unique=1000):
    calc = emt_calculator(atoms)
    
    #dm = atoms.get_all_distances(mic=True)
    #nndist = get_nndist(atoms, dm) + 0.2

    fast_string_maker=get_fast_string_maker(atoms)
    
    mcdft = MCDFT(atoms,fast_string_maker, calc, N_steps, traj,compare=compare,chain_length=chain_length,max_try_for_unique=max_try_for_unique)
    mcdft.build_traj(Temps)
    
 

@click.command()
def main(arg=None):
    with open("mcdft_parameters.json", "r") as jsonfile:
        data = json.load(jsonfile)

    dir = os.getcwd()
    atoms=read(data['structure_file'])
    N_steps=data['N_steps']
    Temps=data['Temps']
    compare=data['compare']
    chain_length=data['chain_length']
    max_try_for_unique=data['max_try_for_unique']
    traj= Trajectory(os.path.join(dir,data['traj_file']), "w")
    run_mcdft(atoms,N_steps,traj,Temps,compare=compare, chain_length=chain_length, max_try_for_unique=max_try_for_unique)
    return 0

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover

from ase import Atoms, Atom
from ase.io import read, write
from itertools import combinations
from ase.ga.utilities import get_nnmat
from mcdft.structure_generator import swap_atoms
from mcdft.calculators import vasp_calculator
import random
import math
import numpy as np
from ase.io.trajectory import Trajectory
from ase.calculators.vasp import Vasp
from mcdft.structure_comparator import StructureComparator
import os


class MCDFT:
    def __init__(self, atoms, fast_string_maker, calculator, N, traj,compare=False, chain_length=10, max_try_for_unique=1000):
        self.atoms = atoms
        self.calculator = calculator
        self.N = N
        self.traj = traj
        self.compare=compare
        if self.compare:
            self.chain_length=chain_length
            self.comparator=StructureComparator(fast_string_maker, self.chain_length)
            self.max_try_for_unique=max_try_for_unique

    def monte_carlo(self, dE,Temp):
        kB = 1.0 / 11604.0
        prob = min(math.exp(-dE / (kB * Temp)), 1.0)
        accept = 1
        if prob < random.random():
            accept = 0
        return accept

    def build_traj(self,Temps):
        print('Ok good luck!')
        E0=self.calculator.calculate_energy(0)
        lowest_struc=self.atoms.copy()
        for Temp in Temps:
            print(f"Temp: {Temp}")
            print(f"Initial Energy : {E0}")
            structures=[lowest_struc]
            energy=[E0]
            for i in range(1, self.N):
                if i%10==0:
                    print(f'{100*i/self.N}%')
                
                atoms = swap_atoms(structures[-1].copy())

                if self.compare:
                    q=0
                    while (not self.comparator.is_unique(atoms.copy())) and q<self.max_try_for_unique:
                        atoms = swap_atoms(structures[-1].copy())
                        
                        q+=1
                        
                self.calculator.structure = atoms
                e = self.calculator.calculate_energy(i)

                if e<E0:
                    E0=e
                    lowest_struc=atoms.copy()

                dE = self.calculator.calculate_dE(energy[-1], e)
                accept = self.monte_carlo(dE,Temp)

                atoms.info['accept'] = accept
                if self.traj is not None:
                    self.traj.write(atoms)
                if accept:
                    structures.append(atoms)
                    energy.append(e)
                
            print(f"Lowest Energy: {E0}")

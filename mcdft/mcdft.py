from ase import Atoms, Atom
from ase.io import read, write
from itertools import combinations
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
    def __init__(self, atoms, calculator, N, traj, chain_length=10):
        self.atoms = atoms
        self.calculator = calculator
        self.N = N
        self.traj = traj
        self.chain_length=chain_length
        self.comparator=StructureComparator(self.chain_length)

    def monte_carlo(self, dE,Temp):
        kB = 1.0 / 11604.0
        prob = min(math.exp(-dE / (kB * Temp)), 1.0)
        accept = 1
        if prob < random.random():
            accept = 0
        return accept

    def build_traj(self,Temps):

        E0=self.calculator.calculate_energy(0)
        lowest_struc=self.atoms.copy()
        for Temp in Temps:
            print(f"Temp: {Temp}")
            print(f"Initial Energy : {E0}")
            structures=[lowest_struc]
            energy=[E0]
            for i in range(1, self.N):
                atoms = swap_atoms(structures[-1].copy())
                self.calculator.structure = atoms
                e = self.calculator.calculate_energy(i)

                if e<E0:
                    E0=e
                    lowest_struc=atoms.copy()

                dE = self.calculator.calculate_dE(energy[-1], e)
                accept = self.monte_carlo(dE,Temp)
                atoms.info["accept"] = accept
                if self.traj is not None:
                    self.traj.write(atoms)
                if accept:
                    structures.append(atoms)
                    energy.append(e)
            print(f"Lowest Energy: {E0}")

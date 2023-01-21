import numpy as np

# from ase.ga.utilities import get_nnmat_string
from ase.calculators.emt import EMT


class StructureComparator:
    def __init__(
        self, fast_string_maker, max_size, pair_cor_cum_diff=0.015, pair_cor_max=0.7
    ):
        self.pair_cor_cum_diff = pair_cor_cum_diff
        self.pair_cor_max = pair_cor_max
        self.max_size = max_size
        self.fast_string_maker = fast_string_maker
        self.count = 0
        self.structures = []
        self.nmstring = []

    def push(self, atoms):
        """Keep max_size number of previous structure and their sort distance
        whenever stack reached max capacity oldest"""
        if len(self.structures) >= self.max_size:
            self.structures.pop()
            self.nmstring.pop()
        atoms_nmstring = self.fast_string_maker.get_nnmat_string(
            atoms * (3, 3, 3)
        )  # get_nnmat_string(atoms*(3,3,3),decimals=1,mic=True)
        self.structures = [atoms] + self.structures
        self.nmstring = [atoms_nmstring] + self.nmstring

    def is_unique(self, atoms):
        """Comapare a given structure with chain_length number of previous structure. return
        from latest to oldest return false as soon as it match with any structure.
        finaly return true if does not match with any other structure in the stack"""

        atoms_nmstring = self.fast_string_maker.get_nnmat_string(
            atoms * (3, 3, 3)
        )  # get_nnmat_string(atoms*(3,3,3),decimals=1,mic=True)
        if len(self.structures) < 1:
            self.push(atoms)
            return True
        else:
            calc = EMT()
            atoms.calc = calc
            atoms_energy = np.round(atoms.get_potential_energy(), 3)
            #print(f"S: {atoms_nmstring} E:{atoms_energy}")
            for structure, s_string in zip(self.structures, self.nmstring):
                # not_unique=self.compare_sorted_dist(atoms,structure,atoms_sort_dist, structure_sort_dist)

                # struc_nmstring=get_nnmat_string(structure)
                not_unique = atoms_nmstring == s_string

                structure.calc = calc

                struc_energy = np.round(structure.get_potential_energy(), 3)

                #print(f"S:{s_string} E:{struc_energy}")

                if not_unique:
                    self.count += 1
                    #print(f"match returning false. stopped {self.count} non unique structure")
                    return False
        self.push(atoms)
        #print("match with no structure returning true")
        return True

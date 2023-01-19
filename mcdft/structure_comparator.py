import numpy as np
from ase.ga.utilities import get_nnmat_string


class StructureComparator:
    def __init__(self, max_size, pair_cor_cum_diff=0.015 , pair_cor_max= 0.7):
        self.pair_cor_cum_diff=pair_cor_cum_diff
        self.pair_cor_max= pair_cor_max
        self.max_size = max_size
        self.count=0
        self.structures = []
        self.sort_dists= []
    
    def push(self, atoms):
        """Keep max_size number of previous structure and their sort distance
        whenever stack reached max capacity oldest"""
        if len(self.structures) >= self.max_size:
            self.structures.pop()
            self.sort_dists.pop()
        sort_dist=self.get_sorted_dist_list(atoms, mic=True)
        self.structures=[atoms]+self.structures
        self.sort_dists=[sort_dist]+self.sort_dists

    def get_sorted_dist_list(self, atoms, mic=True):
        """ Utility method used to calculate the sorted distance list
            describing the cluster in atoms. """
        numbers = atoms.numbers
        unique_types = set(numbers)
        pair_cor = dict()
        for n in unique_types:
            i_un = [i for i in range(len(atoms)) if atoms[i].number == n]
            d = []
            for i, n1 in enumerate(i_un):
                for n2 in i_un[i + 1:]:
                    d.append(atoms.get_distance(n1, n2, mic))
            d.sort()
            pair_cor[n] = np.array(d)
        return pair_cor

    def compare_sorted_dist(self,a1,a2,p1, p2):
        """ Private method for calculating the structural difference. """
        numbers = a1.numbers
        total_cum_diff = 0.
        max_diff = 0
        for n in p1.keys():
            cum_diff = 0.
            c1 = p1[n]
            c2 = p2[n]
            assert len(c1) == len(c2)
            if len(c1) == 0:
                continue
            t_size = np.sum(c1)
            d = np.abs(c1 - c2)
            cum_diff = np.sum(d)
            max_diff = np.max(d)
            ntype = float(sum([i == n for i in numbers]))
            total_cum_diff += cum_diff / t_size * ntype / float(len(numbers))
        
        return (total_cum_diff < self.pair_cor_cum_diff
                and max_diff < self.pair_cor_max)
    
    def is_unique(self,atoms):
        """Comapare a given structure with chain_length number of previous structure. return
        from latest to oldest return false as soon as it match with any structure.
        finaly return true if does not match with any other structure in the stack"""
        
        #atoms_sort_dist=self.get_sorted_dist_list(atoms, mic=True)
        
        if len(self.structures)<1:
            
            return True
        else:
            for structure, structure_sort_dist in zip(self.structures, self.sort_dists):
                #not_unique=self.compare_sorted_dist(atoms,structure,atoms_sort_dist, structure_sort_dist)
                atoms_nmstring=get_nnmat_string(atoms)
                struc_nmstring=get_nnmat_string(structure)
                not_unique = atoms_nmstring==struc_nmstring
                print(atoms_nmstring)
                print(struc_nmstring)
                if not_unique:
                    self.count+=1
                    print(f'match with a structure returning false. stopped {self.count} non unique structure')
                    return False
        print('match with no structure returning true')
        return True
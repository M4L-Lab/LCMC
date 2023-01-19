import numpy as np



class StructureComparator:
    def __init__(self, max_size, pair_cor_cum_diff=0.015 , pair_cor_max= 0.7):
        self.pair_cor_cum_diff=pair_cor_cum_diff
        self.pair_cor_max= pair_cor_max
        self.max_size = max_size
        self.structures = []
        self.sort_dists= []
    
    def push(self, atoms, sort_dist):
        """Keep max_size number of previous structure and their sort distance
        whenever stack reached max capacity oldest"""
        if len(self.structures) >= self.max_size:
            self.structures.pop()
            self.sort_dists.pop()
        self.structures.insert(atoms,0)
        self.sort_dists.insert(sort_dist,0)

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
        
        atoms_sort_dist=self.get_sorted_dist_list(atoms, mic=True)
        
        if len(self.structures)<self.max_size:
            self.push(atoms,atoms_sort_dist)
        else:
            for structure, structure_sort_dist in zip(self.structures, self.sort_dists):
                unique=self.compare_sorted_dist(atoms,structure,atoms_sort_dist, structure_sort_dist)
                if unique:
                    return False
            self.push(atoms,atoms_sort_dist)
        return True
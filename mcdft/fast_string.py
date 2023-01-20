from ase.ga.utilities import get_nndist
import numpy as np

# dm = atoms.get_all_distances(mic=mic)
# nndist = get_nndist(atoms, dm) + 0.2

class FastString:
    def __init__(self,dm,nndist,decimals=2,mic=True):
        self.mic=mic
        self.decimals=decimals
        self.dm = dm #atoms.get_all_distances(mic=self.mic)
        self.nndist = nndist #get_nndist(atoms, self.dm) + 0.2



    def get_nnmat(self,atoms):
        
        if 'data' in atoms.info and 'nnmat' in atoms.info['data']:
            return atoms.info['data']['nnmat']
        elements = sorted(set(atoms.get_chemical_symbols()))
        nnmat = np.zeros((len(elements), len(elements)))
        

        for i in range(len(atoms)):
            row = [j for j in range(len(elements))
                if atoms[i].symbol == elements[j]][0]
            neighbors = [j for j in range(len(self.dm[i])) if self.dm[i][j] < self.nndist]
            for n in neighbors:
                column = [j for j in range(len(elements))
                        if atoms[n].symbol == elements[j]][0]
                nnmat[row][column] += 1
        # divide by the number of that type of atoms in the structure
        for i, el in enumerate(elements):
            nnmat[i] /= len([j for j in range(len(atoms))
                            if atoms[int(j)].symbol == el])
        # makes a single list out of a list of lists
        nnlist = np.reshape(nnmat, (len(nnmat)**2))
        return nnlist
    
    def get_nnmat_string(self, atoms):
        nnmat = self.get_nnmat(atoms)
        s = '-'.join(['{1:2.{0}f}'.format(self.decimals, i)
                    for i in nnmat])
        if len(nnmat) == 1:
            return s + '-'
        return s